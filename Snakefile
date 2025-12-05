# Created By: Sharon Colson
# Creation Date: 12/01/2025
# Last Modified: 12/05/2025

# To run this on TTU HPC:
#     spack load trf snakemake graphviz
# sample run line:
#     snakemake -np --config sample=Guy11_Final_S4
#         OR
#     snakemake --cores 12
# Create DAG: snakemake --dag | dot -Tsvg > test.svg

configfile: "config.yaml"

MATCH     = config["trf_params"]["match"]
MISMATCH  = config["trf_params"]["mismatch"]
DELTA     = config["trf_params"]["delta"]
PM        = config["trf_params"]["pm"]
PI        = config["trf_params"]["pi"]
MINSCORE  = config["trf_params"]["minscore"]
MAXPERIOD = config["trf_params"]["maxperiod"]
OPTIONS   = config["trf_params"]["options"]

SAMPLE_TO_FASTA = config["samples"]
SAMPLES_LIST = list(SAMPLE_TO_FASTA.keys())

def get_trf_fastas(wildcards):
    return SAMPLE_TO_FASTA[wildcards.sample]

# Order for TRF and filename suffix
TRF_NUMERIC_VALUES = [MATCH, MISMATCH, DELTA, PM, PI, MINSCORE, MAXPERIOD]

# Add options for the TRF parameter string
TRF_PARAM_STRING = " ".join(str(num) for num in TRF_NUMERIC_VALUES)
if OPTIONS:
    TRF_PARAM_STRING += f" {OPTIONS}"

# Build the .dat file name suffix, e.g. ".2.7.7.80.10.50.dat"
TRF_SUFFIX = "." + ".".join(str(num) for num in TRF_NUMERIC_VALUES) + ".dat"

rule all:
    input:
        expand("results/{sample}_trf.bed", sample=SAMPLES_LIST)

rule trf1:
    input:
        get_trf_fastas
    output:
        "results/{sample}.fasta" + TRF_SUFFIX
    log:
        "logs/trf1_{sample}.log"
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)

        # Set the work directory for TRF to "results/"
        results_dir = os.path.dirname(output[0])

        # Show TRF where the input will be from "results/"
        input_rel = os.path.relpath(input[0], results_dir)

        trf_command = f"trf {input_rel} {TRF_PARAM_STRING}"

        # This code modified from https://stackoverflow.com/questions/45613881/what-would-be-an-elegant-way-of-preventing-snakemake-from-failing-upon-shell-r-e
        try:
            proc_output = subprocess.check_output(trf_command, shell=True, cwd=results_dir, stderr=subprocess.STDOUT)

            # Log what TRF is doing
            with open(log[0], "wb") as lf:
                lf.write(proc_output)

        # an exception is raised by check_output() for non-zero exit codes (usually returned to indicate failure)
        except subprocess.CalledProcessError as exc:
            # Capture errors in the log as well
            with open(log[0], "wb") as lf:
                if exc.output:
                    lf.write(exc.output)
                lf.write(f"\nTRF exit code: {exc.returncode}\n".encode())

            # Current releases of TRF (up to 4.09.1) return an exit code of the number of TRs processed
            # Anything less than 3 should be treated as an error code. Otherwise, this should pass
            if exc.returncode is not None and exc.returncode >= 3:
                pass
            else:
                raise

        # Make certain that output file is created
        if not os.path.exists(output[0]):
            raise Exception(f"TRF did not produce expected output file: {output[0]}")

rule trf2:
    input:
        "results/{sample}.fasta" + TRF_SUFFIX
    output:
        "results/{sample}_trf.bed"
    log:
        "logs/trf2_{sample}.log"
    shell:
        r"""
        mkdir -p $(dirname {output}) $(dirname {log})
        python3 trf2bed.py \
            --dat {input} \
            --bed {output} \
            --tool repeatseq &> {log}
        """
