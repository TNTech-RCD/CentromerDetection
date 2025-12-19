#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --time=6-00:00:00
#SBATCH --account=phd-ssalimi42
#SBATCH --partition=batch-impulse

set -e  # Stop on error


spack load samtools/f6zgopc py-numpy/5kbd5ic

input_1=$1 #name of genome
#input_2=$2 #long read change the name to genome name like Fo47.fasta (genoem) and Fo47.subreads.bam (long read)
WORKDIR=/home/tntech.edu/ssalimi42/work/Thesis_Finall/sterategy_2/centromere/ncor74a/bed_meth


spack load py-virtualenv
virtualenv kineticTools_env
. kineticTools_env/bin/activate

# 1. call hifi reads with kinetics if needed
# should have added pbccs to $PATH or the used environment
ccsmeth call_hifi --subreads ${WORKDIR}/${input_1}.subreads.bam \
  --threads ${SLURM_CPUS_PER_TASK} \
  --output ${input_1}.hifi.bam


# 2. align hifi reads
# should have added pbmm2 to $PATH or the used environment
ccsmeth align_hifi \
  --hifireads ${input_1}.hifi.bam \
  --ref ${input_1}.fasta\
  --output ${input_1}.hifi.pbmm2.bam \
  --threads ${SLURM_CPUS_PER_TASK}


# 3. call modifications
# output: [--output].modbam.bam
CUDA_VISIBLE_DEVICES=0 ccsmeth call_mods \
  --input ${input_1}.hifi.pbmm2.bam \
  --ref ${input_1}.fasta \
  --model_file /home/tntech.edu/ssalimi42/work/Thesis_Finall/ccsmeth/models/model_ccsmeth_5mCpG_call_mods_attbigru2s_b21.v3.ckpt \
  --output ${input_1}.hifi.pbmm2.call_mods \
  --threads ${SLURM_CPUS_PER_TASK} --threads_call 2 --model_type attbigru2s \
  --mode align


# 4. call modification frequency
# outputs: [--output].[--call_mode].all.bed
# if the input bam file contains haplotags, 
# there will be [--output].[--call_mode].[hp1/hp2].bed in outputs.
# use '--call_mode count':
ccsmeth call_freqb \
  --input_bam ${input_1}.hifi.pbmm2.call_mods.modbam.bam \
  --ref ${input_1}.fasta \
  --output ${input_1}hifi.pbmm2.call_mods.modbam.freq \
  --threads ${SLURM_CPUS_PER_TASK} --sort --bed

# OR, use '--call_mode aggregate':
# NOTE: usually is more accurate than 'count' mode
ccsmeth call_freqb \
  --input_bam ${input_1}.hifi.pbmm2.call_mods.modbam.bam \
  --ref ${input_1}.fasta \
  --output ${input_1}.hifi.pbmm2.call_mods.modbam.freq \
  --threads ${SLURM_CPUS_PER_TASK} --sort --bed \
  --call_mode aggregate \
  --aggre_model /home/tntech.edu/ssalimi42/work/Thesis_Finall/ccsmeth/models/model_ccsmeth_5mCpG_aggregate_attbigru_b11.v2p.ckpt
