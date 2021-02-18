#!/bin/bash -l

# slurm settings if called using sbatch
#SBATCH -J GTDBtk_METABOLIC
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=2-00:00:00
#SBATCH -p batch
#SBATCH --qos=qos-batch

# conda env name or path
CENV="snakemake" # TODO
# config files
SMK_CONFIG="extraction_config.yaml"
# SMK_SLURM="slurm.yaml"
# slurm cluster call
# SMK_CLUSTER="sbatch -p {cluster.partition} -q {cluster.qos} {cluster.explicit} -N {cluster.nodes} -n {cluster.n} -c {threads} -t {cluster.time} --job-name={cluster.job-name}"

conda activate ${CENV} && \
# snakemake --jobs 10 --local-cores 1 --configfile ${SMK_CONFIG} --use-conda --cluster-config ${SMK_SLURM} --cluster "${SMK_CLUSTER}" -rp
snakemake --jobs 50 --configfile ${SMK_CONFIG} --use-conda -s extraction_Snakefile -rp
