#!/bin/bash


source ${HOME}/.bashrc

eval "$(conda shell.bash hook)"

conda activate art

mkdir -p .github/data/fastq

art_illumina \
    --paired \
    --in .github/data/assemblies/NC_000962.3.fa \
    --fcov 12 \
    --len 150 \
    --mflen 400 \
    --sdev 100 \
    --rndSeed 42 \
    --qShift 0 \
    --qShift2 0 \
    --out .github/data/fastq/NC000962.3_R

mv .github/data/fastq/NC000962.3_R1.fq .github/data/fastq/NC000962.3_R1.fastq
mv .github/data/fastq/NC000962.3_R2.fq .github/data/fastq/NC000962.3_R2.fastq

gzip .github/data/fastq/NC000962.3_R1.fastq
gzip .github/data/fastq/NC000962.3_R2.fastq
