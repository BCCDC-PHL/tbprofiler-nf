#!/bin/bash

mkdir -p .github/data/assemblies

curl -o .github/data/assemblies/NC_000962.3.fa "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=text&id=NC_000962.3&db=nucleotide&rettype=fasta"
