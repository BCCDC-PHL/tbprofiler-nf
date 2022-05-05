# TBProfiler-NF: A Nextflow Wrapper for the TBProfiler Mycobacteria Genomics Analysis Tool

This is a Nextflow-based wrapper for the [TBProfiler](https://github.com/jodyphelan/TBProfiler) pipeline for antimicrobial resistance (AMR)
analysis of whole-genome sequence data for Mycobacteria tuberculosis complex samples. This pipeline provides a convenient way to analyze
many samples at once, on an HPC cluster system. It also integrates additional QC analysis of the input data using [fastp](https://github.com/OpenGene/fastp),
and of the generated alignments using [Qualimap](https://github.com/scchess/Qualimap).

## Usage

```
nextflow run BCCDC-PHL/tbprofiler-nf \
  --fastq_input <fastq_input_dir> \
  --outdir <output_dir>
```

### Parameters

| Flag                         | Default Value | Description                                       |
|:-----------------------------|--------------:|:--------------------------------------------------|
| `min_depth`                  |            10 | Minimum depth of coverage used to call a SNP      |
| `min_af_used_for_calling`    |           0.1 | Minimum minor allele fraction used to call a SNP  |
| `min_af_used_for_prediction` |           0.1 | Minimum minor allele fraction used to predict AMR |


## Output

An output directory will be created for each sample under the directory provided for the `--outdir` flag.

The following files will be produced for each sample:

```
.
└── sample-01
    ├── sample-01_fastp.csv
    ├── sample-01_fastp.json
    ├── sample-01_qualimap_alignment_qc.csv
    ├── sample-01_tbprofiler.bam
    ├── sample-01_tbprofiler.bam.bai
    ├── sample-01_tbprofiler_full_report.csv
    ├── sample-01_tbprofiler_full_report.json
    ├── sample-01_tbprofiler_lineage.csv
    ├── sample-01_tbprofiler_resistance.csv
    ├── sample-01_tbprofiler_summary.csv
    └── sample-01_tbprofiler.vcf
```