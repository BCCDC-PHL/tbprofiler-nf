process fastp {

  tag { sample_id }

  publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastp.{json,csv}", mode: 'copy'

  input:
  tuple val(sample_id), path(reads_1), path(reads_2)

  output:
  tuple val(sample_id), path("${sample_id}_fastp.json"), emit: fastp_json
  tuple val(sample_id), path("${sample_id}_fastp.csv"), emit: fastp_csv
  tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: reads
  tuple val(sample_id), path("${sample_id}_fastp_provenance.yml"), emit: provenance

  script:
  """
  printf -- "- process_name: fastp\\n" > ${sample_id}_fastp_provenance.yml
  printf -- "  tool_name: fastp\\n  tool_version: \$(fastp --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_fastp_provenance.yml

  fastp \
    --cut_tail \
    -i ${reads_1} \
    -I ${reads_2} \
    -o ${sample_id}_trimmed_R1.fastq.gz \
    -O ${sample_id}_trimmed_R2.fastq.gz

  mv fastp.json ${sample_id}_fastp.json
  fastp_json_to_csv.py -s ${sample_id} ${sample_id}_fastp.json > ${sample_id}_fastp.csv
  """
}


process tbprofiler {

    tag { sample_id }
    
    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_tbprofiler*.{csv,json}"
    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_tbprofiler*.{bam,bam.bai,vcf}", enabled: !params.rename_ref

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_tbprofiler*.{json,csv}"), emit: reports
    tuple val(sample_id), path("${sample_id}_tbprofiler*.{bam,bam.bai}"), emit: alignment
    tuple val(sample_id), path("${sample_id}_tbprofiler_targets.vcf"), emit: targets_vcf
    tuple val(sample_id), path("${sample_id}_tbprofiler_whole_genome.vcf"), emit: whole_genome_vcf
    tuple val(sample_id), path("${sample_id}_tbprofiler_provenance.yml"), emit: provenance
    
    script:
    """

    tb-profiler profile \
      --threads ${task.cpus} \
      --platform ${params.platform} \
      --mapper ${params.mapper} \
      --caller ${params.caller} \
      --min_depth ${params.min_depth} \
      --af ${params.min_af_used_for_calling} \
      --reporting_af ${params.min_af_used_for_prediction} \
      --read1 ${reads_1} \
      --read2 ${reads_2} \
      --prefix ${sample_id} \
      --csv \
      --call_whole_genome

    mv bam/${sample_id}.bam ./${sample_id}_tbprofiler.bam
    mv bam/${sample_id}.bam.bai ./${sample_id}_tbprofiler.bam.bai

    mv vcf/${sample_id}.targets.csq.vcf.gz ./${sample_id}_tbprofiler_targets.vcf.gz
    gunzip ./${sample_id}_tbprofiler_targets.vcf.gz

    mv vcf/${sample_id}.vcf.gz ./${sample_id}_tbprofiler_whole_genome.vcf.gz
    gunzip ./${sample_id}_tbprofiler_whole_genome.vcf.gz

    cp results/${sample_id}.results.csv ${sample_id}_tbprofiler_full_report.csv
    cp results/${sample_id}.results.json ${sample_id}_tbprofiler_full_report.json

    split_tbprofiler_csv.py -p ${sample_id} -s ${sample_id} ${sample_id}_tbprofiler_full_report.csv

    printf -- "- process_name: tb-profiler\\n" > ${sample_id}_tbprofiler_provenance.yml
    printf -- "  tool_name: tb-profiler\\n  tool_version: \$(tb-profiler profile --version 2>&1 | cut -d ' ' -f 3)\\n" >> ${sample_id}_tbprofiler_provenance.yml
    printf -- "  database_version: \$(grep 'Database version' ${sample_id}_tbprofiler_full_report.csv | cut -d',' -f2)\\n" >> ${sample_id}_tbprofiler_provenance.yml

    """
}

process rename_ref_in_alignment {

  tag { sample_id }
  
  publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_tbprofiler*.{bam,bam.bai}", enabled: params.rename_ref, saveAs: { x -> x.replace("_renamed", "") }
  
  input:
  tuple val(sample_id), path(alignment)

  output:
  tuple val(sample_id), path("${sample_id}_tbprofiler*.{bam,bam.bai}")

  script:
  """
  samtools view -h ${alignment[0]} | sed 's/Chromosome/${params.ref_name}/g' | samtools sort -o ${sample_id}_tbprofiler_renamed.bam
  samtools index ${sample_id}_tbprofiler_renamed.bam
  """
}

process rename_ref_in_variants {

  tag { sample_id }

  publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_tbprofiler*.vcf", enabled: params.rename_ref, saveAs: { x -> x.replace("_renamed", "") }

  input:
  tuple val(sample_id), path(variants)

  output:
  tuple val(sample_id), path("${sample_id}_tbprofiler*.vcf")

  script:
  """
  sed 's/Chromosome/${params.ref_name}/g' ${variants} > ${sample_id}_tbprofiler_renamed.vcf
  """
}

process qualimap_bamqc {

    tag { sample_id }

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_qualimap_alignment_qc.csv"

    input:
    tuple val(sample_id), file(alignment)

    output:
    tuple val(sample_id), path("${sample_id}_qualimap_alignment_qc.csv"), emit: genome_results
    
    script:
    """
    qualimap bamqc -bam ${alignment[0]} --outdir ${sample_id}_bamqc
    qualimap_bamqc_genome_results_to_csv.py -s ${sample_id} ${sample_id}_bamqc/genome_results.txt > ${sample_id}_qualimap_alignment_qc.csv
    """
}


process snpit {

    tag { sample_id }

    conda "$baseDir/environments/snpit.yml"

    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_snpit.tsv"

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}_snpit.tsv")
    
    script:
    """
    snpit --input ${vcf} > ${sample_id}_snpit.tsv
    """
}
