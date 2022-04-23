process fastp {

  tag { sample_id }

  publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastp.{json,csv}", mode: 'copy'

  input:
  tuple val(sample_id), path(reads_1), path(reads_2)

  output:
  tuple val(sample_id), path("${sample_id}_fastp.json"), emit: fastp_json
  tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: reads

  script:
  """
  fastp -i ${reads_1} -I ${reads_2} -o ${sample_id}_trimmed_R1.fastq.gz -O ${sample_id}_trimmed_R2.fastq.gz
  mv fastp.json ${sample_id}_fastp.json
  fastp_json_to_csv.py -s ${sample_id} ${sample_id}_fastp.json > ${sample_id}_fastp.csv
  """
}


process tbprofiler {

    tag { sample_id }
    
    publishDir params.versioned_outdir ? "${params.outdir}/${sample_id}/${params.pipeline_short_name}-v${params.pipeline_minor_version}-output" : "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_tbprofiler*.{csv,json}"

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_tbprofiler*.{json,csv}"), emit: reports
    
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
      --csv

    cp results/${sample_id}.results.csv ${sample_id}_tbprofiler_full_report.csv
    cp results/${sample_id}.results.json ${sample_id}_tbprofiler_full_report.json
    split_tbprofiler_csv.py -p ${sample_id} -s ${sample_id} ${sample_id}_tbprofiler_full_report.csv
    """
}


process snp_it {

    tag { sample_id }
    
    cpus 8

    publishDir "${params.outdir}", mode: 'copy', pattern: "${sample_id}_snpit.txt"	

    input:
    file(vcf)

    output:
    tuple val(sample_id), path("${sample_id}_snpit.txt")
    
    script:
    """
    snpit-run.py --input ${vcf} > ${sample_id}_snpit.txt
    """
}