process fastp {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastp.{json,csv}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}_fastp.csv"), emit: csv
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}_fastp_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: fastp\\n"  >> ${sample_id}_fastp_provenance.yml
    printf -- "  tools:\\n"               >> ${sample_id}_fastp_provenance.yml
    printf -- "    - tool_name: fastp\\n" >> ${sample_id}_fastp_provenance.yml
    printf -- "      tool_version: \$(fastp --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_fastp_provenance.yml
    printf -- "      parameters:\\n"               >> ${sample_id}_fastp_provenance.yml
    printf -- "        - parameter: --cut_tail\\n" >> ${sample_id}_fastp_provenance.yml
    printf -- "          value: null\\n"           >> ${sample_id}_fastp_provenance.yml

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
    
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_tbprofiler*.{csv,json}"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_tbprofiler*.{bam,bam.bai,vcf}", enabled: !params.rename_ref

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_tbprofiler*.json"),                    emit: report_json
    tuple val(sample_id), path("${sample_id}_tbprofiler_full_report.csv"),          emit: full_report_csv
    tuple val(sample_id), path("${sample_id}_tbprofiler_lineage.csv"),              emit: lineage_csv
    tuple val(sample_id), path("${sample_id}_tbprofiler_resistance.csv"),           emit: resistance_csv
    tuple val(sample_id), path("${sample_id}_tbprofiler_resistance_mutations.csv"), emit: resistance_mutations_csv
    tuple val(sample_id), path("${sample_id}_tbprofiler_summary.csv"),              emit: summary_csv
    tuple val(sample_id), path("${sample_id}_tbprofiler*.{bam,bam.bai}"),           emit: alignment
    tuple val(sample_id), path("${sample_id}_tbprofiler_targets.vcf"),              emit: targets_vcf
    tuple val(sample_id), path("${sample_id}_tbprofiler_whole_genome.vcf"),         emit: whole_genome_vcf
    tuple val(sample_id), path("${sample_id}_tbprofiler_provenance.yml"),           emit: provenance
    
    script:
    """

    tb-profiler profile \
      --threads ${task.cpus} \
      --platform ${params.platform} \
      --mapper ${params.mapper} \
      --caller ${params.caller} \
      --depth ${params.min_depth} \
      --af ${params.min_af_used_for_calling} \
      --read1 ${reads_1} \
      --read2 ${reads_2} \
      --prefix ${sample_id} \
      --csv \
      --call_whole_genome

    mv bam/${sample_id}.bam ./${sample_id}_tbprofiler.bam
    mv bam/${sample_id}.bam.bai ./${sample_id}_tbprofiler.bam.bai

    mv vcf/${sample_id}.targets.vcf.gz ./${sample_id}_tbprofiler_targets.vcf.gz
    gunzip ./${sample_id}_tbprofiler_targets.vcf.gz

    mv vcf/${sample_id}.vcf.gz ./${sample_id}_tbprofiler_whole_genome.vcf.gz
    gunzip ./${sample_id}_tbprofiler_whole_genome.vcf.gz

    cp results/${sample_id}.results.csv ${sample_id}_tbprofiler_full_report.csv
    cp results/${sample_id}.results.json ${sample_id}_tbprofiler_full_report.json

    split_tbprofiler_csv.py -p ${sample_id} -s ${sample_id} ${sample_id}_tbprofiler_full_report.csv
    

    printf -- "- process_name: tbprofiler\\n"                  >> ${sample_id}_tbprofiler_provenance.yml
    printf -- "  tools:\\n"                                    >> ${sample_id}_tbprofiler_provenance.yml
    printf -- "    - tool_name: tb-profiler\\n"                >> ${sample_id}_tbprofiler_provenance.yml
    printf -- "      tool_version: \$(tb-profiler profile --version | cut -d ' ' -f 3)\\n" >> ${sample_id}_tbprofiler_provenance.yml
    printf -- "      subcommand: profile\\n"                   >> ${sample_id}_tbprofiler_provenance.yml
    printf -- "      parameters:\\n"                           >> ${sample_id}_tbprofiler_provenance.yml
    printf -- "        - parameter: --platform\\n"             >> ${sample_id}_tbprofiler_provenance.yml
    printf -- "          value: ${params.platform}\\n"         >> ${sample_id}_tbprofiler_provenance.yml

    printf -- "        - parameter: --mapper\\n"               >> ${sample_id}_tbprofiler_provenance.yml
    printf -- "          value: ${params.mapper}\\n"           >> ${sample_id}_tbprofiler_provenance.yml

    printf -- "        - parameter: --caller\\n"               >> ${sample_id}_tbprofiler_provenance.yml
    printf -- "          value: ${params.caller}\\n"           >> ${sample_id}_tbprofiler_provenance.yml

    printf -- "        - parameter: --af\\n"                   >> ${sample_id}_tbprofiler_provenance.yml
    printf -- "          value: ${params.min_af_used_for_calling}\\n" >> ${sample_id}_tbprofiler_provenance.yml

    printf -- "        - parameter: --reporting_af\\n"         >> ${sample_id}_tbprofiler_provenance.yml
    printf -- "          value: ${params.min_af_used_for_prediction}\\n" >> ${sample_id}_tbprofiler_provenance.yml

    printf -- "        - parameter: --prefix\\n"               >> ${sample_id}_tbprofiler_provenance.yml
    printf -- "          value: ${sample_id}\\n"               >> ${sample_id}_tbprofiler_provenance.yml

    printf -- "        - parameter: --csv\\n"                  >> ${sample_id}_tbprofiler_provenance.yml
    printf -- "          value: null\\n"                       >> ${sample_id}_tbprofiler_provenance.yml

    printf -- "        - parameter: --call_whole_genome\\n"    >> ${sample_id}_tbprofiler_provenance.yml
    printf -- "          value: null\\n"                       >> ${sample_id}_tbprofiler_provenance.yml
    """
}

process rename_ref_in_alignment {

  tag { sample_id }
  
  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_tbprofiler*.{bam,bam.bai}", enabled: params.rename_ref, saveAs: { x -> x.replace("_renamed", "") }
  
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

  publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_tbprofiler*.vcf", enabled: params.rename_ref, saveAs: { x -> x.replace("_renamed", "") }

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

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_qualimap_alignment_qc.csv"

    input:
    tuple val(sample_id), file(alignment)

    output:
    tuple val(sample_id), path("${sample_id}_qualimap_alignment_qc.csv"), emit: genome_results
    tuple val(sample_id), path("${sample_id}_qualimap_bamqc_provenance.yml"), emit: provenance
    
    script:
    """
    qualimap bamqc -bam ${alignment[0]} --outdir ${sample_id}_bamqc
    qualimap_bamqc_genome_results_to_csv.py -s ${sample_id} ${sample_id}_bamqc/genome_results.txt > ${sample_id}_qualimap_alignment_qc.csv

    printf -- "- process_name: qualimap\\n"          >> ${sample_id}_qualimap_bamqc_provenance.yml
    printf -- "  tools:\\n"                          >> ${sample_id}_qualimap_bamqc_provenance.yml
    printf -- "    - tool_name: qualimap\\n"         >> ${sample_id}_qualimap_bamqc_provenance.yml
    printf -- "      tool_version: \$(qualimap --version 2>&1 | grep "QualiMap" | awk '{print \$2}')\\n" >> ${sample_id}_qualimap_bamqc_provenance.yml
    printf -- "      subcommand: bamqc\\n"           >> ${sample_id}_qualimap_bamqc_provenance.yml
    printf -- "      parameters:\\n"                 >> ${sample_id}_qualimap_bamqc_provenance.yml
    printf -- "        - parameter: --bam\\n"        >> ${sample_id}_qualimap_bamqc_provenance.yml
    printf -- "          value: ${alignment[0]}\\n"  >> ${sample_id}_qualimap_bamqc_provenance.yml

    """
}


process snpit {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_snpit_unchecked.tsv"

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}_snpit_unchecked.tsv"), emit: snpit_results
    tuple val(sample_id), path("${sample_id}_snpit_provenance.yml"), emit: provenance
    
    script:
    """
    snpit --input ${vcf} > ${sample_id}_snpit_unchecked.tsv

    printf -- "- process_name: snpit\\n"           >> ${sample_id}_snpit_provenance.yml
    printf -- "  tools:\\n"                        >> ${sample_id}_snpit_provenance.yml
    printf -- "    - tool_name: snpit\\n"          >> ${sample_id}_snpit_provenance.yml
    printf -- "      tool_version:  \$(snpit --version 2>&1)\\n" >> ${sample_id}_snpit_provenance.yml
    printf -- "      parameters:\\n"               >> ${sample_id}_snpit_provenance.yml
    printf -- "        - parameter: --input\\n"    >> ${sample_id}_snpit_provenance.yml
    printf -- "          value: ${vcf}\\n"         >> ${sample_id}_snpit_provenance.yml
    """
}


process check_snpit_against_tbprofiler {

    tag { sample_id }

    executor 'local'

    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_snpit.tsv"

    input:
    tuple val(sample_id), path(snpit_results), path(tbprofiler_report_json)

    output:
    tuple val(sample_id), path("${sample_id}_snpit.tsv")
    
    script:
    """
    check_snpit_against_tbprofiler.py \
	--snpit ${snpit_results} \
	--tbprofiler ${tbprofiler_report_json} \
	> ${sample_id}_snpit.tsv
    """
}


process mpileup {

    tag { sample_id }

    input:
    tuple val(sample_id), path(alignment), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}_depths.tsv"), emit: depths
    tuple val(sample_id), path("${sample_id}_mpileup_provenance.yml"), emit: provenance

    script:
    """
    samtools faidx ${ref}
    
    printf "chrom\tpos\tref\tdepth\n" > ${sample_id}_depths.tsv

    samtools mpileup -a \
      --fasta-ref ${ref} \
      --min-BQ 0 \
      --count-orphans \
      ${alignment[0]} | cut -f 1-4 >> ${sample_id}_depths.tsv


    printf -- "- process_name: mpileup\\n"              >> ${sample_id}_mpileup_provenance.yml
    printf -- "  tools:\\n"                             >> ${sample_id}_mpileup_provenance.yml
    printf -- "    - tool_name: samtools\\n"            >> ${sample_id}_mpileup_provenance.yml
    printf -- "      tool_version: \$(samtools --version 2>&1 | awk 'NR==1 {print \$2}')\\n" >> ${sample_id}_mpileup_provenance.yml
    printf -- "      subcommand: mpileup\\n"            >> ${sample_id}_mpileup_provenance.yml
    printf -- "      parameters:\\n"                    >> ${sample_id}_mpileup_provenance.yml
    printf -- "        - parameter: --fasta-ref\\n"     >> ${sample_id}_mpileup_provenance.yml
    printf -- "          value: ${ref}\\n"              >> ${sample_id}_mpileup_provenance.yml

    printf -- "        - parameter: --min-BQ\\n"        >> ${sample_id}_mpileup_provenance.yml
    printf -- "          value: 0\\n"                   >> ${sample_id}_mpileup_provenance.yml

    printf -- "        - parameter: --count-orphans\\n" >> ${sample_id}_mpileup_provenance.yml
    printf -- "          value: null\\n"                >> ${sample_id}_mpileup_provenance.yml
    """
}


process plot_coverage {
    
    tag { sample_id }
    
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_coverage_plot.png"

    input:
    tuple val(sample_id), path(depths), path(resistance_genes_bed)
    
    output:
    tuple val(sample_id), path("${sample_id}_coverage_plot.png")

    script:
    """
    if [ ${params.rename_ref} == "true" ]; then
        sed s/'Chromosome'/'${params.ref_name}'/ ${resistance_genes_bed} > genes.bed
    else
        cp ${resistance_genes_bed} genes.bed
    fi
    plot_coverage.py \
	--input ${depths} \
	--sample-id ${sample_id} \
	--threshold ${params.min_depth} \
	--log-scale \
	--rolling-window 100 \
	--y-limit 500 \
	--genes genes.bed
    """
}


process generate_low_coverage_bed {

    tag { sample_id }
    
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_low_coverage_regions.bed"
    
    input:
    tuple val(sample_id), path(depths)

    output:
    tuple val(sample_id), path("${sample_id}_low_coverage_regions.bed")

    script:
    """
    generate_low_coverage_bed.py \
	--input ${depths} \
	--threshold ${params.min_depth} \
	> ${sample_id}_low_coverage_regions.bed
    """
}


process calculate_gene_coverage {

    tag { sample_id }
    
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_resistance_gene_coverage.csv"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_resistance_drug_coverage.csv"

    input:
    tuple val(sample_id), path(depths), path(resistance_genes_bed), path(resistance_csv)

    output:
    tuple val(sample_id), path("${sample_id}_resistance_gene_coverage.csv")
    tuple val(sample_id), path("${sample_id}_resistance_drug_coverage.csv")

    script:
    """
    calculate_res_gene_depth_qc.py \
	--bed ${resistance_genes_bed} \
	--depth ${depths} \
	--threshold ${params.min_depth} \
	--output ${sample_id}_resistance_gene_coverage.csv

    add_gene_coverage_to_res_csv.py \
    --resistance ${resistance_csv} \
    --coverage ${sample_id}_resistance_gene_coverage.csv \
    --threshold ${params.min_gene_coverage} \
    --output  ${sample_id}_resistance_drug_coverage.csv 


    """
}
