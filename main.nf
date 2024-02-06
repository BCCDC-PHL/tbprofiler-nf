#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { fastp }                          from './modules/tbprofiler.nf'
include { tbprofiler }                     from './modules/tbprofiler.nf'
include { snpit }                          from './modules/tbprofiler.nf'
include { check_snpit_against_tbprofiler } from './modules/tbprofiler.nf'
include { rename_ref_in_alignment }        from './modules/tbprofiler.nf'
include { rename_ref_in_variants as rename_ref_in_targets_variants }       from './modules/tbprofiler.nf'
include { rename_ref_in_variants as rename_ref_in_whole_genome_variants }  from './modules/tbprofiler.nf'
include { qualimap_bamqc }                 from './modules/tbprofiler.nf'
include { mpileup }                        from './modules/tbprofiler.nf'
include { plot_coverage }                  from './modules/tbprofiler.nf'
include { generate_low_coverage_bed }      from './modules/tbprofiler.nf'
include { calculate_gene_coverage }        from './modules/tbprofiler.nf'
include { pipeline_provenance }            from './modules/provenance.nf'
include { collect_provenance }             from './modules/provenance.nf'



workflow {

  ch_start_time = Channel.of(LocalDateTime.now())
  ch_pipeline_name = Channel.of(workflow.manifest.name)
  ch_pipeline_version = Channel.of(workflow.manifest.version)

  ch_pipeline_provenance = pipeline_provenance(ch_pipeline_name.combine(ch_pipeline_version).combine(ch_start_time))

  if (params.samplesheet_input != 'NO_FILE') {
    ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
  } else {
    ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
  }

  ch_resistance_genes_bed = Channel.fromPath("${baseDir}/assets/resistance_genes.bed")

  main:
    fastp(ch_fastq)

    tbprofiler(fastp.out.reads)

    if (params.rename_ref) {
      ch_ref = Channel.fromPath("${baseDir}/assets/NC_000962.3.fa")
      rename_ref_in_alignment(tbprofiler.out.alignment)
      rename_ref_in_targets_variants(tbprofiler.out.targets_vcf)
      rename_ref_in_whole_genome_variants(tbprofiler.out.whole_genome_vcf)
      ch_alignment = rename_ref_in_alignment.out
      ch_whole_genome_variants = rename_ref_in_whole_genome_variants.out
    } else {
      ch_ref = Channel.fromPath("${baseDir}/assets/tbdb_genome.fa")
      ch_alignment = tbprofiler.out.alignment
      ch_whole_genome_variants = tbprofiler.out.whole_genome_vcf
    }

    snpit(ch_whole_genome_variants)

    check_snpit_against_tbprofiler(snpit.out.snpit_results.join(tbprofiler.out.report_json))

    qualimap_bamqc(ch_alignment)

    ch_depths = mpileup(ch_alignment.combine(ch_ref))

    plot_coverage(ch_depths.combine(ch_resistance_genes_bed))

    generate_low_coverage_bed(ch_depths)
    
    calculate_gene_coverage(ch_depths.combine(ch_resistance_genes_bed))

    ch_provenance = fastp.out.provenance
    ch_provenance = ch_provenance.join(tbprofiler.out.provenance).map{ it -> [it[0], [it[1], it[2]]] }
    ch_provenance = ch_provenance.join(snpit.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

    ch_provenance = ch_provenance.join(ch_fastq.map{ it -> it[0] }.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] ] }

    collect_provenance(ch_provenance)
  
}
