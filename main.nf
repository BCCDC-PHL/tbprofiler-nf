#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { fastp }                   from './modules/tbprofiler.nf'
include { tbprofiler }              from './modules/tbprofiler.nf'
include { rename_ref_in_alignment } from './modules/tbprofiler.nf'
include { rename_ref_in_variants as rename_ref_in_targets_variants }       from './modules/tbprofiler.nf'
include { rename_ref_in_variants as rename_ref_in_whole_genome_variants }  from './modules/tbprofiler.nf'
include { qualimap_bamqc }          from './modules/tbprofiler.nf'
include { pipeline_provenance }     from './modules/provenance.nf'
include { collect_provenance }      from './modules/provenance.nf'

include { snpit }                   from './modules/tbprofiler.nf'

include { coverage }                from './modules/tbprofiler.nf'

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

  main:
    fastp(ch_fastq)

    tbprofiler(fastp.out.reads)

    if (params.rename_ref) {
      rename_ref_in_alignment(tbprofiler.out.alignment)
      rename_ref_in_targets_variants(tbprofiler.out.targets_vcf)
      rename_ref_in_whole_genome_variants(tbprofiler.out.whole_genome_vcf)
      qualimap_bamqc(rename_ref_in_alignment.out)
    } else {
      qualimap_bamqc(tbprofiler.out.alignment)
    }

    if (params.rename_ref) {
      snpit(rename_ref_in_whole_genome_variants.out)
    } else {
      snpit(tbprofiler.out.whole_genome_vcf)
    }
    
    ch_coverage(tbprofiler.out.alignment)
    ch_provenance = fastp.out.provenance
    ch_provenance = ch_provenance.join(tbprofiler.out.provenance).map{ it -> [it[0], [it[1], it[2]]] }
    ch_provenance = ch_provenance.join(snpit.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

    ch_provenance = ch_provenance.join(ch_fastq.map{ it -> it[0] }.combine(ch_pipeline_provenance)).map{ it -> [it[0], it[1] ] }

    collect_provenance(ch_provenance)
  
}
