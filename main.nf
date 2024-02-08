#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { hash_files }                     from './modules/hash_files.nf'
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

    ch_workflow_metadata = Channel.value([
	workflow.sessionId,
	workflow.runName,
	workflow.manifest.name,
	workflow.manifest.version,
	workflow.start,
    ])

    if (params.samplesheet_input != 'NO_FILE') {
	ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
    } else {
	ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
    }

    ch_resistance_genes_bed = Channel.fromPath("${baseDir}/assets/resistance_genes.bed")

    main:

    hash_files(ch_fastq.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq-input")))

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

    mpileup(ch_alignment.combine(ch_ref))

    ch_depths = mpileup.out.depths

    plot_coverage(ch_depths.combine(ch_resistance_genes_bed))

    generate_low_coverage_bed(ch_depths)
    
    calculate_gene_coverage(ch_depths.combine(ch_resistance_genes_bed))

    if (params.collect_outputs) {
	fastp.out.csv.map{ it -> it[1] }.collectFile(
	    name: params.collected_outputs_prefix + "_fastp.csv",
	    storeDir: params.outdir,
	    keepHeader: true,
	    sort: { it -> it.readLines()[1].split(',')[0] }
	)
	check_snpit_against_tbprofiler.out.map{ it -> it[1] }.collectFile(
	    name: params.collected_outputs_prefix + "_snpit.tsv",
	    storeDir: params.outdir,
	    keepHeader: true,
	    newLine: true,
	    sort: { it -> it.readLines()[1].split('\\t')[0] }
	)
	tbprofiler.out.resistance_mutations_csv.map{ it -> it[1] }.collectFile(
	    name: params.collected_outputs_prefix + "_tbprofiler_resistance_mutations.csv",
	    storeDir: params.outdir,
	    keepHeader: true,
	    sort: { it -> it.readLines()[1].split(',')[0..1].join('-') }
	)
	tbprofiler.out.summary_csv.map{ it -> it[1] }.collectFile(
	    name: params.collected_outputs_prefix + "_tbprofiler_summary.csv",
	    storeDir: params.outdir,
	    keepHeader: true,
	    sort: { it -> it.readLines()[1].split(',')[0] }
	)
	qualimap_bamqc.out.genome_results.map{ it -> it[1] }.collectFile(
	    name: params.collected_outputs_prefix + "_qualimap_alignment_qc.csv",
	    storeDir: params.outdir,
	    keepHeader: true,
	    sort: { it -> it.readLines()[1].split(',')[0] }
	)
	
    }
    
    // Collect Provenance
    // The basic idea is to build up a channel with the following structure:
    // [sample_id, [provenance_file_1.yml, provenance_file_2.yml, provenance_file_3.yml...]]
    // At each step, we add another provenance file to the list using the << operator...
    // ...and then concatenate them all together in the 'collect_provenance' process.
    ch_sample_ids = ch_fastq.map{ it -> it[0] }
    ch_provenance = ch_sample_ids
    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)
    ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it ->     [it[0], [it[1]]] }
    ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it ->     [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it ->          [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(tbprofiler.out.provenance).map{ it ->     [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(snpit.out.provenance).map{ it ->          [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(mpileup.out.provenance).map{ it ->        [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(qualimap_bamqc.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    

    collect_provenance(ch_provenance)
  

}
