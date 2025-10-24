#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { hash_files as hash_fastq_input }          from './modules/hash_files.nf'
include { hash_files as hash_fastq_output }         from './modules/hash_files.nf'
include { fastp as fastp_input }                    from './modules/downsample_reads.nf'
include { downsample_rasusa ; downsample_bbnorm }   from './modules/downsample_reads.nf'
include { fastp as fastp_output }                   from './modules/downsample_reads.nf'
include { pipeline_provenance }                     from './modules/provenance.nf'
include { collect_provenance }                      from './modules/provenance.nf'

workflow {

    ch_workflow_metadata = Channel.value([
	workflow.sessionId,
	workflow.runName,
	workflow.manifest.name,
	workflow.manifest.version,
	workflow.start,
    ])

    if (params.samplesheet_input != 'NO_FILE') {
	ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], [it['R1'], it['R2']]] }
	ch_coverages = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['COVERAGE'], it['GENOME_SIZE']] }
    } else {
	ch_fastq = Channel.fromFilePairs(params.fastq_search_path, flat: true).map{ it -> [it[0].split('_')[0], [it[1], it[2]]] }.unique{ it -> it[0] }
	ch_coverages = ch_fastq.map{ it -> [it[0], params.coverage, params.genome_size] }
    }

    main:

    hash_fastq_input(ch_fastq.join(ch_coverages).map({ it -> [it[0], it[2], it[1]] }).combine(Channel.of("fastq-input")))
    
    ch_fastp_input = ch_fastq.join(ch_coverages.map({ it -> [it[0], it[2]] }))

    fastp = fastp_input(ch_fastp_input.combine(Channel.of("original")))

    if (params.use_filtered_reads) {
        printf("Using filtered reads for downsampling.\n")
        ch_fastq = fastp.filtered_reads
    }

    if (params.tool == 'rasusa') {
        ch_downsample = downsample_rasusa(ch_fastq.join(ch_coverages))
    } else if (params.tool == 'bbnorm') {
        ch_downsample = downsample_bbnorm(ch_fastq.join(ch_coverages))
    } else {
        error "ERROR: Invalid tool choice. Please choose either 'rasusa' or 'bbnorm'."
    }

    hash_fastq_output(ch_downsample.reads.map{ it -> [it[0], it[3], it[1]] }.combine(Channel.of("fastq-output")))

    fastp_output(ch_downsample.reads)

    fastp_input.out.csv.concat(fastp_output.out.csv).map{ it -> it[1] }.collectFile(name: params.collected_outputs_prefix + "_downsampling_summary.csv", storeDir: params.outdir, keepHeader: true, skip: 1, sort: { it -> it[0] })

    // Collect Provenance
    // The basic idea is to build up a channel with the following structure:
    // [sample_id, coverage, [provenance_file_1.yml, provenance_file_2.yml, provenance_file_3.yml...]]
    // At each step, we add another provenance file to the list using the << operator...
    // ...and then concatenate them all together in the 'collect_provenance' process.
    ch_sample_ids_with_coverages = ch_fastq.map({ it -> it[0] }).join(ch_coverages.map({ it -> [it[0], it[1]] }))
    ch_provenance = ch_sample_ids_with_coverages
    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)
    ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map({ it -> [it[0], it[1], [it[2]]] })
    ch_provenance = ch_provenance.join(hash_fastq_input.out.provenance, by: [0, 1]).map{ it -> [it[0], it[1], it[2] << it[3]] }
    ch_provenance = ch_provenance.join(fastp_input.out.provenance).map{ it -> [it[0], it[1], it[2] << it[4]] }
    ch_provenance = ch_provenance.join(ch_downsample.provenance, by: [0, 1]).map{ it -> [it[0], it[1], it[2] << it[3]] }
    ch_provenance = ch_provenance.join(hash_fastq_output.out.provenance, by: [0, 1]).map{ it -> [it[0], it[1], it[2] << it[3]] }

    collect_provenance(ch_provenance.map{ it -> [it[0], it[1], it[2].minus(null)] })
}
