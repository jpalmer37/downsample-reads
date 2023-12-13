#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { fastp as fastp_input }  from './modules/downsample_reads.nf'
include { downsample }            from './modules/downsample_reads.nf'
include { fastp as fastp_output } from './modules/downsample_reads.nf'

workflow {

    if (params.samplesheet_input != 'NO_FILE') {
	ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }
    } else {
	ch_fastq = Channel.fromFilePairs(params.fastq_search_path, flat: true).map{ it -> [it[0].split('_')[0], [it[1], it[2]]] }.unique{ it -> it[0] }
    }
    if (params.coverages != 'NO_FILE') {
	ch_coverages = Channel.fromPath(params.coverages).splitCsv(header: false)
    } else {
	ch_coverages = Channel.of([params.coverage])
    }

    main:
    fastp_input(ch_fastq.combine(Channel.of("original")))
    downsample(ch_fastq.combine(ch_coverages))
    fastp_output(downsample.out.reads)
    fastp_input.out.csv.concat(fastp_output.out.csv).map{ it -> it[1] }.collectFile(name: "fastp.csv", storeDir: params.outdir, keepHeader: true, skip: 1, sort: { it -> Integer.parseInt(it.readLines()[1].split(',')[1]) })
}
