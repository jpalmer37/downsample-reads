process fastp {

    tag { sample_id }

    input:
    tuple val(sample_id), path(reads), val(target_coverage)

    output:
    tuple val(sample_id), path("${sample_id}_${target_coverage}x_fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}_${target_coverage}x_fastp.csv"), emit: csv

    script:
    """
    fastp \
      -t ${task.cpus} \
      -i ${reads[0]} \
      -I ${reads[1]} \
      --cut_tail \
      -o ${sample_id}_R1.trim.fastq.gz \
      -O ${sample_id}_R2.trim.fastq.gz \
      -j ${sample_id}_${target_coverage}x_fastp.json

    echo "target_coverage"  >> coverage_field.csv
    echo ${target_coverage} >> coverage_field.csv
    fastp_json_to_csv.py -s ${sample_id} ${sample_id}_${target_coverage}x_fastp.json > ${sample_id}_fastp.csv
    paste -d ',' ${sample_id}_fastp.csv coverage_field.csv > ${sample_id}_${target_coverage}x_fastp.csv
    """
}

process downsample {

    tag { sample_id + ' / ' + coverage + 'x' }

    publishDir "${params.outdir}", pattern: "${sample_id}-downsample-*x_R*.fastq.gz", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), val(coverage)

    output:
    tuple val(sample_id), path("${sample_id}-downsample-*x_R*.fastq.gz"), val(coverage), emit: reads

    script:
    """
    rasusa \
      -i ${reads[0]} \
      -i ${reads[1]} \
      --coverage ${coverage} \
      --genome-size ${params.genome_size} \
      -o ${sample_id}-downsample-${coverage}x_R1.fastq.gz \
      -o ${sample_id}-downsample-${coverage}x_R2.fastq.gz
    """
}

