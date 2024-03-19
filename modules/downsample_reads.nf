process fastp {

    tag { sample_id }

    input:
    tuple val(sample_id), path(reads), val(genome_size), val(target_coverage)

    output:
    tuple val(sample_id), path("${sample_id}_${target_coverage}x_fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}_${target_coverage}x_fastp.csv"), emit: csv
    tuple val(sample_id), val(target_coverage), path("${sample_id}_${target_coverage}x_fastp_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: fastp\\n"  >> ${sample_id}_${target_coverage}x_fastp_provenance.yml
    printf -- "  tools:\\n"               >> ${sample_id}_${target_coverage}x_fastp_provenance.yml
    printf -- "    - tool_name: fastp\\n" >> ${sample_id}_${target_coverage}x_fastp_provenance.yml
    printf -- "      tool_version: \$(fastp --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_${target_coverage}x_fastp_provenance.yml
    printf -- "      parameters:\\n"               >> ${sample_id}_${target_coverage}x_fastp_provenance.yml
    printf -- "        - parameter: --cut_tail\\n" >> ${sample_id}_${target_coverage}x_fastp_provenance.yml
    printf -- "          value: null\\n"           >> ${sample_id}_${target_coverage}x_fastp_provenance.yml

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

    echo "genome_size"  >> genome_size_field.csv
    echo ${genome_size} >> genome_size_field.csv

    fastp_json_to_csv.py -s ${sample_id} ${sample_id}_${target_coverage}x_fastp.json > ${sample_id}_fastp.csv
    paste -d ',' ${sample_id}_fastp.csv genome_size_field.csv coverage_field.csv | calculate_estimated_coverage.py > ${sample_id}_${target_coverage}x_fastp.csv
    """
}

process downsample {

    tag { sample_id + ' / ' + genome_size + ' / ' + coverage + 'x' }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}-downsample-*x_R*.fastq.gz", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), val(genome_size), val(coverage)

    output:
    tuple val(sample_id), path("${sample_id}-downsample-*x_R*.fastq.gz"), val(genome_size), val(coverage), emit: reads
    tuple val(sample_id), val(coverage), path("${sample_id}_${coverage}x_downsample_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: downsample\\n"         >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "  tools:\\n"                           >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "    - tool_name: rasusa\\n"            >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "      tool_version: \$(rasusa --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "      parameters:\\n"                  >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "        - parameter: --coverage\\n"    >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "          value: ${coverage}\\n"       >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "        - parameter: --genome-size\\n" >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "          value: ${genome_size}\\n"    >> ${sample_id}_${coverage}x_downsample_provenance.yml
    
    rasusa \
      -i ${reads[0]} \
      -i ${reads[1]} \
      --coverage ${coverage} \
      --genome-size ${genome_size} \
      -o ${sample_id}-downsample-${coverage}x_R1.fastq.gz \
      -o ${sample_id}-downsample-${coverage}x_R2.fastq.gz
    """
}

