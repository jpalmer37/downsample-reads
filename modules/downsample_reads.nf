process fastp {

    tag { sample_id + ' / ' + target_coverage_filename }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_${target_coverage_filename}_downsampling_summary.csv", mode: 'copy'
    // publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_R*.trim.fastq.gz", mode: 'copy', enabled: params.use_filtered_reads

    input:
    tuple val(sample_id), path(reads), val(genome_size), val(target_coverage)

    output:
    tuple val(sample_id), path("${sample_id}_${target_coverage_filename}_fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}_${target_coverage_filename}_downsampling_summary.csv"), emit: csv
    tuple val(sample_id), path("${sample_id}_R*.trim.fastq.gz"), emit: filtered_reads
    tuple val(sample_id), val(target_coverage), path("${sample_id}_original_fastp_provenance.yml"), emit: provenance, optional: true

    script:
    if (target_coverage == 'original') {
	target_coverage_filename = 'original'
    } else {
	target_coverage_filename = target_coverage + 'x'
    }
    if (target_coverage == 'original' && params.enable_quality_trimming) {
	quality_trimming = '--cut_tail'
    } else {
	quality_trimming = ''
    }
    if (target_coverage == 'original' && params.disable_quality_filtering) {
	quality_filtering = '--disable_quality_filtering'
    } else {
	quality_filtering = ''
    }
    """
    printf -- "- process_name: fastp\\n"  >> ${sample_id}_${target_coverage_filename}_fastp_provenance.yml
    printf -- "  tools:\\n"               >> ${sample_id}_${target_coverage_filename}_fastp_provenance.yml
    printf -- "    - tool_name: fastp\\n" >> ${sample_id}_${target_coverage_filename}_fastp_provenance.yml
    printf -- "      tool_version: \$(fastp --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_${target_coverage_filename}_fastp_provenance.yml
    if [[ "${quality_trimming}" != "" || "${quality_filtering}" != "" ]]; then
        printf -- "      parameters:\\n"               >> ${sample_id}_${target_coverage_filename}_fastp_provenance.yml
    fi
    if [[ "${quality_trimming}" != "" ]]; then
        printf -- "        - parameter: --cut_tail\\n" >> ${sample_id}_${target_coverage_filename}_fastp_provenance.yml
        printf -- "          value: null\\n"           >> ${sample_id}_${target_coverage_filename}_fastp_provenance.yml
    fi
    if [[ "${quality_filtering}" != "" ]]; then
        printf -- "        - parameter: --disable_quality_filtering\\n" >> ${sample_id}_${target_coverage_filename}_fastp_provenance.yml
        printf -- "          value: null\\n"           >> ${sample_id}_${target_coverage_filename}_fastp_provenance.yml
    fi

    fastp \
	-t ${task.cpus} \
	-i ${reads[0]} \
	-I ${reads[1]} \
	${quality_trimming} \
	${quality_filtering} \
	-o ${sample_id}_R1.trim.fastq.gz \
	-O ${sample_id}_R2.trim.fastq.gz \
	-j ${sample_id}_${target_coverage_filename}_fastp.json

    echo "target_coverage"  >> coverage_field.csv
    echo ${target_coverage} >> coverage_field.csv

    echo "genome_size"  >> genome_size_field.csv
    echo ${genome_size} >> genome_size_field.csv

    fastp_json_to_csv.py -s ${sample_id} ${sample_id}_${target_coverage_filename}_fastp.json > ${sample_id}_fastp.csv
    paste -d ',' ${sample_id}_fastp.csv genome_size_field.csv coverage_field.csv | calculate_estimated_coverage.py > ${sample_id}_${target_coverage_filename}_downsampling_summary.csv
    """
}

process downsample_rasusa {

    tag { sample_id + ' / ' + genome_size + ' / ' + coverage + 'x' }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}-downsample-rasusa-*x_R*.fastq.gz", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), val(coverage), val(genome_size)

    output:
    tuple val(sample_id), path("${sample_id}-downsample-rasusa-*x_R*.fastq.gz"), val(genome_size), val(coverage), emit: reads
    tuple val(sample_id), val(coverage), path("${sample_id}_${coverage}x_downsample_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: downsample_rasusa\\n"             >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "  tools:\\n"                               >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "    - tool_name: rasusa\\n"                >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "      tool_version: \$(rasusa --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "      parameters:\\n"                      >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "        - parameter: --coverage\\n"        >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "          value: ${coverage}\\n"           >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "        - parameter: --genome-size\\n"     >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "          value: ${genome_size}\\n"        >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "        - parameter: --seed\\n"            >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "          value: ${params.random_seed}\\n" >> ${sample_id}_${coverage}x_downsample_provenance.yml
    
    rasusa reads \
        --seed ${params.random_seed} \
        --coverage ${coverage} \
        --genome-size ${genome_size} \
        ${reads[0]} \
        ${reads[1]} \
        -o ${sample_id}-downsample-rasusa-${coverage}x_R1.fastq.gz \
        -o ${sample_id}-downsample-rasusa-${coverage}x_R2.fastq.gz
    """
}

process downsample_bbnorm {

    tag { sample_id + ' / ' + genome_size + ' / ' + coverage + 'x' }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}-downsample-*x_R*.fastq.gz", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), val(coverage), val(genome_size)

    output:
    tuple val(sample_id), path("${sample_id}-downsample-bbnorm-*x_R*.fastq.gz"), val(genome_size), val(coverage), emit: reads
    tuple val(sample_id), val(coverage), path("${sample_id}_${coverage}x_downsample_provenance.yml"), emit: provenance

    script:
    max_memory_gb = task.memory.toString().split(" ")[0]

    """
    printf -- "- process_name: downsample_bbnorm\\n"             >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "  tools:\\n"                               >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "    - tool_name: bbnorm\\n"                >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "      tool_version: \$(bbnorm.sh --version 2>&1 | head -n 2 | tail -n 1 | cut -d ' ' -f 3)\\n" >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "      parameters:\\n"                      >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "        - parameter: --coverage\\n"        >> ${sample_id}_${coverage}x_downsample_provenance.yml
    printf -- "          value: ${coverage}\\n"           >> ${sample_id}_${coverage}x_downsample_provenance.yml
    
    bbnorm.sh \
        -Xmx${max_memory_gb}g \
        in1=${reads[0]} \
        in2=${reads[1]} \
        out1=${sample_id}-downsample-bbnorm-${coverage}x_R1.fastq \
        out2=${sample_id}-downsample-bbnorm-${coverage}x_R2.fastq \
        target=${params.coverage} \
	

    gzip ${sample_id}-downsample-bbnorm-${coverage}x_R1.fastq
    gzip ${sample_id}-downsample-bbnorm-${coverage}x_R2.fastq

    """
}