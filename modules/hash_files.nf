process hash_files {

    tag { sample_id + " / " + file_type }

    input:
    tuple  val(sample_id), val(coverage), path(files_to_hash), val(file_type)

    output:
    tuple  val(sample_id), val(coverage), path("${sample_id}_${coverage_filename}_${file_type}.sha256.csv"), emit: csv
    tuple  val(sample_id), val(coverage), path("${sample_id}_${coverage_filename}_${file_type}_provenance.yml"), emit: provenance

    script:
    if (coverage == "original") {
	coverage_filename = "original"
    } else {
	coverage_filename = coverage + "x"
    }
    
    """
    shasum -a 256 ${files_to_hash} | tr -s ' ' ',' >> ${sample_id}_${coverage_filename}_${file_type}.sha256.csv
    while IFS=',' read -r hash filename; do
      printf -- "- filename: \$filename\\n"        >> ${sample_id}_${coverage_filename}_${file_type}_provenance.yml;
      printf -- "  file_type: ${file_type}\\n"     >> ${sample_id}_${coverage_filename}_${file_type}_provenance.yml;
      printf -- "  sha256: \$hash\\n"              >> ${sample_id}_${coverage_filename}_${file_type}_provenance.yml;
    done < ${sample_id}_${coverage_filename}_${file_type}.sha256.csv
    """

}
