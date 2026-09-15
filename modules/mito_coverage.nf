process MITO_COVERAGE {
    tag "${id}"
    cpus 2
    container 'docker://mgibio/samtools:v1.21-noble'
    containerOptions "-B ${params.batchDir} -B ${params.rsync}"

    input:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file(bam), file(bai), val(mito_contig)

    output:
      tuple val(id), file("${id}.mito_coverage.txt")

    script:
      """
      samtools coverage \
        --region ${mito_contig} \
        ${bam} \
        -o ${id}.mito_coverage.txt

      DEST_DIR="${params.batchDir}/r04_metrics"
      mkdir -p "\$DEST_DIR"

      TMP_FILE=".${id}.mito_coverage.txt.partial.\$\$"
      attempts=5
      delay=10
      published=0

      for i in \$(seq 1 \$attempts); do
          if ${params.rsync} -a --no-g --checksum --delay-updates --partial \
              "${id}.mito_coverage.txt" "\$DEST_DIR/\$TMP_FILE"; then
              mv -f "\$DEST_DIR/\$TMP_FILE" "\$DEST_DIR/${id}.mito_coverage.txt"
              published=1
              break
          fi
          echo "rsync failed (attempt \$i/\$attempts) — retrying in \${delay}s" >&2
          sleep "\$delay"
      done

      if [[ \$published -eq 0 ]]; then
          echo "Publish failed after \$attempts attempts" >&2
          exit 1
      fi
      """
}
