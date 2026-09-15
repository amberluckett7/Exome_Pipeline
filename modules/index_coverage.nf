process INDEX_COVERAGE {
    tag "batch_${id[0]}"
    cpus 1

    input:
      tuple val(id), val(sex), file(coverage)
      file control
      file target_bed

    output:
      tuple val(id), val(sex), file("Coverage.indexed")

    script:
      """
      touch Coverage.indexed

      DEST_DIR="${params.batchDir}/r04_metrics"
      mkdir -p "\$DEST_DIR"

      ${params.rsync} -a --no-g --checksum --delay-updates --partial \
          "Coverage.indexed" "\$DEST_DIR/Coverage.indexed"
      """
}
