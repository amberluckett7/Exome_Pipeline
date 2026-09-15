process DUPMETRICS {
    cpus 2
    memory '32 GB'
    // publishDir "${params.batchDir}/r04_assembly", mode: 'copy', overwrite: true, failOnError: true
    container 'docker://amazoncorretto:21.0.7'
    containerOptions "-B ${params.resourcesDir} -B ${params.picardJar} -B ${params.batchDir} -B ${params.rsync}" 
    tag "${id}"

    input:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file(bam), file(bai)

    output:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file("${id}.markdup_metrics"), emit: metrics

    script:
      """
      mkdir tmp_dir
      export TMPDIR=\${PWD}/tmp_dir
      export _JAVA_OPTIONS="-Djava.io.tmpdir=\$TMPDIR"
      java -Xmx30g -jar ${params.picardJar} CollectDuplicateMetrics \
      INPUT=${bam} \
      METRICS_FILE=${id}.markdup_metrics \
      VALIDATION_STRINGENCY=SILENT

      DEST_DIR="${params.batchDir}/r04_assembly"
      mkdir -p "\$DEST_DIR"

      SRC_FILE="${id}.markdup_metrics"
      TMP_FILE=".${id}.markdup_metrics.partial.\$\$"

      attempts=5
      delay=10

      for i in \$(seq 1 \$attempts); do
        # copy to temp name first
        if ${params.rsync} -a --no-g --checksum --delay-updates --partial \
            "\$SRC_FILE" "\$DEST_DIR/\$TMP_FILE" ; then
          # then rename into place
          mv -f "\$DEST_DIR/\$TMP_FILE" "\$DEST_DIR/\$SRC_FILE"
          break
        fi
        echo "rsync failed (attempt \$i/\$attempts) — retrying in \${delay}s" >&2
        sleep "\$delay"
      done

      # if temp still exists, we never successfully promoted to final name
      if [[ -e "\$DEST_DIR/\$TMP_FILE" ]]; then
        echo "Publish failed: temp file still present: \$DEST_DIR/\$TMP_FILE" >&2
        exit 1
      fi
      """
}