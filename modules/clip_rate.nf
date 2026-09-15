process CLIP_RATE {
    tag "${id}"
    cpus 2
    
    input:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file(bam), file(bai)

    output:
      tuple val(id), val(sex), file("${id}.clipRate")

    script:
      """
      java -Xmx5g -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2 -cp ${params.javaDir}:${params.gatkJar} BamDoubleClipRate \\
        ${bam} \\
        > ${id}.clipRate

      DEST_DIR="${params.batchDir}/r04_metrics"
      mkdir -p "\$DEST_DIR"

      SRC_FILE="${id}.clipRate"
      TMP_FILE=".${id}.clipRate.partial.\$\$"

      attempts=5
      delay=10

      for i in \$(seq 1 \$attempts); do
        # Swapped 'rsync' command  to run '${params.rsync}' instead
        if ${params.rsync} -a --no-g --checksum --delay-updates --partial \\
            "\$SRC_FILE" "\$DEST_DIR/\$TMP_FILE" ; then
          mv -f "\$DEST_DIR/\$TMP_FILE" "\$DEST_DIR/\$SRC_FILE"
          break
        fi
        echo "rsync failed (attempt \$i/\$attempts) — retrying in \${delay}s" >&2
        sleep "\$delay"
      done

      if [[ -e "\$DEST_DIR/\$TMP_FILE" ]]; then
        echo "Publish failed: temp file still present: \$DEST_DIR/\$TMP_FILE" >&2
        exit 1
      fi
      """
}