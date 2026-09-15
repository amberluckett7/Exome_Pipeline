process INSERT_SIZES {
    tag "${id}"
    cpus 2
    container 'docker://amazoncorretto:21.0.7'
    containerOptions "-B ${params.javaDir} -B ${params.base} -B ${params.batchDir} -B ${params.rsync}"

    input:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file(bam), file(bai)
      file(target_bed)

    output:
      tuple val(id), file("${id}.insertSize.histogram"), file("${id}.insertSize.stats")

    script:
    """
    java -Xmx5g -XX:ParallelGCThreads=${task.cpus} -XX:ConcGCThreads=${task.cpus} \\
      -cp ${params.javaDir}:${params.gatkJar} InsertSizeHistogram \\
      ${bam} > ${id}.insertSize.histogram 2> ${id}.insertSize.stats

    DEST_DIR="${params.batchDir}/r04_metrics"
    mkdir -p "\$DEST_DIR"

    STAGE_DIR="\$DEST_DIR/.stage.${id}.\$\$"
    rm -rf "\$STAGE_DIR"
    mkdir -p "\$STAGE_DIR/.partial"

    FILES=( "${id}.insertSize.histogram" "${id}.insertSize.stats" )

    for f in "\${FILES[@]}"; do
      [[ -s "\$f" ]] || { echo "Missing/empty output: \$f" >&2; exit 1; }
    done

    attempts=5
    delay=10
    ok=0

    for i in \$(seq 1 \$attempts); do
      rm -f "\$STAGE_DIR/"* 2>/dev/null || true
      mkdir -p "\$STAGE_DIR/.partial"

      if ${params.rsync} -a --no-g --checksum --delay-updates \\
          --partial --partial-dir=".partial" \\
          "\${FILES[@]}" "\$STAGE_DIR/"; then
        ok=1
        break
      fi
      sleep "\$delay"
    done

    [[ "\$ok" -eq 1 ]] || { echo "Publish failed: rsync never succeeded" >&2; exit 1; }

    for f in "\${FILES[@]}"; do
      mv -f "\$STAGE_DIR/\$f" "\$DEST_DIR/\$f"
    done
    rm -rf "\$STAGE_DIR"
    """
}
