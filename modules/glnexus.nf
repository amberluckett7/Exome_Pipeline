process GLNEXUS {
    container 'docker://cgrlab/glnexus:v1.4.1'
    containerOptions "-B ${params.batchDir} -B ${params.rsync}"
    tag "${family}"
    cpus 16

    input:
      // Changed file() to path() to unpack the grouped lists of files
      tuple val(ids), val(sexes), val(family), val(famSampleCount), path(gvcfs), path(gvcfcsis)

    output:
      tuple val(ids), val(sexes), val(family), val(famSampleCount), file("${family}.vcf.gz"), file("${family}.vcf.gz.csi")

    script:
    """
    mkdir -p \$PWD/temp
    export TMPDIR=\$PWD/temp
    
    glnexus_cli --config DeepVariantWES ${gvcfs.join(" ")} > ${family}.bcf

    bcftools view ${family}.bcf | bgzip -@ 4 -c > ${family}.vcf.gz
    bcftools index ${family}.vcf.gz

    DEST_DIR="${params.batchDir}/r04_vcfs"
    mkdir -p "\$DEST_DIR"

    STAGE_DIR="\$DEST_DIR/.stage.${family}.\$\$"
    rm -rf "\$STAGE_DIR"
    mkdir -p "\$STAGE_DIR/.partial"

    FILES=( "${family}.vcf.gz" "${family}.vcf.gz.csi" )

    for f in "\${FILES[@]}"; do
      [[ -s "\$f" ]] || { echo "Missing/empty output: \$f" >&2; exit 1; }
    done

    attempts=5
    delay=10
    ok=0

    for i in \$(seq 1 \$attempts); do
      # clean stage between attempts
      rm -f "\$STAGE_DIR/"* 2>/dev/null || true
      mkdir -p "\$STAGE_DIR/.partial"

      if ${params.rsync} -a --no-g --checksum --delay-updates \\
          --partial --partial-dir=".partial" \\
          "\${FILES[@]}" "\$STAGE_DIR/"; then
        ok=1
        break
      fi

      echo "rsync to stage failed (attempt \$i/\$attempts) — retrying in \${delay}s" >&2
      sleep "\$delay"
    done

    [[ "\$ok" -eq 1 ]] || { echo "Publish failed: rsync never succeeded" >&2; exit 1; }

    for f in "\${FILES[@]}"; do
      [[ -s "\$STAGE_DIR/\$f" ]] || { echo "Stage missing/empty: \$STAGE_DIR/\$f" >&2; exit 1; }
    done

    # Promote data files
    for f in "\${FILES[@]}"; do
      mv -f "\$STAGE_DIR/\$f" "\$DEST_DIR/\$f"
    done

    rm -rf "\$STAGE_DIR"

    """
}