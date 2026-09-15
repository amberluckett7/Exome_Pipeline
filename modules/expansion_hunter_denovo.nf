process EXPANSION_HUNTER_DE_NOVO {
    tag "${id}"
    cpus 4
    container 'https://depot.galaxyproject.org/singularity/expansionhunterdenovo:0.9.0--h9ee0642_1'
    containerOptions "-B ${params.batchDir} -B ${params.rsync}"

    input:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file(bam), file(bai)
      file(fasta)
      file(fai)

    output:
      tuple val(id), file("${id}.str_profile.json")

    script:
    """
    ExpansionHunterDenovo profile \\
      --reads ${bam} \\
      --reference ${fasta} \\
      --output-prefix ${id} \\
      --min-anchor-mapq 50

    DEST_DIR="${params.batchDir}/r04_ehdn"
    mkdir -p "\$DEST_DIR"

    STAGE_DIR="\$DEST_DIR/.stage.${id}.\$\$"
    rm -rf "\$STAGE_DIR"
    mkdir -p "\$STAGE_DIR/.partial"

    FILES=( "${id}.str_profile.json" )

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

    mv -f "\$STAGE_DIR/${id}.str_profile.json" "\$DEST_DIR/${id}.str_profile.json"
    rm -rf "\$STAGE_DIR"
    """
}
