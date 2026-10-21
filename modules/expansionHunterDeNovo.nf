process EXPANSION_HUNTER_DE_NOVO {
    tag "${id}"
    cpus 16

    input:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file(bam), file(bai)

    output:
      tuple file("${id}.locus.tsv"), file("${id}.motif.tsv"), file("${id}.str_profile.json")

    script:
      """
      echo "Running expansionHunterDeNovo for sample ${id}"
      ${params.expHunterDeNovo} profile --reads ${bam} --reference ${params.referenceFasta} --output-prefix ${id} --min-anchor-mapq 50 --max-irr-mapq 40 

      if [ ! -f "${id}.locus.tsv" ]; then touch "${id}.locus.tsv"; fi
      if [ ! -f "${id}.motif.tsv" ]; then touch "${id}.motif.tsv"; fi
      if [ ! -f "${id}.str_profile.json" ]; then echo "{}" > "${id}.str_profile.json"; fi

      DEST_DIR="${params.batchDir}/r04_expansionhunterdenovo"
      mkdir -p "\$DEST_DIR"

      STAGE_DIR="\$DEST_DIR/.stage.${id}.\$\$"
      rm -rf "\$STAGE_DIR"
      mkdir -p "\$STAGE_DIR/.partial"

      FILES=( "${id}.locus.tsv" "${id}.motif.tsv" "${id}.str_profile.json" )

      for f in "\${FILES[@]}"; do
        [[ -e "\$f" ]] || { echo "Missing output data file allocation: \$f" >&2; exit 1; }
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
