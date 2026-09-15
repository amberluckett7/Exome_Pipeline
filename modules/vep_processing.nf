#!/usr/bin/env nextflow

process VEP_PROCESSING {
    tag "${family}"
    module 'BCFtools/1.17-GCC-12.2.0'
    
    input:
      tuple val(id), val(sex), val(family), val(famSampleCount), file(vcf)
    
    output:
      tuple val(id), val(sex), val(family), val(famSampleCount), file("${family}_vep_annotated.vcf.gz"), file("${family}_vep_annotated.vcf.gz.csi")
    
    script:
      """
      bcftools view -Oz -o ${family}_vep_annotated.vcf.gz ${vcf}
      bcftools index --csi ${family}_vep_annotated.vcf.gz

      DEST_DIR="${params.batchDir}/r04_vep"
      mkdir -p "\$DEST_DIR"

      STAGE_DIR="\$DEST_DIR/.stage.${family}.\$\$"
      rm -rf "\$STAGE_DIR"
      mkdir -p "\$STAGE_DIR/.partial"

      FILES=( "${family}_vep_annotated.vcf.gz" "${family}_vep_annotated.vcf.gz.csi" )

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

        echo "rsync to stage failed (attempt \$i/\$attempts) — retrying in \${delay}s" >&2
        sleep "\$delay"
      done

      [[ "\$ok" -eq 1 ]] || { echo "Publish failed: rsync never succeeded" >&2; exit 1; }

      for f in "\${FILES[@]}"; do
        [[ -s "\$STAGE_DIR/\$f" ]] || { echo "Stage missing/empty: \$STAGE_DIR/\$f" >&2; exit 1; }
      done

      for f in "\${FILES[@]}"; do
        mv -f "\$STAGE_DIR/\$f" "\$DEST_DIR/\$f"
      done

      rm -rf "\$STAGE_DIR"
      """
}
