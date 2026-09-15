process MITO_VARIANTS {
    tag "${id}"
    cpus 4
    memory '8 GB'
    container 'docker://broadinstitute/gatk:4.5.0.0'
    containerOptions "-B ${params.resourcesDir} -B ${params.batchDir} -B ${params.rsync}"

    input:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file(bam), file(bai), val(mito_contig)
      path ref_fasta
      path ref_fai
      path ref_dict

    output:
      tuple val(id), file("${id}.mito.vcf.gz"), file("${id}.mito.vcf.gz.tbi")

    script:
      """
      gatk Mutect2 \
        -R ${ref_fasta} \
        -I ${bam} \
        --mitochondria-mode \
        -L ${mito_contig} \
        -O ${id}.mito.raw.vcf.gz || true

      gatk FilterMutectCalls \
        -R ${ref_fasta} \
        -V ${id}.mito.raw.vcf.gz \
        --mitochondria-mode \
        -O ${id}.mito.vcf.gz || true

      [[ -s "${id}.mito.vcf.gz" ]] || { echo "Mutect2 failed to produce output" >&2; exit 1; }
      [[ -s "${id}.mito.vcf.gz.tbi" ]] || { echo "Mutect2 index missing" >&2; exit 1; }

      DEST_DIR="${params.batchDir}/r04_mito"
      mkdir -p "\$DEST_DIR"

      for f in ${id}.mito.vcf.gz ${id}.mito.vcf.gz.tbi; do
          TMP_FILE=".\${f}.partial.\$\$"
          attempts=5
          delay=10
          published=0

          for i in \$(seq 1 \$attempts); do
              if ${params.rsync} -a --no-g --checksum --delay-updates --partial \
                  "\$f" "\$DEST_DIR/\$TMP_FILE"; then
                  mv -f "\$DEST_DIR/\$TMP_FILE" "\$DEST_DIR/\$f"
                  published=1
                  break
              fi
              echo "rsync failed (attempt \$i/\$attempts) — retrying in \${delay}s" >&2
              sleep "\$delay"
          done

          if [[ \$published -eq 0 ]]; then
              echo "Publish failed for \$f" >&2
              exit 1
          fi
      done

      exit 0
      """
}
