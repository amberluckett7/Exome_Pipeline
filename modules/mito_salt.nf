process MITO_SALT {
    tag "${id}"
    cpus 8
    memory '32 GB'
    container 'docker://al655/mitosalt:1.1.1'
    containerOptions "-B ${params.resourcesDir} -B ${params.batchDir} -B ${params.rsync}"

    input:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file(bam), file(bai), val(mito_contig)

    output:
      tuple val(id), file("indel/${id}.cluster"), optional: true
      tuple val(id), file("indel/${id}.breakpoint"), optional: true
      tuple val(id), file("indel/${id}.bed"), optional: true
      tuple val(id), file("indel/${id}.count.txt"), optional: true

    script:
      """
      samtools view -@ ${task.cpus} -b ${bam} ${mito_contig} | \\
        samtools sort -@ ${task.cpus} -n | \\
        samtools fastq \\
          -1 ${id}_R1.fastq.gz \\
          -2 ${id}_R2.fastq.gz \\
          -0 /dev/null \\
          -s /dev/null

      perl \$(which MitoSAlt1.1.1.pl) \\
        ${params.mitosaltConfig} \\
        ${id}_R1.fastq.gz \\
        ${id}_R2.fastq.gz \\
        ${id}

      DEST_DIR="${params.batchDir}/r04_mito"
      mkdir -p "\$DEST_DIR"

      for f in indel/${id}.cluster indel/${id}.breakpoint indel/${id}.bed indel/${id}.count.txt; do
          if [[ -s "\$f" ]]; then
              TMP_FILE=".\$(basename \$f).partial.\$\$"
              attempts=5
              delay=10
              published=0

              for i in \$(seq 1 \$attempts); do
                  if ${params.rsync} -a --no-g --checksum --delay-updates --partial \\
                      "\$f" "\$DEST_DIR/\$TMP_FILE"; then
                      mv -f "\$DEST_DIR/\$TMP_FILE" "\$DEST_DIR/\$(basename \$f)"
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
          fi
      done

      exit 0
      """
}
