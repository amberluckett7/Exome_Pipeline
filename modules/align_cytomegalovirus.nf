#!/usr/bin/env nextflow

process ALIGN_CYTOMEGALOVIRUS {
    cpus 16
    tag "${id}"
    
    input:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file(fastq)
      
    output:
      tuple val(id), val(sex), val(family), val(famSampleCount), file("${id}_cytomegalovirus_grepped2.sam")
      
    script:
      """
      ${params.bwaMem2} mem -t 16 -M ${params.cytomegaloFasta} ${fastq} > ${id}_cytomegalovirus.sam

      grep '^@\\|NC_006273\\.2' ${id}_cytomegalovirus.sam > ${id}_cytomegalovirus_grepped2.sam

      DEST_DIR="${params.batchDir}/r04_cytomegalovirus"
      mkdir -p "\$DEST_DIR"

      SRC_FILE="${id}_cytomegalovirus_grepped2.sam"
      TMP_FILE=".${id}_cytomegalovirus_grepped2.sam.partial.\$\$"

      attempts=5
      delay=10

      for i in \$(seq 1 \$attempts); do
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