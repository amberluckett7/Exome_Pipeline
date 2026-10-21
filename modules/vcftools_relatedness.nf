#!/usr/bin/env nextflow

process RELATEDNESS {
    tag "batch_${id[0]}"
    module 'VCFtools/0.1.16-GCC-11.2.0'
    cpus 1
    
    input:
      tuple val(id), val(sex), val(family), val(famSampleCount), file(vcf), file(csi)
      
    output:
      tuple val(id), val(sex), val(family), val(famSampleCount), file("relatedness2.csv")
      
    script:
      """
      declare -a vcfs=( ${vcf.join(" ")} )
      
      # Removed inline tabix write calls on read-only system files
      # Implemented header cleanup to prevent header duplication inside appended CSV streams
      for i in "\${vcfs[@]}"; do
          if [ ! -f relatedness2.csv ]; then
              vcftools --relatedness2 -c --gzvcf "\$i" >> relatedness2.csv
          else
              vcftools --relatedness2 -c --gzvcf "\$i" | sed 1d >> relatedness2.csv
          fi
      done

      DEST_DIR="${params.batchDir}/r04_metrics"
      mkdir -p "\$DEST_DIR"

      SRC_FILE="relatedness2.csv"
      TMP_FILE=".relatedness2.csv.partial.\$\$"

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