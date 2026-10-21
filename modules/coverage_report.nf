process COVERAGE_REPORT {
    tag "batch_metrics"
    cpus 4
    container 'docker://amazoncorretto:21.0.7'
    containerOptions "-B ${params.javaDir} -B ${params.base} -B ${params.batchDir} -B ${params.rsync}"

    input:
      path cov_files
      path target_bed

    output:
      file "coverage_report"

    script:
      """
      echo -e "Sample\\tMean_depth\\t10X\\t20X\\t30X\\t50X\\t100X\\tVDP\\tFold80" > coverage_report

      for cov_file in *.coverage; do
          sample=\${cov_file%.coverage}
          awk -v sample="\${sample}" '
              /^BAIT_SET/ {
                  for (i=1; i<=NF; i++) col[\$i]=i
                  next
              }
              /^targets/ {
                  printf "%s\\t%.2f\\t%.2f\\t%.2f\\t%.2f\\t%.2f\\t%.2f\\t%.4f\\t%.2f\\n",
                      sample,
                      \$col["MEAN_TARGET_COVERAGE"],
                      \$col["PCT_TARGET_BASES_10X"]*100,
                      \$col["PCT_TARGET_BASES_20X"]*100,
                      \$col["PCT_TARGET_BASES_30X"]*100,
                      \$col["PCT_TARGET_BASES_50X"]*100,
                      \$col["PCT_TARGET_BASES_100X"]*100,
                      \$col["HET_SNP_SENSITIVITY"],
                      \$col["FOLD_80_BASE_PENALTY"]
              }
          ' "\${cov_file}"
      done >> coverage_report

      DEST_DIR="${params.batchDir}/r04_metrics"
      mkdir -p "\$DEST_DIR"

      SRC_FILE="coverage_report"
      TMP_FILE=".coverage_report.partial.\$\$"

      attempts=5
      delay=10
      published=0

      for i in \$(seq 1 \$attempts); do
          if ${params.rsync} -a --no-g --checksum --delay-updates --partial \
              "\$SRC_FILE" "\$DEST_DIR/\$TMP_FILE" ; then
              mv -f "\$DEST_DIR/\$TMP_FILE" "\$DEST_DIR/\$SRC_FILE"
              published=1
              break
          fi
          echo "rsync failed (attempt \$i/\$attempts) — retrying in \${delay}s" >&2
          sleep "\$delay"
      done

      if [[ \$published -eq 0 ]]; then
          echo "Publish failed after \$attempts attempts" >&2
          rm -f "\$DEST_DIR/\$TMP_FILE"
          exit 1
      fi
      """
}
