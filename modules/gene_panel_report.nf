process GENE_PANEL_REPORT {
    tag "${family}"
    cpus 2
    memory '8 GB'
    container 'docker://al655/bcftools-python3:1.17'
    containerOptions "-B /usr/lib/locale/:/usr/lib/locale/ -B ${params.batchDir} -B ${params.resourcesDir} -B ${params.rsync} -B ${params.base}"

    input:
      tuple val(id), val(sex), val(family), val(famSampleCount), file(snv_vcf), file(snv_idx)
      tuple val(sv_id), val(sv_sex), val(sv_family), val(sv_famSampleCount), file(sv_vcf), file(sv_tbi)
      path gene_panel_tsv
      path report_script

    output:
      tuple val(family), file("${family}_gene_panel_report.html"), optional: true
      path("no_report_samples.txt"), optional: true

    script:
      """
      # Define destination directory upfront so it's always available
      DEST_DIR="${params.batchDir}/gene_panel_reports"
      mkdir -p "\$DEST_DIR"

      python3 ${report_script} \
          --snv    ${snv_vcf} \
          --sv     ${sv_vcf} \
          --panel  ${gene_panel_tsv} \
          --sample ${family} \
          --out    ${family}_gene_panel_report.html

      # Handles samples where no variants passed filters (no HTML generated)
      if [[ ! -s "${family}_gene_panel_report.html" ]]; then
          echo "WARNING: No variants passed filters for ${family} — logging to missing report list." >&2

          # Safely append to the shared text file using flock
          (
              flock -x 200
              echo "${family}" >> "\$DEST_DIR/no_report_samples.txt"
          ) 200>"\$DEST_DIR/no_report_samples.lock"

          exit 0
      fi

      TMP_FILE=".${family}_gene_panel_report.html.partial.\$\$"
      attempts=5
      delay=10
      published=0

      for i in \$(seq 1 \$attempts); do
          if ${params.rsync} -a --no-g --checksum --delay-updates --partial \
              "${family}_gene_panel_report.html" "\$DEST_DIR/\$TMP_FILE"; then
              mv -f "\$DEST_DIR/\$TMP_FILE" "\$DEST_DIR/${family}_gene_panel_report.html"
              published=1
              break
          fi
          echo "rsync failed (attempt \$i/\$attempts) — retrying in \${delay}s" >&2
          sleep "\$delay"
      done

      if [[ \$published -eq 0 ]]; then
          echo "Publish failed for ${family}_gene_panel_report.html" >&2
          exit 1
      fi
      """
}
