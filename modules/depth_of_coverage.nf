process DEPTH_OF_COVERAGE {
    tag "${id}"
    label 'process_low'
    cpus 4

    input:
      tuple val(id), val(sex), val(family), val(trio), val(count), file(bam), file(bai)
      path target_bed

    output:
      tuple val(id), val(sex), file("${id}.coverage"), emit: coverage
      tuple val(id), file("${id}.coverage_summary.txt"), emit: summary

    script:
      """
      # 1. Convert the BED file into a Picard Interval List
      java -Xmx2g -jar ${params.gatkJar} BedToIntervalList \\
        -I "${target_bed}" \\
        -SD "${bam}" \\
        -O targets.interval_list

      # 2. Run Picard CollectHsMetrics to calculate target distributions
      java -Xmx8g -jar ${params.gatkJar} CollectHsMetrics \\
        -I "${bam}" \\
        -R ${params.referenceFasta} \\
        -BI targets.interval_list \\
        -TI targets.interval_list \\
        -O "${id}.coverage"

      # 3. Compile summary data 
      out_summary="${id}.coverage_summary.txt"
      echo -e "Sample_ID\\tMean_Coverage\\tPct_10x\\tPct_20x\\tPct_30x" > "\${out_summary}"

      awk '
        /^BAIT_SET/ { for(i=1;i<=NF;i++) col[\$i]=i; next }
        /^targets/ && col["MEAN_TARGET_COVERAGE"] {
            printf "%s\\t%.3fx\\t%.2f%%\\t%.2f%%\\t%.2f%%\\n",
              id,
              \$col["MEAN_TARGET_COVERAGE"],
              \$col["PCT_TARGET_BASES_10X"]*100,
              \$col["PCT_TARGET_BASES_20X"]*100,
              \$col["PCT_TARGET_BASES_30X"]*100
        }
      ' id="${id}" "${id}.coverage" >> "\${out_summary}"

      # 4. Export a copy directly to your designated results directory
      metrics_dir="${params.batchDir}/r04_metrics"
      if [ -d "\${metrics_dir}" ]; then
          cp -a "\${out_summary}" "\${metrics_dir}/"
      fi

      rm -f targets.interval_list
      """
}
