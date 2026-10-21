process XY_COVERAGE {
    tag "batch_${ids[0]}-${ids[-1]}"
    cpus 1
    publishDir("${params.batchDir}/r04_metrics", mode: 'copy', overwrite: true, failOnError: true)
    container 'docker://amazoncorretto:21.0.7'
    containerOptions "-B ${params.javaDir} -B ${params.gatkJar} -B ${params.xy_bitmap} -B ${params.batchDir} -B ${params.rsync}"

    input:
      tuple val(ids), val(sex), path(coverageBinner, stageAs: "?/*")

    output:
      tuple val(ids), val(sex), file("XY_coverage")

    script:
      """
      id_flags="${ ids.withIndex()
                    .findAll { id, idx ->
                        def s = sex[idx]?.toString()?.toUpperCase()
                        s == 'MALE' || s == 'FEMALE'
                    }
                    .collect { id, idx ->
                        sex[idx]?.toString()?.toUpperCase() == 'MALE' ? "-m ${id}" : "-f ${id}"
                    }
                    .join(' ') }"

      # Find staged files across subdirectories created by stageAs: "?/*"
      staged_files=\$(find . -name "*.coverageBinner" | sort | tr '\n' ' ')

      if [[ -z "\${staged_files}" ]]; then
          echo "ERROR: No coverageBinner files found in work directory" >&2
          exit 1
      fi

      java -Xmx1g -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -cp ${params.javaDir}:${params.javaDir}/jars/htsjdk-4.2.0.jar XYCoverageBinned \\
        -bitmap ${params.xy_bitmap} \\
        -v exome \\
        \$id_flags \\
        \${staged_files} \\
        > XY_coverage
      """
}
