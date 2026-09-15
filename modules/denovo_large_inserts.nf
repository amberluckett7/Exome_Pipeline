#!/usr/bin/env nextflow

process DENOVO_LARGE_INSERTS {
    tag { trio_ids.join('-') }
    cpus 4
    container 'docker://amazoncorretto:21.0.7'
    containerOptions "-B ${params.javaDir} -B ${params.gatkJar} -B ${params.batchDir} -B ${params.rsync}"    
    
    input:
      // Swapped out val(id) and val(sex) for array-aware names (ids, sexes) 
      // to map cleanly to the list arrays from trio workflow logic.
      tuple val(ids), val(sexes), val(trio_ids), file(bam), file(bai)
      
    output:
      tuple val(ids), val(sexes), val(trio_ids), file("*.csv")
      
    script:
      """
      # Edited javaclass typo to keep mother and father in pipeline workflow
      java -Xmx9g -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 \\
        -cp ${params.javaDir}:${params.gatkJar}:${params.javaDir}/jars/htsjdk-4.2.0.jar FindLargeInsertSizesNF \\
        ${bam[0]} \\
        ${bam[1]} \\
        ${bam[2]} \\
        > ${trio_ids.join("-")}.csv

      DEST_DIR="${params.batchDir}/r04_denovolargeinserts"
      mkdir -p "\$DEST_DIR"

      SRC_FILE="${trio_ids.join("-")}.csv"
      TMP_FILE=".${trio_ids.join("-")}.csv.partial.\$\$"

      attempts=5
      delay=10

      for i in \$(seq 1 \$attempts); do
        # copy to temp name first
        if ${params.rsync} -a --no-g --checksum --delay-updates --partial \\
            "\$SRC_FILE" "\$DEST_DIR/\$TMP_FILE" ; then
          # then rename into place
          mv -f "\$DEST_DIR/\$TMP_FILE" "\$DEST_DIR/\$SRC_FILE"
          break
        fi
        echo "rsync failed (attempt \$i/\$attempts) — retrying in \${delay}s" >&2
        sleep "\$delay"
      done

      # if temp still exists, we never successfully promoted to final name
      if [[ -e "\$DEST_DIR/\$TMP_FILE" ]]; then
        echo "Publish failed: temp file still present: \$DEST_DIR/\$TMP_FILE" >&2
        exit 1
      fi
      """
}