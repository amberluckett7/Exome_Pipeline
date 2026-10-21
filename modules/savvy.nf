#!/usr/bin/env nextflow

process HOMOZYGOSITY {
    tag "${family[0]}"
    cpus 4
    memory '9 GB' 
    module 'BCFtools/1.17-GCC-12.2.0'
    
    input:
      tuple val(id), val(sex), val(family), val(famSampleCount), file(vcf), file(csi)
      
    output:
      tuple val(id), val(sex), val(family), val(famSampleCount), file("*_homozygosity.csv")
      
    script:
      """
      # Build index inside the container workspace
      tabix -p vcf ${vcf}

      # Single-array loop matching the genome pipeline
      declare -a ids=( ${id.join(" ")} )

      for i in "\${ids[@]}"; do
          # Optimized: Kept at -Xmx9g. Iterating sample-by-sample drops the RAM footprint dramatically
          java -Xmx9g -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 \\
              -cp ${params.javaDir}:${params.gatkJar} VcfToHomozygosity8 \\
              "${vcf}" "\$i" > "\${i}_homozygosity.csv"
      done

      DEST_DIR="${params.batchDir}/r04_metrics"
      mkdir -p "\$DEST_DIR"

      STAGE_DIR="\$DEST_DIR/.stage.\${ids[0]}.sampleHZ.\$\$"
      rm -rf "\$STAGE_DIR"
      mkdir -p "\$STAGE_DIR/.partial"

      FILES=( *_homozygosity.csv )

      for f in "\${FILES[@]}"; do
        [[ -s "\$f" ]] || { echo "Missing/empty output: \$f" >&2; exit 1; }
      done

      attempts=5
      delay=10
      ok=0

      for i in \$(seq 1 \$attempts); do
        rm -f "\$STAGE_DIR/"* 2>/dev/null || true
        mkdir -p "\$STAGE_DIR/.partial"

        if ${params.rsync} -a --checksum --delay-updates \\
            --partial --partial-dir=".partial" \\
            "\${FILES[@]}" "\$STAGE_DIR/"; then
          ok=1
          break
        fi
        sleep "\$delay"
      done

      [[ "\$ok" -eq 1 ]] || { echo "Publish failed: rsync never succeeded" >&2; exit 1; }

      for f in "\${FILES[@]}"; do
        mv -f "\$STAGE_DIR/\$f" "\$DEST_DIR/\$f"
      done
      rm -rf "\$STAGE_DIR"      
      """
}

process SHARED_HAPLOTYPES {
    tag "${family[0]}"
    module 'BCFtools/1.17-GCC-12.2.0'
    cpus 4
    memory '9 GB'
    
    input:
      tuple val(id), val(sex), val(family), val(famSampleCount), file(vcf), file(csi)
      
    output:
      tuple val(id), val(sex), val(family), val(famSampleCount), file("*_homozygosity.csv")
      
    script:
      """
      tabix -p vcf ${vcf}

      declare -a ids=( ${id.join(" ")} )

      for (( i=0; i < \${#ids[@]}; i++ )); do       
        for (( j=i+1; j < \${#ids[@]}; j++ )); do 
          java -Xmx9g -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 \\
              -cp "${params.javaDir}:${params.javaDir}/jars/*" VcfToHomozygosity8 \\
              "${vcf}" -p "\${ids[i]}" "\${ids[j]}" > "\${ids[i]}-\${ids[j]}_homozygosity.csv"
        done
      done
          
      DEST_DIR="${params.batchDir}/r04_metrics"
      mkdir -p "\$DEST_DIR"

      STAGE_DIR="\$DEST_DIR/.stage.\${ids[0]}.shared_haps.\$\$"
      rm -rf "\$STAGE_DIR"
      mkdir -p "\$STAGE_DIR/.partial"

      FILES=( *_homozygosity.csv )

      for f in "\${FILES[@]}"; do
        [[ -s "\$f" ]] || { echo "Missing/empty output: \$f" >&2; exit 1; }
      done

      attempts=5
      delay=10
      ok=0

      for i in \$(seq 1 \$attempts); do
        rm -f "\$STAGE_DIR/"* 2>/dev/null || true
        mkdir -p "\$STAGE_DIR/.partial"

        if ${params.rsync} -a --checksum --delay-updates \\
            --partial --partial-dir=".partial" \\
            "\text{\${FILES[@]}}" "\$STAGE_DIR/"; then
          ok=1
          break
        fi
        sleep "\$delay"
      done

      [[ "\$ok" -eq 1 ]] || { echo "Publish failed: rsync never succeeded" >&2; exit 1; }

      for f in "\text{\text{\${FILES[@]}}}"; do
        mv -f "\$STAGE_DIR/\$f" "\$DEST_DIR/\$f"
      done
      rm -rf "\$STAGE_DIR"    
      """
}

process UPD {
    tag "${family[0]}"
    module 'BCFtools/1.17-GCC-12.2.0'
    cpus 4
    memory '9 GB'
    
    input:
      tuple val(id), val(sex), val(family), val(famSampleCount), file(vcf), file(csi)
      
    output:
      tuple val(id), val(sex), val(family), val(famSampleCount), file("*_UPD.csv")
      
    script:
      """
      tabix -p vcf ${vcf}

      declare -a ids=( ${id.join(" ")} )

      for (( i=0; i < \${#ids[@]}; i++ )); do
        for (( j=0; j < \${#ids[@]}; j++ )); do
          [[ \$i -eq \$j ]] && continue

          java -Xmx9g -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 \\
              -cp ${params.javaDir}:${params.gatkJar} VcfToUPD2 \\
              "${vcf}" "\${ids[\${i}]}" "\${ids[\${j}]}" \\
              > "\${ids[\${i}]}-\${ids[\${j}]}_UPD.csv"
        done
      done

      DEST_DIR="${params.batchDir}/r04_metrics"
      mkdir -p "\$DEST_DIR"

      STAGE_DIR="\$DEST_DIR/.stage.\text{\${ids[0]}}.UPD.\$\$"
      rm -rf "\$STAGE_DIR"
      mkdir -p "\$STAGE_DIR/.partial"

      FILES=( *_UPD.csv )

      for f in "\${FILES[@]}"; do
        [[ -s "\$f" ]] || { echo "Missing/empty output: \$f" >&2; exit 1; }
      done

      attempts=5
      delay=10
      ok=0

      for i in \$(seq 1 \$attempts); do
        rm -f "\$STAGE_DIR/"* 2>/dev/null || true
        mkdir -p "\$STAGE_DIR/.partial"

        if ${params.rsync} -a --checksum --delay-updates \\
            --partial --partial-dir=".partial" \\
            "\text{\${FILES[@]}}" "\$STAGE_DIR/"; then
          ok=1
          break
        fi
        sleep "\$delay"
      done

      [[ "\$ok" -eq 1 ]] || { echo "Publish failed: rsync never succeeded" >&2; exit 1; }

      for f in "\text{\${FILES[@]}}"; do
        mv -f "\$STAGE_DIR/\$f" "\$DEST_DIR/\$f"
      done
      rm -rf "\$STAGE_DIR"
      """
}

process BATCH_HOMOZYGOSITY {
    tag "batch_${id[0]}"
    module 'BCFtools/1.17-GCC-12.2.0'
    cpus 4
    memory '9 GB'
    
    input:
      tuple val(id), val(sex), val(family), val(famSampleCount), file(vcf), file(csi)
      
    output:
      tuple val(id), val(sex), val(family), val(famSampleCount), file("homozygosity.csv")
      
    script:
      """
      vcfs=( ${vcf.join(" ")} )
      for i in \${vcfs[@]}; do
          tabix -p vcf \${i}
      done

      java -Xmx9g -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 \\
          -cp "${params.javaDir}:${params.javaDir}/jars/*" VcfToHomozygosity8 \\
          ${vcf.join(" ")} -com -dp 80 > homozygosity.csv

      DEST_DIR="${params.batchDir}/r04_metrics"
      mkdir -p "\$DEST_DIR"

      SRC_FILE="homozygosity.csv"
      TMP_FILE=".homozygosity.csv.partial.\$\$"

      attempts=5
      delay=10

      for i in \$(seq 1 \$attempts); do
        if ${params.rsync} -a --no-g --checksum --delay-updates --partial \\
            "\$SRC_FILE" "\$DEST_DIR/\$TMP_FILE" ; then
          mv -f "\$DEST_DIR/\$TMP_FILE" "\$DEST_DIR/\$SRC_FILE"
          break
        fi
        sleep "\$delay"
      done

      if [[ -e "\$DEST_DIR/\$TMP_FILE" ]]; then
        echo "Publish failed: temp file still present: \$DEST_DIR/\$TMP_FILE" >&2
        exit 1
      fi
      """
}
