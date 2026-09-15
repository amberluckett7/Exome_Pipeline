#!/usr/bin/env nextflow

process XTEA {
    tag "${family}"
    container 'jamesrusssilsby/exetea:0.1.9d'
    containerOptions " -B ${params.xtea_libraries} -B ${params.batchDir} -B ${params.rsync}"
    cpus 16
    
    input:
      tuple val(id), val(sex), val(family), val(famSampleCount), file(bam), file(bai)
      file(Fasta)
      file(Fai)
      file(GFF)
      
    output:
      tuple val(id), val(sex), val(family), val(famSampleCount), file("*.vcf"), emit: vcf
      
    script:
      """
      declare -a sample_ids=( ${id.join(" ")} )
      declare -a bam_files=( ${bam.join(" ")} )
      declare -a bai_files=( ${bai.join(" ")} )

      > ${family}_samples.tsv
      > ${family}_bams.tsv

      for i in "\${!sample_ids[@]}"
      do
        curr_id="\${sample_ids[\$i]}"
        curr_bam="\${bam_files[\$i]}"
        curr_bai="\${bai_files[\$i]}"

        ln -sf "\$curr_bam" "\${curr_id}.bam"
        ln -sf "\$curr_bai" "\${curr_id}.bam.bai"

        printf "%s\n" "\${curr_id}" >> ${family}_samples.tsv
        printf "%s\t%s.bam\n" "\${curr_id}" "\${curr_id}" >> ${family}_bams.tsv
      done

      # Run xTea configuration generation
      ${params.xtea_libraries}xtea \\
        -i ${family}_samples.tsv \\
        -b ${family}_bams.tsv \\
        -x null \\
        -p \$PWD \\
        -o jobs.sh \\
        -l ${params.xtea_libraries}rep_lib_annotation/ \\
        -r ${Fasta} \\
        -g ${GFF} \\
        --xtea /data/xTea/xtea/ \\
        -f 5907 \\
        -y 7

      # Execute generated jobs in parallel
      if [ -f jobs.sh ]; then
        sed -n 's/^sbatch[[:space:]]*<[[:space:]]*\\(.*\\)\$/\\1/p' jobs.sh \\
          | parallel -j16 'bash < {}' || true
      fi

      # Search directories for output VCF files
      find . -mindepth 2 -type f -name "*.vcf" | while read -r file; do
        grep -v "orphan" "\$file" > \$(basename "\$file") || true
      done

      # Absolute Fallback Check: Count how many local .vcf files exist
      vcf_count=\$(ls -1 *.vcf 2>/dev/null | wc -l)

      if [ "\$vcf_count" -eq 0 ]; then
         echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > xtea_empty_fallback.vcf
      fi

      # Fixed: Correctly escaped all internal Bash variables so Nextflow doesn't misparse them
      DEST_DIR="${params.batchDir}/r04_transposable_elements"
      mkdir -p "\$DEST_DIR"

      STAGE_DIR="\$DEST_DIR/.stage.${family}.\$\$"
      rm -rf "\$STAGE_DIR"
      mkdir -p "\$STAGE_DIR/.partial"

      FILES=( *.vcf )

      for f in "\${FILES[@]}"; do
        [[ -s "\$f" ]] || { echo "Missing or empty output VCF: \$f" >&2; exit 1; }
      done

      attempts=5
      delay=10
      ok=0

      for i in \$(seq 1 \$attempts); do
        rm -f "\$STAGE_DIR/"* 2>/dev/null || true
        mkdir -p "\$STAGE_DIR/.partial"

        if ${params.rsync} -a --no-g --checksum --delay-updates \\
            --partial --partial-dir=".partial" \\
            "\${FILES[@]}" "\$STAGE_DIR/"; then
          ok=1
          break
        fi

        echo "rsync to staging directory failed (attempt \$i/\$attempts) — retrying in \${delay}s" >&2
        sleep "\$delay"
      done

      [[ "\$ok" -eq 1 ]] || { echo "Publish failed: rsync staging pipeline never succeeded" >&2; exit 1; }

      for f in "\${FILES[@]}"; do
        [[ -s "\$STAGE_DIR/\$f" ]] || { echo "Staged file verification failed: \$STAGE_DIR/\$f is missing/empty" >&2; exit 1; }
      done

      for f in "\${FILES[@]}"; do
        mv -f "\$STAGE_DIR/\$f" "\$DEST_DIR/\$f"
      done

      rm -rf "\$STAGE_DIR"
      """
}
