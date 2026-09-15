#!/usr/bin/env nextflow

process MANTA {
    cpus 16
    memory '100G'
    tag { family } 
    container 'https://depot.galaxyproject.org/singularity/manta:1.6.0--h9ee0642_1'
    containerOptions("-B /usr/lib/locale/:/usr/lib/locale/ -B ${params.batchDir} -B ${params.rsync} -B ${params.base}")

    input:
      tuple val(id), val(sex), val(family), val(famSampleCount), file(bam), file(bai)
      file(fasta)
      file(fai)
      file(mantaBed)
      file(mantaBedTbi)

    output:
      tuple val(id), val(sex), val(family), val(famSampleCount), file("*.SV.vcf.gz"), file("*.SV.vcf.gz.tbi")

    script:
      def inBamStr = bam.collect { "--bam ${it}" }.join(' ')
      def prefix = id instanceof List ? id.join("-") : id
      

      """
      mkdir manta_rundir

      configManta.py \\
        --reference ${fasta} \\
        --runDir manta_rundir \\
        --exome \\
        --callRegions ${mantaBed} \\
        ${inBamStr}

      # Mitigate known Manta local smtplib timeout hang bug
      sed -i "s/smtplib\\.SMTP('localhost')/smtplib.SMTP('localhost', timeout=2)/" manta_rundir/runWorkflow.py
      
      python manta_rundir/runWorkflow.py -j ${task.cpus} -m local -g 100
      
      mv manta_rundir/results/variants/diploidSV.vcf.gz ${prefix}.SV.vcf.gz
      mv manta_rundir/results/variants/diploidSV.vcf.gz.tbi ${prefix}.SV.vcf.gz.tbi

      # Publish output
      DEST_DIR="${params.batchDir}/r04_manta"
      mkdir -p "\$DEST_DIR"

      STAGE_DIR="\$DEST_DIR/.stage.${prefix}.\$\$"
      rm -rf "\$STAGE_DIR"
      mkdir -p "\$STAGE_DIR/.partial"

      FILES=( ${prefix}.SV.vcf.gz ${prefix}.SV.vcf.gz.tbi )

      for f in "\${FILES[@]}"; do
        [[ -s "\$f" ]] || { echo "Missing/empty output: \$f" >&2; exit 1; }
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

        echo "rsync to stage failed (attempt \$i/\$attempts) — retrying in \${delay}s" >&2
        sleep "\$delay"
      done

      [[ "\$ok" -eq 1 ]] || { echo "Publish failed: rsync never succeeded" >&2; exit 1; }

      for f in "\${FILES[@]}"; do
        [[ -s "\$STAGE_DIR/\$f" ]] || { echo "Stage missing/empty: \$STAGE_DIR/\$f" >&2; exit 1; }
      done

      for f in "\${FILES[@]}"; do
        mv -f "\$STAGE_DIR/\$f" "\$DEST_DIR/\$f"
      done

      rm -rf "\$STAGE_DIR"
      """
}
