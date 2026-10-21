#!/usr/bin/env nextflow

process SURVINDEL2 {
    cpus 16
    tag { id }
    container 'jamescraufurd/survindel2:latest'
    containerOptions("-B /usr/lib/locale/:/usr/lib/locale/ -B ${params.survindel2model} -B ${params.batchDir} -B ${params.rsync}")

    input:
      tuple val(id), val(sex), val(trio), val(family), val(famSampleCount), file(bam), file(bai)
      file(fasta)
      file(fai)
      file(target_bed) 

    output:
      tuple val(id), val(sex), val(family), val(famSampleCount), file("${id}.survindel.cnv.vcf.gz")
      tuple val(id), val(sex), val(family), val(famSampleCount), file("${id}.survindel.pass-ml.cnv.vcf.gz")
      tuple val(id), val(sex), val(family), val(famSampleCount), file("${id}.genotyped.cnv.vcf.gz")
      tuple val(id), val(sex), val(family), val(famSampleCount), file("${id}.genotyped.pass-ml.cnv.vcf.gz")

    script:
      """
      mkdir survindel_rundir
      mkdir survtyper_base_rundir
      mkdir survtyper_ml_rundir

      python3 /survindel2/survindel2.py \\
        --threads ${task.cpus} \\
        ${bam} \\
        survindel_rundir \\
        ${fasta}

      # Run filtering 
      python3 /survindel2/run_classifier.py \\
        survindel_rundir/out.vcf.gz \\
        survindel_rundir/out.pass-ml.vcf.gz \\
        survindel_rundir/stats.txt \\
        ALL \\
        ${params.survindel2model}

      mv survindel_rundir/out.pass.vcf.gz ${id}.survindel.cnv.vcf.gz
      mv survindel_rundir/out.pass-ml.vcf.gz ${id}.survindel.pass-ml.cnv.vcf.gz

      # Genotype discovered sites back against the BAM
      python3 /survtyper/survtyper.py \\
        --threads ${task.cpus} \\
        ${id}.survindel.cnv.vcf.gz \\
        ${bam} \\
        survtyper_base_rundir \\
        ${fasta}

      python3 /survtyper/survtyper.py \\
        --threads ${task.cpus} \\
        ${id}.survindel.pass-ml.cnv.vcf.gz \\
        ${bam} \\
        survtyper_ml_rundir \\
        ${fasta}

      mv survtyper_base_rundir/genotyped.vcf.gz ${id}.genotyped.cnv.vcf.gz
      mv survtyper_ml_rundir/genotyped.vcf.gz ${id}.genotyped.pass-ml.cnv.vcf.gz

      DEST_DIR="${params.batchDir}/r04_survindel2"
      mkdir -p "\$DEST_DIR"

      STAGE_DIR="\$DEST_DIR/.stage.${id}.\$\$"
      rm -rf "\$STAGE_DIR"
      mkdir -p "\$STAGE_DIR/.partial"

      FILES=( "${id}.survindel.cnv.vcf.gz" "${id}.survindel.pass-ml.cnv.vcf.gz" "${id}.genotyped.cnv.vcf.gz" "${id}.genotyped.pass-ml.cnv.vcf.gz" )

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

      # Promote files to production directory
      for f in "\${FILES[@]}"; do
        mv -f "\$STAGE_DIR/\$f" "\$DEST_DIR/\$f"
      done

      rm -rf "\$STAGE_DIR"
      """
}
