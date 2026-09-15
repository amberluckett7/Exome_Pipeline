process CNVKIT {
    tag "${id}"
    cpus 4
    container 'docker://etal/cnvkit:0.9.10'
    containerOptions "-B ${params.batchDir} -B ${params.rsync} -B /lustre"

    input:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file(bam), file(bai)
      path ref_fasta
      path target_bed

    output:
      tuple val(id), file("${id}.cns"), file("${id}.call.cns"), file("${id}.bintest.cns"), file("${id}.diagram.pdf"), file("${id}.scatter.png"), optional: true

    script:
    """
    cnvkit.py batch ${bam} \\
        -r ${params.cnvref} \\
        --output-dir . \\
        --processes ${task.cpus} \\
        --diagram \\
        --scatter

    cns_file=\$(ls *.cns 2>/dev/null | grep -v call | grep -v bintest | head -1)
    call_file=\$(ls *.call.cns 2>/dev/null | head -1)
    bintest_file=\$(ls *.bintest.cns 2>/dev/null | head -1)
    diagram_file=\$(ls *-diagram.pdf 2>/dev/null | head -1)
    scatter_file=\$(ls *-scatter.png *-scatter.pdf 2>/dev/null | head -1)

    [ -n "\$cns_file" ]     && [ "\$cns_file" != "${id}.cns" ]             && mv "\$cns_file" "${id}.cns"
    [ -n "\$call_file" ]    && [ "\$call_file" != "${id}.call.cns" ]       && mv "\$call_file" "${id}.call.cns"
    [ -n "\$bintest_file" ] && [ "\$bintest_file" != "${id}.bintest.cns" ] && mv "\$bintest_file" "${id}.bintest.cns"
    [ -n "\$diagram_file" ] && [ "\$diagram_file" != "${id}.diagram.pdf" ] && mv "\$diagram_file" "${id}.diagram.pdf"
    [ -n "\$scatter_file" ] && [ "\$scatter_file" != "${id}.scatter.png" ] && mv "\$scatter_file" "${id}.scatter.png"

    DEST_DIR="${params.batchDir}/r04_cnvs"
    mkdir -p "\$DEST_DIR"
    cp -f ${id}.cns ${id}.call.cns ${id}.bintest.cns ${id}.diagram.pdf ${id}.scatter.png "\$DEST_DIR/" 2>/dev/null || true
    """
}
