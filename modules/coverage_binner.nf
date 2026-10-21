process COVERAGE_BINNER {
    tag "${id}"
    cpus 1
    container 'docker://amazoncorretto:21.0.7'
    containerOptions "-B ${params.javaDir} -B ${params.base} -B ${params.batchDir} -B ${params.rsync} -B /lustre/projects/Research_Project-MRC147594/software"

    input:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file(bam), file(bai)
      file(target_bed)

    output:
      tuple val(id), val(sex), file("${id}.coverageBinner")

    script:
    """
    java -cp ${params.javaDir}:${params.gatkJar} CoverageBinner ${bam} > ${id}.coverageBinner

    DEST_DIR="${params.batchDir}/r04_metrics"
    mkdir -p "\$DEST_DIR"
    cp -f ${id}.coverageBinner "\$DEST_DIR/"
    """
}
