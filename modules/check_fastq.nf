#!/usr/bin/env nextflow
process CHECK_FASTQ {
    tag "${meta.id}"
    container 'docker://amazoncorretto:21.0.7'
    // containerOptions edit to bind mount the correct raw fastq directory path
    containerOptions "-B ${params.javaDir} -B ${params.base} -B ${params.batchDir} -B ${params.batchDir} -B ${params.rsync}"

    input:
      tuple val(meta), file(fastq1), file(fastq2)

    output:
      tuple val(meta.id), val(meta.sex), val(meta.family), val(meta.famSampleCount), file("${meta.id}_checkFastq.txt")

    script:
    """
    mkdir -p tmp

    java -Xmx4g -XX:ParallelGCThreads=4 -XX:ConcGCThreads=4 -Djava.io.tmpdir=\$PWD/tmp -cp ${params.javaDir} CheckFastq2 ${ fastq1.withIndex().collect { f1, idx -> "${f1} ${fastq2[idx]}" }.join(' ') } \\
    > ${meta.id}_checkFastq.txt

    DEST_DIR="${params.batchDir}/r04_metrics"
    mkdir -p "\$DEST_DIR"

    SRC_FILE="${meta.id}_checkFastq.txt"
    TMP_FILE=".${meta.id}_checkFastq.txt.partial.\$\$"

    attempts=5
    delay=10

    for i in \$(seq 1 \$attempts); do
        if ${params.rsync} -a --no-g --checksum --delay-updates --partial "\$SRC_FILE" "\$DEST_DIR/\$TMP_FILE" ; then
            mv -f "\$DEST_DIR/\$TMP_FILE" "\$DEST_DIR/\$SRC_FILE"
            break
        fi
        echo "rsync failed (attempt \$i/\$attempts) — retrying in \${delay}s" >&2
        sleep "\$delay"
    done

    if [[ -e "\$DEST_DIR/\$TMP_FILE" ]]; then
        echo "Publish failed: temp file still present: \$DEST_DIR/\$TMP_FILE" >&2
        exit 1
    fi 
    """
}
