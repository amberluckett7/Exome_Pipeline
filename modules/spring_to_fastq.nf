#!/usr/bin/env nextflow
process SPRING2FQ {
    cpus 8
    tag "${id}"
    stageInMode 'copy'
    container 'docker://jamesrusssilsby/spring:latest'

    input:
      tuple val(id), val(platform), val(sex), val(family), val(trio),
            val(flowcell), val(laneCount), val(famSampleCount),
            file(spring)

    output:
      tuple val(id), val(platform), val(sex), val(family), val(trio),
            val(flowcell), val(laneCount), val(famSampleCount),
            file("*_1.fq.gz"), file("*_2.fq.gz")

    script:
      """
      spring -d -i ${spring} -g -o \$(basename ${spring} .spring)_1.fq.gz \$(basename ${spring} .spring)_2.fq.gz
      """
      
}