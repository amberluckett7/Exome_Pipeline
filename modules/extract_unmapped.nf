#!/usr/bin/env nextflow

process EXTRACT_UNMAPPED {
    tag "${id}"
    module 'SAMtools/1.17-GCC-12.2.0'
    
    input:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file(bam), file(bai)
      
    output:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file("${id}_cytomegalovirus_1_samtools_fastq")
      
    script:
      """
      samtools fastq -f 0x4 -F 0xF00 ${bam} > ${id}_cytomegalovirus_1_samtools_fastq 
      """
}