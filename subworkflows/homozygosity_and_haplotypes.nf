#!/usr/bin/env nextflow

include { HOMOZYGOSITY } from '../modules/savvy.nf'
include { SHARED_HAPLOTYPES } from '../modules/savvy.nf'
include { UPD } from '../modules/savvy.nf'

workflow HOMOZYGOSITY_AND_HAPLOTYPES {
    take:
    ch_mixed_vcf
    statusById

    main:

    // Filter and normalize for individual sample homozygosity
    ch_mixed_vcf
        .filter { ids, sex, family, n, vcf, vcfcsi ->
            // Resolve key whether ids is a single string value or an array list
            def rep = (ids instanceof List ? ids[0] : ids)
            return statusById[rep]?.sample_homoz_needed
        }
        .map { ids, sex, family, n, vcf, vcfcsi ->
            // Ensure elements are wrapped into lists to satisfy savvy.nf loop mechanics
            def safeIds  = (ids instanceof List)  ? ids  : [ids]
            def safeSex  = (sex instanceof List)  ? sex  : [sex]
            def safeFam  = (family instanceof List) ? family : [family]
            def safeVcf  = (vcf instanceof List)  ? vcf  : [vcf]
            def safeCsi  = (vcfcsi instanceof List) ? vcfcsi : [vcfcsi]
            tuple(safeIds, safeSex, safeFam, n, safeVcf, safeCsi)
        }
        .set { ch_HZ_vcf }   

    HOMOZYGOSITY(ch_HZ_vcf)

    // Filter out singletons (only process families with more than 1 member)
    ch_mixed_vcf
        .filter { row -> 
            // Handle both scalar integers and lists safely
            def count = (row[3] instanceof List) ? row[3][0] : row[3]
            return count > 1 
        }
        .set { ch_multisample_vcf }

    // Filter and normalize for shared haplotypes
    ch_multisample_vcf
        .filter { ids, sex, family, n, vcf, vcfcsi ->
            def rep = (ids instanceof List ? ids[0] : ids)
            return statusById[rep]?.shared_haps_needed
        }
        .map { ids, sex, family, n, vcf, vcfcsi ->
            def safeIds  = (ids instanceof List)  ? ids  : [ids]
            def safeSex  = (sex instanceof List)  ? sex  : [sex]
            def safeFam  = (family instanceof List) ? family : [family]
            def safeVcf  = (vcf instanceof List)  ? vcf  : [vcf]
            def safeCsi  = (vcfcsi instanceof List) ? vcfcsi : [vcfcsi]
            tuple(safeIds, safeSex, safeFam, n, safeVcf, safeCsi)
        }
        .set { ch_haps_vcf }   

    SHARED_HAPLOTYPES(ch_haps_vcf)

    // Filter and normalize for Uniparental Disomy (UPD)
    ch_multisample_vcf
        .filter { ids, sex, family, n, vcf, vcfcsi ->
            def rep = (ids instanceof List ? ids[0] : ids)
            return statusById[rep]?.upd_needed
        }
        .map { ids, sex, family, n, vcf, vcfcsi ->
            def safeIds  = (ids instanceof List)  ? ids  : [ids]
            def safeSex  = (sex instanceof List)  ? sex  : [sex]
            def safeFam  = (family instanceof List) ? family : [family]
            def safeVcf  = (vcf instanceof List)  ? vcf  : [vcf]
            def safeCsi  = (vcfcsi instanceof List) ? vcfcsi : [vcfcsi]
            tuple(safeIds, safeSex, safeFam, n, safeVcf, safeCsi)
        }
        .set { ch_UPD_vcf }   

    UPD(ch_UPD_vcf)
}