#!/usr/bin/env nextflow

process VEP {
    tag "${family}"
    cpus 16
    container 'docker://ensemblorg/ensembl-vep:release_114.1'
    containerOptions "-B ${params.vepData}:/data"
    
    input:
      // Modified input signature to accept output from VCF_FILTER
      tuple val(id), val(sex), val(family), val(famSampleCount), file(vcf), file(csi)
    
    output:
      tuple val(id), val(sex), val(family), val(famSampleCount), file("${family}_vep_annotated.vcf")
    
    script:
      """
      vep --dir /data --fork 16 --cache --offline --assembly GRCh38 --regulatory --fasta /data/homo_sapiens/114_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz --format vcf --vcf --force_overwrite \\
        --input_file ${vcf} \\
        --output_file ${family}_vep_annotated.vcf \\
        --dir_plugins /data/plugins/ \\
        --plugin AlphaMissense,file=/data/AlphaMissense_hg38.tsv.gz \\
        --plugin REVEL,file=/data/new_tabbed_revel_grch38.tsv.gz \\
        --plugin NMD \\
        --plugin UTRAnnotator,file=/data/uORF_5UTR_GRCh38_PUBLIC.txt \\
        --plugin LoF,loftee_path:/data/plugins/ \\
        --plugin LOEUF,file=/data/loeuf_dataset_grch38.tsv.gz,match_by=transcript \\
        --plugin SpliceAI,snv=/data/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/data/spliceai_scores.raw.indel.hg38.vcf.gz \\
        --plugin Enformer,file=/data/enformer_grch38.vcf.gz \\
        --plugin RiboseqORFs,file=/data/Ribo-seq_ORFs.bed.gz \\
        --custom file=/data/primate_ai_3d.vcf.gz,short_name=PrimateAI_3D,format=vcf,type=exact,coords=0,fields=score%percentile%prediction \\
        --custom file=/data/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN \\
        --custom file=/data/gnomad_v4_1_anno.vcf.gz,short_name=GnomAD_v4_1,format=vcf,type=exact,coords=0,fields=AF_afr%AC_afr%N_Hom_afr%AN_afr%AF_nfe%AC_nfe%N_Hom_nfe%AN_nfe%AF_eas%AC_eas%N_Hom_eas%AN_eas%AF_sas%AC_sas%N_Hom_sas%AN_sas%AF_amr%AC_amr%N_Hom_amr%AN_amr%AF_fin%AC_fin%N_Hom_fin%AN_fin%AF_asj%AC_asj%N_Hom_asj%AN_asj%AF_mid%AC_mid%N_Hom_mid%AN_mid%AF_ami%AC_ami%N_Hom_ami%AN_ami%AF_remaining%AC_remaining%N_Hom_remaining%AN_remaining%AF_popmax%AC_popmax%N_Hom_popmax%AN_popmax%AF_all%AC_all%N_Hom_all%AN_all \\
        --custom file=/data/allofus_250k_genomes.vcf.gz,short_name=AllofUs250k,format=vcf,type=exact,coords=0,fields=gvs_all_ac%gvs_all_an%gvs_all_af%gvs_max_af%gvs_max_ac%gvs_max_an%gvs_max_subpop%gvs_afr_ac%gvs_afr_an%gvs_afr_af%gvs_amr_ac%gvs_amr_an%gvs_amr_af%gvs_eas_ac%gvs_eas_an%gvs_eas_af%gvs_eur_ac%gvs_eur_an%gvs_eur_af%gvs_mid_ac%gvs_mid_an%gvs_mid_af%gvs_oth_ac%gvs_oth_an%gvs_oth_af%gvs_sas_ac%gvs_sas_an%gvs_sas_af \\
        --hgvs \\
        --mane  
      """
}