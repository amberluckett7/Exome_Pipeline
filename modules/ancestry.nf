process ANCESTRY {
    tag "${id}"
    cpus 16
    container 'jamesrusssilsby/gnomadtools:latest'
    containerOptions "-B ${params.resourcesDir}/gnomad_pca -B ${params.batchDir} -B ${params.rsync}"

    input:
      tuple val(id), val(sex), val(family), val(trio), val(famSampleCount), file(gvcf), file(csi)

    output:
      file("${id}.pca_scores.tsv")
      file("${id}.ancestry_assignment.tsv")

    script:
    """
    mkdir -p \$PWD/temp
    mkdir -p \$PWD/ivy_cache
    export TMPDIR=\$PWD/temp

    python3 <<EOF
import os
import hail as hl
import onnx
from hail.vds.combiner import transform_gvcf
from gnomad.sample_qc.ancestry import (
    apply_onnx_classification_model,
    assign_population_pcs,
)

read_if_exists = True
v3_num_pcs = 16
v3_min_prob = 0.75

work = os.getcwd()
tmp_dir = f"{work}/temp"
ivy_dir = f"{work}/ivy_cache"
os.makedirs(tmp_dir, exist_ok=True)
os.makedirs(ivy_dir, exist_ok=True)
os.environ["TMPDIR"] = tmp_dir

hl.init(
    tmp_dir=tmp_dir,
    spark_conf={
        "spark.driver.memory": "8g",
        "spark.executor.memory": "8g",
        "spark.local.dir": tmp_dir,
        "spark.jars.ivy": ivy_dir,
        "spark.driver.maxResultSize": "0",
        "spark.sql.shuffle.partitions": "64",
        "spark.ui.enabled": "false",
        "spark.eventLog.enabled": "false",
    },
)

hl.default_reference('GRCh38')

gnomad_v3_loadings = "${params.resourcesDir}/gnomad_pca/gnomad.v3.1.pca_loadings.ht"
gnomad_v3_onnx_rf  = "${params.resourcesDir}/gnomad_pca/gnomad.v3.1.RF_fit.onnx"

with open(gnomad_v3_onnx_rf, "rb") as f:
    v3_onx_fit = onnx.load(f)

mt_output_path         = os.path.join(work, "${id}_gnomad_v3.1_ancestry_rf.mt")
scores_output_path     = os.path.join(work, "${id}_gnomad_v3.1_ancestry_rf.scores.ht")
gnomad_assignment_path = os.path.join(work, "${id}_gnomad_v3.1_ancestry_rf.assignment.ht")

v3_loading_ht = hl.read_table(gnomad_v3_loadings)

# Our BAMs/gVCFs use no-'chr' contig naming (e.g. '1', 'X', 'M'), but Hail's
# built-in GRCh38 reference genome and the gnomAD v3.1 PCA loadings both
# expect 'chr'-prefixed contigs (e.g. 'chr1', 'chrX', 'chrM'). Recode on
# import so loci match the reference/loadings without altering the gvcf file.
contig_recoding = {str(i): f"chr{i}" for i in range(1, 23)}
contig_recoding.update({"X": "chrX", "Y": "chrY", "M": "chrM", "MT": "chrM"})

mt = hl.import_vcf(
    "${gvcf}",
    reference_genome='GRCh38',
    array_elements_required=False,
    force_bgz=True,
    contig_recoding=contig_recoding
)

sample_vds = transform_gvcf(mt, reference_entry_fields_to_keep=["LA", "LGT", "GQ", "DP", "LAD"])
sample_vds = hl.vds.split_multi(sample_vds, filter_changed_loci=True)
sample_vds = hl.vds.filter_variants(sample_vds, v3_loading_ht)
mt = hl.vds.to_dense_mt(sample_vds)

mt = mt.checkpoint(
    mt_output_path,
    overwrite=not read_if_exists,
    _read_if_exists=read_if_exists
)

v3_pcs_ht = hl.experimental.pc_project(
    mt.GT,
    v3_loading_ht.loadings,
    v3_loading_ht.pca_af,
)

v3_pcs_ht = v3_pcs_ht.checkpoint(
    scores_output_path,
    overwrite=not read_if_exists,
    _read_if_exists=read_if_exists,
)

ht, model = assign_population_pcs(
    v3_pcs_ht,
    pc_cols=v3_pcs_ht.scores[:v3_num_pcs],
    fit=v3_onx_fit,
    min_prob=v3_min_prob,
    apply_model_func=apply_onnx_classification_model,
)

ht = ht.checkpoint(
    gnomad_assignment_path,
    overwrite=not read_if_exists,
    _read_if_exists=read_if_exists,
)

v3_pcs_ht.export("${id}.pca_scores.tsv", header=True, delimiter="\t")
ht.export("${id}.ancestry_assignment.tsv", header=True, delimiter="\t")
EOF

    DEST_DIR="${params.batchDir}/r04_metrics"
    mkdir -p "\$DEST_DIR"

    STAGE_DIR="\$DEST_DIR/.stage.${id}.\$\$"
    rm -rf "\$STAGE_DIR"
    mkdir -p "\$STAGE_DIR/.partial"

    FILES=( "${id}.pca_scores.tsv" "${id}.ancestry_assignment.tsv" )

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
