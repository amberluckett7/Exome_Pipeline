#!/usr/bin/env python3
"""
gene_panel_report.py

Filters VEP-annotated SNV and SV VCFs against a gene panel with 
inheritance-aware thresholds, then produces a per-sample HTML report
with separate sections for Dominant and Recessive findings.

Usage example:
    python3 gene_panel_report.py \
        --snv   Sample1234_vep_annotated.vcf.gz \
        --sv    Sample1234_vep_annotated.SV.vcf.gz \
        --panel monogenic_gene_panel.tsv \
        --sample Sample1234 \
        --out   Sample1234_gene_panel_report.html

Dependencies:
    bcftools >= 1.17 (must be on PATH)
    python3 with no non-stdlib requirements
"""

import argparse
import subprocess
import sys
import os
import re
import html
from collections import defaultdict

# ---------------------------------------------------------------------------
# Filtering thresholds
# ---------------------------------------------------------------------------
DOM_AC_MAX   = 3     # GnomAD4 AC < 3  for dominant genes
REC_AC_MAX   = 50    # GnomAD4 AC < 50 for recessive genes
SPLICEAI_MIN = 0.5   # SpliceAI delta score threshold
SPLICE_UPSTREAM   = 50   # bases upstream of splice site
SPLICE_DOWNSTREAM = 10   # bases downstream of splice site

# Consequence terms that count as coding/non-synonymous
CODING_NONSYN = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_region_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",   # included only if SpliceAI passes
}

# Consequences to always include regardless of synonymous status
ALWAYS_INCLUDE = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
}

# Rank used to pick the "worst"/most severe transcript row when no
# MANE Select transcript is available for a variant
IMPACT_RANK = {'HIGH': 4, 'MODERATE': 3, 'LOW': 2, 'MODIFIER': 1}

# Display order + labels for HTML sections
INHERITANCE_SECTIONS = [
    ('dominant',  'Dominant'),
    ('recessive', 'Recessive'),
]

# ---------------------------------------------------------------------------
# Parse gene panel TSV
# ---------------------------------------------------------------------------
def load_panel(panel_tsv):
    """Returns dict: gene_symbol -> {'inheritance': str, 'phenotype': str, ...}"""
    panel = {}
    with open(panel_tsv) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            if not parts[0]:
                continue
            gene = parts[0].strip()
            row = dict(zip(header, parts))
            inh_raw = row.get('Inheritance', '').strip().lower()
            # Normalise inheritance
            if 'dominant or recessive' in inh_raw or 'dominant and recessive' in inh_raw:
                inh = 'both'
            elif 'x-linked' in inh_raw or 'x linked' in inh_raw:
                inh = 'recessive'   # treat X-linked as recessive for AC/biallelic logic
            elif 'dominant' in inh_raw:
                inh = 'dominant'
            else:
                inh = 'recessive'
            panel[gene] = {
                'inheritance': inh,
                'phenotype':   row.get('Phenotype', ''),
                'transcript':  row.get('MANE_RefSeq_Transcript', ''),
                'notes':       row.get('Notes', ''),
            }
    return panel

# ---------------------------------------------------------------------------
# VCF field indices (from your bcftools +split-vep -l output)
# ---------------------------------------------------------------------------
# SNV VCF CSQ field indices
SNV_IDX = {
    'Consequence':        1,
    'IMPACT':             2,
    'SYMBOL':             3,
    'Gene':               4,
    'Feature':            6,
    'BIOTYPE':            7,
    'EXON':               8,
    'INTRON':             9,
    'HGVSc':             10,
    'HGVSp':             11,
    'MANE_SELECT':       24,
    'SpliceAI_DS_AG':    51,
    'SpliceAI_DS_AL':    52,
    'SpliceAI_DS_DG':    53,
    'SpliceAI_DS_DL':    54,
    'SpliceAI_DP_AG':    47,
    'SpliceAI_DP_AL':    48,
    'SpliceAI_DP_DG':    49,
    'SpliceAI_DP_DL':    50,
    'ClinVar_CLNSIG':    72,
    'ClinVar_CLNDN':     74,
    'GnomAD_AC_all':    121,
    'GnomAD_AC_popmax': 117,
    'GnomAD_AF_all':    120,
    'LoF':               42,
    'REVEL':             35,
    'am_class':          33,
}

# SV VCF CSQ field indices (fewer fields — no gnomAD/SpliceAI)
SV_IDX = {
    'Consequence':  1,
    'IMPACT':       2,
    'SYMBOL':       3,
    'Gene':         4,
    'Feature':      6,
    'BIOTYPE':      7,
    'HGVSc':       10,
    'HGVSp':       11,
    'MANE_SELECT': 24,
}

# ---------------------------------------------------------------------------
# Transcript-row deduplication
# ---------------------------------------------------------------------------
# bcftools +split-vep -d emits one row PER TRANSCRIPT annotation for a given
# variant, so the same variant (same CHROM/POS/REF/ALT, same AC/AF/GT) shows
# up multiple times with only the transcript-level fields (Consequence,
# HGVSc, HGVSp, MANE_SELECT, IMPACT) differing. We collapse these down to a
# single representative row per variant before filtering, so AC/AF and the
# variant itself are only counted/reported once.
def pick_best_row(rows, idx):
    """
    Given multiple parsed rows (same CHROM/POS/REF/ALT, different
    transcripts), pick the single best representative row:
      1. Prefer a row where MANE_SELECT is non-empty (canonical transcript)
      2. Otherwise, prefer the row with the highest-impact consequence
      3. Fall back to the first row if nothing distinguishes them
    """
    if len(rows) == 1:
        return rows[0]

    mane_idx = idx['MANE_SELECT']
    impact_idx = idx['IMPACT']

    def has_mane(r):
        csq = r[-1]
        return len(csq) > mane_idx and csq[mane_idx].strip() not in ('', '.')

    mane_rows = [r for r in rows if has_mane(r)]
    if mane_rows:
        return mane_rows[0]

    def impact_rank(r):
        csq = r[-1]
        if len(csq) <= impact_idx:
            return 0
        return IMPACT_RANK.get(csq[impact_idx].strip().upper(), 0)

    return max(rows, key=impact_rank)


def dedupe_by_variant(raw_rows, idx, key_len):
    """
    raw_rows: list of tuples, where the first key_len elements are the
    variant-identifying fields (e.g. chrom, pos, ref, alt for SNVs, or
    chrom, pos, end, svtype for SVs) and the last element is the csq list.

    Returns one row per unique variant key, using pick_best_row to choose
    the representative transcript annotation.
    """
    groups = defaultdict(list)
    for row in raw_rows:
        key = row[:key_len]
        groups[key].append(row)

    return [pick_best_row(rows, idx) for rows in groups.values()]

# ---------------------------------------------------------------------------
# Is this variant already ClinVar pathogenic / likely pathogenic?
# ---------------------------------------------------------------------------
# Population AC thresholds (DOM_AC_MAX / REC_AC_MAX) are arbitrary cutoffs meant to
# screen out common variation. We bypass the AC filter (but
# NOT the consequence/SpliceAI filter) for variants already called in ClinVar.
#
# NOTE: this is a simple substring-based check on the raw ClinVar CLNSIG
# string and is intentionally conservative, It requires "pathogenic" to be
# present without "benign", and explicitly excludes "conflicting" calls
# (e.g. "Conflicting_interpretations_of_pathogenicity" contains the
# substring "pathogenic" but is NOT a pathogenic call and must not bypass
# the filter). If CLNSIG formatting changes upstream (e.g. multiple
# ClinVar submissions concatenated with different delimiters), re-check
# this logic against real output before trusting it blindly.
def is_clinvar_pathogenic(csq, idx):
    i = idx.get('ClinVar_CLNSIG')
    if i is None or len(csq) <= i:
        return False
    val = csq[i].strip()
    if not val or val == '.':
        return False
    val_lower = val.lower()
    if 'conflicting' in val_lower:
        return False
    return 'pathogenic' in val_lower and 'benign' not in val_lower

# ---------------------------------------------------------------------------
# Process SNV VCF
# ---------------------------------------------------------------------------
def process_snv_vcf(vcf_path, panel, sample_id):
    """
    Returns list of dicts, one per passing variant row (deduplicated across
    transcript annotations first).
    For recessive genes, also returns raw het dict for compound-het detection.
    """
    fmt = '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%GQ]\t[%DP]\t%CSQ\n'
    idx = SNV_IDX

    cmd = [
        'bcftools', '+split-vep',
        '-d', '-A', 'tab',
        '-f', fmt,
        vcf_path
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL, text=True)

    raw_rows = []
    for line in proc.stdout:
        line = line.rstrip('\n')
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) < 8:
            continue

        chrom, pos, ref, alt, gt, gq, dp = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6]
        csq = parts[7:]

        if len(csq) <= max(idx.values()):
            continue

        # key: chrom, pos, ref, alt  ->  rest carried through for later use
        raw_rows.append((chrom, pos, ref, alt, gt, gq, dp, csq))
    proc.wait()

    # Collapse multiple transcript rows per variant down to one representative row
    deduped_rows = dedupe_by_variant(raw_rows, idx, key_len=4)

    passing  = []
    raw_hets = defaultdict(list)   # gene -> list of variant dicts (for comp-het)

    for chrom, pos, ref, alt, gt, gq, dp, csq in deduped_rows:
        symbol = csq[idx['SYMBOL']].strip()
        if symbol not in panel:
            continue

        gene_info = panel[symbol]
        inh = gene_info['inheritance']

        csq_terms   = csq[idx['Consequence']]
        spliceai    = max_spliceai(csq, idx)
        ac          = get_ac(csq, idx)
        is_path     = is_clinvar_pathogenic(csq, idx)

        if not passes_consequence(csq_terms, spliceai):
            continue

        # Determine which AC threshold applies
        thresholds = []
        if inh in ('dominant', 'both'):
            thresholds.append(('dominant', DOM_AC_MAX))
        if inh in ('recessive', 'both'):
            thresholds.append(('recessive', REC_AC_MAX))

        for mode, ac_max in thresholds:
            # AC filter: if AC is known and exceeds threshold, skip —
            # UNLESS this variant already has a ClinVar pathogenic/likely
            # pathogenic call, in which case population frequency alone
            # should not override a pathogenic classification.
            if ac is not None and ac >= ac_max and not is_path:
                continue

            # Genotype check
            gt_clean = gt.replace('|', '/').strip()
            is_hom   = gt_clean in ('1/1', '1|1')
            is_het   = gt_clean in ('0/1', '1/0', '0|1', '1|0')

            var = {
                'type':        'SNV/Indel',
                'chrom':       chrom,
                'pos':         pos,
                'ref':         ref,
                'alt':         alt,
                'gene':        symbol,
                'inheritance': mode,
                'phenotype':   gene_info['phenotype'],
                'consequence': csq_terms,
                'impact':      csq[idx['IMPACT']],
                'hgvsc':       csq[idx['HGVSc']],
                'hgvsp':       csq[idx['HGVSp']],
                'biotype':     csq[idx['BIOTYPE']],
                'mane':        csq[idx['MANE_SELECT']],
                'gt':          gt_clean,
                'gq':          gq,
                'dp':          dp,
                'spliceai':    f'{spliceai:.3f}' if spliceai is not None else '.',
                'ac_gnomad':   str(ac) if ac is not None else '.',
                'af_gnomad':   csq[idx['GnomAD_AF_all']] if len(csq) > idx['GnomAD_AF_all'] else '.',
                'clinvar':     csq[idx['ClinVar_CLNSIG']] if len(csq) > idx['ClinVar_CLNSIG'] else '.',
                'clinvar_dn':  csq[idx['ClinVar_CLNDN']] if len(csq) > idx['ClinVar_CLNDN'] else '.',
                'lof':         csq[idx['LoF']] if len(csq) > idx['LoF'] else '.',
                'revel':       csq[idx['REVEL']] if len(csq) > idx['REVEL'] else '.',
                'am_class':    csq[idx['am_class']] if len(csq) > idx['am_class'] else '.',
            }
            if mode == 'dominant':
                if is_het or is_hom:
                    passing.append(var)

            elif mode == 'recessive':
                if is_hom:
                    var['zygosity'] = 'Homozygous'
                    passing.append(var)
                elif is_het:
                    var['zygosity'] = 'Heterozygous'
                    raw_hets[symbol].append(var)

    # Compound het: 2+ het variants in the same recessive gene
    for gene, hets in raw_hets.items():
        if len(hets) >= 2:
            for v in hets:
                v['zygosity'] = 'Compound het'
                passing.append(v)

    return passing

# ---------------------------------------------------------------------------
# Process SV VCF
# ---------------------------------------------------------------------------
def process_sv_vcf(vcf_path, panel, sample_id):
    """Returns list of dicts for SVs in panel genes (deduplicated across transcripts)."""
    fmt = '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t[%GT]\t[%GQ]\t%CSQ\n'
    idx = SV_IDX

    cmd = [
        'bcftools', '+split-vep',
        '-d', '-A', 'tab',
        '-f', fmt,
        vcf_path
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL, text=True)

    raw_rows = []
    for line in proc.stdout:
        line = line.rstrip('\n')
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) < 8:
            continue

        chrom, pos, end, svtype, svlen, gt, gq = (
            parts[0], parts[1], parts[2], parts[3], parts[4], parts[5], parts[6])
        csq = parts[7:]

        if len(csq) <= max(idx.values()):
            continue

        # key: chrom, pos, end, svtype -> rest carried through
        raw_rows.append((chrom, pos, end, svtype, svlen, gt, gq, csq))
    proc.wait()

    deduped_rows = dedupe_by_variant(raw_rows, idx, key_len=4)

    passing  = []
    raw_hets = defaultdict(list)

    for chrom, pos, end, svtype, svlen, gt, gq, csq in deduped_rows:
        symbol = csq[idx['SYMBOL']].strip()
        if symbol not in panel:
            continue

        gene_info = panel[symbol]
        inh = gene_info['inheritance']
        csq_terms = csq[idx['Consequence']]

        thresholds = []
        if inh in ('dominant', 'both'):
            thresholds.append(('dominant', DOM_AC_MAX))
        if inh in ('recessive', 'both'):
            thresholds.append(('recessive', REC_AC_MAX))

        for mode, _ in thresholds:
            gt_clean = gt.replace('|', '/').strip()
            is_hom = gt_clean in ('1/1', '1|1')
            is_het = gt_clean in ('0/1', '1/0', '0|1', '1|0')
            var = {
                'type':        f'SV ({svtype})',
                'chrom':       chrom,
                'pos':         pos,
                'end':         end,
                'svtype':      svtype,
                'svlen':       svlen,
                'gene':        symbol,
                'inheritance': mode,
                'phenotype':   gene_info['phenotype'],
                'consequence': csq_terms,
                'impact':      csq[idx['IMPACT']],
                'hgvsc':       csq[idx['HGVSc']],
                'hgvsp':       csq[idx['HGVSp']],
                'biotype':     csq[idx['BIOTYPE']],
                'mane':        csq[idx['MANE_SELECT']],
                'gt':          gt_clean,
                'gq':          gq,
                'dp':          '.',
                'spliceai':    '.',
                'ac_gnomad':   '.',
                'af_gnomad':   '.',
                'clinvar':     '.',
                'clinvar_dn':  '.',
                'lof':         '.',
                'revel':       '.',
                'am_class':    '.',
            }

            if mode == 'dominant':
                if is_het or is_hom:
                    passing.append(var)
            elif mode == 'recessive':
                if is_hom:
                    var['zygosity'] = 'Homozygous'
                    passing.append(var)
                elif is_het:
                    var['zygosity'] = 'Heterozygous'
                    raw_hets[symbol].append(var)

    for gene, hets in raw_hets.items():
        if len(hets) >= 2:
            for v in hets:
                v['zygosity'] = 'Compound het'
                passing.append(v)

    return passing

# ---------------------------------------------------------------------------
# Helper: max SpliceAI delta score across AG/AL/DG/DL
# ---------------------------------------------------------------------------
def max_spliceai(csq, idx):
    scores = []
    for key in ('SpliceAI_DS_AG', 'SpliceAI_DS_AL', 'SpliceAI_DS_DG', 'SpliceAI_DS_DL'):
        i = idx.get(key)
        if i is not None and len(csq) > i:
            val = csq[i].strip()
            if val and val != '.':
                try:
                    scores.append(float(val))
                except ValueError:
                    pass
    return max(scores) if scores else None

# ---------------------------------------------------------------------------
# Helper: GnomAD AC (max of AC_all / AC_popmax)
# ---------------------------------------------------------------------------
def get_ac(csq, idx):
    acs = []
    for key in ('GnomAD_AC_all', 'GnomAD_AC_popmax'):
        i = idx.get(key)
        if i is not None and len(csq) > i:
            val = csq[i].strip()
            if val and val != '.':
                try:
                    acs.append(int(float(val)))
                except ValueError:
                    pass
    return max(acs) if acs else None

# ---------------------------------------------------------------------------
# Helper: does this pass the coding/splice filter?
# ---------------------------------------------------------------------------
def passes_consequence(csq_terms, spliceai_score):
    terms = set(t.strip() for t in csq_terms.split('&') if t.strip())

    if terms & ALWAYS_INCLUDE:
        return True

    if 'synonymous_variant' in terms:
        return spliceai_score is not None and spliceai_score >= SPLICEAI_MIN

    if terms & CODING_NONSYN:
        return True

    # Not in our coding/non-synonymous list at all — still allow through if
    # SpliceAI flags it as a likely splice-disrupting variant regardless of
    # the VEP-assigned consequence term.
    if spliceai_score is not None and spliceai_score >= SPLICEAI_MIN:
        return True

    return False

# ---------------------------------------------------------------------------
# HTML report generation
# ---------------------------------------------------------------------------
CONSEQUENCE_COLOURS = {
    'HIGH':     '#ffcccc',
    'MODERATE': '#fff3cc',
    'LOW':      '#e8f4e8',
    'MODIFIER': '#f0f0f0',
}

def badge(text, colour):
    return f'<span style="background:{colour};padding:2px 6px;border-radius:3px;font-size:0.85em;font-weight:bold">{html.escape(str(text))}</span>'

def impact_badge(impact):
    col = CONSEQUENCE_COLOURS.get(impact.upper(), '#f0f0f0')
    return badge(impact, col)

def render_variant_rows(variants):
    """Render <tr> rows for a list of variants, grouped by gene within the section."""
    by_gene = defaultdict(list)
    for v in variants:
        by_gene[v['gene']].append(v)

    rows_html = []
    for gene in sorted(by_gene):
        for v in by_gene[gene]:
            row_class = 'sv-row' if 'SV' in v['type'] else 'snv-row'
            zygosity = v.get('zygosity', ('Het' if '0/1' in v['gt'] or '1/0' in v['gt'] else ('Hom' if '1/1' in v['gt'] else v['gt'])))

            coord = f"{v['chrom']}:{v['pos']}"
            if 'end' in v:
                coord += f"-{v['end']}"
            rows_html.append(f"""
            <tr class="{row_class}">
                <td><strong>{html.escape(gene)}</strong><br>
                    <small style="color:#666">{html.escape(v['phenotype'])}</small></td>
                <td>{html.escape(v['type'])}</td>
                <td style="font-family:monospace;font-size:0.85em">{html.escape(coord)}<br>
                    <small>{html.escape(v['ref'] if 'ref' in v else '')} → {html.escape(v['alt'] if 'alt' in v else v.get('svtype',''))}</small></td>
                <td style="font-size:0.85em">{html.escape(v['consequence'])}</td>
                <td>{impact_badge(v['impact'])}</td>
                <td style="font-family:monospace;font-size:0.85em">
                    {html.escape(v['hgvsc'])}<br>
                    <small>{html.escape(v['hgvsp'])}</small></td>
                <td>{badge(zygosity, '#ddeeff')}<br>
                    <small>GT: {html.escape(v['gt'])} | GQ: {html.escape(v['gq'])} | DP: {html.escape(v['dp'])}</small></td>
                <td>{html.escape(v['ac_gnomad'])}<br>
                    <small>AF: {html.escape(v['af_gnomad'])}</small></td>
                <td>{html.escape(v['spliceai'])}</td>
                <td style="font-size:0.85em">{html.escape(v['clinvar'])}<br>
                    <small>{html.escape(v['clinvar_dn'])}</small></td>
                <td style="font-size:0.85em">
                    LoF: {html.escape(v['lof'])}<br>
                    REVEL: {html.escape(v['revel'])}<br>
                    AM: {html.escape(v['am_class'])}</td>
            </tr>""")
    return ''.join(rows_html)


def render_section(label, variants):
    if not variants:
        return ''
    n_genes = len(set(v['gene'] for v in variants))
    return f"""
<h2 style="color:#2c3e50;border-bottom:2px solid #3498db;padding-bottom:4px;margin-top:32px">
  {html.escape(label)} <span style="font-size:0.6em;color:#888;font-weight:normal">
  ({len(variants)} variant{'s' if len(variants) != 1 else ''} in {n_genes} gene{'s' if n_genes != 1 else ''})</span>
</h2>
<table>
<thead>
<tr>
  <th>Gene / Phenotype</th>
  <th>Type</th>
  <th>Coordinates</th>
  <th>Consequence</th>
  <th>Impact</th>
  <th>HGVS</th>
  <th>Zygosity / QC</th>
  <th>GnomAD4 AC / AF</th>
  <th>SpliceAI</th>
  <th>ClinVar</th>
  <th>Pathogenicity scores</th>
</tr>
</thead>
<tbody>
{render_variant_rows(variants)}
</tbody>
</table>"""


def render_html(variants, sample_id):
    if not variants:
        return None   # caller handles "no output" case

    n_snv = sum(1 for v in variants if 'SV' not in v['type'])
    n_sv  = sum(1 for v in variants if 'SV' in v['type'])
    n_genes_total = len(set(v['gene'] for v in variants))

    # Split into sections by inheritance mode. A gene classified as 'both'
    # in the panel will have produced separate dominant- and
    # recessive-tagged variant entries upstream, so it appears in
    # both sections here without any special-casing.
    sections_html = []
    for mode, label in INHERITANCE_SECTIONS:
        mode_variants = [v for v in variants if v['inheritance'] == mode]
        sections_html.append(render_section(label, mode_variants))

    return f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Gene Panel Report — {html.escape(sample_id)}</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 20px; background: #fafafa; color: #222; }}
  h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 8px; }}
  .summary {{ display:flex; gap:16px; margin: 16px 0; flex-wrap:wrap; }}
  .card {{ background:white; border:1px solid #ddd; border-radius:6px; padding:12px 20px;
           box-shadow: 0 1px 3px rgba(0,0,0,0.1); min-width:140px; text-align:center; }}
  .card .num {{ font-size:2em; font-weight:bold; color:#2c3e50; }}
  .card .lbl {{ font-size:0.85em; color:#666; }}
  table {{ border-collapse:collapse; width:100%; background:white; border-radius:6px;
           box-shadow: 0 1px 3px rgba(0,0,0,0.1); overflow:hidden; margin-top:16px; }}
  th {{ background:#2c3e50; color:white; padding:10px 8px; text-align:left; font-size:0.9em; }}
  td {{ padding:8px; border-bottom:1px solid #eee; vertical-align:top; font-size:0.9em; }}
  .snv-row:hover {{ background:#f0f7ff; }}
  .sv-row  {{ background:#fffbf0; }}
  .sv-row:hover {{ background:#fff3cc; }}
  .filter-info {{ background:#e8f4fd; border:1px solid #bee3f8; border-radius:6px;
                  padding:12px 16px; margin:12px 0; font-size:0.9em; }}
  .no-results {{ text-align:center; padding:40px; color:#888; font-size:1.2em; }}
  th, td {{ white-space: nowrap; }}
  td:nth-child(4), td:nth-child(6), td:nth-child(10) {{ white-space: normal; max-width:200px; }}
</style>
</head>
<body>
<h1>Gene Panel Variant Report</h1>
<p><strong>Sample:</strong> {html.escape(sample_id)} &nbsp;|&nbsp;
   <strong>Generated:</strong> {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M')}</p>

<div class="filter-info">
  <strong>Filtering criteria applied:</strong><br>
  <b>Dominant:</b> GnomAD4 AC &lt; {DOM_AC_MAX} (bypassed for ClinVar Pathogenic/Likely_pathogenic calls) | coding/non-synonymous or SpliceAI ≥ {SPLICEAI_MIN} (±50/+10 bp) | heterozygous or homozygous<br>
  <b>Recessive:</b> GnomAD4 AC &lt; {REC_AC_MAX} (bypassed for ClinVar Pathogenic/Likely_pathogenic calls) | coding/non-synonymous or SpliceAI ≥ {SPLICEAI_MIN} (±50/+10 bp) | homozygous OR compound het (≥2 het variants in same gene)<br>
  <b>Genes with both inheritance patterns</b> are evaluated independently against each threshold and may appear in both sections below.
</div>
<div class="summary">
  <div class="card"><div class="num">{len(variants)}</div><div class="lbl">Total variant entries</div></div>
  <div class="card"><div class="num">{n_genes_total}</div><div class="lbl">Genes affected</div></div>
  <div class="card"><div class="num">{n_snv}</div><div class="lbl">SNV/Indel</div></div>
  <div class="card"><div class="num">{n_sv}</div><div class="lbl">Structural</div></div>
</div>

{''.join(sections_html)}

<p style="margin-top:24px;color:#aaa;font-size:0.8em">
  SNV rows shown in white; SV rows shown in yellow. Compound het pairs shown together.
  GnomAD AC = max(AC_all, AC_popmax). Each variant is deduplicated across
  transcript annotations, preferring the MANE Select transcript.
</p>
</body>
</html>"""

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description='Gene panel variant filter and HTML report')
    ap.add_argument('--snv',    required=True, help='VEP-annotated SNV VCF (.vcf.gz)')
    ap.add_argument('--sv',     required=False, default=None, help='VEP-annotated SV VCF (.vcf.gz)')
    ap.add_argument('--panel',  required=True, help='Gene panel TSV (monogenic_gene_panel.tsv)')
    ap.add_argument('--sample', required=True, help='Sample ID for report title')
    ap.add_argument('--out',    required=True, help='Output HTML file path')
    args = ap.parse_args()

    print(f'Loading gene panel from {args.panel}...')
    panel = load_panel(args.panel)
    print(f'  {len(panel)} genes loaded.')

    all_variants = []

    print(f'Processing SNV VCF: {args.snv}')
    snv_vars = process_snv_vcf(args.snv, panel, args.sample)
    print(f'  {len(snv_vars)} SNV/indel variant entries passed filters (post-dedup).')
    all_variants.extend(snv_vars)

    if args.sv:
        print(f'Processing SV VCF: {args.sv}')
        sv_vars = process_sv_vcf(args.sv, panel, args.sample)
        print(f'  {len(sv_vars)} SV variant entries passed filters (post-dedup).')
        all_variants.extend(sv_vars)

    if not all_variants:
        print(f'No variants passed filters for {args.sample} — no HTML generated.')
        sys.exit(0)

    print(f'Generating HTML report: {args.out}')
    report_html = render_html(all_variants, args.sample)
    with open(args.out, 'w') as f:
        f.write(report_html)
    print(f'Done. {len(all_variants)} variant entries in report.')

if __name__ == '__main__':
    main()

