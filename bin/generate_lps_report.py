#!/usr/bin/env python3
"""Generate a single self-contained HTML report for the LPS typing pipeline.

Reads the aggregated TSVs produced in ``10_report/`` plus the LPS database
assets (GenBank files, phenotype lookup and images, optional gene colours) and
renders one self-contained ``LPS_typing_report.html`` (all figures and images
embedded; no external assets, no network access required to view).

The script is intentionally tolerant of missing inputs: any optional step that
was skipped simply omits its section. The only hard requirement is the subtype
report (``10_*_subtype_report.tsv``).
"""

import argparse
import base64
import io
import json
import re
import sys
import xml.etree.ElementTree as ET
from collections import OrderedDict
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgb  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402
import pandas as pd  # noqa: E402
from jinja2 import Template  # noqa: E402

# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #

# Known aggregated filenames in 10_report/ (the "Illumina" infix is fixed in the
# pipeline). Each is optional except the subtype report.
F_SUBTYPE = "10_Illumina_subtype_report.tsv"
F_PHENOTYPE = "10_Illumina_phenotype_report.tsv"  # one row per sample (lps_phenotype.py)
F_KAPTIVE = "7_Illumina_kaptive_results.tsv"
F_QUAST = "4_Illumina_quast_report.tsv"
F_SHOVILL = "3_Illumina_shovill_stats.tsv"
F_CHECKM = "5_Illumina_checkm_lineage_wf_results.tsv"
F_SYLPH = "6_Illumina_sylph_summary.tsv"
F_AMR = "12_Illumina_amrfinder.tsv"
F_MULTIQC = "2_Illumina_multiqc_report.html"

# QC thresholds for warning flags.
COMPLETENESS_MIN = 90.0
CONTAMINATION_MAX = 5.0
EXPECTED_SPECIES = "pasteurella multocida"

NA_VALUES = {"", "NA", "na", "N/A", ".", "-"}

# Deterministic fallback gene-track palette (used when gene_colors.tsv absent or
# a gene is missing from it). Ordered to roughly echo the collaborator's figures.
FALLBACK_PALETTE = [
    "#D85A30", "#378ADD", "#7F77DD", "#639922", "#BA7517",
    "#E24B4A", "#1D9E75", "#D4537E", "#888780", "#4682B4",
    "#8B4513", "#DA70D6",
]


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #

def is_na(value):
    return value is None or (isinstance(value, float) and pd.isna(value)) or \
        (str(value).strip() in NA_VALUES)


def clean(value):
    """Return a stripped string, or '' for NA-like values."""
    return "" if is_na(value) else str(value).strip()


def read_tsv(path, **kwargs):
    """Read a TSV if it exists and is non-empty, else return None."""
    if path is None or not path.exists() or path.stat().st_size == 0:
        return None
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, **kwargs)
    except Exception as exc:  # pragma: no cover - defensive
        sys.stderr.write(f"warning: could not read {path}: {exc}\n")
        return None
    return df if not df.empty else None


def fig_to_data_uri(fig):
    """Render a matplotlib figure to a base64 SVG data URI and close it."""
    buf = io.StringIO()
    fig.savefig(buf, format="svg", bbox_inches="tight")
    plt.close(fig)
    encoded = base64.b64encode(buf.getvalue().encode("utf-8")).decode("ascii")
    return f"data:image/svg+xml;base64,{encoded}"


_SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", _SVG_NS)
ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")


def _inject_titles(svg_str, titles):
    """Add a <title> child to every <g id=...> whose id is in `titles` so the
    browser shows a native tooltip on hover (requires the SVG to be inlined)."""
    body = svg_str[svg_str.index("<svg"):]
    root = ET.fromstring(body)
    for g in root.iter(f"{{{_SVG_NS}}}g"):
        gid = g.get("id")
        if gid in titles:
            t = ET.Element(f"{{{_SVG_NS}}}title")
            t.text = titles[gid]
            g.insert(0, t)
    return ET.tostring(root, encoding="unicode")


def _namespace_ids(svg_str, prefix):
    """Prefix all internal ids + their references so multiple inlined SVGs in
    one HTML document don't collide (clip-paths, glyph defs, gids)."""
    ids = sorted(set(re.findall(r'id="([^"]+)"', svg_str)), key=len, reverse=True)
    for i in ids:
        svg_str = svg_str.replace(f'id="{i}"', f'id="{prefix}{i}"')
        svg_str = svg_str.replace(f'url(#{i})', f'url(#{prefix}{i})')
        svg_str = svg_str.replace(f'href="#{i}"', f'href="#{prefix}{i}"')
    return svg_str


def fig_to_inline_svg(fig, titles, prefix):
    """Render a figure to inline SVG markup (not an <img> data URI) with native
    <title> tooltips attached and all ids namespaced."""
    buf = io.StringIO()
    fig.savefig(buf, format="svg", bbox_inches="tight")
    plt.close(fig)
    return _namespace_ids(_inject_titles(buf.getvalue(), titles), prefix)


def short_vartype(vartype):
    """Shorten a VARTYPE label for plotting, e.g. 'snp S92*' -> 'S92*',
    '19bp deletion' -> '19bp del'. The database is inconsistent about the space
    after 'snp' (e.g. 'snpY226*'), so both forms are stripped."""
    v = clean(vartype)
    if not v:
        return ""
    if v.lower().startswith("snp"):
        stripped = v[3:].strip()
        if stripped:
            return stripped
    v = re.sub(r"\binsertion\b", "ins", v, flags=re.I)
    v = re.sub(r"\bdeletion\b", "del", v, flags=re.I)
    return v


# --------------------------------------------------------------------------- #
# GenBank gene-feature parsing (lightweight, no Biopython)
# --------------------------------------------------------------------------- #

_LOCUS_LEN_RE = re.compile(r"^LOCUS\s+\S+\s+(\d+)\s+bp", re.I)
_GENE_LOC_RE = re.compile(
    r"^\s{1,8}gene\s+(complement\()?\s*<?(\d+)\.\.>?(\d+)\)?\s*$"
)
_GENE_NAME_RE = re.compile(r'/gene="([^"]+)"')


def parse_genbank_genes(gb_path):
    """Return (locus_length, [{'name','start','end','strand'}...]) from a .gb.

    Only top-level ``gene`` features are used. ``start``/``end`` are 1-based
    inclusive genomic coordinates; ``strand`` is '+' or '-'.
    """
    locus_len = 0
    genes = []
    pending = None  # (start, end, strand) awaiting its /gene= name
    with open(gb_path) as fh:
        for line in fh:
            if not locus_len:
                m = _LOCUS_LEN_RE.match(line)
                if m:
                    locus_len = int(m.group(1))
            m = _GENE_LOC_RE.match(line)
            if m:
                strand = "-" if m.group(1) else "+"
                pending = (int(m.group(2)), int(m.group(3)), strand)
                continue
            if pending is not None:
                nm = _GENE_NAME_RE.search(line)
                if nm:
                    start, end, strand = pending
                    genes.append({
                        "name": nm.group(1).strip(),
                        "start": start, "end": end, "strand": strand,
                    })
                    pending = None
    # De-duplicate genes that appear twice (gene + CDS both carry /gene=),
    # keeping first occurrence by (name,start,end).
    seen = set()
    uniq = []
    for g in genes:
        key = (g["name"], g["start"], g["end"])
        if key not in seen:
            seen.add(key)
            uniq.append(g)
    if not locus_len and uniq:
        locus_len = max(g["end"] for g in uniq)
    return locus_len, uniq


def load_reference_map(lps_db_dir):
    """Parse reference_LPS.txt -> {type: gb_filename}."""
    path = lps_db_dir / "reference_LPS.txt"
    mapping = {}
    if not path.exists():
        return mapping
    for line in path.read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2 and parts[0].strip():
            mapping[parts[0].strip()] = parts[1].strip()
    return mapping


def load_gene_colors(lps_db_dir):
    """Optional gene_colors.tsv -> {gene: hex}."""
    path = lps_db_dir / "gene_colors.tsv"
    colors = {}
    if not path.exists():
        return colors
    for i, line in enumerate(path.read_text().splitlines()):
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2 and parts[0].strip() and not (
                i == 0 and parts[1].strip().lower() in {"hex", "color", "colour"}):
            colors[parts[0].strip()] = parts[1].strip()
    return colors


def color_for_genes(genes, gene_colors):
    """Assign a colour to every gene, using gene_colors.tsv where available and
    the fallback palette otherwise (deterministic by gene order)."""
    out = {}
    pi = 0
    for g in genes:
        name = g["name"]
        if name in gene_colors:
            out[name] = gene_colors[name]
        else:
            out[name] = FALLBACK_PALETTE[pi % len(FALLBACK_PALETTE)]
            pi += 1
    return out


# --------------------------------------------------------------------------- #
# Phenotype images
# --------------------------------------------------------------------------- #

_IMG_MIME = {".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg",
             ".gif": "image/gif", ".svg": "image/svg+xml", ".webp": "image/webp"}


def load_phenotype_images(lps_db_dir):
    """Map phenotype id (file stem, e.g. 'L3_P4') -> base64 data URI."""
    images = {}
    img_dir = lps_db_dir / "phenotype_images"
    if not img_dir.is_dir():
        return images
    for path in sorted(img_dir.iterdir()):
        if path.stem == "template":
            continue
        mime = _IMG_MIME.get(path.suffix.lower())
        if not mime:
            continue
        data = base64.b64encode(path.read_bytes()).decode("ascii")
        images[path.stem] = f"data:{mime};base64,{data}"
    return images


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #

def _gene_arrow_xy(start, end, y, half_h, strand, head):
    """Five-point directional gene arrow, pointing left on the minus strand."""
    hl = min(head, 0.5 * (end - start))
    if strand == "-":
        return [(end, y - half_h), (start + hl, y - half_h), (start, y),
                (start + hl, y + half_h), (end, y + half_h)]
    return [(start, y - half_h), (end - hl, y - half_h), (end, y),
            (end - hl, y + half_h), (start, y + half_h)]


def _contrast_text(hex_color):
    """Near-black or white, whichever reads better on `hex_color`."""
    lin = [c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
           for c in to_rgb(hex_color)]
    lum = 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]
    return "#1A1A1A" if lum > 0.45 else "white"


def make_lollipop(lps_type, chrom, genes, locus_len, variants, gcolors):
    """Build a lollipop figure for one LPS type and return inline SVG markup.

    ``variants`` is a list of dicts: {'pos','label','gene','count'}; height = number
    of genomes carrying the mutation. Mutation labels are placed in leader-line
    stacked tiers to avoid overlap; the gene track is drawn as directional arrows
    labelled only where the text fits (the rest are covered by the gene-colour
    legend). Markers and gene arrows carry <title> tooltips. Width scales with
    locus length and variant count.
    """
    titles = {}
    max_count = max((v["count"] for v in variants), default=1)
    xmax = max(locus_len, max(g["end"] for g in genes) if genes else locus_len, 1)

    fig_w = min(16.0, max(9.0, xmax / 900.0 + len(variants) * 0.18))
    fig, ax = plt.subplots(figsize=(fig_w, 4.0))
    # panel width, not figure width: the y axis label and ticks take ~1.1 in,
    # and both the label-overlap test and the gene-name fit test measure in
    # panel pixels
    fig_px = max(2.0, fig_w - 1.1) * fig.dpi

    band = max_count * 0.18       # gene-arrow height, in count units
    half_h = band / 2.0
    track_y = -half_h - max_count * 0.05   # arrow centre line, below the baseline

    # Gene track: directional arrows on a thin backbone. Strand comes straight
    # from the complement() wrapper in the .gb, so the arrows cannot disagree
    # with the reference. A name is drawn inside its arrow only where the
    # straight body has room; the rest are covered by the gene-colour legend.
    ax.plot([1, locus_len], [track_y, track_y], color="#888780",
            linewidth=0.8, zorder=1)
    for g in genes:
        gid = f"gene_{g['name']}"
        head = min(0.35 * (g["end"] - g["start"]), 0.02 * xmax)
        ax.add_patch(Polygon(
            _gene_arrow_xy(g["start"], g["end"], track_y, half_h, g["strand"], head),
            closed=True, facecolor=gcolors[g["name"]], edgecolor="#333333",
            linewidth=0.4, zorder=2, gid=gid))
        titles[gid] = f"{g['name']} · {g['start']}–{g['end']} ({g['strand']})"
        body_px = (g["end"] - g["start"] - head) / xmax * fig_px
        if body_px >= len(g["name"]) * 8.5:
            cx = ((g["start"] + head + g["end"]) / 2.0 if g["strand"] == "-"
                  else (g["start"] + g["end"] - head) / 2.0)
            ax.text(cx, track_y, g["name"], ha="center", va="center", fontsize=8,
                    color=_contrast_text(gcolors[g["name"]]), zorder=3)

    ax.plot([0, xmax], [0, 0], color="#888780", linewidth=0.8, zorder=1)

    # Lollipops (each marker its own gid -> hover tooltip).
    for i, v in enumerate(sorted(variants, key=lambda d: d["pos"])):
        col = gcolors.get(v["gene"], "#5F5E5A")
        ax.plot([v["pos"], v["pos"]], [0, v["count"]], color="#C9C7BE",
                linewidth=1.1, zorder=2)
        gid = f"m{i}"
        ax.plot([v["pos"]], [v["count"]], marker="o", markersize=9, color=col,
                markeredgecolor="white", markeredgewidth=0.8, zorder=4, gid=gid)
        titles[gid] = (f"{v['label']} · {v['gene'] or 'intergenic'} · "
                       f"pos {v['pos']} · {v['count']} genome(s)")

    # Leader-line stacked labels: greedy tier assignment by estimated label width.
    label_base = max_count * 1.18
    dy = max_count * 0.16
    px_per_data = fig_px / xmax
    tiers_last = []  # rightmost placed-label edge (data units) per tier
    for v in sorted(variants, key=lambda d: d["pos"]):
        if not v["label"]:
            continue
        half = (len(v["label"]) * 6.6) / px_per_data / 2 + xmax * 0.004
        tier = 0
        while tier < len(tiers_last) and v["pos"] - half < tiers_last[tier]:
            tier += 1
        if tier == len(tiers_last):
            tiers_last.append(v["pos"] + half)
        else:
            tiers_last[tier] = v["pos"] + half
        y_lab = label_base + tier * dy
        ax.plot([v["pos"], v["pos"]], [v["count"], y_lab - dy * 0.25],
                color="#D8D6CD", linewidth=0.7, zorder=3)
        ax.text(v["pos"], y_lab, v["label"], ha="center", va="bottom",
                fontsize=8, color="#2C2C2A", zorder=5)
    n_tiers = len(tiers_last) if tiers_last else 1

    ax.set_xlim(-xmax * 0.02, xmax * 1.02)
    ax.set_ylim(track_y - half_h - band * 0.35,
                label_base + n_tiers * dy + max_count * 0.18)
    ax.set_yticks(range(0, max_count + 1, max(1, max_count // 5)))
    ax.set_ylabel("genomes with mutation", fontsize=9)
    ax.set_xlabel(f"{chrom}  ·  position (bp)", fontsize=9)
    ax.set_title(f"{lps_type} — observed mutations", fontsize=11, loc="left")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(labelsize=8)
    return fig_to_inline_svg(fig, titles, f"ll{lps_type}_")


def make_bar(counts, title, xlabel):
    """Horizontal bar chart from an ordered dict {label: count}."""
    labels = list(counts.keys())
    values = [counts[k] for k in labels]
    height = max(1.6, 0.42 * len(labels) + 0.8)
    fig, ax = plt.subplots(figsize=(6.2, height))
    ypos = range(len(labels))
    ax.barh(list(ypos), values, color="#378ADD", edgecolor="white")
    ax.set_yticks(list(ypos))
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_title(title, fontsize=11, loc="left")
    for i, val in enumerate(values):
        ax.text(val, i, f" {val}", va="center", fontsize=8.5, color="#2C2C2A")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(labelsize=8)
    ax.margins(x=0.12)
    return fig_to_data_uri(fig)


# --------------------------------------------------------------------------- #
# Data assembly
# --------------------------------------------------------------------------- #

def index_by_sample(df, sample_col):
    """Return {sample: first-row-dict} for a dataframe keyed on a sample column."""
    out = {}
    if df is None or sample_col not in df.columns:
        return out
    for _, row in df.iterrows():
        s = clean(row[sample_col])
        if s and s not in out:
            out[s] = row.to_dict()
    return out


def parse_sylph(df):
    """Sylph taxonomy summary -> {sample: {'species': <top species>}}.

    Keyed on ``sample`` (not ``sampleID``); the genome's species for the QC
    table and the off-target flag is ``top_species_by_taxonomic_abundance``.
    """
    out = {}
    if df is None or "sample" not in df.columns:
        return out
    for _, row in df.iterrows():
        s = clean(row["sample"])
        if s and s not in out:
            out[s] = {"species": clean(row.get("top_species_by_taxonomic_abundance"))}
    return out


def parse_quast(df):
    """QUAST matrix: metric rows x sample columns (index col = 'Assembly')."""
    out = {}
    if df is None:
        return out
    metric_col = df.columns[0]
    sample_cols = [c for c in df.columns[1:]]
    metrics = {clean(r[metric_col]): r for _, r in df.iterrows()}
    for s in sample_cols:
        rec = {}
        for label, key in (("Total length", "total_length"),
                           ("# contigs", "n_contigs"),
                           ("N50", "n50"),
                           ("GC (%)", "gc")):
            if label in metrics:
                rec[key] = clean(metrics[label][s])
        out[clean(s)] = rec
    return out


PETG_CLASS = {"functional": "ok", "nonfunctional": "bad", "absent": "bad", "uncertain": "flag"}


def basis_categories(basis):
    """Short labels for the PHENOTYPE_BASIS text, used for badges and filtering."""
    basis = clean(basis)
    if not basis:
        return []
    labels = []
    for text, label in (("uncertain", "uncertain"), ("uncurated", "uncurated"),
                        ("low-confidence call", "low-confidence"), ("mixed call", "mixed call"),
                        ("untypeable", "untypeable"), ("no phenotype rules", "no rules")):
        if text in basis:
            labels.append(label)
    if not labels or basis.startswith("curated"):
        labels.insert(0, "curated")
    return labels


def build_records(report_df, qc, pheno_df=None):
    """Collapse the multi-row subtype report into one record per genome.

    Phenotype, petG status and the phenotype basis come from the per-sample
    phenotype report when it is available.
    """
    pheno = index_by_sample(pheno_df, "SAMPLE") if pheno_df is not None else {}
    records = OrderedDict()
    samples = {clean(v) for v in report_df["SAMPLE"] if clean(v)} | set(pheno)
    for s in sorted(samples):
        rows = report_df[report_df["SAMPLE"] == s]
        types = sorted({clean(r["TYPE"]) for _, r in rows.iterrows() if clean(r["TYPE"])})
        subtypes = sorted({clean(r["SUBTYPE"]) for _, r in rows.iterrows() if clean(r["SUBTYPE"])})
        phenos = OrderedDict()
        variants = []
        petg = ""
        mlst = ""
        for _, r in rows.iterrows():
            ph = clean(r.get("PREDICTED_PHENOTYPE"))
            if ph:
                phenos[ph] = clean(r.get("PREDICTED_PHENOTYPE_DESCRIPTION"))
            if clean(r.get("PETG_PRESENT")).lower() == "yes":
                petg = "yes"
            if not mlst:
                mlst = clean(r.get("MLST"))
            pos = clean(r.get("POS"))
            if pos and pos not in NA_VALUES:
                variants.append({
                    "gene": clean(r.get("GENE")),
                    "pos": pos,
                    "ref": clean(r.get("REF")),
                    "alt": clean(r.get("ALT")),
                    "vartype": clean(r.get("VARTYPE")),
                    "subtype": clean(r.get("SUBTYPE")),
                    "note": clean(r.get("NOTE")),
                })

        lps_type = types[0] if types else ""
        ph_row = pheno.get(s)
        basis = ""
        gene_states = ""
        petg_detail = ""
        if ph_row is not None:
            phenos = OrderedDict()
            names = [p.strip() for p in clean(ph_row.get("PREDICTED_PHENOTYPE")).split("|") if p.strip()]
            descs = [d.strip() for d in clean(ph_row.get("PREDICTED_PHENOTYPE_DESCRIPTION")).split("|")]
            for i, name in enumerate(names):
                phenos[name] = descs[i] if i < len(descs) else ""
            petg = clean(ph_row.get("PETG_STATUS"))
            petg_detail = clean(ph_row.get("PETG_DETAIL"))
            basis = clean(ph_row.get("PHENOTYPE_BASIS"))
            gene_states = clean(ph_row.get("LPS_GENE_STATES"))
            lps_type = clean(ph_row.get("TYPE")) or lps_type
            mlst = mlst or clean(ph_row.get("MLST"))
            subtypes = [t for t in clean(ph_row.get("SUBTYPES")).split(";") if t] or subtypes
        # QC joins
        ck = qc["checkm"].get(s, {})
        completeness = clean(ck.get("Completeness"))
        contamination = clean(ck.get("Contamination"))
        tx = qc["taxonomy"].get(s, {})
        species = clean(tx.get("species"))
        kp = qc["kaptive"].get(s, {})
        qd = qc["quast"].get(s, {})
        sh = qc["shovill"].get(s, {})

        flags = []
        if lps_type == "untypeable" or not lps_type:
            flags.append("untypeable")
        if species and EXPECTED_SPECIES not in species.lower():
            flags.append("off-target species")
        try:
            if completeness and float(completeness) < COMPLETENESS_MIN:
                flags.append("low completeness")
        except ValueError:
            pass
        try:
            if contamination and float(contamination) > CONTAMINATION_MAX:
                flags.append("high contamination")
        except ValueError:
            pass

        records[s] = {
            "sample": s,
            "type": lps_type,
            "subtypes": subtypes,
            "phenotypes": phenos,
            "variants": variants,
            "petg": petg or "no",
            "petg_class": PETG_CLASS.get(petg, ""),
            "petg_detail": petg_detail,
            "basis": basis,
            "basis_labels": basis_categories(basis),
            "gene_states": gene_states,
            "mlst": mlst,
            "species": species,
            "completeness": completeness,
            "contamination": contamination,
            "kaptive_confidence": clean(kp.get("Match confidence")),
            "kaptive_identity": clean(kp.get("Identity")),
            "kaptive_coverage": clean(kp.get("Coverage")),
            "kaptive_problems": clean(kp.get("Problems")),
            "assembly_size": clean(qd.get("total_length")) or clean(sh.get("assembly_size")),
            "n_contigs": clean(qd.get("n_contigs")) or clean(sh.get("nb_contigs")),
            "n50": clean(qd.get("n50")),
            "coverage": clean(sh.get("asssembly_coverage")),
            "amr": [],
            "flags": flags,
        }
    return records


def attach_amr(records, amr_df):
    if amr_df is None or "Name" not in amr_df.columns:
        return
    for _, r in amr_df.iterrows():
        s = clean(r.get("Name"))
        if s in records:
            records[s]["amr"].append({
                "symbol": clean(r.get("Element symbol")),
                "name": clean(r.get("Element name")),
                "cls": clean(r.get("Class")),
                "subclass": clean(r.get("Subclass")),
                "identity": clean(r.get("% Identity to reference")),
            })


def build_lollipops(report_df, lps_db_dir, ref_map, gene_colors):
    """One lollipop per LPS type present in the run (with observed variants)."""
    figs = []
    # observed variants per type: key (type) -> {(chrom,pos,ref,alt): {...,'samples':set}}
    have_pos = report_df[report_df["POS"].apply(lambda x: clean(x) not in NA_VALUES)]
    for lps_type in sorted({clean(t) for t in have_pos["TYPE"] if clean(t) and clean(t) != "untypeable"}):
        sub = have_pos[have_pos["TYPE"] == lps_type]
        gb_name = ref_map.get(lps_type)
        if not gb_name:
            continue
        gb_path = lps_db_dir / gb_name
        if not gb_path.exists():
            continue
        locus_len, genes = parse_genbank_genes(gb_path)
        if not genes:
            continue
        chrom = clean(sub.iloc[0]["CHROM"]) or gb_path.stem
        agg = {}
        for _, r in sub.iterrows():
            try:
                pos = int(clean(r["POS"]))
            except ValueError:
                continue
            key = (clean(r["CHROM"]), pos, clean(r["REF"]), clean(r["ALT"]))
            rec = agg.setdefault(key, {
                "pos": pos, "gene": clean(r["GENE"]),
                "label": short_vartype(r["VARTYPE"]), "samples": set()})
            rec["samples"].add(clean(r["SAMPLE"]))
        gcolors = color_for_genes(genes, gene_colors)
        variants = [{"pos": v["pos"], "gene": v["gene"], "label": v["label"],
                     "count": len(v["samples"]),
                     "color": gcolors.get(v["gene"], "#5F5E5A")} for v in agg.values()]
        if not variants:
            continue
        legend = [{"name": g["name"], "color": gcolors[g["name"]],
                   "start": g["start"], "end": g["end"]} for g in genes]
        figs.append({
            "type": lps_type,
            "n_variants": len(variants),
            "svg": make_lollipop(lps_type, chrom, genes, locus_len, variants, gcolors),
            "legend": legend,
            "variants": sorted(variants, key=lambda d: d["pos"]),
        })
    return figs


# --------------------------------------------------------------------------- #
# HTML template
# --------------------------------------------------------------------------- #

TEMPLATE = r"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>LPS typing report</title>
<style>
  :root { --fg:#1f2328; --muted:#656d76; --line:#d8dee4; --bg:#fff;
          --accent:#0969da; --warn-bg:#fff8c5; --warn-fg:#7d4e00; --danger:#cf222e; }
  * { box-sizing: border-box; }
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
         color: var(--fg); background: var(--bg); margin: 0; line-height: 1.5; }
  .wrap { max-width: 1100px; margin: 0 auto; padding: 32px 24px 80px; }
  h1 { font-size: 26px; margin: 0 0 4px; }
  h2 { font-size: 20px; margin: 40px 0 12px; padding-bottom: 6px; border-bottom: 1px solid var(--line); }
  h3 { font-size: 15px; margin: 18px 0 8px; }
  .sub { color: var(--muted); font-size: 14px; }
  .cards { display: flex; flex-wrap: wrap; gap: 12px; margin: 16px 0; }
  .metric { background: #f6f8fa; border-radius: 8px; padding: 12px 16px; min-width: 130px; }
  .metric .v { font-size: 24px; font-weight: 600; }
  .metric .l { font-size: 12px; color: var(--muted); text-transform: uppercase; letter-spacing: .04em; }
  table { border-collapse: collapse; width: 100%; font-size: 13.5px; }
  th, td { text-align: left; padding: 7px 10px; border-bottom: 1px solid var(--line); vertical-align: top; }
  th { background: #f6f8fa; user-select: none; white-space: nowrap; }
  .data-table th { cursor: pointer; position: sticky; top: 0; z-index: 2; }
  .data-table th[data-sort="none"] { cursor: default; }
  .data-table th .arrow { color: var(--muted); font-size: 11px; margin-left: 3px; }
  .data-table tr.filter-row th { position: sticky; top: 33px; z-index: 1; padding: 4px 8px; cursor: auto; }
  .data-table tr.filter-row input, .data-table tr.filter-row select {
    width: 100%; font-size: 12px; padding: 3px 5px; border: 1px solid var(--line);
    border-radius: 4px; background: #fff; color: var(--fg); }
  .data-table tr.filter-row .range { display: flex; gap: 4px; }
  .data-table tr.filter-row .range input { width: 50%; }
  tr:hover td { background: #fafbfc; }
  .pill { display: inline-block; padding: 1px 8px; border-radius: 999px; font-size: 12px; background: #eaeef2; margin: 1px 2px; white-space: nowrap; }
  .nowrap { white-space: nowrap; }
  .flag { background: var(--warn-bg); color: var(--warn-fg); }
  .flag.untypeable { background: #ffd8d3; color: var(--danger); }
  .ok { color: #1a7f37; }
  .bad { color: var(--danger); }
  .pill.basis-curated { background: #dafbe1; color: #1a7f37; }
  .pill.basis-uncertain, .pill.basis-mixed, .pill.basis-low-confidence { background: var(--warn-bg); color: var(--warn-fg); }
  .pill.basis-uncurated { background: #ddf4ff; color: #0550ae; }
  .tablebar { display: flex; flex-wrap: wrap; align-items: center; gap: 12px; margin: 8px 0 12px; }
  .tablebar input[type=text] { font-size: 14px; padding: 7px 10px; width: 280px; border: 1px solid var(--line); border-radius: 6px; }
  .count { color: var(--muted); font-size: 13px; }
  .toggle { font-size: 13px; color: var(--fg); display: inline-flex; align-items: center; gap: 5px; cursor: pointer; }
  .btn { font-size: 13px; padding: 6px 11px; border: 1px solid var(--line); border-radius: 6px; background: #f6f8fa; cursor: pointer; color: var(--fg); }
  .btn:hover { background: #eef1f4; }
  img.fig { max-width: 100%; height: auto; }
  .gallery { display: flex; flex-wrap: wrap; gap: 16px; }
  .pheno { border: 1px solid var(--line); border-radius: 8px; padding: 12px; width: 220px; }
  .pheno img { width: 100%; height: 140px; object-fit: contain; background: #f6f8fa; border-radius: 6px; }
  .pheno .noimg { display: flex; align-items: center; justify-content: center; height: 140px;
                  background: #f6f8fa; border-radius: 6px; color: var(--muted); font-size: 13px; }
  .pheno .cap { font-size: 13px; margin-top: 8px; }
  .pheno .cap b { display: block; }
  details { border: 1px solid var(--line); border-radius: 8px; margin: 8px 0; padding: 0 14px; }
  details > summary { cursor: pointer; padding: 12px 0; font-weight: 600; list-style: none; }
  details > summary::-webkit-details-marker { display: none; }
  details > summary:before { content: "\25B8 "; color: var(--muted); }
  details[open] > summary:before { content: "\25BE "; }
  .kv { display: grid; grid-template-columns: max-content 1fr; gap: 2px 16px; font-size: 13.5px; margin: 6px 0 12px; }
  .kv .k { color: var(--muted); }
  .muted { color: var(--muted); }
  .figblock { margin: 18px 0 26px; }
  .lolli-svg { overflow-x: auto; }
  .lolli-svg svg { max-width: 100%; height: auto; }
  .lolli-svg g[id*="_m"], .lolli-svg g[id*="_gene_"] { cursor: help; }
  .legend { margin: 4px 0 12px; }
  .leg { display: inline-block; font-size: 13px; margin: 2px 14px 2px 0; white-space: nowrap; }
  .sw { display: inline-block; width: 12px; height: 12px; border-radius: 2px; vertical-align: -1px; margin-right: 5px; }
  .coord { color: var(--muted); font-size: 11px; }
  .vtbl { width: auto; min-width: 360px; max-width: 560px; font-size: 13px; margin: 4px 0 10px; }
  .vtbl th { cursor: default; }
  .vtbl th, .vtbl td { padding: 5px 10px; }
  footer { margin-top: 48px; padding-top: 16px; border-top: 1px solid var(--line); color: var(--muted); font-size: 13px; }
  a { color: var(--accent); }
</style></head>
<body><div class="wrap">

<h1>LPS typing report</h1>
<div class="sub">{{ meta.pipeline_name }}{% if meta.pipeline_version %} v{{ meta.pipeline_version }}{% endif %} · generated {{ meta.date }}</div>

<h2>Run overview</h2>
<div class="cards">
  <div class="metric"><div class="v">{{ records|length }}</div><div class="l">genomes</div></div>
  <div class="metric"><div class="v">{{ type_counts|length }}</div><div class="l">LPS types</div></div>
  <div class="metric"><div class="v">{{ pheno_counts|reject("equalto", "uncertain")|list|length }}</div><div class="l">phenotypes</div></div>
  <div class="metric"><div class="v">{{ n_flagged }}</div><div class="l">flagged genomes</div></div>
</div>
{% if meta.steps %}<p class="sub"><b>Steps run:</b> {{ meta.steps }}</p>{% endif %}
{% if meta.skipped %}<p class="sub"><b>Skipped:</b> {{ meta.skipped }}</p>{% endif %}
{% if meta.params %}<p class="sub">{{ meta.params }}</p>{% endif %}

<h2>Genome summary</h2>
<div class="tablebar" data-for="summary">
  <input class="tsearch" type="text" placeholder="Search genomes…">
  <label class="toggle"><input type="checkbox" class="flagged-only"> Flagged only</label>
  <span class="count"></span>
</div>
<table id="summary" class="data-table">
  <thead><tr>
    <th data-filter="text">Sample</th><th data-filter="select">Species</th>
    <th data-filter="select">Type</th><th data-filter="multi">Subtype(s)</th><th data-filter="multi">Phenotype</th>
    <th data-filter="multi">Basis</th>
    <th data-filter="select">MLST</th><th data-filter="select">petG</th>
    <th data-type="num" data-filter="range">Compl.</th><th data-type="num" data-filter="range">Contam.</th><th data-filter="multi">Flags</th>
  </tr></thead>
  <tbody>
  {% for r in records.values() %}
    <tr>
      <td><a href="#g-{{ r.sample }}">{{ r.sample }}</a></td>
      <td>{% if r.species %}<i>{{ r.species }}</i>{% else %}<span class="muted">—</span>{% endif %}</td>
      <td>{{ r.type or "—" }}</td>
      <td data-sort="{{ r.subtypes|join(', ') }}">{% for s in r.subtypes %}<span class="pill">{{ s }}</span>{% else %}<span class="muted">—</span>{% endfor %}</td>
      <td data-sort="{{ r.phenotypes.keys()|join(', ') }}">{% for p, d in r.phenotypes.items() %}<span class="pill" title="{{ d }}">{{ p }}</span>{% else %}<span class="muted">—</span>{% endfor %}</td>
      <td data-sort="{{ r.basis_labels|join(', ') }}">{% for b in r.basis_labels %}<span class="pill basis-{{ b.split(' ')[0] }}" title="{{ r.basis }}">{{ b }}</span>{% else %}<span class="muted">—</span>{% endfor %}</td>
      <td>{{ r.mlst or "—" }}</td>
      <td data-sort="{{ r.petg }}"><span class="nowrap {{ r.petg_class }}" title="{{ r.petg_detail }}">{{ r.petg }}</span></td>
      <td>{{ r.completeness or "—" }}</td>
      <td>{{ r.contamination or "—" }}</td>
      <td data-sort="{{ r.flags|join(', ') }}">{% for f in r.flags %}<span class="pill flag {{ 'untypeable' if f == 'untypeable' else '' }}">{{ f }}</span>{% else %}<span class="muted">—</span>{% endfor %}</td>
    </tr>
  {% endfor %}
  </tbody>
</table>

{% if distributions %}
<h2>Distributions</h2>
{% for d in distributions %}<div class="figblock"><img class="fig" src="{{ d }}" alt="distribution"></div>{% endfor %}
{% endif %}

{% if pheno_gallery %}
<h2>Phenotypes observed</h2>
<div class="gallery">
{% for p in pheno_gallery %}
  <div class="pheno">
    {% if p.img %}<img src="{{ p.img }}" alt="{{ p.id }}">{% else %}<div class="noimg">no image available</div>{% endif %}
    <div class="cap"><b>{{ p.id }}</b>{{ p.desc }}<div class="muted">{{ p.count }} genome(s)</div></div>
  </div>
{% endfor %}
</div>
{% endif %}

{% if lollipops %}
<h2>Mutations per LPS type</h2>
<p class="sub">Lollipop height = number of genomes in this run carrying each mutation. Only mutations matching the LPS subtype database are shown. Hover a lollipop or gene box for details.</p>
{% for l in lollipops %}
<div class="figblock">
  <div class="lolli-svg">{{ l.svg|safe }}</div>
  <div class="legend">
    {% for g in l.legend %}<span class="leg"><span class="sw" style="background: {{ g.color }}"></span>{{ g.name }} <span class="coord">{{ g.start }}–{{ g.end }}</span></span>{% endfor %}
  </div>
  <details><summary>{{ l.type }} mutation details ({{ l.n_variants }})</summary>
    <table class="vtbl">
      <thead><tr><th>Gene</th><th>Pos</th><th>Change</th><th>Genomes</th></tr></thead>
      <tbody>
      {% for v in l.variants %}<tr>
        <td><span class="sw" style="background: {{ v.color }}"></span>{{ v.gene or "intergenic" }}</td>
        <td>{{ v.pos }}</td><td>{{ v.label }}</td><td>{{ v.count }}</td>
      </tr>{% endfor %}
      </tbody>
    </table>
  </details>
</div>
{% endfor %}
{% endif %}

<h2>Quality control</h2>
{% if multiqc %}<p>Read-level QC: <a href="{{ multiqc }}">{{ multiqc }}</a> (MultiQC).</p>{% endif %}
<div class="tablebar" data-for="qc">
  <input class="tsearch" type="text" placeholder="Search…">
  <span class="count"></span>
</div>
<table id="qc" class="data-table">
  <thead><tr><th data-filter="text">Sample</th><th data-filter="select">Species</th>
  <th data-type="num" data-filter="range">Completeness</th><th data-type="num" data-filter="range">Contamination</th>
  <th data-type="num" data-filter="range">Assembly size</th><th data-type="num" data-filter="range">Contigs</th>
  <th data-type="num" data-filter="range">N50</th><th data-type="num" data-filter="range">Coverage</th>
  <th data-filter="select">Kaptive confidence</th></tr></thead>
  <tbody>
  {% for r in records.values() %}
    <tr>
      <td>{{ r.sample }}</td><td><i>{{ r.species or "—" }}</i></td>
      <td>{{ r.completeness or "—" }}</td><td>{{ r.contamination or "—" }}</td>
      <td>{{ r.assembly_size or "—" }}</td><td>{{ r.n_contigs or "—" }}</td>
      <td>{{ r.n50 or "—" }}</td><td>{{ r.coverage or "—" }}</td>
      <td>{{ r.kaptive_confidence or "—" }}</td>
    </tr>
  {% endfor %}
  </tbody>
</table>

{% if has_amr %}
<h2>Antimicrobial resistance</h2>
<div class="tablebar" data-for="amr">
  <input class="tsearch" type="text" placeholder="Search…">
  <span class="count"></span>
</div>
<table id="amr" class="data-table">
  <thead><tr><th data-filter="select">Sample</th><th data-filter="select">Gene</th><th data-filter="text">Element name</th><th data-filter="select">Class</th><th data-filter="select">Subclass</th><th data-type="num" data-filter="range">% identity</th></tr></thead>
  <tbody>
  {% for r in records.values() %}{% for a in r.amr %}
    <tr><td>{{ r.sample }}</td><td>{{ a.symbol }}</td><td>{{ a.name }}</td>
        <td>{{ a.cls }}</td><td>{{ a.subclass }}</td><td>{{ a.identity }}</td></tr>
  {% endfor %}{% endfor %}
  </tbody>
</table>
{% endif %}

<h2>Per-genome detail</h2>
{% if records|length > 1 %}
<div class="tablebar">
  <input id="detail-q" type="text" placeholder="Filter by sample or type…">
  <button type="button" class="btn" id="expand-all">Expand all</button>
  <button type="button" class="btn" id="collapse-all">Collapse all</button>
  <span class="count" id="detail-count"></span>
</div>
{% endif %}
{% for r in records.values() %}
<details id="g-{{ r.sample }}" {% if records|length == 1 %}open{% endif %}>
  <summary>{{ r.sample }} — {{ r.type or "untypeable" }}{% if r.phenotypes %} · {{ r.phenotypes.keys()|join(", ") }}{% endif %}</summary>
  <div style="display:flex; flex-wrap:wrap; gap:24px;">
    <div style="flex:1; min-width:320px;">
      <div class="kv">
        <span class="k">LPS type</span><span>{{ r.type or "—" }}</span>
        <span class="k">Subtype(s)</span><span>{{ r.subtypes|join(", ") or "—" }}</span>
        <span class="k">Phenotype</span><span>{% for p,d in r.phenotypes.items() %}{{ p }}{% if d %} ({{ d }}){% endif %}{% if not loop.last %}; {% endif %}{% else %}—{% endfor %}</span>
        <span class="k">MLST</span><span>{{ r.mlst or "—" }}</span>
        <span class="k">Phenotype basis</span><span>{{ r.basis or "—" }}</span>
        <span class="k">LPS gene states</span><span>{{ r.gene_states or "all rule genes functional" }}</span>
        <span class="k">petG</span><span><span class="{{ r.petg_class }}">{{ r.petg }}</span>{% if r.petg_detail %} <span class="muted">({{ r.petg_detail }})</span>{% endif %}</span>
        <span class="k">Species</span><span><i>{{ r.species or "—" }}</i></span>
        <span class="k">Completeness / contam.</span><span>{{ r.completeness or "—" }} / {{ r.contamination or "—" }}</span>
        <span class="k">Kaptive</span><span>{{ r.kaptive_confidence or "—" }} (id {{ r.kaptive_identity or "—" }}, cov {{ r.kaptive_coverage or "—" }})</span>
        {% if r.flags %}<span class="k">Flags</span><span>{% for f in r.flags %}<span class="pill flag">{{ f }}</span>{% endfor %}</span>{% endif %}
      </div>
      {% if r.variants %}
      <h3>Matched mutations ({{ r.variants|length }})</h3>
      <table><thead><tr><th>Gene</th><th>Pos</th><th>Ref&rarr;Alt</th><th>Type</th><th>Subtype</th></tr></thead>
      <tbody>{% for v in r.variants %}<tr><td>{{ v.gene or "—" }}</td><td>{{ v.pos }}</td>
        <td>{{ v.ref }}&rarr;{{ v.alt }}</td><td>{{ v.vartype }}</td><td>{{ v.subtype }}</td></tr>{% endfor %}</tbody></table>
      {% else %}<p class="muted">No subtype-database mutations matched.</p>{% endif %}
    </div>
    {% if r.pheno_img %}<div style="width:240px;"><img src="{{ r.pheno_img }}" alt="phenotype" style="width:100%; background:#f6f8fa; border-radius:8px;"></div>{% endif %}
  </div>
</details>
{% endfor %}

<footer>
  <p>Database: {{ meta.db_files }}</p>
  <p>{{ meta.pipeline_name }}{% if meta.pipeline_version %} v{{ meta.pipeline_version }}{% endif %}. This report is self-contained — all figures and images are embedded.</p>
</footer>

</div>
<script>
// Sortable + filterable tables. Each table.data-table is paired with a
// preceding .tablebar[data-for="<table id>"] holding its search box / counter.
function enhanceTable(table) {
  var thead = table.tHead, tbody = table.tBodies[0];
  if (!thead || !tbody) return;
  var ths = Array.prototype.slice.call(thead.rows[0].cells);
  var bar = document.querySelector('.tablebar[data-for="' + table.id + '"]');
  var searchInput = bar ? bar.querySelector('.tsearch') : null;
  var countEl = bar ? bar.querySelector('.count') : null;
  var flaggedOnly = bar ? bar.querySelector('.flagged-only') : null;
  var controls = [];

  function sortVal(tr, i) {
    var td = tr.cells[i];
    if (!td) return '';
    return td.dataset.sort != null ? td.dataset.sort : td.textContent.trim();
  }

  // --- sorting (column index is the real cell index, fixing the old offset bug) ---
  ths.forEach(function (th, i) {
    if (th.dataset.sort === 'none') return;
    var arrow = document.createElement('span');
    arrow.className = 'arrow'; arrow.textContent = '↕';
    th.appendChild(arrow);
    th.addEventListener('click', function () {
      var asc = th.dataset.dir !== 'asc';
      ths.forEach(function (o) {
        if (o !== th) { delete o.dataset.dir; var a = o.querySelector('.arrow'); if (a) a.textContent = '↕'; }
      });
      th.dataset.dir = asc ? 'asc' : 'desc';
      arrow.textContent = asc ? '▲' : '▼';
      var num = th.dataset.type === 'num';
      var rows = Array.prototype.slice.call(tbody.rows);
      rows.sort(function (a, b) {
        var x = sortVal(a, i), y = sortVal(b, i);
        if (num) {
          var nx = parseFloat(x), ny = parseFloat(y);
          if (isNaN(nx) && isNaN(ny)) return 0;
          if (isNaN(nx)) return 1;            // missing values always last
          if (isNaN(ny)) return -1;
          return asc ? nx - ny : ny - nx;
        }
        return asc ? x.localeCompare(y) : y.localeCompare(x);
      });
      rows.forEach(function (r) { tbody.appendChild(r); });
    });
  });

  // --- per-column filter row ---
  var filterRow = document.createElement('tr');
  filterRow.className = 'filter-row';
  ths.forEach(function (th, i) {
    var cell = document.createElement('th');
    var kind = th.dataset.filter || 'none';
    if (kind === 'select' || kind === 'multi') {
      var seen = {};
      Array.prototype.slice.call(tbody.rows).forEach(function (tr) {
        var raw = sortVal(tr, i);
        (kind === 'multi' ? raw.split(',') : [raw]).forEach(function (t) {
          t = t.trim(); if (t && t !== '—') seen[t] = 1;
        });
      });
      var sel = document.createElement('select');
      sel.innerHTML = '<option value="">All</option>';
      Object.keys(seen).sort(function (a, b) { return a.localeCompare(b); }).forEach(function (o) {
        var op = document.createElement('option'); op.value = o; op.textContent = o; sel.appendChild(op);
      });
      sel.addEventListener('change', apply);
      cell.appendChild(sel);
      controls.push({ i: i, kind: kind, el: sel });
    } else if (kind === 'range') {
      var wrap = document.createElement('div'); wrap.className = 'range';
      var lo = document.createElement('input'), hi = document.createElement('input');
      lo.type = hi.type = 'number'; lo.placeholder = 'min'; hi.placeholder = 'max';
      lo.addEventListener('input', apply); hi.addEventListener('input', apply);
      wrap.appendChild(lo); wrap.appendChild(hi); cell.appendChild(wrap);
      controls.push({ i: i, kind: kind, lo: lo, hi: hi });
    } else if (kind === 'text') {
      var inp = document.createElement('input'); inp.type = 'text'; inp.placeholder = 'filter';
      inp.addEventListener('input', apply);
      cell.appendChild(inp);
      controls.push({ i: i, kind: kind, el: inp });
    }
    filterRow.appendChild(cell);
  });
  thead.appendChild(filterRow);

  function matches(tr) {
    if (searchInput && searchInput.value &&
        tr.textContent.toLowerCase().indexOf(searchInput.value.toLowerCase()) === -1) return false;
    if (flaggedOnly && flaggedOnly.checked && !tr.querySelector('.flag')) return false;
    for (var c = 0; c < controls.length; c++) {
      var ctl = controls[c], val = sortVal(tr, ctl.i);
      if (ctl.kind === 'select') {
        if (ctl.el.value && val.trim() !== ctl.el.value) return false;
      } else if (ctl.kind === 'multi') {
        if (ctl.el.value && val.split(',').map(function (t) { return t.trim(); }).indexOf(ctl.el.value) === -1) return false;
      } else if (ctl.kind === 'text') {
        if (ctl.el.value && val.toLowerCase().indexOf(ctl.el.value.toLowerCase()) === -1) return false;
      } else if (ctl.kind === 'range') {
        var n = parseFloat(val);
        if (ctl.lo.value !== '' && (isNaN(n) || n < parseFloat(ctl.lo.value))) return false;
        if (ctl.hi.value !== '' && (isNaN(n) || n > parseFloat(ctl.hi.value))) return false;
      }
    }
    return true;
  }

  function apply() {
    var shown = 0, total = tbody.rows.length;
    Array.prototype.slice.call(tbody.rows).forEach(function (tr) {
      var ok = matches(tr); tr.style.display = ok ? '' : 'none'; if (ok) shown++;
    });
    if (countEl) countEl.textContent = 'Showing ' + shown + ' of ' + total;
  }

  if (searchInput) searchInput.addEventListener('input', apply);
  if (flaggedOnly) flaggedOnly.addEventListener('change', apply);
  apply();
}
document.querySelectorAll('table.data-table').forEach(enhanceTable);

// Per-genome detail: text filter + expand/collapse all.
(function () {
  var dq = document.getElementById('detail-q');
  var details = Array.prototype.slice.call(document.querySelectorAll('details[id^="g-"]'));
  var dc = document.getElementById('detail-count');
  function applyDetail() {
    var q = dq ? dq.value.toLowerCase() : '', shown = 0;
    details.forEach(function (d) {
      var s = d.querySelector('summary');
      var ok = !q || (s && s.textContent.toLowerCase().indexOf(q) > -1);
      d.style.display = ok ? '' : 'none'; if (ok) shown++;
    });
    if (dc) dc.textContent = 'Showing ' + shown + ' of ' + details.length;
  }
  if (dq) dq.addEventListener('input', applyDetail);
  var ea = document.getElementById('expand-all'), ca = document.getElementById('collapse-all');
  if (ea) ea.addEventListener('click', function () { details.forEach(function (d) { d.open = true; }); });
  if (ca) ca.addEventListener('click', function () { details.forEach(function (d) { d.open = false; }); });
  applyDetail();
})();
</script>
</body></html>
"""


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--report-dir", type=Path, default=Path("."),
                    help="Directory containing the aggregated 10_report TSVs.")
    ap.add_argument("--lps-db-dir", type=Path, required=True,
                    help="LPS database directory (reference_LPS.txt, .gb files, "
                         "phenotype_lookup.tsv, phenotype_images/, gene_colors.tsv).")
    ap.add_argument("--out", type=Path, default=Path("LPS_typing_report.html"))
    ap.add_argument("--pipeline-name", default="Pasteurella multocida LPS typing pipeline")
    ap.add_argument("--pipeline-version", default="")
    ap.add_argument("--steps", default="", help="Comma-separated steps that ran.")
    ap.add_argument("--skipped", default="", help="Comma-separated steps skipped.")
    ap.add_argument("--params", default="", help="Free-text note of key parameters.")
    ap.add_argument("--run-date", default="", help="Override report date (YYYY-MM-DD).")
    args = ap.parse_args(argv)

    rdir = args.report_dir
    lps = args.lps_db_dir

    report_df = read_tsv(rdir / F_SUBTYPE)
    if report_df is None:
        sys.exit(f"error: required subtype report not found/empty: {rdir / F_SUBTYPE}")

    qc = {
        "checkm": index_by_sample(read_tsv(rdir / F_CHECKM), "sampleID"),
        "taxonomy": parse_sylph(read_tsv(rdir / F_SYLPH)),
        "kaptive": index_by_sample(read_tsv(rdir / F_KAPTIVE), "sampleID"),
        "quast": parse_quast(read_tsv(rdir / F_QUAST)),
        "shovill": index_by_sample(read_tsv(rdir / F_SHOVILL), "sample"),
    }

    ref_map = load_reference_map(lps)
    gene_colors = load_gene_colors(lps)
    pheno_images = load_phenotype_images(lps)

    records = build_records(report_df, qc, read_tsv(rdir / F_PHENOTYPE))
    attach_amr(records, read_tsv(rdir / F_AMR))

    # attach per-genome phenotype image (first phenotype with an image)
    for r in records.values():
        r["pheno_img"] = None
        for p in r["phenotypes"]:
            if p in pheno_images:
                r["pheno_img"] = pheno_images[p]
                break

    # distributions
    type_counts = OrderedDict()
    pheno_counts = OrderedDict()
    mlst_counts = OrderedDict()
    pheno_desc = {}
    for r in records.values():
        t = r["type"] or "untypeable"
        type_counts[t] = type_counts.get(t, 0) + 1
        if r["mlst"]:
            mlst_counts[r["mlst"]] = mlst_counts.get(r["mlst"], 0) + 1
        # An uncertain call lists several alternatives; count it once, as "uncertain".
        if len(r["phenotypes"]) > 1:
            pheno_counts["uncertain"] = pheno_counts.get("uncertain", 0) + 1
        for p, d in r["phenotypes"].items():
            if len(r["phenotypes"]) == 1:
                pheno_counts[p] = pheno_counts.get(p, 0) + 1
            pheno_desc[p] = d

    distributions = []
    if len(records) > 1:
        distributions.append(make_bar(OrderedDict(sorted(type_counts.items())),
                                      "Genomes per LPS type", "genomes"))
        if pheno_counts:
            distributions.append(make_bar(OrderedDict(sorted(pheno_counts.items())),
                                          "Genomes per phenotype", "genomes"))
        if mlst_counts:
            distributions.append(make_bar(OrderedDict(sorted(mlst_counts.items())),
                                          "Genomes per MLST sequence type", "genomes"))

    # phenotype gallery (observed, sorted)
    pheno_gallery = []
    for p in sorted(k for k in pheno_counts if k != "uncertain"):
        pheno_gallery.append({
            "id": p, "desc": (": " + pheno_desc[p]) if pheno_desc.get(p) else "",
            "count": pheno_counts[p], "img": pheno_images.get(p),
        })

    lollipops = build_lollipops(report_df, lps, ref_map, gene_colors)

    n_flagged = sum(1 for r in records.values() if r["flags"])
    has_amr = any(r["amr"] for r in records.values())

    db_files = ", ".join(sorted(
        p.name for p in [lps / "LPS_subtype_database_v3.txt",
                         lps / "phenotype_rules.tsv",
                         lps / "phenotype_lookup.tsv"] if p.exists()))

    meta = {
        "pipeline_name": args.pipeline_name,
        "pipeline_version": args.pipeline_version,
        "date": args.run_date or date.today().isoformat(),
        "steps": args.steps.replace(",", ", "),
        "skipped": args.skipped.replace(",", ", "),
        "params": args.params,
        "db_files": db_files or "(LPS database)",
    }

    html = Template(TEMPLATE).render(
        meta=meta, records=records, distributions=distributions,
        pheno_gallery=pheno_gallery, lollipops=lollipops,
        type_counts=type_counts, pheno_counts=pheno_counts,
        n_flagged=n_flagged, has_amr=has_amr,
        multiqc=F_MULTIQC if (rdir / F_MULTIQC).exists() else None,
    )
    args.out.write_text(html, encoding="utf-8")
    sys.stderr.write(f"wrote {args.out} ({len(records)} genomes, "
                     f"{len(lollipops)} lollipop plots)\n")


if __name__ == "__main__":
    main()
