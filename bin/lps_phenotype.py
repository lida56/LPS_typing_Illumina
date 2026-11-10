#!/usr/bin/env python3
"""Rule-based LPS phenotype assignment shared by the Illumina and ONT pipelines.

This file is kept byte-identical in both pipeline repositories (a test checks it).
It uses the Python standard library only, so it runs in any container with python3.

Overview
--------
1. Every called variant is assigned to the LPS genes it overlaps, using the CDS
   coordinates of the type's reference GenBank file.
2. Each gene of the type's phenotype rules gets one state:
   ``inactivating`` > ``uncertain`` > ``modifying`` > ``functional``.
   The state comes from curated subtype database EFFECT values, snpEff effects and
   the net indel length in the gene (no translation is done).
3. petG is classified from its blastn alignment against the X73 reference as
   ``functional``, ``nonfunctional``, ``uncertain`` or ``absent``.
4. The phenotype is the rule row of the lowest inactivated biosynthetic step that
   matches the petG state (``phenotype_rules.tsv``). Uncertain genes or an
   uncertain petG produce every distinct alternative phenotype.

Outputs
-------
* ``<prefix>_phenotype_report.tsv``: one row per sample.
* ``<prefix>_subtype_report.tsv``: one row per matched subtype database variant
  (plus one fallback row per sample without matches).
* ``<prefix>_petG_status.tsv``: petG status per sample.
"""

import argparse
import csv
import itertools
import re
import sys
from collections import OrderedDict
from pathlib import Path

STATES_ORDER = ("inactivating", "uncertain", "modifying", "functional")
EFFECT_VALUES = {"", "none", "modifying", "inactivating"}
PETG_VALUES = {"functional", "nonfunctional", "any"}
STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODONS = {"ATG", "GTG", "TTG"}
NA_VALUES = {"", "NA", "na", "N/A", ".", "-"}
DEFAULT_MIN_COVERAGE = 0.90
# Allele fraction below which a call is reported as mixed (same cut-off as the
# Illumina snippy_typing_calls.sh "heterogeneous" flag).
MIXED_AF = 0.90

# snpEff / Sequence Ontology terms that inactivate a gene on their own.
INACTIVATING_TERMS = {"stop_gained", "start_lost"}
# Terms that leave the reading frame open-ended; their effect is not predictable.
UNCERTAIN_TERMS = {"stop_lost"}


def is_na(value):
    return value is None or str(value).strip() in NA_VALUES


def clean(value):
    return "" if is_na(value) else str(value).strip()


def read_tsv(path):
    """Read a TSV with a header row into a list of dicts (empty list if no file)."""
    if path is None:
        return []
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, columns, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for row in rows:
            writer.writerow([row.get(column, "") for column in columns])


# --------------------------------------------------------------------------- #
# Reference annotation
# --------------------------------------------------------------------------- #

class Gene:
    """A CDS on an LPS reference record, 1-based inclusive coordinates."""

    def __init__(self, name, start, end, strand):
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand

    @property
    def length(self):
        return self.end - self.start + 1

    def start_codon_span(self):
        if self.strand == "+":
            return self.start, self.start + 2
        return self.end - 2, self.end

    def stop_codon_span(self):
        if self.strand == "+":
            return self.end - 2, self.end
        return self.start, self.start + 2


def parse_genbank(path):
    """Return (record name, sequence, [Gene]) from a single-record GenBank file."""
    text = Path(path).read_text()
    locus = re.search(r"^LOCUS\s+(\S+)", text, re.M)
    record = locus.group(1) if locus else Path(path).stem
    genes = []
    feature_re = re.compile(
        r"^ {5}CDS {13}(complement\()?<?(\d+)\.\.>?(\d+)\)?[^\n]*\n((?: {21}[^\n]*\n)*)",
        re.M,
    )
    for match in feature_re.finditer(text):
        qualifiers = match.group(4)
        gene = re.search(r'/gene="([^"]+)"', qualifiers)
        if not gene:
            continue
        strand = "-" if match.group(1) else "+"
        genes.append(Gene(gene.group(1), int(match.group(2)), int(match.group(3)), strand))
    sequence = ""
    if "ORIGIN" in text:
        sequence = re.sub(r"[^acgtnACGTN]", "", text.split("ORIGIN", 1)[1].split("//")[0]).upper()
    return record, sequence, genes


class Reference:
    def __init__(self, lps_type, gb_path):
        self.type = lps_type
        self.gb_path = Path(gb_path)
        self.record, self.sequence, gene_list = parse_genbank(gb_path)
        self.genes = OrderedDict((gene.name, gene) for gene in gene_list)


def load_references(lps_db_dir):
    """Map LPS type -> Reference using reference_LPS.txt (TYPE, GenBank, FASTA)."""
    lps_db_dir = Path(lps_db_dir)
    references = {}
    table = lps_db_dir / "reference_LPS.txt"
    for line in table.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 2 or not fields[0].strip():
            continue
        gb_path = lps_db_dir / fields[1].strip()
        if gb_path.is_file():
            references[fields[0].strip()] = Reference(fields[0].strip(), gb_path)
    return references


# --------------------------------------------------------------------------- #
# Phenotype rules
# --------------------------------------------------------------------------- #

class Rule:
    def __init__(self, lps_type, step, genes, petg, requires, phenotype):
        self.type = lps_type
        self.step = step  # int, or None for "full"
        self.genes = genes
        self.petg = petg
        self.requires = requires  # list of (gene, state)
        self.phenotype = phenotype


class RuleSet:
    """All rule rows of one LPS type."""

    def __init__(self, lps_type):
        self.type = lps_type
        self.rules = []

    @property
    def step_genes(self):
        """Ordered mapping step -> genes (steps ascending)."""
        steps = OrderedDict()
        for rule in sorted((r for r in self.rules if r.step is not None), key=lambda r: r.step):
            steps.setdefault(rule.step, [])
            for gene in rule.genes:
                if gene not in steps[rule.step]:
                    steps[rule.step].append(gene)
        return steps

    @property
    def genes(self):
        names = []
        for genes in self.step_genes.values():
            names.extend(g for g in genes if g not in names)
        for rule in self.rules:
            for gene, _ in rule.requires:
                if gene not in names:
                    names.append(gene)
        return names

    @property
    def petg_relevant(self):
        return any(rule.petg != "any" for rule in self.rules)


def parse_requires(value):
    value = clean(value)
    if not value:
        return []
    conditions = []
    for item in value.split(";"):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(f"REQUIRES entry '{item}' must look like gene:state")
        gene, state = (part.strip() for part in item.split(":", 1))
        if state not in STATES_ORDER:
            raise ValueError(f"REQUIRES state '{state}' must be one of {', '.join(STATES_ORDER)}")
        conditions.append((gene, state))
    return conditions


def load_rules(path, references=None):
    """Load and validate phenotype_rules.tsv into {type: RuleSet}."""
    rulesets = OrderedDict()
    errors = []
    seen = set()
    for line_no, row in enumerate(read_tsv(path), start=2):
        lps_type = clean(row.get("TYPE"))
        step_value = clean(row.get("STEP"))
        petg = clean(row.get("PETG")) or "any"
        phenotype = clean(row.get("PHENOTYPE"))
        genes = [g.strip() for g in clean(row.get("GENES")).split(",") if g.strip()]
        where = f"{path}:{line_no}"
        if not lps_type or not phenotype:
            errors.append(f"{where}: TYPE and PHENOTYPE are required")
            continue
        if step_value == "full":
            step = None
            if genes:
                errors.append(f"{where}: a 'full' row must not list GENES")
        else:
            try:
                step = int(step_value)
            except ValueError:
                errors.append(f"{where}: STEP must be an integer or 'full', got '{step_value}'")
                continue
            if not genes:
                errors.append(f"{where}: STEP {step} needs at least one gene in GENES")
        if petg not in PETG_VALUES:
            errors.append(f"{where}: PETG must be one of {', '.join(sorted(PETG_VALUES))}")
        try:
            requires = parse_requires(row.get("REQUIRES"))
        except ValueError as exc:
            errors.append(f"{where}: {exc}")
            continue
        key = (lps_type, step_value, petg, tuple(requires))
        if key in seen:
            errors.append(f"{where}: duplicate rule for TYPE/STEP/PETG/REQUIRES {key}")
        seen.add(key)
        ruleset = rulesets.setdefault(lps_type, RuleSet(lps_type))
        ruleset.rules.append(Rule(lps_type, step, genes, petg, requires, phenotype))

    for lps_type, ruleset in rulesets.items():
        steps = {}
        for rule in ruleset.rules:
            if rule.step is not None:
                for gene in rule.genes:
                    if steps.setdefault(gene, rule.step) != rule.step:
                        errors.append(f"{lps_type}: gene {gene} is listed at more than one STEP")
        if not any(rule.step is None for rule in ruleset.rules):
            errors.append(f"{lps_type}: no 'full' rule row")
        if references is not None:
            reference = references.get(lps_type)
            if reference is None:
                errors.append(f"{lps_type}: no reference GenBank in reference_LPS.txt")
                continue
            for gene in ruleset.genes:
                if gene not in reference.genes:
                    errors.append(f"{lps_type}: gene {gene} is not a CDS in {reference.gb_path.name}")
    if errors:
        raise ValueError("Invalid phenotype rules:\n  " + "\n  ".join(errors))
    return rulesets


# --------------------------------------------------------------------------- #
# Subtype database
# --------------------------------------------------------------------------- #

class DbVariant:
    def __init__(self, row):
        self.type = clean(row.get("TYPE"))
        self.subtype = clean(row.get("SUBTYPE"))
        self.isolate = clean(row.get("ISOLATE"))
        self.chrom = clean(row.get("CHROM"))
        self.pos = clean(row.get("POS"))
        self.vartype = clean(row.get("VARTYPE"))
        self.ref = clean(row.get("REF"))
        self.alt = clean(row.get("ALT"))
        self.gene = clean(row.get("GENE"))
        self.effect = clean(row.get("EFFECT")).lower()
        self.note = clean(row.get("NOTE"))

    @property
    def key(self):
        return (self.chrom, self.pos, self.ref, self.alt)

    @property
    def matchable(self):
        return all(self.key)


def load_subtype_db(path):
    entries = [DbVariant(row) for row in read_tsv(path)]
    bad = [e.subtype for e in entries if e.effect not in EFFECT_VALUES]
    if bad:
        raise ValueError(
            "Invalid EFFECT in subtype database for "
            + ", ".join(bad)
            + "; allowed: blank, none, modifying, inactivating"
        )
    index = {}
    for entry in entries:
        if entry.matchable:
            index.setdefault(entry.key, []).append(entry)
    return entries, index


# --------------------------------------------------------------------------- #
# Variant calls
# --------------------------------------------------------------------------- #

class Variant:
    def __init__(self, sample, chrom, pos, ref, alt, effects=None, af=None, mixed=False):
        self.sample = sample
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref.upper()
        self.alt = alt.upper()
        self.effects = set(effects or ())
        self.af = af
        self.mixed = mixed
        self.db = []  # matched DbVariant entries
        # Why an uncurated call is not trusted on its own (None = confident call)
        self.low_confidence = None

    @property
    def key(self):
        return (self.chrom, str(self.pos), self.ref, self.alt)

    @property
    def length_change(self):
        return len(self.alt) - len(self.ref)

    @property
    def is_indel(self):
        return self.length_change != 0

    def affected_span(self):
        """Return (start, end) of reference bases changed, or (pos, None) for an insertion point.

        For an insertion the returned start is the base after which bases are inserted.
        """
        prefix = 0
        while prefix < min(len(self.ref), len(self.alt)) and self.ref[prefix] == self.alt[prefix]:
            prefix += 1
        ref_rest = self.ref[prefix:]
        if not ref_rest:
            return self.pos + prefix - 1, None
        return self.pos + prefix, self.pos + len(self.ref) - 1

    def overlaps(self, start, end):
        span_start, span_end = self.affected_span()
        if span_end is None:
            # insertion between span_start and span_start + 1
            return start <= span_start < end
        return span_start <= end and span_end >= start

    def label(self):
        change = self.length_change
        if change > 0:
            kind = f"{change}bp ins"
        elif change < 0:
            kind = f"{-change}bp del"
        else:
            kind = f"{self.ref}>{self.alt}"
        return f"{kind} at {self.pos}"


def homopolymer_run(sequence, variant):
    """Length of the reference run of the inserted/deleted base at an indel, else 0.

    For an insertion the run is counted on both sides of the insertion point; for a
    deletion it includes the deleted bases. Uses the plain reference sequence only.
    """
    if not variant.is_indel or not sequence:
        return 0
    prefix = 0
    while prefix < min(len(variant.ref), len(variant.alt)) and variant.ref[prefix] == variant.alt[prefix]:
        prefix += 1
    changed = variant.alt[prefix:] if variant.length_change > 0 else variant.ref[prefix:]
    if not changed or len(set(changed)) != 1:
        return 0
    base = changed[0]
    if variant.length_change > 0:
        left = variant.pos + prefix - 1  # inserted after this 1-based reference position
        right = left + 1
    else:
        left = variant.pos + prefix - 1
        right = variant.pos + prefix + len(changed)
    run = 0
    i = left
    while i >= 1 and sequence[i - 1] == base:
        run += 1
        i -= 1
    j = right
    while j <= len(sequence) and sequence[j - 1] == base:
        run += 1
        j += 1
    if variant.length_change < 0:
        run += len(changed)
    return run


def format_af(af):
    """Allele fraction with 3 decimals, so a value just below a cut-off (0.8988) never reads as 0.90."""
    return f"{af:.3f}"


def can_change_state(variant):
    """Only indels and stop/start effects can change a gene state; plain substitutions cannot."""
    return variant.is_indel or bool(variant.effects & (INACTIVATING_TERMS | UNCERTAIN_TERMS))


def assess_confidence(variant, sequence, min_af=MIXED_AF, homopolymer_min_run=0):
    """Set variant.low_confidence for an uncurated call that is not trusted on its own.

    Curated variants always count. An uncurated call is low-confidence when its allele
    fraction is below ``min_af``, or (when ``homopolymer_min_run`` > 0, used for ONT)
    when it is an indel of a base that sits in a reference run of at least that length.
    """
    if variant.db:
        variant.low_confidence = None
        return
    reasons = []
    if variant.af is not None and variant.af < min_af:
        reasons.append(f"AF={format_af(variant.af)}")
    if homopolymer_min_run:
        run = homopolymer_run(sequence, variant)
        if run >= homopolymer_min_run:
            base = (variant.alt if variant.length_change > 0 else variant.ref)[-1]
            reasons.append(f"{'extends' if variant.length_change > 0 else 'in'} a {base * run} run")
    variant.low_confidence = ", ".join(reasons) or None


def split_effect_terms(text):
    """Return the set of SO terms in a snippy EFFECT or snpEff annotation string."""
    text = clean(text)
    if not text:
        return set()
    first = text.split(" ", 1)[0]
    return {term for term in first.split("&") if term}


def read_snippy_variants(path):
    """Read the aggregated snippy table (sampleID + snps.tab columns)."""
    variants = []
    for row in read_tsv(path):
        sample = clean(row.get("sampleID") or row.get("SAMPLE"))
        if not sample or is_na(row.get("POS")):
            continue
        variants.append(Variant(
            sample, clean(row.get("CHROM")), row["POS"], clean(row.get("REF")), clean(row.get("ALT")),
            split_effect_terms(row.get("EFFECT")),
        ))
    return variants


def parse_vcf_info(info):
    result = {}
    for item in info.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        elif item:
            result[item] = True
    return result


def read_vcf_variants(path, sample_from_file=None, pass_only=True):
    """Read variants from a (possibly multi-sample-aggregated) snpEff-annotated VCF.

    The aggregated ONT table prefixes each record with the sample ID (first column);
    a plain VCF is read with ``sample_from_file`` as the sample name.
    """
    variants = []
    path = Path(path)
    if not path.is_file():
        return variants
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if sample_from_file is None:
                sample, fields = fields[0], fields[1:]
            else:
                sample = sample_from_file
            if len(fields) < 8 or fields[0] in ("CHROM", "#CHROM"):
                continue
            chrom, pos, _id, ref, alts, _qual, filt, info = fields[:8]
            if pass_only and filt not in ("PASS", "."):
                continue
            genotype = ""
            if len(fields) >= 10:
                fmt = fields[8].split(":")
                values = fields[9].split(":")
                genotype = dict(zip(fmt, values)).get("GT", "")
            alleles = re.split(r"[/|]", genotype) if genotype else []
            info_map = parse_vcf_info(info)
            annotations = info_map.get("ANN", "") if isinstance(info_map.get("ANN", ""), str) else ""
            for index, alt in enumerate(alts.split(","), start=1):
                if alt in ("*", ".") or alt.startswith("<"):
                    continue
                if alleles and str(index) not in alleles:
                    continue
                effects = set()
                for annotation in annotations.split(","):
                    parts = annotation.split("|")
                    if len(parts) > 1 and parts[0] == alt:
                        effects.update(t for t in parts[1].split("&") if t)
                af = None
                if len(fields) >= 10:
                    fmt_map = dict(zip(fields[8].split(":"), fields[9].split(":")))
                    try:
                        af = float(fmt_map.get("AF", "").split(",")[index - 1])
                    except (ValueError, IndexError):
                        af = None
                variants.append(Variant(sample, chrom, pos, ref, alt, effects, af=af,
                                        mixed=af is not None and af < MIXED_AF))
    return variants


def read_rescued(paths):
    """Return {(sample, chrom, pos, ref, alt): af} for heterogeneous rescued 0/1 calls."""
    mixed = {}
    for path in paths:
        path = Path(path)
        sample = path.name.split("_snps.rescued.tsv")[0]
        for row in read_tsv(path):
            if clean(row.get("heterogeneous")).lower() not in ("yes", "true", "1"):
                continue
            try:
                af = float(row.get("AF"))
            except (TypeError, ValueError):
                af = None
            mixed[(sample, clean(row.get("CHROM")), clean(row.get("POS")),
                   clean(row.get("REF")).upper(), clean(row.get("ALT")).upper())] = af
    return mixed


# --------------------------------------------------------------------------- #
# Coverage
# --------------------------------------------------------------------------- #

def read_coverage(path, min_depth=1):
    """Return {chrom: set(covered 1-based positions)} or None if unavailable.

    Accepts a snippy ``snps.aligned.fa`` (``-`` = no coverage, ``N`` = below
    --mincov) or a ``samtools depth -a`` table (chrom, pos, depth).
    """
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return None
    covered = {}
    with path.open() as handle:
        first = handle.readline()
        handle.seek(0)
        if first.startswith(">"):
            name, chunks = None, []
            for line in itertools.chain(handle, [">"]):
                line = line.strip()
                if line.startswith(">"):
                    if name is not None:
                        sequence = "".join(chunks)
                        covered[name] = {i + 1 for i, base in enumerate(sequence) if base not in "-Nn"}
                    name, chunks = (line[1:].split()[0] if len(line) > 1 else None), []
                else:
                    chunks.append(line)
        else:
            for line in handle:
                fields = line.split("\t")
                if len(fields) < 3:
                    continue
                try:
                    if int(fields[2]) >= min_depth:
                        covered.setdefault(fields[0], set()).add(int(fields[1]))
                except ValueError:
                    continue
    return covered


def gene_coverage(covered, chrom, gene):
    positions = covered.get(chrom)
    if positions is None:
        return 0.0
    hits = sum(1 for p in range(gene.start, gene.end + 1) if p in positions)
    return hits / gene.length


# --------------------------------------------------------------------------- #
# Gene states
# --------------------------------------------------------------------------- #

class GeneCall:
    def __init__(self, gene, state, reason="", variants=None, uncurated=False, mixed=None):
        self.gene = gene
        self.state = state
        self.reason = reason
        self.variants = variants or []
        self.uncurated = uncurated
        self.mixed = mixed or []

    def describe(self):
        return f"{self.gene}:{self.state}"


def subtypes_of(variants):
    names = []
    for variant in variants:
        for entry in variant.db:
            if entry.subtype not in names:
                names.append(entry.subtype)
    return names


def classify_gene(gene, variants, coverage_fraction=None, min_coverage=DEFAULT_MIN_COVERAGE):
    """Decide the state of one gene from the variants that overlap it.

    Low-confidence uncurated calls (see assess_confidence) are not trusted either way:
    if the state changes when they are left out, the gene is ``uncertain``.
    """
    call = _classify(gene, variants, coverage_fraction, min_coverage)
    weak = [v for v in variants if v.low_confidence and can_change_state(v)]
    if not weak:
        return call
    without = _classify(gene, [v for v in variants if v not in weak], coverage_fraction, min_coverage)
    if without.state == call.state:
        return call
    detail = "; ".join(f"{v.label()}: {v.low_confidence}" for v in weak)
    used = list(dict.fromkeys(call.variants + weak))
    return GeneCall(gene.name, "uncertain", f"depends on a low-confidence call ({detail})", used, True,
                    [v for v in used if v.mixed and v.db])


def _classify(gene, variants, coverage_fraction=None, min_coverage=DEFAULT_MIN_COVERAGE):
    """Gene state from its variants, taking every call at face value."""
    curated = [v for v in variants if v.db]
    curated_effects = {e.effect for v in curated for e in v.db}

    def build(state, reason, used):
        uncurated = any(not v.db for v in used) and state in ("inactivating", "uncertain")
        # Uncurated calls are explained by the confidence gate; the suffix is for curated ones.
        mixed = [v for v in used if v.mixed and v.db]
        return GeneCall(gene.name, state, reason, used, uncurated, mixed)

    curated_inactivating = [v for v in curated if any(e.effect == "inactivating" for e in v.db)]
    if curated_inactivating:
        names = ", ".join(subtypes_of(curated_inactivating))
        return build("inactivating", f"curated inactivating ({names})", curated_inactivating)

    # Variants curated as having no effect are left out of the arithmetic.
    active = [v for v in variants if not (v.db and all(e.effect == "none" for e in v.db))]

    for variant in active:
        terms = variant.effects & INACTIVATING_TERMS
        if terms:
            return build("inactivating", f"{sorted(terms)[0]} {variant.label()}", [variant])

    start_span = gene.start_codon_span()
    stop_span = gene.stop_codon_span()
    for variant in active:
        if variant.is_indel and (variant.overlaps(*start_span) or variant.overlaps(*stop_span)):
            where = "start" if variant.overlaps(*start_span) else "stop"
            return build("inactivating", f"indel in {where} codon {variant.label()}", [variant])

    indels = [v for v in active if v.is_indel]
    net = sum(v.length_change for v in indels)
    if net % 3 != 0:
        sign = "+" if net > 0 else ""
        return build("inactivating", f"frameshift, net {sign}{net} bp ({'; '.join(v.label() for v in indels)})", indels)
    if indels:
        sign = "+" if net > 0 else ""
        if any(v.length_change % 3 != 0 for v in indels):
            reason = f"compensating frameshifts, net {sign}{net} bp"
        else:
            reason = f"in-frame indel, net {sign}{net} bp"
        names = subtypes_of(indels)
        if names:
            reason += f" ({'+'.join(names)})"
        return build("uncertain", reason, indels)

    stop_lost = [v for v in active if v.effects & UNCERTAIN_TERMS]
    if stop_lost:
        return build("uncertain", f"stop_lost {stop_lost[0].label()}", stop_lost)

    if coverage_fraction is not None and coverage_fraction < min_coverage:
        return GeneCall(gene.name, "uncertain", f"only {coverage_fraction:.0%} covered", [], False)

    modifying = [v for v in curated if any(e.effect == "modifying" for e in v.db)]
    if modifying:
        names = ", ".join(subtypes_of(modifying))
        return build("modifying", f"curated modifying ({names})", modifying)

    return GeneCall(gene.name, "functional", "", variants)


# --------------------------------------------------------------------------- #
# petG
# --------------------------------------------------------------------------- #

PETG_BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore", "qlen", "qseq", "sseq",
]


def classify_petg_alignment(qseq, sseq, qstart, qlen):
    """Classify a petG blastn alignment (query = X73 reference CDS, plus strand).

    Returns (status, detail). No translation: frameshifts come from gap arithmetic,
    nonsense and start loss from looking up the sample codon at changed codons.
    """
    qseq = qseq.upper()
    sseq = sseq.upper()
    # Build the sample bases aligned to each reference position (1-based) and
    # collect insertions (relative to the reference) after a reference position.
    sample_at = {}
    inserted_after = {}
    ref_pos = qstart - 1
    for q_base, s_base in zip(qseq, sseq):
        if q_base == "-":
            inserted_after[ref_pos] = inserted_after.get(ref_pos, 0) + 1
            continue
        ref_pos += 1
        sample_at[ref_pos] = s_base
    deletions = [p for p, base in sample_at.items() if base == "-"]
    events = []
    for pos, count in sorted(inserted_after.items()):
        events.append((pos, count))
    run_start, run_len = None, 0
    for pos in sorted(deletions) + [None]:
        if pos is not None and run_start is not None and pos == run_start + run_len:
            run_len += 1
            continue
        if run_start is not None:
            events.append((run_start, -run_len))
        run_start, run_len = (pos, 1) if pos is not None else (None, 0)
    events.sort()
    net = sum(change for _, change in events)

    def fmt(change):
        return f"+{change}" if change > 0 else str(change)

    if net % 3 != 0:
        first = next((e for e in events if e[1] % 3 != 0), events[0])
        return "nonfunctional", f"frameshift: net {fmt(net)} bp (first {fmt(first[1])} bp at ref {first[0]})"

    # Codon checks only make sense where the frame is the reference frame, i.e.
    # before the first indel event.
    first_event = events[0][0] if events else qlen + 1
    ref_stop_codon = qlen - 2
    start_codon = "".join(sample_at.get(p, "") for p in range(1, 4))
    if qstart == 1 and len(start_codon) == 3 and "-" not in start_codon and start_codon not in START_CODONS:
        return "nonfunctional", f"start lost: {start_codon}"
    for codon_start in range(1, min(first_event, ref_stop_codon), 3):
        codon = "".join(sample_at.get(p, "") for p in range(codon_start, codon_start + 3))
        if len(codon) == 3 and "-" not in codon and codon in STOP_CODONS:
            return "nonfunctional", f"premature stop {codon} at codon {(codon_start + 2) // 3}"
    if events:
        if any(change % 3 != 0 for _, change in events):
            return "uncertain", f"compensating frameshifts, net {fmt(net)} bp"
        return "uncertain", f"in-frame indel, net {fmt(net)} bp at ref {events[0][0]}"
    return "functional", ""


def classify_petg_table(path, min_identity, min_length, max_end_offset=3):
    """Classify one sample's raw petG blast table. Returns (status, detail)."""
    rows = []
    for line in Path(path).read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < len(PETG_BLAST_COLUMNS) or fields[0] == "qseqid":
            continue
        rows.append(dict(zip(PETG_BLAST_COLUMNS, fields)))
    passing = []
    for row in rows:
        qstart, qend, qlen = int(row["qstart"]), int(row["qend"]), int(row["qlen"])
        sstart, send = int(row["sstart"]), int(row["send"])
        span = abs(send - sstart) + 1
        if (float(row["pident"]) >= min_identity and span > min_length
                and qstart <= 1 + max_end_offset and qend >= qlen - max_end_offset):
            passing.append(row)
    if not passing:
        return "absent", "no X73-like petG hit passing identity, length and coverage thresholds"
    best = max(passing, key=lambda r: float(r["bitscore"]))
    status, detail = classify_petg_alignment(best["qseq"], best["sseq"], int(best["qstart"]), int(best["qlen"]))
    location = f"{best['sseqid']}:{min(int(best['sstart']), int(best['send']))}-{max(int(best['sstart']), int(best['send']))}"
    detail = f"{detail}; {location}" if detail else location
    if len(passing) > 1:
        detail += f"; {len(passing)} passing hits, best bitscore used"
    return status, detail


# --------------------------------------------------------------------------- #
# Phenotype
# --------------------------------------------------------------------------- #

def requires_met(rule, states):
    return all(states.get(gene) == state for gene, state in rule.requires)


def evaluate(ruleset, states, petg_state):
    """Return (phenotype, determining gene or '') for definite gene/petG states."""
    determining_step, determining_gene = None, ""
    for step, genes in ruleset.step_genes.items():
        inactive = [g for g in genes if states.get(g) == "inactivating"]
        if inactive:
            determining_step, determining_gene = step, inactive[0]
            break
    candidates = [r for r in ruleset.rules if r.step == determining_step]
    candidates = [r for r in candidates if r.petg in ("any", petg_state)]
    for rule in sorted(candidates, key=lambda r: 0 if r.requires else 1):
        if requires_met(rule, states):
            return rule.phenotype, determining_gene
    return "", determining_gene


def petg_for_rules(status):
    """petG states the rules see; 'uncertain' expands to both alternatives."""
    if status == "functional":
        return ["functional"]
    if status in ("uncertain", "not assessed"):
        return ["functional", "nonfunctional"]
    return ["nonfunctional"]  # nonfunctional or absent


def predict(ruleset, gene_calls, petg_status):
    """Return (phenotypes in order, determining genes, alternatives detail)."""
    base_states = {name: call.state for name, call in gene_calls.items()}
    uncertain_genes = [name for name, call in gene_calls.items() if call.state == "uncertain"]
    petg_options = petg_for_rules(petg_status) if ruleset.petg_relevant else ["functional"]
    outcomes = OrderedDict()
    for combo in itertools.product(("functional", "inactivating"), repeat=len(uncertain_genes)):
        states = dict(base_states)
        states.update(zip(uncertain_genes, combo))
        for petg_state in petg_options:
            phenotype, gene = evaluate(ruleset, states, petg_state)
            outcomes.setdefault(phenotype, []).append((dict(zip(uncertain_genes, combo)), petg_state, gene))
    return outcomes


# --------------------------------------------------------------------------- #
# Kaptive / MLST
# --------------------------------------------------------------------------- #

def read_kaptive(path):
    """Return {sample: type or 'untypeable'} from the aggregated Kaptive table."""
    types = OrderedDict()
    for row in read_tsv(path):
        sample = clean(row.get("sampleID") or row.get("Assembly"))
        if not sample:
            continue
        sample = re.sub(r"_contigs$", "", sample)
        if clean(row.get("Match confidence")) == "Typeable":
            locus = clean(row.get("Best match locus")).split("-")[0]
            types[sample] = locus.replace("LPS", "L")
        else:
            types[sample] = "untypeable"
    return types


def read_mlst(paths):
    """Return {sample: ST} from per-sample <sample>_mlst.csv files (ST is column 3)."""
    result = {}
    for path in paths:
        path = Path(path)
        if not path.is_file() or not path.name.endswith("_mlst.csv"):
            continue
        lines = path.read_text().splitlines()
        fields = lines[0].split(",") if lines else []
        if len(fields) >= 3:
            result[path.name[: -len("_mlst.csv")]] = fields[2]
    return result


# --------------------------------------------------------------------------- #
# Main assembly of the reports
# --------------------------------------------------------------------------- #

PHENOTYPE_COLUMNS = [
    "SAMPLE", "MLST", "TYPE", "SUBTYPES", "PETG_STATUS", "PETG_DETAIL", "LPS_GENE_STATES",
    "DETERMINING_GENE", "PREDICTED_PHENOTYPE", "PREDICTED_PHENOTYPE_DESCRIPTION", "PHENOTYPE_BASIS",
]
SUBTYPE_COLUMNS = [
    "SAMPLE", "MLST", "TYPE", "SUBTYPE", "VARTYPE", "ISOLATE_DATABASE", "CHROM", "POS", "REF", "ALT",
    "GENE", "EFFECT", "PETG_STATUS", "NOTE",
]
PETG_COLUMNS = ["SAMPLE", "PETG_STATUS", "PETG_DETAIL"]


def phenotype_for_sample(sample, lps_type, ruleset, reference, variants, petg_status, coverage,
                         min_coverage):
    """Return a dict of phenotype report fields for one sample."""
    result = {"LPS_GENE_STATES": "", "DETERMINING_GENE": "", "PREDICTED_PHENOTYPE": "", "PHENOTYPE_BASIS": ""}
    if lps_type == "untypeable" or not lps_type:
        result["PHENOTYPE_BASIS"] = "untypeable"
        return result, []
    if ruleset is None or reference is None:
        result["PHENOTYPE_BASIS"] = "no phenotype rules for type"
        return result, []

    on_reference = [v for v in variants if v.chrom == reference.record]
    gene_calls = OrderedDict()
    for name in ruleset.genes:
        gene = reference.genes[name]
        overlapping = [v for v in on_reference if v.overlaps(gene.start, gene.end)]
        fraction = gene_coverage(coverage, reference.record, gene) if coverage is not None else None
        gene_calls[name] = classify_gene(gene, overlapping, fraction, min_coverage)

    outcomes = predict(ruleset, gene_calls, petg_status)
    phenotypes = [p for p in outcomes if p]
    result["PREDICTED_PHENOTYPE"] = " | ".join(phenotypes)

    determining = []
    for phenotype in phenotypes:
        for _, _, gene in outcomes[phenotype]:
            if gene and gene not in determining:
                determining.append(gene)
    result["DETERMINING_GENE"] = ";".join(determining) if determining else ("none" if phenotypes else "")
    # Only genes that are not plainly functional; substitution-only genes are left out.
    result["LPS_GENE_STATES"] = ";".join(
        call.describe() for call in gene_calls.values() if call.state != "functional")

    # Basis text
    notes = []
    uncertain = [c for c in gene_calls.values() if c.state == "uncertain"]
    for call in uncertain:
        options = []
        for phenotype, cases in outcomes.items():
            for combo, _petg, _gene in cases:
                if call.gene in combo:
                    options.append((combo[call.gene], phenotype))
        functional = sorted({p for s, p in options if s == "functional" and p})
        inactive = sorted({p for s, p in options if s == "inactivating" and p})
        if functional != inactive:
            notes.append(f"uncertain: {call.gene} {call.reason}; "
                         f"{'/'.join(functional)} if functional, {'/'.join(inactive)} if not")
        else:
            notes.append(f"uncertain: {call.gene} {call.reason} (no effect on phenotype)")
    if petg_status == "uncertain" and ruleset.petg_relevant and len(phenotypes) > 1:
        notes.append("uncertain: petG reading frame altered; both petG alternatives shown")
    if petg_status == "not assessed" and ruleset.petg_relevant and len(phenotypes) > 1:
        notes.append("uncertain: petG not assessed; both petG alternatives shown")
    uncurated_calls = [c for c in gene_calls.values() if c.uncurated]
    for call in uncurated_calls:
        novel = [v.label() for v in call.variants if not v.db]
        if call.state == "uncertain":
            # The reason is already given in the uncertain note
            notes.append(f"uncurated: {call.gene} (not in subtype DB: {', '.join(novel)})")
        else:
            notes.append(f"uncurated: {call.gene} {call.state} by {call.reason} (not in subtype DB: {', '.join(novel)})")
    if not notes:
        notes.append("curated")
    mixed = [v for c in gene_calls.values() if c.state != "functional" for v in c.mixed]
    if mixed:
        notes.append("mixed call: " + ", ".join(
            f"{v.label()}" + (f" AF={format_af(v.af)}" if v.af is not None else "") for v in mixed))
    if not phenotypes:
        notes.insert(0, "no matching phenotype rule")
    result["PHENOTYPE_BASIS"] = "; ".join(notes)
    return result, uncurated_calls


def build_reports(args):
    references = load_references(args.lps_db_dir)
    rulesets = load_rules(args.rules, references)
    db_entries, db_index = load_subtype_db(args.subtype_db)
    lookup = {clean(r.get("PREDICTED_PHENOTYPE")): clean(r.get("PREDICTED_PHENOTYPE_DESCRIPTION"))
              for r in read_tsv(args.phenotype_lookup)} if args.phenotype_lookup else {}

    kaptive = read_kaptive(args.kaptive)
    mlst = read_mlst(args.mlst or [])

    variants = []
    if args.snippy_table:
        variants.extend(read_snippy_variants(args.snippy_table))
    for vcf in args.vcf or []:
        sample = Path(vcf).name.split(args.vcf_suffix)[0] if args.vcf_suffix else None
        variants.extend(read_vcf_variants(vcf, sample_from_file=sample, pass_only=not args.include_nonpass))
    mixed = read_rescued(args.rescued or [])
    for variant in variants:
        variant.db = list(db_index.get(variant.key, []))
        af = mixed.get((variant.sample,) + variant.key)
        if (variant.sample,) + variant.key in mixed:
            variant.mixed = True
            variant.af = af

    sequences = {ref.record: ref.sequence for ref in references.values()}
    for variant in variants:
        assess_confidence(variant, sequences.get(variant.chrom, ""), args.min_confident_af,
                          args.homopolymer_min_run)

    by_sample = OrderedDict()
    for variant in variants:
        by_sample.setdefault(variant.sample, []).append(variant)

    petg = {}
    if args.petg_skipped:
        pass
    else:
        for path in args.petg_blast or []:
            sample = Path(path).name.split("_petG_blast")[0]
            petg[sample] = classify_petg_table(path, args.petg_min_identity, args.petg_min_length)

    coverage_paths = {}
    for path in args.coverage or []:
        sample = Path(path).name.split(args.coverage_suffix)[0]
        coverage_paths[sample] = path

    samples = list(OrderedDict.fromkeys(list(kaptive) + list(by_sample)))
    phenotype_rows, subtype_rows, petg_rows = [], [], []
    for sample in samples:
        sample_variants = by_sample.get(sample, [])
        lps_type = kaptive.get(sample, "")
        if not lps_type:
            types = {e.type for v in sample_variants for e in v.db if e.type}
            lps_type = types.pop() if len(types) == 1 else ""
        if args.petg_skipped:
            petg_status, petg_detail = "not assessed", "petG step skipped"
        else:
            petg_status, petg_detail = petg.get(sample, ("not assessed", "no petG blast result"))
        petg_rows.append({"SAMPLE": sample, "PETG_STATUS": petg_status, "PETG_DETAIL": petg_detail})

        coverage = None
        if sample in coverage_paths:
            coverage = read_coverage(coverage_paths[sample], args.min_depth)
        ruleset = rulesets.get(lps_type)
        reference = references.get(lps_type)
        fields, uncurated_calls = phenotype_for_sample(
            sample, lps_type, ruleset, reference, sample_variants, petg_status, coverage, args.min_coverage)
        if coverage is None and ruleset is not None and args.coverage is not None:
            fields["PHENOTYPE_BASIS"] += "; gene coverage not assessed"

        matched = [(v, e) for v in sample_variants for e in v.db]
        subtypes = []
        for _, entry in matched:
            if entry.subtype not in subtypes:
                subtypes.append(entry.subtype)
        if uncurated_calls:
            subtypes.append("novel")
        if not subtypes and ruleset is not None and reference is not None:
            # Reference-like subtypes (no variant) are not matchable; report none.
            subtypes = []
        descriptions = [lookup.get(p.strip(), "") for p in fields["PREDICTED_PHENOTYPE"].split("|") if p.strip()]
        phenotype_rows.append({
            "SAMPLE": sample,
            "MLST": mlst.get(sample, ""),
            "TYPE": lps_type,
            "SUBTYPES": ";".join(subtypes),
            "PETG_STATUS": petg_status,
            "PETG_DETAIL": petg_detail,
            "PREDICTED_PHENOTYPE_DESCRIPTION": " | ".join(descriptions) if any(descriptions) else "",
            **fields,
        })

        if matched:
            for variant, entry in matched:
                subtype_rows.append({
                    "SAMPLE": sample, "MLST": mlst.get(sample, ""),
                    "TYPE": lps_type if lps_type not in ("", "untypeable") else entry.type,
                    "SUBTYPE": entry.subtype, "VARTYPE": entry.vartype, "ISOLATE_DATABASE": entry.isolate,
                    "CHROM": entry.chrom, "POS": entry.pos, "REF": entry.ref, "ALT": entry.alt,
                    "GENE": entry.gene, "EFFECT": entry.effect, "PETG_STATUS": petg_status, "NOTE": entry.note,
                })
        else:
            subtype_rows.append({
                "SAMPLE": sample, "MLST": mlst.get(sample, ""), "TYPE": lps_type or "untypeable",
                "SUBTYPE": "NA", "VARTYPE": "NA", "ISOLATE_DATABASE": "NA", "CHROM": "NA", "POS": "NA",
                "REF": "NA", "ALT": "NA", "GENE": "NA", "EFFECT": "", "PETG_STATUS": petg_status, "NOTE": "",
            })

    prefix = args.prefix
    write_tsv(f"{prefix}_phenotype_report.tsv", PHENOTYPE_COLUMNS, phenotype_rows)
    write_tsv(f"{prefix}_subtype_report.tsv", SUBTYPE_COLUMNS, subtype_rows)
    write_tsv(f"{prefix}_petG_status.tsv", PETG_COLUMNS, petg_rows)
    return phenotype_rows, subtype_rows, petg_rows


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--lps-db-dir", required=True, help="LPS database directory (reference_LPS.txt, GenBank files)")
    parser.add_argument("--subtype-db", required=True, help="Subtype database (v3 schema with EFFECT)")
    parser.add_argument("--rules", required=True, help="phenotype_rules.tsv")
    parser.add_argument("--phenotype-lookup", help="phenotype_lookup.tsv with descriptions")
    parser.add_argument("--kaptive", required=True, help="Aggregated Kaptive results table")
    parser.add_argument("--mlst", nargs="*", help="Per-sample mlst CSV files")
    parser.add_argument("--snippy-table", help="Aggregated snippy snps table (Illumina)")
    parser.add_argument("--vcf", nargs="*", help="Per-sample snpEff-annotated VCFs (ONT)")
    parser.add_argument("--vcf-suffix", default="_snpeff.vcf", help="Suffix removed from VCF file names to get the sample")
    parser.add_argument("--include-nonpass", action="store_true", help="Also use VCF records whose FILTER is not PASS")
    parser.add_argument("--rescued", nargs="*", help="Per-sample <sample>_snps.rescued.tsv files (Illumina)")
    parser.add_argument("--coverage", nargs="*", help="Per-sample coverage files (snippy aligned FASTA or samtools depth -a)")
    parser.add_argument("--coverage-suffix", default="_snps.aligned.fa", help="Suffix removed from coverage file names")
    parser.add_argument("--min-coverage", type=float, default=DEFAULT_MIN_COVERAGE,
                        help="Minimum covered fraction of a CDS before it is called uncertain (default 0.9)")
    parser.add_argument("--min-depth", type=int, default=1, help="Minimum depth counted as covered for depth tables")
    parser.add_argument("--min-confident-af", type=float, default=MIXED_AF,
                        help="Uncurated calls below this allele fraction are low-confidence (default 0.9)")
    parser.add_argument("--homopolymer-min-run", type=int, default=0,
                        help="Uncurated indels of a base in a reference run of at least this length are "
                             "low-confidence (0 = off; the ONT pipeline uses 3)")
    parser.add_argument("--petg-blast", nargs="*", help="Per-sample <sample>_petG_blast.tsv (with qlen qseq sseq)")
    parser.add_argument("--petg-skipped", action="store_true", help="petG step was skipped")
    parser.add_argument("--petg-min-identity", type=float, default=95.0)
    parser.add_argument("--petg-min-length", type=int, default=1570)
    parser.add_argument("--prefix", required=True, help="Output prefix, e.g. 10_Illumina")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        build_reports(args)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
