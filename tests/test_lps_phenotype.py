#!/usr/bin/env python3
"""Tests for bin/lps_phenotype.py (stdlib unittest; also runs under pytest).

Run: python3 -m unittest discover -s tests -p 'test_*.py'
"""

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DB = REPO / "databases" / "LPS"
spec = importlib.util.spec_from_file_location("lps_phenotype", REPO / "bin" / "lps_phenotype.py")
lp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lp)

REFERENCES = lp.load_references(DB)
RULES = lp.load_rules(DB / "phenotype_rules.tsv", REFERENCES)
DB_ENTRIES, DB_INDEX = lp.load_subtype_db(DB / "LPS_subtype_database_v3.txt")


def db_variant(subtype, sample="S"):
    """A called variant identical to a subtype database row.

    Nonsense rows get the stop_gained term that snpEff adds in the pipelines.
    """
    entry = next(e for e in DB_ENTRIES if e.subtype == subtype)
    effects = {"stop_gained"} if "*" in entry.vartype else set()
    variant = lp.Variant(sample, entry.chrom, entry.pos, entry.ref, entry.alt, effects)
    variant.db = list(DB_INDEX[variant.key])
    return variant


def novel(chrom, pos, ref, alt, effects=(), sample="S"):
    return lp.Variant(sample, chrom, pos, ref, alt, set(effects))


def full_coverage(lps_type):
    reference = REFERENCES[lps_type]
    return {reference.record: set(range(1, len(reference.sequence) + 1))}


def call(lps_type, variants, petg="functional", coverage="full", called=True):
    """Phenotype fields for one sample; every position is covered unless a coverage map is given."""
    if isinstance(coverage, str):
        coverage = full_coverage(lps_type) if lps_type in REFERENCES else None
    fields, _ = lp.phenotype_for_sample(
        "S", lps_type, RULES.get(lps_type), REFERENCES.get(lps_type), variants, petg, coverage, 0.9, called)
    return fields


L1 = "L1_X73_H1"
L3 = "L3_NC002663_LPS"


class WorkedExamples(unittest.TestCase):
    def test_a_pcga_deletion_petg_absent_is_p4(self):
        fields = call("L1", [db_variant("L1-2")], petg="absent")
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L1_P4")
        self.assertEqual(fields["DETERMINING_GENE"], "pcgA")
        self.assertEqual(fields["PHENOTYPE_BASIS"], "curated")

    def test_a_pcga_deletion_petg_functional_is_p3(self):
        self.assertEqual(call("L1", [db_variant("L1-2")])["PREDICTED_PHENOTYPE"], "L1_P3")

    def test_b_earliest_step_wins_l3(self):
        fields = call("L3", [db_variant("L3-30"), db_variant("L3-46")])
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L3_P4")
        self.assertEqual(fields["DETERMINING_GENE"], "gatG")

    def test_c_novel_gata_stop_is_p5_uncurated(self):
        gata = REFERENCES["L1"].genes["gatA"]
        variant = novel(L1, gata.start + 300, "C", "T", {"stop_gained"})
        fields = call("L1", [variant], petg="nonfunctional")
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L1_P5")
        self.assertTrue(fields["PHENOTYPE_BASIS"].startswith("uncurated: gatA inactivating"))

    def test_d_no_lps_variants_petg_nonfunctional_is_p2(self):
        fields = call("L1", [], petg="nonfunctional")
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L1_P2")
        self.assertEqual(fields["DETERMINING_GENE"], "none")

    def test_wild_type_petg_functional_is_p1(self):
        self.assertEqual(call("L1", [])["PREDICTED_PHENOTYPE"], "L1_P1")

    def test_hpte_ignores_petg(self):
        for petg in ("functional", "absent", "uncertain"):
            self.assertEqual(call("L1", [db_variant("L1-7")], petg=petg)["PREDICTED_PHENOTYPE"], "L1_P6")


class CompensatingFrameshifts(unittest.TestCase):
    def test_l1_3_plus_l1_4_is_uncertain(self):
        fields = call("L1", [db_variant("L1-3"), db_variant("L1-4")])
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L1_P1 | L1_P3")
        self.assertIn("pcgB:uncertain", fields["LPS_GENE_STATES"])
        self.assertIn("compensating frameshifts, net +9 bp (L1-3+L1-4)", fields["PHENOTYPE_BASIS"])
        self.assertIn("L1_P1 if functional, L1_P3 if not", fields["PHENOTYPE_BASIS"])

    def test_l1_3_plus_l1_4_petg_absent(self):
        fields = call("L1", [db_variant("L1-3"), db_variant("L1-4")], petg="absent")
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L1_P2 | L1_P4")

    def test_merged_complex_record_behaves_the_same(self):
        # snippy/freebayes may merge L1-3 and L1-4 into one +9 bp record
        ref = REFERENCES["L1"].sequence[2542:2548]  # positions 2543..2548
        alt = ref[0] + "AAAATAA" + ref[1:5] + "A" + "AT"[1:] + ref[5:]
        variant = novel(L1, 2543, ref, alt)
        self.assertEqual(variant.length_change, 9)
        fields = call("L1", [variant])
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L1_P1 | L1_P3")

    def test_each_alone_is_inactivating(self):
        for subtype in ("L1-3", "L1-4"):
            self.assertEqual(call("L1", [db_variant(subtype)])["PREDICTED_PHENOTYPE"], "L1_P3")

    def test_curated_none_resolves_uncertain(self):
        a, b = db_variant("L1-3"), db_variant("L1-4")
        for v in (a, b):
            for entry in v.db:
                entry_copy = lp.DbVariant({})
                entry_copy.__dict__.update(entry.__dict__)
                entry_copy.effect = "none"
                v.db = [entry_copy]
        fields = call("L1", [a, b], petg="absent")
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L1_P2")
        self.assertEqual(fields["PHENOTYPE_BASIS"], "curated")


class GeneStateRules(unittest.TestCase):
    def test_curated_inactivating_missense(self):
        variant = db_variant("L1-2")
        gene = REFERENCES["L1"].genes["gatA"]
        missense = novel(L1, gene.start + 30, "A", "G", {"missense_variant"})
        entry = lp.DbVariant({"TYPE": "L1", "SUBTYPE": "L1-x", "CHROM": L1, "POS": str(missense.pos),
                              "REF": "A", "ALT": "G", "GENE": "gatA", "EFFECT": "inactivating"})
        missense.db = [entry]
        fields = call("L1", [variant, missense])
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L1_P5")

    def test_plain_substitutions_are_functional(self):
        gene = REFERENCES["L1"].genes["gatA"]
        fields = call("L1", [novel(L1, gene.start + 10, "A", "G", {"missense_variant"})])
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L1_P1")
        self.assertEqual(fields["SUBTYPE_STATUS"], "uncurated_variants")
        self.assertEqual(fields["LPS_GENE_STATES"], "")
        self.assertIn("uncurated variants: A>G at", fields["PHENOTYPE_BASIS"])
        self.assertIn("(gatA)", fields["PHENOTYPE_BASIS"])

    def test_in_frame_deletion_is_uncertain(self):
        gene = REFERENCES["L3"].genes["gatG"]
        pos = gene.start + 100
        ref = REFERENCES["L3"].sequence[pos - 1:pos + 3]
        fields = call("L3", [novel(L3, pos, ref, ref[0])])
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "Not predicted")
        self.assertIn("the rules alone would give L3_P1 | L3_P4", fields["PHENOTYPE_BASIS"])
        self.assertIn("in-frame indel", fields["PHENOTYPE_BASIS"])

    def test_l3_29_large_deletion_hits_gatf(self):
        fields = call("L3", [db_variant("L3-29")])
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L3_P5")
        self.assertEqual(fields["DETERMINING_GENE"], "gatF")

    def test_intergenic_indel_ignored(self):
        fields = call("L1", [novel(L1, 10, "A", "AT")])
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L1_P1")
        self.assertIn("(intergenic)", fields["PHENOTYPE_BASIS"])

    def test_low_coverage_gene_is_uncertain(self):
        reference = REFERENCES["L1"]
        gene = reference.genes["gatA"]
        covered = {reference.record: set(range(1, len(reference.sequence) + 1)) - set(range(gene.start, gene.end))}
        fields = call("L1", [db_variant("L1-2")], coverage=covered)
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L1_P3 | L1_P5")
        self.assertIn("gatA only 0% covered", fields["PHENOTYPE_BASIS"])

    def test_requires_row(self):
        ruleset = lp.RuleSet("LX")
        ruleset.rules = [
            lp.Rule("LX", 1, ["hptE"], "any", [], "LX_P9"),
            lp.Rule("LX", None, [], "any", [("gctD", "modifying")], "LX_P2"),
            lp.Rule("LX", None, [], "any", [], "LX_P1"),
        ]
        self.assertEqual(lp.evaluate(ruleset, {"gctD": "modifying"}, "functional")[0], "LX_P2")
        self.assertEqual(lp.evaluate(ruleset, {"gctD": "functional"}, "functional")[0], "LX_P1")
        self.assertEqual(lp.evaluate(ruleset, {"hptE": "inactivating", "gctD": "modifying"}, "functional")[0], "LX_P9")

    def test_type_without_rules(self):
        fields = call("L4", [])
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "Not predicted")
        self.assertEqual(fields["PHENOTYPE_BASIS"], "no phenotype rules for type")
        self.assertEqual(fields["SUBTYPE_STATUS"], "reference_identical")

    def test_untypeable(self):
        self.assertEqual(call("untypeable", [])["PHENOTYPE_BASIS"], "untypeable")

    def test_petg_not_assessed_gives_both(self):
        fields = call("L1", [], petg="not assessed")
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L1_P1 | L1_P2")
        self.assertIn("petG not assessed", fields["PHENOTYPE_BASIS"])


class VcfInput(unittest.TestCase):
    def test_pass_only_and_mixed_af(self):
        vcf = "\n".join([
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
            f"{L3}\t2041\t.\tC\tA\t24.5\tPASS\tF;ANN=A|stop_gained|HIGH|gatG|x\tGT:GQ:DP:AD:AF\t1:24:38:14,22:0.5789",
            f"{L3}\t5349\t.\tT\tA\t0\tLowQual\tF;ANN=A|stop_gained|HIGH|hptE|x\tGT:GQ:DP:AD:AF\t1:0:10:4,6:0.6",
            f"{L3}\t100\t.\tA\tG\t30\tPASS\tP\tGT:GQ:DP:AD:AF\t0:30:40:40,0:0.0",
        ]) + "\n"
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "S1_clair3.snpeff.vcf"
            path.write_text(vcf)
            variants = lp.read_vcf_variants(path, sample_from_file="S1")
        self.assertEqual([v.pos for v in variants], [2041])
        self.assertEqual(variants[0].effects, {"stop_gained"})
        self.assertTrue(variants[0].mixed)
        # Uncurated stop at AF 0.58 is low-confidence: gatG becomes uncertain, not inactivating
        lp.assess_confidence(variants[0], REFERENCES["L3"].sequence)
        fields = call("L3", variants)
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "Not predicted")
        self.assertIn("the rules alone would give L3_P1 | L3_P4", fields["PHENOTYPE_BASIS"])
        self.assertIn("gatG depends on a low-confidence call (C>A at 2041: AF=0.579)", fields["PHENOTYPE_BASIS"])


class ConfidenceGate(unittest.TestCase):
    """Option 2: uncurated low-confidence calls can only widen a call to uncertain."""

    def hpte_pair(self, ins_af, del_af):
        ins = novel(L1, 5155, "C", "CG")
        ins.af = ins_af
        dele = novel(L1, 5166, "GT", "G")
        dele.af = del_af
        return ins, dele

    def assess(self, variants, homopolymer_min_run=3):
        for v in variants:
            lp.assess_confidence(v, REFERENCES["L1"].sequence, 0.9, homopolymer_min_run)

    def test_homopolymer_run(self):
        ins, dele = self.hpte_pair(0.63, 0.92)
        seq = REFERENCES["L1"].sequence
        self.assertEqual(lp.homopolymer_run(seq, ins), 3)
        self.assertEqual(lp.homopolymer_run(seq, dele), 1)
        # L1-4 inserts "AT" next to an A run: not a single-base run length change
        self.assertEqual(lp.homopolymer_run(seq, db_variant("L1-4")), 0)
        # a 1 bp A insertion at the same place would extend the AAAAAA run
        self.assertEqual(lp.homopolymer_run(seq, novel(L1, 2548, "A", "AA")), 6)

    def test_pm3289_hpte_stays_uncertain_with_reason(self):
        ins, dele = self.hpte_pair(0.63, 0.92)
        self.assess([ins, dele])
        self.assertIn("AF=0.630", ins.low_confidence)
        self.assertIn("extends a GGG run", ins.low_confidence)
        self.assertIsNone(dele.low_confidence)
        fields = call("L1", [ins, dele, db_variant("L1-3"), db_variant("L1-4")])
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L1_P1 | L1_P3 | L1_P6")
        self.assertIn("uncertain: hptE depends on a low-confidence call (1bp ins at 5155: AF=0.630, extends a GGG run)",
                      fields["PHENOTYPE_BASIS"])
        self.assertIn("uncertain: pcgB compensating frameshifts, net +9 bp (L1-3+L1-4)", fields["PHENOTYPE_BASIS"])
        self.assertIn("uncurated: hptE (not in subtype DB: 1bp ins at 5155, 1bp del at 5166)", fields["PHENOTYPE_BASIS"])

    def test_weak_substitutions_are_not_listed_and_af_truncated(self):
        ins, dele = self.hpte_pair(0.63, 0.92)
        snv = novel(L1, 5159, "C", "G", {"missense_variant"})
        snv.af = 0.8988
        self.assess([ins, dele, snv])
        self.assertEqual(snv.low_confidence, "AF=0.899")
        fields = call("L1", [ins, dele, snv])
        self.assertNotIn("5159", fields["PHENOTYPE_BASIS"].split("uncurated")[0])

    def test_weak_frameshift_alone_is_uncertain_not_inactivating(self):
        ins, _ = self.hpte_pair(0.63, 0.92)
        self.assess([ins])
        self.assertEqual(call("L1", [ins])["PREDICTED_PHENOTYPE"], "L1_P1 | L1_P6")

    def test_confident_novel_stop_still_drives(self):
        stop = novel(L3, 5349, "T", "A", {"stop_gained"})
        stop.af = 0.98
        lp.assess_confidence(stop, REFERENCES["L3"].sequence, 0.9, 3)
        fields = call("L3", [stop])
        # L3 gives no phenotype without a database subtype, but the basis still says what the rules see
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "Not predicted")
        self.assertIn("the rules alone would give L3_P7", fields["PHENOTYPE_BASIS"])
        self.assertIn("uncurated: hptE inactivating", fields["PHENOTYPE_BASIS"])
        # L1 has no such restriction
        hpte = REFERENCES["L1"].genes["hptE"]
        l1_stop = novel(L1, hpte.start + 300, "C", "T", {"stop_gained"})
        self.assertEqual(call("L1", [l1_stop])["PREDICTED_PHENOTYPE"], "L1_P6")

    def test_curated_variant_in_homopolymer_counts(self):
        variants = [db_variant("L1-4")]
        variants[0].af = 0.7
        variants[0].mixed = True
        self.assess(variants)
        self.assertIsNone(variants[0].low_confidence)
        fields = call("L1", variants)
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L1_P3")
        self.assertIn("mixed call: 2bp ins at 2548 AF=0.700", fields["PHENOTYPE_BASIS"])

    def test_illumina_has_no_homopolymer_check(self):
        ins, dele = self.hpte_pair(None, None)
        self.assess([ins, dele], homopolymer_min_run=0)
        self.assertIsNone(ins.low_confidence)
        self.assertEqual(call("L1", [ins, dele])["PREDICTED_PHENOTYPE"], "L1_P1 | L1_P6")


def non_rule_gene(lps_type):
    """A CDS of the reference that no phenotype rule uses."""
    return next(g for name, g in REFERENCES[lps_type].genes.items() if name not in RULES[lps_type].genes)


class SubtypeStatus(unittest.TestCase):
    """SUBTYPE_STATUS and the L3 rule that a phenotype needs a subtype in the database."""

    def test_l3_without_variants_is_reference_identical(self):
        fields = call("L3", [])
        self.assertEqual(fields["SUBTYPE_STATUS"], "reference_identical")
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L3_P1")
        self.assertEqual(fields["PHENOTYPE_BASIS"], "no variants in the locus (identical to the reference)")

    def test_l3_with_unlisted_variant_is_not_predicted(self):
        gene = REFERENCES["L3"].genes["natC"]
        fields = call("L3", [novel(L3, gene.start + 10, "A", "G", {"missense_variant"})])
        self.assertEqual(fields["SUBTYPE_STATUS"], "uncurated_variants")
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "Not predicted")
        self.assertTrue(fields["PHENOTYPE_BASIS"].startswith(
            "not predicted: no subtype in the database matches this isolate; the rules alone would give L3_P1"))
        self.assertIn("uncurated variants: A>G at", fields["PHENOTYPE_BASIS"])
        self.assertIn("(natC)", fields["PHENOTYPE_BASIS"])

    def test_l3_matched_subtype_is_predicted(self):
        fields = call("L3", [db_variant("L3-6")])
        self.assertEqual(fields["SUBTYPE_STATUS"], "assigned")
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L3_P2")
        self.assertEqual(fields["PHENOTYPE_BASIS"], "curated")

    def test_matched_subtype_plus_uncurated_variant_that_changes_the_result(self):
        # L3-6 (natC, step 6) alone is L3_P2, but a new hptE stop (step 1) would make it L3_P7
        stop = novel(L3, 5349, "T", "A", {"stop_gained"})
        fields = call("L3", [db_variant("L3-6"), stop])
        self.assertEqual(fields["SUBTYPE_STATUS"], "assigned")
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "Not predicted")
        self.assertEqual(fields["DETERMINING_GENE"], "")
        self.assertIn("an uncurated variant changes the phenotype; the rules alone would give L3_P7",
                      fields["PHENOTYPE_BASIS"])
        self.assertIn("natC:inactivating", fields["LPS_GENE_STATES"])
        self.assertIn("hptE:inactivating", fields["LPS_GENE_STATES"])

    def test_matched_subtype_plus_uncurated_variant_that_changes_nothing(self):
        # L3-12 (hptE, step 1) wins whatever happens at step 6, so L3_P7 is safe to report
        gene = REFERENCES["L3"].genes["natC"]
        stop = novel(L3, gene.start + 300, "C", "A", {"stop_gained"})
        fields = call("L3", [db_variant("L3-12"), stop])
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L3_P7")
        self.assertEqual(fields["DETERMINING_GENE"], "hptE")
        self.assertIn("uncurated: natC inactivating", fields["PHENOTYPE_BASIS"])

    def test_l1_keeps_its_rules_without_a_subtype(self):
        gene = REFERENCES["L1"].genes["hptE"]
        fields = call("L1", [novel(L1, gene.start + 300, "C", "T", {"stop_gained"})])
        self.assertEqual(fields["SUBTYPE_STATUS"], "uncurated_variants")
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L1_P6")

    def test_calling_that_did_not_run_is_not_assessed(self):
        for lps_type in ("L1", "L3", "L4"):
            fields = call(lps_type, [], called=False)
            self.assertEqual(fields["SUBTYPE_STATUS"], "not_assessed")
            self.assertEqual(fields["PREDICTED_PHENOTYPE"], "Not predicted")
            self.assertTrue(fields["PHENOTYPE_BASIS"].startswith("not assessed: no variant calling output"))

    def test_no_variants_needs_coverage_data(self):
        fields = call("L3", [], coverage=None)
        self.assertEqual(fields["SUBTYPE_STATUS"], "not_assessed")
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "Not predicted")
        self.assertIn("no coverage data", fields["PHENOTYPE_BASIS"])

    def test_no_variants_but_a_gap_anywhere_in_the_locus_is_not_assessed(self):
        # L3 has only rule genes; L1 also has genes that no rule uses (a gap there must count too)
        cases = [("L3", REFERENCES["L3"].genes["natB"]), ("L1", non_rule_gene("L1"))]
        for lps_type, gene in cases:
            with self.subTest(type=lps_type, gene=gene.name):
                reference = REFERENCES[lps_type]
                covered = {reference.record: set(range(1, len(reference.sequence) + 1)) - set(range(gene.start, gene.end + 1))}
                fields = call(lps_type, [], coverage=covered)
                self.assertEqual(fields["SUBTYPE_STATUS"], "not_assessed")
                self.assertEqual(fields["PREDICTED_PHENOTYPE"], "Not predicted")
                self.assertIn(f"{gene.name} 0%", fields["PHENOTYPE_BASIS"])

    def test_untypeable_is_not_applicable(self):
        fields = call("untypeable", [])
        self.assertEqual(fields["SUBTYPE_STATUS"], "not_applicable")
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "Not predicted")

    def test_reference_subtypes_are_marked_in_the_database(self):
        marked = {e.type: e.subtype for e in DB_ENTRIES if e.is_reference}
        self.assertEqual(marked["L1"], "L1-1")
        self.assertEqual(marked["L3"], "L3-1")
        # L3-3 and L3-4 have no coordinates either, but they are not the reference sequence
        self.assertFalse(next(e for e in DB_ENTRIES if e.subtype == "L3-3").is_reference)
        self.assertFalse(next(e for e in DB_ENTRIES if e.subtype == "L3-4").is_reference)
        self.assertFalse(any(e.matchable for e in DB_ENTRIES if e.is_reference))

    def test_l1_6_is_labelled_as_the_deletion_it_is(self):
        entry = next(e for e in DB_ENTRIES if e.subtype == "L1-6")
        self.assertEqual((entry.ref, entry.alt, entry.vartype), ("TA", "T", "1bp deletion"))

    def test_vartype_matches_allele_lengths(self):
        for entry in DB_ENTRIES:
            if not entry.matchable:
                continue
            change = len(entry.alt) - len(entry.ref)
            with self.subTest(subtype=entry.subtype):
                if "insertion" in entry.vartype:
                    self.assertGreater(change, 0)
                elif "deletion" in entry.vartype:
                    self.assertLess(change, 0)
                elif entry.vartype.startswith("snp"):
                    self.assertEqual(change, 0)

    def test_l3_6_with_l3_12_is_p7_by_the_step_rule(self):
        # The old v2 override "L3-6 + L3-32 = L3_P7" is now just the lowest-step rule
        fields = call("L3", [db_variant("L3-6"), db_variant("L3-12")])
        self.assertEqual(fields["PREDICTED_PHENOTYPE"], "L3_P7")

    def test_two_reference_subtypes_for_one_type_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "db.txt"
            path.write_text("TYPE\tSUBTYPE\tISOLATE\tCHROM\tPOS\tVARTYPE\tREF\tALT\tGENE\tEFFECT\tNOTE\n"
                            "L3\tL3-1\tA\tc\tNA\treference\tNA\tNA\tNA\t\t\n"
                            "L3\tL3-2\tB\tc\tNA\treference\tNA\tNA\tNA\t\t\n")
            with self.assertRaisesRegex(ValueError, "More than one reference subtype"):
                lp.load_subtype_db(path)


class ReportFiles(unittest.TestCase):
    """The whole script on small input files."""

    def run_reports(self, tmp, called_samples):
        tmp = Path(tmp)
        reference = REFERENCES["L3"]
        kaptive = tmp / "kaptive.tsv"
        header = "sampleID\tBest match locus\tMatch confidence"
        kaptive.write_text("\n".join([
            header,
            "ASSIGNED\tL3-x\tTypeable",
            "REFERENCE\tL3-x\tTypeable",
            "NOVEL\tL3-x\tTypeable",
            "FAILED\tL3-x\tTypeable",
            "UNTYPED\tnone\tNot typeable",
        ]) + "\n")
        entry = next(e for e in DB_ENTRIES if e.subtype == "L3-6")
        gene = reference.genes["natC"]
        snippy = tmp / "snps.tsv"
        snippy.write_text("\n".join([
            "sampleID\tCHROM\tPOS\tREF\tALT\tEFFECT",
            f"ASSIGNED\t{entry.chrom}\t{entry.pos}\t{entry.ref}\t{entry.alt}\tstop_gained",
            f"NOVEL\t{reference.record}\t{gene.start + 10}\tA\tG\tmissense_variant",
        ]) + "\n")
        depth_files = []
        for sample in ("ASSIGNED", "REFERENCE", "NOVEL", "FAILED"):
            path = tmp / f"{sample}_depth.tsv"
            path.write_text("".join(f"{reference.record}\t{i}\t30\n" for i in range(1, len(reference.sequence) + 1)))
            depth_files.append(str(path))
        called = [f for f in depth_files if Path(f).name.split("_depth.tsv")[0] in called_samples]
        prefix = str(tmp / "10_T")
        code = lp.main([
            "--lps-db-dir", str(DB), "--subtype-db", str(DB / "LPS_subtype_database_v3.txt"),
            "--rules", str(DB / "phenotype_rules.tsv"), "--phenotype-lookup", str(DB / "phenotype_lookup.tsv"),
            "--kaptive", str(kaptive), "--snippy-table", str(snippy),
            "--coverage", *depth_files, "--coverage-suffix", "_depth.tsv", "--min-depth", "10",
            "--called", *called, "--called-suffix", "_depth.tsv", "--petg-skipped", "--prefix", prefix,
        ])
        self.assertEqual(code, 0)
        def read(name):
            with open(f"{prefix}_{name}.tsv", newline="") as handle:
                return list(csv.DictReader(handle, delimiter="\t"))
        return {r["SAMPLE"]: r for r in read("phenotype_report")}, read("subtype_report")

    def test_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            pheno, subtype_rows = self.run_reports(tmp, {"ASSIGNED", "REFERENCE", "NOVEL"})
        expected = {
            "ASSIGNED": ("L3-6", "assigned", "L3_P2"),
            "REFERENCE": ("L3-1", "reference_identical", "L3_P1"),
            "NOVEL": ("Not assigned", "uncurated_variants", "Not predicted"),
            "FAILED": ("Not assigned", "not_assessed", "Not predicted"),
            "UNTYPED": ("Not assigned", "not_applicable", "Not predicted"),
        }
        for sample, (subtypes, status, phenotype) in expected.items():
            with self.subTest(sample=sample):
                self.assertEqual(pheno[sample]["SUBTYPES"], subtypes)
                self.assertEqual(pheno[sample]["SUBTYPE_STATUS"], status)
                self.assertEqual(pheno[sample]["PREDICTED_PHENOTYPE"], phenotype)
                rows = [r for r in subtype_rows if r["SAMPLE"] == sample]
                self.assertEqual({r["SUBTYPE"] for r in rows}, {subtypes})
                self.assertEqual({r["SUBTYPE_STATUS"] for r in rows}, {status})
        self.assertEqual(pheno["NOVEL"]["PREDICTED_PHENOTYPE_DESCRIPTION"], "Not predicted")
        reference_row = next(r for r in subtype_rows if r["SAMPLE"] == "REFERENCE")
        self.assertEqual((reference_row["VARTYPE"], reference_row["ISOLATE_DATABASE"]), ("reference", "Pm70"))
        untyped_row = next(r for r in subtype_rows if r["SAMPLE"] == "UNTYPED")
        self.assertEqual(untyped_row["TYPE"], "untypeable")

    def test_sample_without_calling_output_is_never_reference_identical(self):
        with tempfile.TemporaryDirectory() as tmp:
            pheno, _ = self.run_reports(tmp, {"ASSIGNED", "NOVEL"})
        # REFERENCE has a depth file but no calling output, exactly like a failed calling run
        self.assertEqual(pheno["REFERENCE"]["SUBTYPE_STATUS"], "not_assessed")
        self.assertEqual(pheno["REFERENCE"]["PREDICTED_PHENOTYPE"], "Not predicted")


class RulesValidation(unittest.TestCase):
    def write_rules(self, text):
        handle = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
        handle.write("TYPE\tSTEP\tGENES\tPETG\tREQUIRES\tPHENOTYPE\n" + text)
        handle.close()
        return handle.name

    def test_unknown_gene_rejected(self):
        path = self.write_rules("L1\t1\tnotAGene\tany\t-\tL1_P6\nL1\tfull\t-\tany\t-\tL1_P1\n")
        with self.assertRaisesRegex(ValueError, "notAGene is not a CDS"):
            lp.load_rules(path, REFERENCES)

    def test_missing_full_rejected(self):
        path = self.write_rules("L1\t1\thptE\tany\t-\tL1_P6\n")
        with self.assertRaisesRegex(ValueError, "no 'full' rule row"):
            lp.load_rules(path, REFERENCES)

    def test_curated_only_is_set_for_l3_only(self):
        self.assertTrue(RULES["L3"].curated_only)
        self.assertFalse(RULES["L1"].curated_only)

    def test_curated_only_must_be_the_same_for_a_type(self):
        handle = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
        handle.write("TYPE\tSTEP\tGENES\tPETG\tREQUIRES\tPHENOTYPE\tCURATED_ONLY\n"
                     "L1\t1\thptE\tany\t-\tL1_P6\tyes\nL1\tfull\t-\tany\t-\tL1_P1\tno\n")
        handle.close()
        with self.assertRaisesRegex(ValueError, "CURATED_ONLY must be the same"):
            lp.load_rules(handle.name, REFERENCES)

    def test_db_effect_values(self):
        allowed = {"", "none", "modifying", "inactivating"}
        self.assertTrue(all(e.effect in allowed for e in DB_ENTRIES))

    def test_l3_32_blank_and_unmatchable(self):
        entry = next(e for e in DB_ENTRIES if e.subtype == "L3-32")
        self.assertFalse(entry.matchable)
        self.assertEqual([e.subtype for e in DB_INDEX[(L3, "5983", "C", "T")]], ["L3-12"])


class ReferenceIntegrity(unittest.TestCase):
    CODON = dict(zip(
        [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"],
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))

    def test_rule_genes_are_intact_orfs(self):
        complement = str.maketrans("ACGT", "TGCA")
        for lps_type, ruleset in RULES.items():
            reference = REFERENCES[lps_type]
            for name in ruleset.genes:
                gene = reference.genes[name]
                seq = reference.sequence[gene.start - 1:gene.end]
                if gene.strand == "-":
                    seq = seq.translate(complement)[::-1]
                with self.subTest(type=lps_type, gene=name):
                    self.assertEqual(len(seq) % 3, 0)
                    self.assertIn(seq[:3], lp.START_CODONS)
                    protein = "".join(self.CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq), 3))
                    self.assertNotIn("*", protein[:-1])
                    self.assertTrue(protein.endswith("*"))


class PetgAlignment(unittest.TestCase):
    REF = (DB / "petG_X73_NZ_CM001580.fasta").read_text().split("\n", 1)[1].replace("\n", "").upper()

    def classify(self, qseq, sseq, qstart=1):
        return lp.classify_petg_alignment(qseq, sseq, qstart, len(self.REF))

    def test_intact(self):
        self.assertEqual(self.classify(self.REF, self.REF)[0], "functional")

    def test_one_bp_deletion(self):
        sseq = self.REF[:400] + "-" + self.REF[401:]
        status, detail = self.classify(self.REF, sseq)
        self.assertEqual(status, "nonfunctional")
        self.assertIn("frameshift", detail)

    def test_one_bp_insertion(self):
        qseq = self.REF[:400] + "-" + self.REF[400:]
        sseq = self.REF[:400] + "A" + self.REF[400:]
        self.assertEqual(self.classify(qseq, sseq)[0], "nonfunctional")

    def test_nonsense(self):
        sseq = self.REF[:300] + "TAA" + self.REF[303:]
        status, detail = self.classify(self.REF, sseq)
        self.assertEqual(status, "nonfunctional")
        self.assertIn("premature stop TAA at codon 101", detail)

    def test_start_lost(self):
        sseq = "ATA" + self.REF[3:]
        self.assertEqual(self.classify(self.REF, sseq)[0], "nonfunctional")

    def test_compensating_frameshifts(self):
        qseq = self.REF[:400] + "-" + self.REF[400:]
        sseq = self.REF[:400] + "A" + self.REF[400:450] + "-" + self.REF[451:]
        status, detail = self.classify(qseq, sseq)
        self.assertEqual(status, "uncertain")
        self.assertIn("compensating", detail)

    def test_blast_table_gates(self):
        ref = self.REF
        qlen = len(ref)
        rows = [
            # split across contigs: two partial hits, neither covers the gene
            ["petG", "c1", "100.0", "900", "0", "0", "1", "900", "1", "900", "0", "1600", str(qlen), ref[:900], ref[:900]],
            ["petG", "c2", "100.0", "681", "0", "0", "901", str(qlen), "1", "681", "0", "1200", str(qlen), ref[900:], ref[900:]],
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "S_petG_blast.tsv"
            path.write_text("\n".join("\t".join(r) for r in rows) + "\n")
            self.assertEqual(lp.classify_petg_table(path, 95, 1570)[0], "absent")
            # reverse-strand full hit (sseq is reported in query orientation by blastn)
            full = ["petG", "c3", "99.9", str(qlen), "1", "0", "1", str(qlen), str(qlen + 10), "11", "0", "2900",
                    str(qlen), ref, ref]
            path.write_text("\t".join(full) + "\n")
            status, detail = lp.classify_petg_table(path, 95, 1570)
            self.assertEqual(status, "functional")
            self.assertIn("c3:11-", detail)


class SharedCopy(unittest.TestCase):
    def test_identical_in_ont_repo(self):
        other = REPO.parent / ("LPS_typing_ONT" if REPO.name == "LPS_typing_Illumina" else "LPS_typing_Illumina")
        mirror = other / "bin" / "lps_phenotype.py"
        if not mirror.is_file():
            self.skipTest("sibling pipeline repository not checked out next to this one")
        self.assertEqual(mirror.read_bytes(), (REPO / "bin" / "lps_phenotype.py").read_bytes())
        for name in ("LPS_subtype_database_v3.txt", "phenotype_rules.tsv", "phenotype_lookup.tsv"):
            with self.subTest(file=name):
                self.assertEqual((other / "databases" / "LPS" / name).read_bytes(), (DB / name).read_bytes())


if __name__ == "__main__":
    unittest.main()
