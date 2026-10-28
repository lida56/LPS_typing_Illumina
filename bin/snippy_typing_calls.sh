#!/usr/bin/env bash
#
# Rebuild snippy's typing call set from snps.raw.vcf, keeping 0/1 calls whose
# alt fraction is a majority.
#
# Snippy 4.6.0 runs freebayes in diploid mode (-p 2) and keeps only GT=1/1, so a
# real haploid allele with a minority of reference reads (e.g. a repeat-length
# indel at AF 0.78, GT 0/1) never reaches snps.tab (tseemann/snippy#343). This script
# re-applies snippy's own rule, adds 0/1 records with AF >= --het-min-af, and
# then runs the same normalise / annotate / tab steps snippy does, so every
# record snippy keeps comes out byte-identical.
#
# Runs inside the snippy container, from snippy's --outdir.
#
# Usage:
#   snippy_typing_calls.sh --raw snps.raw.vcf --filt snps.filt.vcf --refdir reference \
#       --minqual 100 --mincov 10 --minfrac 0 --het-min-af 0.5 [--prefix snps.typing]
#
# Outputs:
#   <prefix>.filt.vcf  selected, normalised, tag-stripped (cf. snps.filt.vcf)
#   <prefix>.vcf       snpEff-annotated (cf. snps.vcf)
#   <prefix>.tab       snippy table (cf. snps.tab)
#   snps.rescued.tsv   records in <prefix>.filt.vcf but not in --filt, with AF

set -euo pipefail

usage() { sed -n '15,17p' "$0" >&2; exit 2; }

raw="" filt="" refdir="" minqual="" mincov="" minfrac="" het_min_af=""
prefix="snps.typing"
audit="snps.rescued.tsv"
while [[ $# -gt 0 ]]; do
	case "$1" in
		--raw) raw="$2"; shift 2 ;;
		--filt) filt="$2"; shift 2 ;;
		--refdir) refdir="$2"; shift 2 ;;
		--minqual) minqual="$2"; shift 2 ;;
		--mincov) mincov="$2"; shift 2 ;;
		--minfrac) minfrac="$2"; shift 2 ;;
		--het-min-af) het_min_af="$2"; shift 2 ;;
		--prefix) prefix="$2"; shift 2 ;;
		*) echo "Unknown argument: $1" >&2; usage ;;
	esac
done
for v in raw filt refdir minqual mincov minfrac het_min_af; do
	[[ -n "${!v}" ]] || { echo "Missing --${v//_/-}" >&2; usage; }
done
for f in "$raw" "$filt" "$refdir/ref.fa" "$refdir/ref.gff"; do
	[[ -s "$f" ]] || { echo "Missing or empty input: $f" >&2; exit 1; }
done

# Tags snippy keeps (snippy 4.6.0, $KEEP_VCF_TAGS).
KEEP_VCF_TAGS='^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL'

# Select on the raw VCF, so this stage does not depend on representation.
# First clause is snippy's $bcf_filter verbatim; the second admits biallelic 0/1
# calls where the alt allele is the majority, i.e. the haploid isolate's genotype.
select_calls() {
	local snippy_rule="FMT/GT=\"1/1\" && QUAL>=$minqual && FMT/DP>=$mincov && (FMT/AO)/(FMT/DP)>=$minfrac"
	local het_rule="FMT/GT=\"0/1\" && N_ALT=1 && QUAL>=$minqual && FMT/DP>=$mincov && (FMT/AO)/(FMT/DP)>=$het_min_af"
	bcftools view --include "($snippy_rule) || ($het_rule)" "$raw"
}

normalise() {
	vt normalize -r "$refdir/ref.fa" - | bcftools annotate --remove "$KEEP_VCF_TAGS"
}

# Snippy skips snpEff when the reference has no features; follow its lead by
# checking for the database snippy built.
annotate() {
	if [[ -s "$refdir/ref/snpEffectPredictor.bin" ]]; then
		snpEff ann -noLog -noStats -no-downstream -no-upstream -no-utr \
			-c "$refdir/snpeff.config" -dataDir . ref "$1"
	else
		cat "$1"
	fi
}

to_tab() {
	snippy-vcf_to_tab --gff "$refdir/ref.gff" --ref "$refdir/ref.fa" --vcf "$1"
}

# Records gained over snippy, keyed on the normalised CHROM/POS/REF/ALT.
# heterogeneous = AF < 0.9 flags calls with a sizeable reference fraction.
audit_rescued() {
	local fmt='%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%QUAL\t[%DP]\t[%AO]\n'
	{
		printf 'CHROM\tPOS\tREF\tALT\tGT\tQUAL\tDP\tAO\tAF\theterogeneous\n'
		awk -F'\t' -v OFS='\t' '
			NR == FNR { seen[$1 FS $2 FS $3 FS $4] = 1; next }
			!(($1 FS $2 FS $3 FS $4) in seen) {
				af = $7 > 0 ? $8 / $7 : 0
				print $0, sprintf("%.4f", af), (af < 0.9 ? "yes" : "no")
			}' <(bcftools query -f "$fmt" "$filt") <(bcftools query -f "$fmt" "$1")
	} > "$audit"
}

select_calls | normalise > "$prefix.filt.vcf"
annotate "$prefix.filt.vcf" > "$prefix.vcf"
to_tab "$prefix.vcf" > "$prefix.tab"
audit_rescued "$prefix.filt.vcf"

# The process ignores errors, so fail loudly here rather than emit a headerless table.
head -1 "$prefix.tab" | grep -q '^CHROM' || { echo "No header in $prefix.tab" >&2; exit 1; }
echo "snippy_typing_calls: $(($(wc -l < "$audit") - 1)) record(s) added over snippy's filter, see $audit" >&2
