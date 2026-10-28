#!/usr/bin/env bash
#
# Fixture test for bin/snippy_typing_calls.sh on PM1533's snippy 4.6.0 output.
# Run from the repository root inside the snippy container, e.g.
#   docker run --rm -v "$PWD":/p -w /p quay.io/staphb/snippy:4.6.0-SC2 tests/test_snippy_typing_calls.sh
#   apptainer exec docker://quay.io/staphb/snippy:4.6.0-SC2 tests/test_snippy_typing_calls.sh
#
# The raw VCF holds four records: two 1/1 SNPs snippy keeps, the real 0/1 natC
# insertion at 859 (AF 0.78) and a QUAL ~0 complex record at 870 that must stay out.

set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
fixture="$root/tests/snippy_typing_calls"
script="$root/bin/snippy_typing_calls.sh"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

fail() { echo "FAIL: $*" >&2; exit 1; }

# Build the reference dir (ref.fa, ref.gff, snpEff db) the same way snippy does
snippy --cpus 1 --ram 2 --force --quiet --outdir "$work/ref" \
	--ref "$root/databases/LPS/L3_NC002663_LPS.gb" --ctgs "$root/databases/LPS/L3_NC002663_LPS.fasta" \
	> "$work/ref.log" 2>&1 || { cat "$work/ref.log" >&2; fail "reference build failed"; }

run() {
	local dir="$work/$1"
	mkdir -p "$dir"
	cp "$fixture/PM1533_snps.raw.vcf" "$dir/snps.raw.vcf"
	cp "$fixture/PM1533_snps.filt.vcf" "$dir/snps.filt.vcf"
	ln -s "$work/ref/reference" "$dir/reference"
	(cd "$dir" && "$script" --raw snps.raw.vcf --filt snps.filt.vcf --refdir reference \
		--minqual 100 --mincov 10 --minfrac 0 --het-min-af "$2" 2> run.log)
}

run het05 0.5
tab="$work/het05/snps.typing.tab"
awk -F'\t' '$2 == 859 && $3 == "ins" && $4 == "A" && $5 == "ATATTTAT" && $11 ~ /^frameshift_variant/ && $13 == "natC"' "$tab" | grep -q . \
	|| fail "859 A>ATATTTAT frameshift in natC missing from snps.typing.tab"
awk -F'\t' '$2 == 870' "$tab" | grep -q . && fail "QUAL ~0 record at 870 was kept"
# Every row snippy wrote must come through unchanged, plus exactly the rescued one
diff <(grep -v -P '\t859\t' "$tab") "$fixture/PM1533_snps.snippy.tab" > /dev/null \
	|| fail "rows other than 859 differ from snippy's snps.tab"
[[ $(tail -n +2 "$work/het05/snps.rescued.tsv" | cut -f2,10) == $'859\tyes' ]] \
	|| fail "snps.rescued.tsv should list only 859, flagged heterogeneous"

# Above the call's AF (0.78) the output must collapse to snippy's own table
run het08 0.8
diff "$work/het08/snps.typing.tab" "$fixture/PM1533_snps.snippy.tab" > /dev/null \
	|| fail "at --het-min-af 0.8 snps.typing.tab differs from snippy's snps.tab"
[[ $(wc -l < "$work/het08/snps.rescued.tsv") -eq 1 ]] || fail "at 0.8 snps.rescued.tsv should be header only"

echo "PASS: snippy_typing_calls.sh"
