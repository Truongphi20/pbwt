#!/usr/bin/env bash
set -euo pipefail

target="tests/data/OMNI.vcf"
output="tests/data/query_OMNI.vcf.gz"

samples_file=$(mktemp)
tmp_subset=$(mktemp --suffix=.vcf.gz)

echo ">>> Compressing and indexing input"

bgzip -c "$target" > "${target}.gz"
tabix -f -p vcf "${target}.gz"

echo ">>> Extract first 10 samples"

bcftools query -l "${target}.gz" | head -n 10 > "$samples_file" || true

echo ">>> Subset 10 samples"

bcftools view \
    -S "$samples_file" \
    -Oz \
    -o "$tmp_subset" \
    "${target}.gz"

echo ">>> Keep first 5 variants"

bcftools view \
    -H "$tmp_subset" | head -n 5 > variants.tmp

bcftools view -h "$tmp_subset" > header.tmp

cat header.tmp variants.tmp | bgzip -c > "$output"

tabix -f -p vcf "$output"

echo ">>> Cleaning up"
rm -f "$samples_file" "$tmp_subset" variants.tmp header.tmp

echo "Done: $output"