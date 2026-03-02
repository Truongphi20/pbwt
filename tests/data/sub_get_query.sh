#!/usr/bin/env bash
set -euo pipefail

target="tests/data/OMNI.vcf"
output="tests/data/query_OMNI.vcf.gz"

samples_file=$(mktemp)

echo ">>> Compressing and indexing input"

bgzip -c "$target" > "${target}.gz"
tabix -f -p vcf "${target}.gz"

echo ">>> Extract first 10 samples"

bcftools query -l "${target}.gz" | head -n 10 > "$samples_file" || true

bcftools view \
    -S "$samples_file" \
    -Oz \
    -o "$output" \
    "${target}.gz"

tabix -f -p vcf "$output"

echo ">>> Cleaning up"
rm -f "$samples_file"

echo "Done: $output"

## Generate pbwt for query
## ./build/pbwt -readVcfGT tests/data/query_OMNI.vcf.gz -write tests/data/query_OMNI.pbwt