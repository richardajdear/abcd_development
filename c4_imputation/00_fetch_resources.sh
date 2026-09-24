#!/bin/bash
# Step 0: download the public reference panel and the two Java tools into
# c4_imputation/resources/ (gitignored), and check them against the checksums
# recorded when the pipeline was built (2026-09-24).  Run on the laptop and on
# a CSD3 login node (compute nodes have no internet).
#
#   bash c4_imputation/00_fetch_resources.sh
#
# Panel: Sekar et al. 2016 (Nature 530:177) C4 reference haplotypes, 111 HapMap3
# CEU individuals (222 haplotypes), 7,752 markers incl. the multi-allelic C4
# marker, GRCh37 VCF as distributed by the imputec4 protocol
# (github.com/freeseek/imputec4).  European only -- see README.md, "Panel".
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p resources; cd resources
fetch() {  # url sha256
  local f; f=$(basename "$1")
  [[ -s "$f" ]] || curl -sSfL --max-time 600 -o "$f" "$1"
  local got; got=$( (sha256sum "$f" 2>/dev/null || shasum -a 256 "$f") | cut -d' ' -f1)
  [[ "$got" == "$2" ]] || { echo "CHECKSUM MISMATCH $f: $got (expected $2)" >&2; exit 3; }
  echo "ok $f"
}
fetch https://personal.broadinstitute.org/giulio/panels/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh37.vcf.gz \
  7ec66e1b4225a46e4d15520620f09055f66c5c4539a1c5f8cc89e82952ff513a
fetch https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar \
  7319f4af9638be05c18dcc1bfb8fb41a58a09293507ebf0d54617d0e40df5a70
fetch https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar \
  87728763085e198fbb274f319d3855b106e3939b3d1a5e3fb8373a50465371a6
