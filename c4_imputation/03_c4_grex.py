"""Step 3: imputed C4 structural alleles -> per-child C4 copy numbers and GREx.

Usage:
  python 03_c4_grex.py --vcf <beagle out>.vcf.gz --out-calls work/c4_calls.tsv \
      --out-summary results/c4_imputation_summary.tsv [--keep-ids eur.keep] \
      [--panel resources/MHC_haplotypes_CEU_HapMap3_ref_panel.GRCh37.vcf.gz]

Reads the multi-allelic C4 marker (GRCh37 6:31948000, ID "C4") from the Beagle
output, which must have been run with ap=true so that every haplotype carries a
posterior over the 25 structural alleles of the Sekar et al. 2016 panel.

Per haplotype, the expected number of each C4 gene segment is
    E[seg] = sum_k AP_k * count(seg in allele k),  seg in {AL, AS, BL, BS}
(AL = C4A long (HERV-carrying), AS = C4A short, BL/BS the same for C4B; the
numeric suffix in allele names such as "AL-BS-3" marks the SNP-haplotype group
and is ignored for copy number).  Per child, the two haplotypes are summed.

GREx (genetically regulated expression in brain) uses the weights printed in
Hernandez et al. 2023 (Genome Biology 24:42, Methods), who cite Sekar 2016:
    C4A_GREx = 0.47*AL + 0.47*AS + 0.20*BL
    C4B_GREx = 1.03*BL + 0.88*BS
The weights are a module constant (GREX_WEIGHTS) so a correction is one edit.

QC fields written per child (not applied here; 04_c4_assoc.R applies them):
  post_struct_h1/h2  posterior of the most likely STRUCTURE (suffix groups
                     collapsed) on each haplotype
  post_mean          mean of the two -- Hernandez et al. excluded children with
                     average posterior < 0.7
  common5            both hard-called structures in {AL, AL-AL, AL-BL, AL-BS, BS}

Per-subject output (--out-calls) is individual-level genotype-derived data:
never commit it (genetic_analysis/README_HPC.md rule 16).  The summary table is
aggregate and is the file to commit.
"""
from __future__ import annotations

import argparse
import gzip
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

C4_POS_GRCH37 = 31948000
SEGMENTS = ("AL", "AS", "BL", "BS")
GREX_WEIGHTS = {"C4A_GREx": {"AL": 0.47, "AS": 0.47, "BL": 0.20},
                "C4B_GREx": {"BL": 1.03, "BS": 0.88}}
COMMON5 = {"AL", "AL-AL", "AL-BL", "AL-BS", "BS"}


def structure(allele: str) -> str:
    """'AL-BS-3' -> 'AL-BS'; 'AL-BL-other' -> 'AL-BL'; 'BS' -> 'BS'."""
    return "-".join(p for p in allele.split("-") if p in SEGMENTS)


def seg_counts(allele: str) -> dict:
    c = Counter(p for p in allele.split("-") if p in SEGMENTS)
    return {s: c.get(s, 0) for s in SEGMENTS}


def token(x: str) -> str:
    """The 8-character ABCD subject token (README_HPC.md rule 1)."""
    m = re.search(r"([A-Z0-9]{8})$", str(x))
    return m.group(1) if m else str(x)


def read_c4_record(vcf: str, bcftools: str = "bcftools"):
    """Return (alt alleles, samples, per-sample (GT, AP1, AP2), INFO dict)."""
    rec = subprocess.run([bcftools, "view", "-H", "-i", 'ID=="C4"', vcf], check=True,
                         capture_output=True, text=True).stdout.strip().splitlines()
    if len(rec) != 1:
        sys.exit(f"expected exactly one C4 record in {vcf}, found {len(rec)}")
    samples = subprocess.run([bcftools, "query", "-l", vcf], check=True,
                             capture_output=True, text=True).stdout.split()
    f = rec[0].split("\t")
    if int(f[1]) != C4_POS_GRCH37:
        sys.exit(f"C4 record at {f[0]}:{f[1]}, expected GRCh37 {C4_POS_GRCH37}")
    alts = [x.strip("<>") for x in f[4].split(",")]
    fmt = f[8].split(":")
    if "AP1" not in fmt or "AP2" not in fmt:
        sys.exit("C4 record has no AP1/AP2 -- rerun Beagle with ap=true")
    info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
    gi, a1, a2 = fmt.index("GT"), fmt.index("AP1"), fmt.index("AP2")
    per = []
    for s in f[9:]:
        x = s.split(":")
        per.append((x[gi], np.array(x[a1].split(","), float), np.array(x[a2].split(","), float)))
    return alts, samples, per, info


def panel_structure_freqs(panel: str) -> pd.Series:
    """Structure frequencies among the reference haplotypes (for the summary)."""
    with gzip.open(panel, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] == "C4":
                alleles = ["REF"] + [x.strip("<>") for x in f[4].split(",")]
                haps = [alleles[int(h)] for g in f[9:] for h in re.split(r"[|/]", g)]
                return pd.Series([structure(h) for h in haps]).value_counts(normalize=True)
    sys.exit(f"no C4 record in {panel}")


def calls_table(alts, samples, per) -> pd.DataFrame:
    counts = np.array([[seg_counts(al)[s] for s in SEGMENTS] for al in alts], float)  # K x 4
    structs = [structure(al) for al in alts]
    ustruct = sorted(set(structs))
    S = np.array([[1.0 if s == u else 0.0 for u in ustruct] for s in structs])      # K x U
    rows = []
    for iid, (gt, p1, p2) in zip(samples, per):
        if len(p1) != len(alts) or len(p2) != len(alts):
            sys.exit(f"AP length mismatch for {iid}")
        s1, s2 = p1.sum(), p2.sum()
        # Beagle prints AP to 2 decimals over 25 alleles, so sums drift by a few
        # hundredths; renormalise before taking expectations.
        p1, p2 = p1 / s1, p2 / s2
        e = (p1 + p2) @ counts                  # expected segment counts, both haplotypes
        sp1, sp2 = p1 @ S, p2 @ S               # structure posteriors per haplotype
        h = re.split(r"[|/]", gt)
        call = [structs[int(x) - 1] if x not in (".", "0") else "NA" for x in h]
        r = dict(IID=iid, tok=token(iid), hap1=call[0], hap2=call[1],
                 post_struct_h1=sp1.max(), post_struct_h2=sp2.max(),
                 ap_sum_h1=s1, ap_sum_h2=s2,
                 **{f"E_{s}": v for s, v in zip(SEGMENTS, e)})
        r["post_mean"] = (r["post_struct_h1"] + r["post_struct_h2"]) / 2
        r["C4A_copies"] = r["E_AL"] + r["E_AS"]
        r["C4B_copies"] = r["E_BL"] + r["E_BS"]
        r["HERV_copies"] = r["E_AL"] + r["E_BL"]
        for g, w in GREX_WEIGHTS.items():
            r[g] = sum(wt * r[f"E_{s}"] for s, wt in w.items())
        r["common5"] = int(call[0] in COMMON5 and call[1] in COMMON5)
        rows.append(r)
    d = pd.DataFrame(rows)
    # Beagle writes AP over the ALT alleles only.  REF ("G") is a placeholder
    # that no reference haplotype carries, so each haplotype's AP must sum to 1
    # up to print rounding (<= 0.125 for 25 alleles at 2 dp; > 0.15 means real
    # mass on REF, i.e. a panel problem).
    bad = (d[["ap_sum_h1", "ap_sum_h2"]].sub(1).abs() > 0.15).any(axis=1)
    if bad.mean() > 0.001:
        sys.exit(f"{bad.sum()} samples have AP not summing to 1 -- REF carries mass; check the panel")
    return d


def summary_table(d: pd.DataFrame, keep_ids, panel, alts, info) -> pd.DataFrame:
    groups = {"all": d}
    if keep_ids:
        keep = {token(x.split()[-1]) for x in open(keep_ids) if x.strip()}
        groups["keep"] = d[d.tok.isin(keep)]
    srows = []
    for gname, g in groups.items():
        srows += [dict(group=gname, metric="n_samples", key="", value=len(g)),
                  dict(group=gname, metric="post_mean_median", key="", value=g.post_mean.median()),
                  dict(group=gname, metric="frac_post_mean_ge_0.7", key="", value=(g.post_mean >= 0.7).mean()),
                  dict(group=gname, metric="frac_common5", key="", value=g.common5.mean())]
        for col in ("C4A_copies", "C4B_copies", "HERV_copies", "C4A_GREx", "C4B_GREx"):
            srows += [dict(group=gname, metric=f"{col}_mean", key="", value=g[col].mean()),
                      dict(group=gname, metric=f"{col}_sd", key="", value=g[col].std())]
        srows.append(dict(group=gname, metric="corr_C4A_C4B_GREx", key="", value=g.C4A_GREx.corr(g.C4B_GREx)))
        hf = pd.concat([g.hap1, g.hap2]).value_counts(normalize=True)
        srows += [dict(group=gname, metric="structure_freq_hardcall", key=k, value=v) for k, v in hf.items()]
    if panel:
        srows += [dict(group="reference_panel", metric="structure_freq", key=k, value=v)
                  for k, v in panel_structure_freqs(panel).items()]
    srows.append(dict(group="all", metric="beagle_DR2_per_allele", key=",".join(alts), value=info.get("DR2", "")))
    return pd.DataFrame(srows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out-calls", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--keep-ids", help="one ID per line (any ABCD spelling, last column); summary also split by it")
    ap.add_argument("--panel", help="reference panel VCF, for panel-vs-target frequency rows")
    ap.add_argument("--bcftools", default="bcftools")
    a = ap.parse_args()

    alts, samples, per, info = read_c4_record(a.vcf, a.bcftools)
    d = calls_table(alts, samples, per)
    Path(a.out_calls).parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(a.out_calls, sep="\t", index=False, float_format="%.5g")
    s = summary_table(d, a.keep_ids, a.panel, alts, info)
    Path(a.out_summary).parent.mkdir(parents=True, exist_ok=True)
    s.to_csv(a.out_summary, sep="\t", index=False, float_format="%.4g")
    print(f"{len(d)} samples; median structure posterior {d.post_mean.median():.3f}; "
          f"{(d.post_mean >= 0.7).mean():.1%} >= 0.7; common5 {d.common5.mean():.1%}")


if __name__ == "__main__":
    main()
