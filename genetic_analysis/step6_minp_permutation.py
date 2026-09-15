#!/usr/bin/env python
"""Family-block min-p permutation over the 8 nested C+T thresholds, 7.0 export.

README_HPC.md 4 rule 6: multiplicity over nested C+T thresholds is a
family-block permutation of min-p, not Bonferroni (8 thresholds are ~4.1
effective tests).  It applies to C+T ONLY -- PRS-CS / SBayesR / SBayesRC each
give a single score, so there is nothing to correct.

Reuses the statistics in setup/prs_paired_delta.py unchanged (its
`minp_permutation`, `Stratum`, `load`); only the inputs differ:
  * the 7.0 phenotype export (FID = family id, which the family blocks need),
  * the settled five phenotypes instead of the retired v3 four,
  * TWO designs, because the 5 benchmark (p_perm 0.007, SCZ EUR) was NOT
    ancestry-matched -- prs_ct_v3's SCZ scores are v1's, built from the
    multi-ancestry PGC3 *primary* file and read on the EUR target:
      - matched   : prs_final/CT/<arm> scores in the rule-4 stratum
                    (European GWAS -> EUR, multi-ancestry GWAS -> full)
      - legacy    : prs_ct_v3 scores, both strata, for a like-for-like
                    comparison with the 6.0 benchmark.
Covariates are the module's "primary arm" (sex + 10 PCs, no age), which is
what the benchmark used.

    python genetic_analysis/step6_minp_permutation.py [--n-perm 2000]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "genetic_analysis" / "setup"))
import prs_paired_delta as m  # noqa: E402

PHENOS = ["global_slope", "baseline_thickness", "slope_PC1", "slope_PC2", "slope_PC3"]
FINAL_LEGACY = REPO / "legacy/hpc_v2/work/results_v2/prs_final/CT"
CT_V3 = REPO / "legacy/hpc_v2/work/results_v2/prs_ct_v3"
# trait-arm dir -> rule-4 stratum
MATCHED = {"SCZ_pooled": "full", "MDD_pooled": "full",
           "SCZ_eur": "EUR", "MDD_eur": "EUR", "ASD": "EUR", "ALZ": "EUR",
           "ALZ_noAPOE": "EUR", "ALZ_IGAP": "EUR", "ALZ_IGAP_noAPOE": "EUR", "EA": "EUR"}


def token_in(d: Path) -> str:
    toks = {re.match(r"score_(.+)_[^_]+\.profile$", p.name).group(1)
            for p in d.glob("score_*.profile")}
    if len(toks) != 1:
        raise SystemExit(f"{d}: expected one score token, found {sorted(toks)}")
    return toks.pop()


def run(pheno_dir: Path, eur_ids: Path, prs_dir: Path, token: str,
        strata_names: list[str], n_perm: int, seed: int, design: str,
        trait_arm: str) -> pd.DataFrame:
    m.DISORDERS = [token]
    m.NEW = [p for p in PHENOS if p != m.REF]      # all five, REF last
    d, cov, scores = m.load(pheno_dir, prs_dir, eur_ids, with_age=False)
    if not scores:
        raise SystemExit(f"no score_{token}_<thr>.profile under {prs_dir}")
    strata = [m.Stratum(d, s, cov) for s in strata_names]
    out = m.minp_permutation(strata, scores, n_perm, seed)
    out.insert(0, "design", design)
    out.insert(1, "trait_arm", trait_arm)
    out["matched"] = out.stratum.map(lambda s: "yes" if design == "legacy" and False
                                     else ("yes" if MATCHED.get(trait_arm) == s else "no"))
    print(f"  {design:<8} {trait_arm:<16} {token:<14} strata={strata_names} "
          f"{len(scores)} scores -> {len(out)} rows", flush=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pheno-dir", default=str(REPO / "genetic_analysis/work/pheno_70tab"))
    ap.add_argument("--eur-ids", default=str(REPO / "legacy/hpc/work/results/ancestry/eur_anchor.keep"))
    ap.add_argument("--out", default=str(REPO / "genetic_analysis/work/results_70tab/prs_final/table_minp_permutation.tsv"))
    ap.add_argument("--n-perm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260908 + 1)   # the module's min-p seed
    a = ap.parse_args()
    pheno, eur = Path(a.pheno_dir), Path(a.eur_ids)
    parts = []
    print("=== matched design (rule 4) ===", flush=True)
    for arm, stratum in MATCHED.items():
        d = FINAL_LEGACY / arm
        parts.append(run(pheno, eur, d, token_in(d), [stratum], a.n_perm, a.seed, "matched", arm))
    print("=== legacy-like design (prs_ct_v3, both strata; benchmark comparison) ===", flush=True)
    for tok in ["SCZ", "MDD", "ASD", "ALZ", "ALZnoAPOE"]:
        parts.append(run(pheno, eur, CT_V3, tok, ["EUR", "full"], a.n_perm, a.seed, "legacy", tok))
    res = pd.concat(parts, ignore_index=True)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(a.out, sep="\t", index=False)
    print(f"\nwrote {a.out} ({len(res)} rows)")
    show = res[(res.phenotype == "global_slope")]
    print(show[["design", "trait_arm", "stratum", "n", "min_p", "p_perm", "p_bonferroni"]]
          .to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
