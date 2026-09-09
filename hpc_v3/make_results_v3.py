"""Assemble the Experiment A results table (hpc_v3/results_v3_summary.csv).

Run as ``python hpc_v3/make_results_v3.py`` from the repo root.

One long table, one row per (panel, phenotype, pipeline, disorder, stratum),
drawn from the committed result files listed in SOURCES below.  A value the
source table does not carry gets a clearly-marked ``status=pending`` row with
dummy values, so the slide layout survives a missing input; **as of the CSD3
run of 2026-09-08 every row is observed.**  No hand-typed numbers (the project
rule: a figure reads statistics from a committed table, never from prose).

Two things changed after that run landed:

**1. The PRS source was wrong.**  It read ``hpc/work/results/prs/``, the
ARRAY-genotype arm (n = 3,725 EUR / 4,126 pooled) — not the imputed arm every
published PRS number in this project comes from.  It now reads
``hpc_v3/prs_tables/prs_association_v3.tsv``: imputed genotypes, all nine
phenotypes in one table, one EUR definition.  The difference is not cosmetic —
``global_slope`` × SCZ × EUR was −0.0385 (p = 0.015) under the array arm and is
−0.0468 (p = 0.0022) under the imputed one.

**2. Two v1 comparison rows were added.**  ``baseline_thickness`` and
``global_slope`` are also carried as ``pipeline="v1"``, read from v1's own
result files, so the slide shows v1 against v2 for the two anchors.

READ THE PIPELINE COMPARISON WITH THESE TWO CAVEATS — they are recorded in the
table's own ``caveat`` column so the figure can print them:

* **h² is not the same estimator.**  v1 is single-GRM REML with relatives
  DROPPED (n = 5,649); v2 is Zaitlen two-GRM with relatives KEPT (n = 8,082).
  The gap for ``baseline_thickness`` (0.247 → 0.539) is mostly that change of
  model, not a better measurement.
* **"EUR" is not the same subset.**  v1's PRS EUR arm is a PC-distance cut
  (n = 5,361); the v3 table uses the 5,656-member anchor set (n = 4,116
  phenotyped), which is what the EUR GWAS arm and the control-disorder runs
  use.  **The pooled ("full") stratum IS the same 8,082 in both**, so pooled is
  the clean pipeline comparison and EUR is not.

Choices fixed here, so the slide cannot drift from them:

* h²: Zaitlen two-GRM ``h2_snp`` on the PC-Relate GRM, relatives kept
  (n = 8,082) — one estimate per phenotype as agreed. The kinship-threshold
  sensitivity (0.183 → 0.115 for global_slope) is a caveat for text, not a
  second error bar.
* PRS: population association at C+T p < 0.5 (the project's headline
  threshold), pooled ("full") and EUR strata.  Threshold multiplicity is
  corrected with the effective-test count, p_corr = min(1, p × M_EFF) with
  M_EFF = 4.6 — the most conservative of the three estimators in
  hpc/README_HPC.md (range 2.7–4.6) — not Bonferroni × 8.
* MAGMA: the full locus pools only (SCZ_locus_pool, MDD_pool), marginal
  model, EUR.  Unsigned test — say so wherever these are shown.
* rg: LDSC vs SCZ/MDD, EUR, with the source table's ``underpowered`` flag
  (either trait h² z < 4) carried through.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parent

PHENOS = ["baseline_thickness", "global_slope", "slope_topDelta",
          "slope_topC3", "slope_projDelta", "slope_projC3"]
#: the two anchors are additionally carried through v1's own pipeline
V1_PHENOS = ["baseline_thickness", "global_slope"]

M_EFF = 4.6          # conservative end of hpc/README_HPC.md's 2.7-4.6
PRS_THRESHOLD = "0p5"

SOURCES = {
    "h2": "hpc_v2/work/results_v2/reml_zaitlen_pcrel/reml_zaitlen_summary.tsv",
    # imputed arm, anchor-set EUR -- NOT hpc/work/results/prs/ (array arm)
    "prs": "hpc_v3/prs_tables/prs_association_v3.tsv",
    # two files each: the settled five (magma_prio_eur, from v2) and the four
    # new ones (magma_prio_eur_v3, from hpc_v3/05_prio_gsa.sbatch).  Same
    # harness, same gene sets, separate output dirs so neither rewrites the
    # other's <TAG>_gsa_summary.tsv -- see 05_prio_gsa.sbatch.
    "magma_scz": ["hpc_v2/work/results_v2/magma_prio_eur/sczprio_gsa_summary.tsv",
                  "hpc_v2/work/results_v2/magma_prio_eur_v3/sczprio_gsa_summary.tsv"],
    "magma_mdd": ["hpc_v2/work/results_v2/magma_prio_eur/mddhc_gsa_summary.tsv",
                  "hpc_v2/work/results_v2/magma_prio_eur_v3/mddhc_gsa_summary.tsv"],
    "rg": "hpc_v2/work/results_v2/ldsc_eur_v3/ldsc_rg_summary.tsv",
}

#: v1's equivalents, for the two anchor phenotypes only.
SOURCES_V1 = {
    "h2": "hpc/work/results/reml_imp_pooled/reml_summary.tsv",
    "prs": "hpc/work/results/prs_imp/prs_association.tsv",
    "magma_scz": "hpc/work/results/magma_prio_imp_eur/sczprio_gsa_summary.tsv",
    "magma_mdd": "hpc/work/results/magma_prio_imp_eur/mddhc_gsa_summary.tsv",
    "rg": "hpc/work/results/ldsc_imp_eur/ldsc_rg_summary.tsv",
}

CAVEAT = {
    "h2": "v1 single-GRM REML, relatives dropped (n=5,649); "
          "v2 Zaitlen two-GRM, relatives kept (n=8,082)",
    "prs": "same scores, so pooled (n=8,082) is the IDENTICAL estimate in both "
           "and is drawn once; v1's row shows only EUR, a PC-distance cut "
           "(n=5,361) against v3's anchor set (n=4,116).",
    "magma": "v1 vs v2 GWAS, same gene sets and MAGMA code",
    "rg": "v1 vs v2 GWAS, same LDSC code and LD reference",
}

#: dummy values for pending rows -- chosen to be visually plausible so the
#: layout is realistic, and flagged status=pending so the slide greys them out.
DUMMY = {
    "h2": dict(estimate=0.15, se=0.10),
    "prs": dict(estimate=0.0, se=0.016),
    "magma": dict(estimate=0.0, se=0.05),
    "rg": dict(estimate=0.0, se=0.11),
}


def read(spec: str | list[str]) -> pd.DataFrame:
    """Read one source, or several concatenated (a later file wins on a tie).

    Several exist because a v3 output directory is always separate from the v2
    one it extends -- the collectors in this project rebuild their summary
    tables from a glob and would otherwise overwrite published rows.
    """
    paths = [spec] if isinstance(spec, str) else spec
    frames = [pd.read_csv(REPO / q, sep="\t") for q in paths
              if (REPO / q).exists()]
    if not frames:
        raise FileNotFoundError(paths[0])
    return pd.concat(frames, ignore_index=True)


def rows() -> list[dict]:
    out: list[dict] = []

    def add(panel, pheno, pipeline, disorder, stratum, hit, source,
            status="pending", **kw):
        """A row.  `hit=None` means the source table has no value; `status`
        then says WHY -- "pending" (the input has not been produced) or
        "not_estimable" (it was produced and the method declined to return a
        number).  The figure must not draw those two the same way: one is a
        gap in this run, the other is a result."""
        base = dict(panel=panel, phenotype=pheno, pipeline=pipeline,
                    disorder=disorder, stratum=stratum, source=source,
                    caveat=CAVEAT[panel])
        if hit is None:
            out.append(base | DUMMY[panel if panel in DUMMY else "prs"]
                       | dict(p=np.nan, n=np.nan, status=status) | kw)
        else:
            out.append(base | hit | dict(status="observed") | kw)

    def build(pipeline: str, src: dict, phenos: list[str]) -> None:
        # ---- h2 -------------------------------------------------------------
        h2 = read(src["h2"])
        # v1's collector names the columns h2/se; v2's two-GRM one h2_snp/se_snp
        est, sd = (("h2_snp", "se_snp") if "h2_snp" in h2.columns
                   else ("h2", "se"))
        for ph in phenos:
            m = h2[h2.phenotype == ph]
            hit = (dict(estimate=float(m[est].iloc[0]),
                        se=float(m[sd].iloc[0]), p=float(m.pval.iloc[0]),
                        n=int(m.n.iloc[0])) if len(m) else None)
            # stratum="full": both REML arms are the pooled multi-ancestry
            # sample, so this must carry the pooled marker, not the EUR one.
            add("h2", ph, pipeline, "", "full", hit, src["h2"])

        # ---- PRS ------------------------------------------------------------
        # A polygenic score applies consortium weights to ABCD genotypes and
        # never touches our GWAS, so for the POOLED stratum v1 and v2 are the
        # same estimate on the same 8,082 subjects -- drawing it twice would
        # look like a replication and is only a duplicate.  We assert the
        # equality instead (a real check that the two tables agree) and give
        # v1 only its EUR row, which is what genuinely differs: a PC-distance
        # cut (n = 5,361) against the anchor set (n = 4,116).
        prs = read(src["prs"])
        prs = prs[prs.threshold == PRS_THRESHOLD]
        strata = ("EUR", "full") if pipeline == "v2" else ("EUR",)
        for ph in phenos:
            for dis in ("SCZ", "MDD"):
                if pipeline == "v1":
                    a = prs[(prs.phenotype == ph) & (prs.disorder == dis)
                            & (prs.stratum == "full")]
                    b = [r for r in out if r["panel"] == "prs"
                         and r["phenotype"] == ph and r["pipeline"] == "v2"
                         and r["disorder"] == dis and r["stratum"] == "full"]
                    if len(a) and b:
                        assert abs(float(a.beta.iloc[0]) - b[0]["estimate"]) < 1e-9, (
                            f"v1 and v2 pooled PRS disagree for {ph} x {dis}: "
                            f"{a.beta.iloc[0]} vs {b[0]['estimate']}")
                for stratum in strata:
                    m = prs[(prs.phenotype == ph) & (prs.disorder == dis)
                            & (prs.stratum == stratum)]
                    hit = (dict(estimate=float(m.beta.iloc[0]),
                                se=float(m.se.iloc[0]), p=float(m.p.iloc[0]),
                                n=int(m.n.iloc[0])) if len(m) else None)
                    add("prs", ph, pipeline, dis, stratum, hit, src["prs"],
                        m_eff=M_EFF, threshold=PRS_THRESHOLD)

        # ---- MAGMA locus pools ----------------------------------------------
        for dis, key, gset in (("SCZ", "magma_scz", "SCZ_locus_pool"),
                               ("MDD", "magma_mdd", "MDD_pool")):
            g = read(src[key])
            g = g[(g.model == "marginal") & (g["set"] == gset)].drop_duplicates(
                subset=["trait", "set"])
            for ph in phenos:
                m = g[g.trait == ph]
                hit = (dict(estimate=float(m.beta.iloc[0]),
                            se=float(m.se.iloc[0]), p=float(m.p.iloc[0]),
                            n=int(m.n_genes.iloc[0])) if len(m) else None)
                add("magma", ph, pipeline, dis, "EUR", hit,
                    src[key] if isinstance(src[key], str) else src[key][-1],
                    gene_set=gset)

        # ---- LDSC rg --------------------------------------------------------
        rg = read(src["rg"])
        for ph in phenos:
            for dis in ("SCZ", "MDD"):
                m = rg[(rg.phenotype == ph) & (rg.disorder == dis)]
                if len(m) and pd.notna(m.rg.iloc[0]):
                    add("rg", ph, pipeline, dis, "EUR",
                        dict(estimate=float(m.rg.iloc[0]),
                             se=float(m.se.iloc[0]), p=float(m.p.iloc[0]),
                             n=np.nan), src["rg"],
                        h2_z=float(m.h2_z.iloc[0]),
                        underpowered=str(m.underpowered.iloc[0]))
                else:
                    # LDSC returns no rg when it cannot normalise by
                    # sqrt(h2_1 * h2_2): at h2 z ~ 0.3 the jackknife h2 goes
                    # negative in some blocks and the denominator has no
                    # square root.  That is a FINDING about power, not a
                    # missing input -- see hpc_v3/README_HPC.md section 7.
                    add("rg", ph, pipeline, dis, "EUR", None, src["rg"],
                        status="not_estimable",
                        h2_z=np.nan, underpowered="not estimable")

    build("v2", SOURCES, PHENOS)
    build("v1", SOURCES_V1, V1_PHENOS)
    return out


def main() -> int:
    df = pd.DataFrame(rows())
    n_obs = int((df.status == "observed").sum())
    out = HERE / "results_v3_summary.csv"
    df.to_csv(out, index=False)
    print(f"{out}: {len(df)} rows, {n_obs} observed / "
          f"{len(df) - n_obs} pending")
    print(df[df.status == "observed"]
          .groupby("panel").size().rename("observed_rows").to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
