"""Step 3 (CSD3, ~10-30 min per arm x outcome set): partition tests.

    $D7PY code/03_partition_test.py --arm SCZ25_META --cell pooled --outcomes imaging_hcp \
        --work $WORK --pheno-dir $PHENO_HCP --strata $STRATA --eur-keep $EUR_KEEP \
        --out $WORK/partition/SCZ25_META__imaging_hcp.tsv

--outcomes: imaging_hcp | imaging_dk (--pheno-dir = a prs_final_1lmm/pheno export)
            symptoms    (--symptoms d7_symptom_outcomes.tsv; PCs from --pheno-dir)
            c3axis      (--pheno-dir = results_70tab_hcp/c3axis/pheno; optional)

Model (per outcome, matching tools/prs_assoc.R except OLS for lmer, gated below):
    imaging    scale(y) ~ scale(score) + sex + age_c + PC1..PC10
    symptoms   y_LATE ~ scale(score) + y_BASE + age_LATE + sex + site + PC1..PC10
               (Figure-1 panel h definition plus PCs; y scaled to SD 1)
Family-clustered SEs throughout.  Pooled cell: every gene column gets the linear
form of the `_zanc` transform (centre within ancestry cluster, divide by the
cluster SD of the full score), so parts still add to the Figure-1 score.  EUR
cell: raw score, EUR-arm children only.

Per set T (partition_core docstring): beta_T, s_T, f_T, ER_T = s_T / f_T,
family-bootstrap CI of ER, and two matched random-set nulls from T's universe:
  null A  bins = decile of n weighted SNPs per gene          -> p_f   (disorder side:
          does T carry more of the disorder score's variance than its SNP count predicts?)
  null B  bins = quintile of n SNPs x quintile of d_j        -> p_ER  (thinning side:
          given its share of score variance, does T carry more or less of the
          association with y than expected?)
Read ER against its matched null, not against 1: intergenic SNPs are part of the
full score, so a genic set's ER exceeds 1 whenever genic SNPs carry more than
their share.  ER_rel = ER / median(null B) is the enrichment beyond a random set
of matched genes.  ER is only interpretable where the full score is associated with y; rows carry
`full_p` and `er_interpretable` (full_p < 0.05).

GATE: beta of the full score (with MHC) in the OLS here vs the lmer value in the
step-9 association table, if found, |diff| < 0.005.  Recorded, not fatal.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import partition_core as pc  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
CBCL = ["pfactor", "internal", "external", "depress", "thought"]
PAIRS = [("K2", "K2_neuronal", "K2_glial"), ("K3", "K3_syngo", "K3_oligodendrocyte"),
         ("S1", "S1_snrna_pc1_top", "S1_snrna_pc1_bottom"), ("S2", "S2_c3_top", "S2_c3_bottom"),
         ("S3", "S3_c1_top", "S3_c1_bottom"), ("S3", "S3_c2_top", "S3_c2_bottom"),
         ("S4", "S4_pls2_top", "S4_pls2_bottom")]


def token(x) -> str:
    """Rule 1: join on the 8-character id token."""
    x = re.sub(r"^sub-", "", str(x))
    x = re.sub(r"^NDAR_?INV", "", x)
    return x.replace("_", "").upper()


def load_outcomes(kind, pheno_dir, symptoms):
    pdir = Path(pheno_dir)
    q = pd.read_csv(pdir / "covar_quant.txt", sep=r"\s+")
    q["token"] = q.IID.map(token)
    pcs = [c for c in q.columns if re.fullmatch(r"PC\d+", c)]
    if kind == "symptoms":
        s = pd.read_csv(symptoms, sep="\t", dtype={"token": str})
        d = s.merge(q[["token", *pcs]], on="token", how="inner")
        out = []
        for o in CBCL:
            out.append(dict(name=o, y=f"{o}_LATE", num=["age_LATE", f"{o}_BASE", *pcs],
                            cat=["sex", "site"], family="family"))
        return d, out
    ph = pd.read_csv(pdir / "phenotypes_gcta.txt", sep=r"\s+")
    cc = pd.read_csv(pdir / "covar_categorical.txt", sep=r"\s+")
    man = pd.read_csv(pdir / "phenotype_manifest.tsv", sep="\t")
    for t in (ph, cc):
        t["token"] = t.IID.map(token)
    d = ph.merge(q.drop(columns=["FID", "IID"]), on="token").merge(
        cc[["token", "sex"]], on="token")
    d["family"] = d.FID
    d["age_c"] = d.baseline_age - d.baseline_age.mean()
    names = [n for n in man.name if n in ph.columns]
    return d, [dict(name=n, y=n, num=["age_c", *pcs], cat=["sex"], family="family") for n in names]


def design(d, num, cat):
    X = [np.ones(len(d))]
    for c in num:
        X.append(d[c].astype(float).values)
    for c in cat:
        X.append(pd.get_dummies(d[c].astype(str), drop_first=True).astype(float).values)
    return np.column_stack(X)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True)
    ap.add_argument("--cell", choices=["pooled", "EUR"], required=True)
    ap.add_argument("--outcomes", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--pheno-dir", required=True)
    ap.add_argument("--symptoms")
    ap.add_argument("--strata", required=True)
    ap.add_argument("--eur-keep", required=True)
    ap.add_argument("--assoc-ref", help="step-9 association tsv for the beta gate (optional)")
    ap.add_argument("--n-null", type=int, default=5000)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    work = Path(a.work)

    G = np.load(work / "gene_scores" / f"{a.arm}.npy").astype(np.float64)
    cols = pd.read_csv(work / "gene_scores" / f"{a.arm}_cols.tsv", sep="\t")
    iid = pd.read_csv(work / "gene_scores" / "iids.txt", header=None)[0].map(token).values
    gcols = np.flatnonzero((cols.n_snp > 0).values & ~cols.entrez.isin(["INTERGENIC", "MHC"]).values)
    inter = int(cols.index[cols.entrez == "INTERGENIC"][0])
    mhc = int(cols.index[cols.entrez == "MHC"][0])
    ent2col = pd.Series(gcols, index=cols.entrez.values[gcols].astype(int))
    nsnp = cols.n_snp.values.astype(float)

    full_with_mhc = G.sum(1)
    if a.cell == "pooled":
        st = pd.read_csv(a.strata, sep="\t")
        cl = pd.Series(st.stratum.values, index=st.IID.map(token)).reindex(iid).fillna("UNASSIGNED").values
        G, full_with_mhc = pc.within_cluster_transform(G, full_with_mhc, cl)
    keep_eur = set(pd.read_csv(a.eur_keep, sep=r"\s+", header=None)[0].map(token))

    sets = pd.read_csv(HERE / "gene_sets" / "d7_gene_sets.tsv", sep="\t")
    meta = pd.read_csv(HERE / "gene_sets" / "d7_sets_meta.tsv", sep="\t").set_index("set")
    uni_tab = pd.read_csv(HERE / "gene_sets" / "d7_universes.tsv", sep="\t")
    universes = {"genome": gcols}
    for u, dd in uni_tab.groupby("universe"):
        universes[u] = np.intersect1d(ent2col.reindex(dd.entrez).dropna().astype(int).values, gcols)
    set_cols = {s: np.intersect1d(ent2col.reindex(dd.entrez).dropna().astype(int).values,
                                  universes[meta.loc[s, "universe"]])
                for s, dd in sets.groupby("set")}
    # K1 is ~90 % of the genome; its informative side is the complement
    set_cols["K1c_not_brain_expressed"] = np.setdiff1d(gcols, set_cols["K1_brain_expressed"])
    meta.loc["K1c_not_brain_expressed"] = dict(tier="primary", contrast="K1", universe="genome",
                                               n_genes=len(set_cols["K1c_not_brain_expressed"]),
                                               source="complement of K1 in the genome universe")

    d, outcomes = load_outcomes(a.outcomes, a.pheno_dir, a.symptoms)
    row_of = pd.Series(np.arange(len(iid)), index=iid)
    d = d[d.token.isin(row_of.index)].copy()
    if a.cell == "EUR":
        d = d[d.token.isin(keep_eur)].copy()
    ref = pd.read_csv(a.assoc_ref, sep="\t") if a.assoc_ref and Path(a.assoc_ref).exists() else None

    rows, pair_rows = [], []
    draws_cache: dict = {}
    for oc in outcomes:
        dd = d.dropna(subset=[oc["y"], *oc["num"], *oc["cat"]])
        r = row_of.reindex(dd.token).values
        ok = np.isfinite(full_with_mhc[r])
        dd, r = dd[ok], r[ok]
        n = len(r)
        if n < 500:
            continue
        C = design(dd, oc["num"], oc["cat"])
        fam = dd[oc["family"]].astype(str).values
        y = dd[oc["y"]].astype(float).values
        y = (y - y.mean()) / y.std(ddof=1)
        Gr = G[r]
        S_nomhc = Gr[:, gcols].sum(1) + Gr[:, inter]
        sdS = S_nomhc.std(ddof=1)
        Gt = pc.residualise(Gr / sdS, C)
        St = Gt[:, gcols].sum(1) + Gt[:, inter]
        yt = pc.residualise(y[:, None], C)[:, 0]
        Ft = pc.residualise(((full_with_mhc[r] - full_with_mhc[r].mean()) / full_with_mhc[r].std(ddof=1))[:, None], C)[:, 0]
        b_full, se_full, p_full = pc.ols_cluster(yt, Ft, fam)
        b_nm, se_nm, p_nm = pc.ols_cluster(yt, St, fam)
        gate = np.nan
        if ref is not None:
            m = ref[(ref.phenotype == oc["name"]) & (ref.stratum == ("full" if a.cell == "pooled" else "EUR"))]
            if len(m):
                gate = float(b_full - m.beta.iloc[0])
        c, dcov, cyS, vS = pc.gene_moments(yt, Gt, St)
        f_mhc_excluded = float(np.var(Gt[:, mhc], ddof=1))

        parts, est = {}, {}
        for sname, idx in set_cols.items():
            if len(idx) < 10:
                continue
            uni = universes[meta.loc[sname, "universe"]]
            s, f, er = pc.share(idx, c, dcov, cyS, vS)
            ST = Gt[:, idx].sum(1)
            b_T, se_T, p_T = pc.ols_cluster(yt, ST / ST.std(ddof=1), fam)
            key = (sname, "A")
            if sname == "K1_brain_expressed":            # nearly the whole universe: no null
                pf = fe = per = lo_n = hi_n = med_n = np.nan
                nA = nB = np.full(1, np.nan)
            else:
                if key not in draws_cache:           # score-only bins: reuse across outcomes
                    binsA = pc.bin_genes(uni, nsnp, q=(10,))
                    binsB = pc.bin_genes(uni, nsnp, dcov, q=(5, 5))
                    draws_cache[key] = pc.matched_draws(idx, uni, binsA, a.n_null, rng)
                    draws_cache[(sname, "B")] = pc.matched_draws(idx, uni, binsB, a.n_null, rng)
                _, nfA, _ = pc.null_distribution(draws_cache[key], c, dcov, cyS, vS)
                _, _, nB = pc.null_distribution(draws_cache[(sname, "B")], c, dcov, cyS, vS)
                nA = nfA
                pf, fe = pc.p_upper(f, nfA), f / np.median(nfA)
                per = pc.p_two_sided(er, nB)
                lo_n, med_n, hi_n = np.percentile(nB, [2.5, 50, 97.5])
            parts[sname] = ST
            est[sname] = dict(er=er, nullB=nB)
            rows.append(dict(arm=a.arm, cell=a.cell, outcome_set=a.outcomes, outcome=oc["name"], n=n,
                             set=sname, tier=meta.loc[sname, "tier"], contrast=meta.loc[sname, "contrast"],
                             universe=meta.loc[sname, "universe"], n_genes=len(idx),
                             n_snp=int(nsnp[idx].sum()),
                             full_beta=b_full, full_se=se_full, full_p=p_full,
                             full_nomhc_beta=b_nm, full_nomhc_p=p_nm, gate_dbeta_vs_lmer=gate,
                             mhc_var_share=f_mhc_excluded / (f_mhc_excluded + vS),
                             beta=b_T, se=se_T, p=p_T, s=s, f=f, ER=er,
                             f_enrich=fe, p_f=pf, ER_null_lo=lo_n, ER_null_med=med_n, ER_null_hi=hi_n,
                             ER_rel=er / med_n, p_ER=per,
                             er_interpretable=bool(p_nm < 0.05)))
        boot = pc.family_bootstrap(yt, St, parts, fam, a.n_boot, rng)
        for rr in rows[-len(parts):]:
            bs = boot[rr["set"]][:, 2]
            rr["ER_lo"], rr["ER_hi"] = np.percentile(bs, [2.5, 97.5])
        for contrast, sa, sb in PAIRS:
            if sa not in boot or sb not in boot:
                continue
            dl = boot[sa][:, 2] - boot[sb][:, 2]
            obs = est[sa]["er"] - est[sb]["er"]
            nd = est[sa]["nullB"] - est[sb]["nullB"]
            pair_rows.append(dict(arm=a.arm, cell=a.cell, outcome_set=a.outcomes, outcome=oc["name"],
                                  contrast=contrast, set_a=sa, set_b=sb, dER=obs,
                                  dER_lo=np.percentile(dl, 2.5), dER_hi=np.percentile(dl, 97.5),
                                  p_null=pc.p_two_sided(obs, nd), er_interpretable=bool(p_nm < 0.05)))
        print(f"{a.arm} {a.cell} {oc['name']}: n={n} full beta={b_full:+.4f} (p {p_full:.2g}); "
              f"gate dbeta={gate:+.4f}" if np.isfinite(gate) else
              f"{a.arm} {a.cell} {oc['name']}: n={n} full beta={b_full:+.4f} (p {p_full:.2g})", flush=True)

    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(pair_rows).to_csv(out.with_name(out.stem + "__pairs.tsv"), sep="\t", index=False,
                                   float_format="%.6g")


if __name__ == "__main__":
    main()
