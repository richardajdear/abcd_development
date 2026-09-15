"""Order-of-operations check: per-region LMMs then average the slope BLUPs (the
pipeline's global_slope) vs one LMM on the per-scan whole-cortex mean
(global_slope_1lmm).  Matched SCZ/MDD cells, both atlases."""
import pandas as pd, glob, re, sys
M={"SCZ_pooled":"full","MDD_pooled":"full","SCZ_eur":"EUR","MDD_eur":"EUR"}
rows=[]
for atlas,root in (("dk","genetic_analysis/work/results_70tab"),("hcp","genetic_analysis/work/results_70tab_hcp")):
    base=pd.read_csv(f"{root}/prs_final/table_main.tsv",sep="\t")
    base=base[(base.matched=="yes")&(base.phenotype=="global_slope")&(base.trait_arm.isin(M))]
    for f in glob.glob(f"{root}/prs_final_1lmm/assoc_*.tsv"):
        m=re.match(r".*/assoc_(CT|PRSCS|SBayesR|SBayesRC)_(SCZ_pooled|SCZ_eur|MDD_pooled|MDD_eur)(_zanc)?\.tsv$",f)
        if not m: continue
        meth,arm,z=m.group(1),m.group(2),("zanc" if m.group(3) else "raw")
        want=("zanc" if M[arm]=="full" else "raw")
        if z!=want: continue
        t=pd.read_csv(f,sep="\t"); t=t[(t.phenotype=="global_slope_1lmm")&(t.stratum==M[arm])]
        if t.empty: continue
        one=t.sort_values("p").iloc[0]
        b=base[(base.trait_arm==arm)&(base.method==meth)&(base.score==z)&(base.target_stratum==M[arm])]
        if b.empty: continue
        b=b.sort_values("p").iloc[0]
        rows.append(dict(atlas=atlas,trait_arm=arm,stratum=M[arm],method=meth,score=z,
            beta_perregion=b.beta,se_perregion=b.se,p_adj_perregion=b.p_adj,
            beta_1lmm=one.beta,se_1lmm=one.se,p_adj_1lmm=one.p_adj,n=int(one.n)))
out=pd.DataFrame(rows).sort_values(["atlas","trait_arm","method"])
for atlas,root in (("dk","genetic_analysis/work/results_70tab"),("hcp","genetic_analysis/work/results_70tab_hcp")):
    out[out.atlas==atlas].to_csv(f"{root}/prs_final_1lmm/table_order_of_operations.tsv",sep="\t",index=False)
pd.set_option("display.width",220)
f=out.copy()
for c in ["beta_perregion","se_perregion","beta_1lmm","se_1lmm"]: f[c]=f[c].round(4)
for c in ["p_adj_perregion","p_adj_1lmm"]: f[c]=f[c].map(lambda v:f"{v:.3g}")
print(f.drop(columns=["score","n"]).to_string(index=False))
d=(out.beta_1lmm-out.beta_perregion); print(f"\nbeta_1lmm - beta_perregion: mean {d.mean():+.4f}, max |diff| {d.abs().max():.4f} (max |beta| {out.beta_perregion.abs().max():.4f}); SE ratio 1lmm/perregion mean {(out.se_1lmm/out.se_perregion).mean():.3f}")
