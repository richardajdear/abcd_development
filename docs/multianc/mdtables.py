"""Markdown tables without a tabulate dependency (not installed in work/envs/abcd)."""
def md(df, floatfmt="{:.1f}"):
    def f(v):
        if isinstance(v,float): return floatfmt.format(v)
        return str(v)
    head = "| " + " | ".join([df.index.name or ""] + [str(c) for c in df.columns]) + " |"
    sep  = "|" + "|".join(["---"]*(len(df.columns)+1)) + "|"
    rows = ["| " + " | ".join([str(i)] + [f(v) for v in df.loc[i]]) + " |" for i in df.index]
    return "\n".join([head, sep] + rows)
