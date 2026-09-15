import pandas as pd, sys
b=pd.read_parquet(f"{sys.argv[1]}/fits/blups.parquet"); print(b.shape, b.label.unique(), "subjects", b.subject.nunique())
