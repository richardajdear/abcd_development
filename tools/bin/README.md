# tools/bin — local binaries (gitignored)

| path | what | source |
|:--|:--|:--|
| `magma` (symlink) | MAGMA v1.10, **Linux static** build used on CSD3 | `~/Git/AHBA/magma/magma` (also `magma_v1.10_static.zip`) — will not run on macOS |
| `magma_mac/magma` | MAGMA v1.10, macOS x86_64 build (runs under Rosetta) | https://cncr.nl/research/magma/ → "Mac OS (64 bits)"; the surfsara link redirects to `vu.data.surf.nl` |
| `NCBI37.3.gene.loc` (symlink) | gene locations, build 37 | `~/Git/AHBA/magma/gene_locations/` |

`ahba_pls/code/07_magma_gene_property.py` uses `magma_mac/magma` if it exists, else `magma`.
