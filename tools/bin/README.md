# tools/bin — local binaries (gitignored)

| path | what | source |
|:--|:--|:--|
| `magma` (symlink) | MAGMA v1.10, **Linux static** build used on CSD3 | `~/Git/AHBA/magma/magma` (also `magma_v1.10_static.zip`) — will not run on macOS |
| **`magma_src/magma`** | **MAGMA v1.10, native macOS arm64, built from source — used first** | source zip from https://cncr.nl/research/magma/ ("Source code", share `1OOi7bxLWef0GwY`); build recipe below. Reproduces the x86 build's `magma_all_marginal.tsv` exactly (18 cells, Δβ = Δp = 0) |
| `magma_mac/magma` | MAGMA v1.10, macOS x86_64 build (runs under Rosetta — **not installed on this machine as of 2026-09-24**, so it fails with "Bad CPU type") | https://cncr.nl/research/magma/ → "Mac OS (64 bits)"; the surfsara link redirects to `vu.data.surf.nl` |
| `NCBI37.3.gene.loc` (symlink) | gene locations, build 37 | `~/Git/AHBA/magma/gene_locations/` |

The `ahba_pls/code/*` MAGMA scripts use the first of `magma_src/magma`, `magma_mac/magma`, `magma` that exists.

Build recipe (Apple clang 21; the default toolchain cannot find the libc++ headers, and the bundled
Eigen's NEON path is 32-bit ARM assembly, hence the three flags):

```sh
cd tools/bin/magma_src && unzip magma_src.zip
S=/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk
make -j8 CXX="clang++ -isysroot $S -nostdinc++ -isystem $S/usr/include/c++/v1 -DEIGEN_DONT_VECTORIZE -DEIGEN_DONT_ALIGN"
```
