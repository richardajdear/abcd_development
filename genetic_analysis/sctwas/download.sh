#!/bin/bash
# Sequential download of SingleBrain eQTL files (Zenodo 14908182) in priority order.
# Resumes partial files; verifies size + md5 against the Zenodo API manifest (raw/zenodo_manifest.txt).
set -u
cd "$(dirname "$0")/raw"
CTS="${CTS:-Ext1 Ext OPC IN OD Ast MG Ext2 Ext6 Ext3 Ext4 Ext5 Ext7 Ext8}"
for ct in $CTS; do
  for kind in top full; do
    f=${ct}_eqtl_${kind}_assoc.tsv.gz
    [[ -e $f.done ]] && continue
    read -r _ size md5 < <(awk -v f=$f '$1==f' zenodo_manifest.txt)
    for try in 1 2 3 4 5 6; do
      have=$(stat -c %s $f 2>/dev/null || echo 0)
      [[ $have -ge $size ]] && break
      curl -sSfL -C - --retry 5 -o $f https://zenodo.org/records/14908182/files/$f || sleep 30
    done
    if [[ "md5:$(md5sum $f | cut -d' ' -f1)" == "$md5" ]]; then touch $f.done; echo "$(date) OK $f"; else echo "$(date) FAIL $f"; fi
  done
done
echo "$(date) ALL DONE"
