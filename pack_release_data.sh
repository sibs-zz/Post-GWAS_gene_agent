#!/usr/bin/env bash
# Build a standalone SoySAGE release bundle (scripts + all data). For maintainers only.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
VERSION="${1:-v1.0.0}"
DIR="SoySAGE_full_${VERSION}"
OUT="$ROOT/${DIR}.tar.gz"
STAGING="$ROOT/$DIR"

RELEASE_FILES=(
  post_gwas_gene_agent_v1.py
  post_gwas_gene_agent.py
  calc_priority.py
  gene_ai_analysis.py
  soysage_common.py
  README.md
  requirements.txt
  requirements-agent.txt
  LICENSE
  TPM_class.txt
  tst.list
  2898zhugao.pheno.clean.mlma.ma.gene.fastbat
  stringtie_gene_314_TPM.txt
  ZH13.gene.iprscn.mod.txt
  ZH13.gene.Pfam.mod.txt
  ZH13.gene.uniprot.plants.simple.txt
  ZH13.GN.txt
  ZH13.GO.txt
  ZH13.KEGG.txt
  STRING_ZH13v2.1_physical_links.tsv
  Soybean_Gene_Embeddings_Full.npy
  Soybean_Gene_Semantic_Profiles_Full.csv
  Gene_Phenotype_Associations_Readable_clean.csv
  ZH13v3_53403_314tissue_gene_TPM.txt
  ZHV3.0_mRNA_53403_20250402_DE.txt
  ZH13v3.GN.txt
  ZH13v3_ZH13v2_Wm82_Ath.txt
  ZHV3.0_mRNA_53403_20250402_gene_GO.txt
  ZHV3.0_mRNA_53403_20250402_gene_KEGG_pathway.txt
  STRING_ZH13v3_physical_links.tsv
)

cd "$ROOT"
missing=()
for f in "${RELEASE_FILES[@]}"; do
  [[ -f "$f" ]] || missing+=("$f")
done
if ((${#missing[@]})); then
  echo "Missing files:" >&2
  printf '  %s\n' "${missing[@]}" >&2
  exit 1
fi

rm -rf "$STAGING"
mkdir "$STAGING"
cp -f "${RELEASE_FILES[@]}" "$STAGING/"
tar -czvf "$OUT" "$DIR"
rm -rf "$STAGING"
ls -lah "$OUT"
echo "Created $OUT ($(tar -tzf "$OUT" | wc -l) entries)"
