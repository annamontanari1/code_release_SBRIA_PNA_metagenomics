#!/usr/bin/env bash
# Seção 19 do README (23/10/2026): árvore de HzsA v2, com grupo externo Scalindua. Não altera a árvore antiga (12_anammox_contigs/tree_hzsA).
# Conjunto: 12 referências validadas (19_hzsA_v2/hzsA_reference_candidates.tsv, coluna proposed = yes) + 5 HzsA deste estudo
# (mesmas sequências de 12_anammox_contigs/tree_hzsA/in.faa). Pipeline idêntico ao da árvore antiga (classify_anammox_contigs.py):
#   mafft --auto (v7.526) -> trimal -automated1 (1.5.1) -> iqtree3 -m MFP -B 1000 -alrt 1000 (3.1.4, ambiente phylo)
# Únicas diferenças: -o com as quatro Scalindua (só enraíza a árvore de saída; não muda a inferência) e -seed 12345 (reprodutível).
# Duas árvores: main (17 sequências) e sens (sem NODE_17695 e NODE_29475, fragmentos de 115 e 408 aa).
set -euo pipefail
RC=~/sbria_anammox_shotgun/results/review_checks; O=$RC/19_hzsA_v2; E=~/miniconda3/envs/phylo/bin; SK=/opt/homebrew/bin/seqkit   # seqkit 2.x (Homebrew)
OUTG="Scalindua_japonica_WP_096896276v1,Scalindua_rubra_ODS33869v1,Scalindua_brodae_KHE92665v1,Scalindua_sp15L_AEW50032v1"
# rótulos (sem "." nem espaços, no estilo da árvore antiga)
LAB="WP_096896276.1 Scalindua_japonica_WP_096896276v1
ODS33869.1 Scalindua_rubra_ODS33869v1
KHE92665.1 Scalindua_brodae_KHE92665v1
AEW50032.1 Scalindua_sp15L_AEW50032v1
WP_099325812.1 Kuenenia_stuttgartiensis_WP_099325812v1
TLD39960.1 Jettenia_ecosi_TLD39960v1
AEW50031.1 Jettenia_asiatica_AEW50031v1
WP_080325085.1 Brocadia_sapporoensis_WP_080325085v1
GAN34140.1 Brocadia_sinica_GAN34140v1
KKO20072.1 Brocadia_fulgida_KKO20072v1
BBO18370.1 Brocadia_pituitae_BBO18370v1
WP_285442451.1 Brocadia_spAMX2_WP_285442451v1"   # bash 3.2 do macOS: sem arrays associativos
build() {   # $1 = nome; $2 = regex de contigs a excluir (ou vazio)
  d=$O/tree_$1; mkdir -p $d; : > $d/in.faa
  echo "$LAB" | while read a lab; do $SK seq -i $O/hzsA_reference_candidates.faa 2>/dev/null | $SK grep -p "$a" 2>/dev/null \
      | sed "1s/.*/>$lab/" >> $d/in.faa; done
  $SK grep -r -p "^THISSTUDY_" $RC/12_anammox_contigs/tree_hzsA/in.faa 2>/dev/null | { if [ -n "$2" ]; then $SK grep -v -r -p "$2" 2>/dev/null; else cat; fi; } >> $d/in.faa
  echo "$1: $(grep -c '>' $d/in.faa) sequências"
  $E/mafft --auto --thread 8 $d/in.faa > $d/aln.faa 2> $d/mafft.log
  $E/trimal -in $d/aln.faa -out $d/aln_trim.faa -automated1
  $E/iqtree3 -s $d/aln_trim.faa -m MFP -B 1000 -alrt 1000 -T AUTO -seed 12345 -o "$OUTG" --prefix $d/hzsA_v2_$1 -redo -quiet
}
build main ""
build sens "NODE_(17695|29475)_"
