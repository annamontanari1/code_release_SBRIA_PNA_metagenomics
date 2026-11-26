#!/usr/bin/env bash
# Secao 17 do README (22/10/2026): barrnap e tRNAscan-SE sobre TODOS os contigs de cada um dos 25 bins (sem o corte de 2.500 pb do DRAM).
# Ambiente: conda "rrna" criado em 21/10/2026 com --platform osx-arm64 (barrnap 1.10.5, tRNAscan-SE 2.0.13, Infernal 1.1.5);
# o ~/miniconda3/.condarc geral tem subdir: osx-64 e nao foi alterado.
# Saidas: results/review_checks/16_rrna_trna/raw/<bin>/ (barrnap.gff, *.rRNA.tbl do cmsearch, trnascan.txt, logs).
# bin.003 (Nanoarchaeota) -> barrnap --kingdom arc e tRNAscan-SE -A; demais -> bac e -B.
set -euo pipefail
ROOT=~/sbria_anammox_shotgun
BINS="$ROOT/das_tool_result (1)/das_tool_output_dir_DASTool_bins"
OUT=$ROOT/results/review_checks/16_rrna_trna/raw
E=~/miniconda3/envs/rrna/bin
export PATH=$E:$PATH
JOBS=${JOBS:-4}; THR=${THR:-2}
mkdir -p "$OUT"
{ echo "barrnap: $($E/barrnap --version 2>&1 | tail -1)"; echo "tRNAscan-SE: $($E/tRNAscan-SE -h 2>&1 | grep -m1 'tRNAscan-SE 2')";
  echo "infernal: $($E/cmsearch -h | sed -n 2p)"; date; } > "$OUT/../versions.txt"

run_bin() {
  b=$1; d="$OUT/$b"; mkdir -p "$d"
  if [ "$b" = "bin.003" ]; then K=arc; T=-A; else K=bac; T=-B; fi
  cp "$BINS/$b.fasta" "$d/$b.fasta"                      # copia local: o barrnap --debug escreve ao lado do FASTA
  ( cd "$d" && barrnap --debug --kingdom $K --threads $THR "$b.fasta" > barrnap.gff 2> barrnap.log )
  tRNAscan-SE $T --thread $THR -Q -o "$d/trnascan.txt" -m "$d/trnascan_stats.txt" "$d/$b.fasta" > "$d/trnascan.log" 2>&1
  rm -f "$d/$b.fasta" "$d/$b.fasta.fai"
  echo "$b done"
}
export -f run_bin; export OUT BINS E THR
ls "$BINS"/bin.*.fasta | sed -E 's|.*/(bin\.[0-9]+)\.fasta|\1|' | xargs -P "$JOBS" -I{} bash -c 'run_bin {}'
