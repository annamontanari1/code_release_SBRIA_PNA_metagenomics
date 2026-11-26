# =============================================================================
# Figura 4 — Detecção padrão (DRAM/KOfam) vs. BLAST/tBLASTn dirigido
# Atualização 19/10/2026 (cópia: *_v_antiga6.R): valores lidos da Table S2 via make_FigS2_values.py; 22 barras
# (novas: amoA/amoB comammox de NODE_13511 atribuídos ao bin.007 com "\u2021", amoC comammox de NODE_2578 e nosZ
# do bin.017); altura 145 -> 177 mm na mesma proporção das barras (layout e fontes inalterados).
# Projeto: SBRIA PN/A (anammox) — amostra NGS1140_6S
# Fonte: logs/reproducibility_log_NGS1140_6S.md (Seções 10, 11, 18, 20, 23,
# 31) + compilado_SBRIA_anammox_NGS1140_6S.md (quadro metabólico, linha 131).
#
# NOTA DE CORREÇÃO (já aplicada): a linha "hao — Nitrosomonas (bin.014)" foi
# removida. bin.014 foi reidentificado como Methylocystis, não Nitrosomonas;
# aquele hit é cross-reatividade pmoCAB/hao-family, não confirmação real de
# Nitrosomonas (ver Seção 31 do log, "achado novo: norB de bin.014 é
# genuíno... nunca havia sido verificada sua classificação taxonômica própria").
#
# NOTA DE RASTREABILIDADE (consolidação aplicada nesta versão): checagem
# linha a linha contra o log encontrou valores sem fonte específica
# rastreável, já corrigidos aqui:
#   - hzsA(81,2%)/hzsB(80,9%)/hzsC(86,2%)/hdh(87,5%) individuais: os três
#     primeiros pertencem ao contig NODE_1 da Seção 16 (montagem exploratória
#     e NÃO CONCLUSIVA dos reads "Kuenenia" não mapeados), não à confirmação
#     principal de bin.009 (Seção 10, que reporta só a faixa agregada
#     78,4-94,6%, sem quebra por gene); hdh=87,5% não tem fonte localizada.
#     -> consolidado em 1 linha ("hzsA/B/C, hdh, operon"). Valor usado = PONTO
#     MÉDIO da faixa agregada da Seção 10 (78,4-94,6% -> 86,5%), não a média
#     dos 4 valores da Seção 16 (que são de outro bin/investigação e dariam
#     83,95% — não usar, não é a mesma fonte). Mesmo critério da Figura 2
#     ("Recovery evidence") e da Figura 5 (célula "Anammox").
#   - nxrA/nxrB de bin.007: a narrativa do log/compilado só menciona uma faixa
#     combinada com bin.004 ("84-99%", compilado linha 131), sem valor
#     específico de bin.007. Buscando no dado BRUTO (não na narrativa) —
#     results/targeted_search/tblastn_nxr_full_assembly.tsv, filtrado pelos
#     contigs do próprio bin.007.fasta — os dois hits de nxrA de bin.007 são
#     84,875% (NODE_4141, 1124/1124 aa) e 87,100% (NODE_348, 1124/1124 aa),
#     e-value 0,0 em ambos. Média = 86,0% — valor real e específico de
#     bin.007, não mais uma faixa "não específica de bin". Isso substitui a
#     versão anterior desta figura, que desenhava bin.007 como barra
#     "flutuante" (84-99%) por acreditar que não havia valor próprio do bin.
#
# CORREÇÃO DE DIRECIONALIDADE (dsrA/dsrB, bin.011) — marcada com "^":
# investigação filogenética para a Figura 2 (árvore dsrA, 7 táxons: os 4
# MAGs deste projeto com dsrA + Desulfovibrio/Archaeoglobus/Allochromatium)
# mostrou bin.011 filogeneticamente muito mais próximo do dsrA OXIDATIVO
# (Allochromatium, distância ML 0,52-0,54) do que de qualquer referência
# redutora (Archaeoglobus/Desulfovibrio, 1,03-1,30). Confirmado de forma
# independente na literatura primária de descrição da espécie (Okubo et al.
# 2021, MicrobiologyOpen, já citado no projeto a propósito do nirS): o
# dsrAB de Ca. Desulfobacillus denitrificans é classificado pelos autores
# como do tipo OXIDATIVO (via reversa do DSR, rDSR: sulfeto -> S0 -> sulfito
# -> sulfato), não redução dissimilatória clássica — e a espécie pertence à
# ordem Rhodocyclales, fora do grupo taxonômico clássico de redutores de
# sulfato (concentrado em Desulfobacterota/Deltaproteobacteria). Duas linhas
# de evidência independentes (filogenia local + literatura da espécie)
# convergem na mesma direção. Os 100% de identidade permanecem válidos como
# confirmação de PRESENÇA do gene e identidade taxonômica; o que muda é a
# interpretação funcional (oxidativo, não redutor) — refletida também na
# Figura 5 (bin.011 migrou de "Sulfate reduction" para "Sulfur oxidation
# (rDSR)", coluna própria, distinta de "Sulfur oxidation (Sox)": rDSR e Sox
# são vias bioquímicas diferentes, mesmo erro de mistura já evitado no
# nxrA/nxrB desta figura).
# =============================================================================

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(forcats)
  library(ggplot2)
})

# ---- 1. Dados ----------------------------------------------------------------
# Versão final (FigS2): valores lidos de results/tables_text/FigS2_values.tsv,
# gerado por results/tables_text/make_FigS2_values.py a partir dos arquivos
# brutos de results/targeted_search/ (melhor hit = maior bit score; arquivo e
# linha de origem de cada valor registrados na tabela). Status DRAM: Yes = KO
# atribuído ao gene do hit; Partial = só Pfam PF00384, sem KO; No = sem
# anotação ou contig fora de bin consolidado.
dram_vs_blast <- read.delim("results/tables_text/FigS2_values.tsv",
                            check.names = FALSE) |>
  # 22/10/2026 (cópia *_v_antiga7.R): rótulo das barras hzs/hdh alinhado ao texto (contigs afiliados a Ca. Brocadia);
  # só o rótulo muda, a tabela de valores não
  mutate(organism = sub("^Anammox contig \\(", "Brocadia-affiliated contig (", organism)) |>
  mutate(
    identidade = identity_pct,
    detectado_dram = factor(dram, levels = c("No", "Partial", "Yes")),
    id = fct_reorder(paste(axis_label, organism, sep = " | "), identidade),
    # Hit parcial: alinhamento < 70% do comprimento real da referência
    # (ref_len lido dos fastas de results/targeted_search/). Mesmo critério
    # para todos os genes; o comprimento do alinhamento vai ao lado da barra.
    parcial = as.logical(partial),
    # Anotações extras (16/10/2026): segunda cópia de comprimento total
    # (barra = cópia de maior bit score), gene partido em ORFs adjacentes
    # (cobertura somada) e hit abaixo de 40% de identidade (não suportado).
    copia = sub("^NODE_\\d+ ", "", other_full_length_copies),
    extra = paste0(
      ifelse(parcial, sprintf("; %d aa", aln_len), ""),
      ifelse(!is.na(copia) & copia != "", paste0("; second copy ", copia), ""),
      ifelse(grepl("split gene", note), "; 2 ORFs", ""),
      ifelse(blast_status == "not supported" & grepl("identity", note),
             "; identity < 40%", "")),
    rotulo = sprintf("%s (%.1f%%%s)", organism, identidade, extra)
  )
stopifnot(nrow(dram_vs_blast) == 22, !anyNA(dram_vs_blast$identidade),
          !anyNA(dram_vs_blast$parcial))
eixo_lab <- setNames(dram_vs_blast$axis_label, as.character(dram_vs_blast$id))

# ---- Rótulos em plotmath (v_antiga5) ------------------------------------------
# Gene: primeira palavra em itálico, descritor em romano. Organismo: gênero em
# itálico (sufixo GTDB e códigos em romano); "Ca." em itálico e o nome em
# romano; códigos (VBCG01, Anammox contig, NODE_..., bin.0xx) em romano. Os
# valores entre parênteses ficam exatamente como antes, dentro de aspas.
pm_q <- function(x) paste0('"', gsub('"', '\\"', x), '"')
pm_gene <- function(lbl) {
  g <- sub(" .*$", "", lbl); rest <- sub("^\\S+ ?", "", lbl)
  ifelse(nzchar(rest), paste0("italic(", pm_q(g), ")~", pm_q(rest)),
         paste0("italic(", pm_q(g), ")"))
}
generos <- c("Nitrosomonas", "Nitrospira", "Sphaerotilus", "Thiothrix")
pm_org <- function(o) vapply(o, function(x) {
  if (startsWith(x, "Ca. ")) return(paste0('italic("Ca.")~', pm_q(sub("^Ca\\. ", "", x))))
  g <- sub("[_ ].*$", "", x)
  if (g %in% generos) {
    rest <- substring(x, nchar(g) + 1)
    if (startsWith(rest, "_")) return(paste0("italic(", pm_q(g), ")*", pm_q(rest)))
    return(paste0("italic(", pm_q(g), ")~", pm_q(sub("^ ", "", rest))))
  }
  pm_q(x)
}, character(1), USE.NAMES = FALSE)
dram_vs_blast <- dram_vs_blast |>
  mutate(valores = sprintf("(%.1f%%%s)", identidade, extra),
         rotulo_pm = paste0(pm_org(organism), "~", pm_q(valores)))
stopifnot(all(paste(dram_vs_blast$organism, dram_vs_blast$valores) == dram_vs_blast$rotulo))
eixo_pm <- setNames(pm_gene(dram_vs_blast$axis_label), as.character(dram_vs_blast$id))
invisible(parse(text = c(eixo_pm, dram_vs_blast$rotulo_pm)))  # falha cedo se inválido

# Hachura desenhada à mão (ggpattern não instalado): diagonais dentro de cada
# barra parcial, recortadas a [0, identidade]; a cor de preenchimento continua
# sendo a do status DRAM.
larg <- 0.62; passo <- 3; sobe <- 4
hachura <- dram_vs_blast |>
  filter(parcial) |>
  mutate(pos = as.numeric(id)) |>
  rowwise() |>
  reframe(pos = pos, top = identidade, y0 = seq(-sobe, identidade, by = passo)) |>
  mutate(t0 = pmax(0, (0 - y0) / sobe), t1 = pmin(1, (top - y0) / sobe)) |>
  filter(t1 > t0) |>
  mutate(x  = pos - larg / 2 + larg * t0, xend = pos - larg / 2 + larg * t1,
         y  = y0 + sobe * t0,             yend = y0 + sobe * t1,
         tipo = "Partial alignment (<70% of reference length)")

# ---- 2. Gráfico ---------------------------------------------------------------
base_size <- 8
font_family <- Sys.getenv("FIG_FONT",
  if (Sys.info()[["sysname"]] == "Darwin") "Helvetica" else "Liberation Sans")

# Mesma paleta colorblind-safe (Okabe-Ito) das Figuras 1-2: vermillion/dourado/
# verde-azulado, progressão "não detectado -> parcial -> confirmado".
detect_pal <- c("No" = "#D55E00", "Partial" = "#E69F00", "Yes" = "#009E73")

p <- ggplot(dram_vs_blast, aes(x = id, y = identidade, fill = detectado_dram)) +
  geom_col(width = 0.62, colour = "black", linewidth = 0.3) +
  geom_segment(data = hachura, inherit.aes = FALSE,
               aes(x = x, xend = xend, y = y, yend = yend, linetype = tipo),
               colour = "black", linewidth = 0.25, key_glyph = "abline") +
  scale_linetype_manual(values = "solid", name = NULL) +
  geom_text(aes(label = rotulo_pm), parse = TRUE, hjust = -0.03,
            size = (base_size - 1) / .pt, family = font_family) +
  coord_flip(clip = "off") +
  scale_x_discrete(labels = function(b) parse(text = eixo_pm[b])) +
  scale_fill_manual(values = detect_pal, name = "Detected by DRAM/KOfam?",
                    drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0.78)),
                     name = "Amino acid identity, targeted BLAST/tBLASTn (%)") +
  labs(x = NULL) +
  theme_classic(base_size = base_size, base_family = font_family) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    axis.line         = element_line(linewidth = 0.35),
    axis.ticks        = element_line(linewidth = 0.35),
    axis.text         = element_text(colour = "black"),
    axis.text.y       = element_text(size = base_size),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.title      = element_text(size = base_size, face = "bold"),
    legend.text       = element_text(size = base_size - 1),
    legend.key.size   = unit(3.2, "mm"),
    legend.spacing.x  = unit(1.5, "mm"),
    legend.margin     = margin(0, 0, 0, 0),
    plot.caption      = element_text(hjust = 0, size = base_size - 2,
                                     colour = "grey30", lineheight = 1.15),
    legend.justification = "left",
    legend.box.just   = "left",
    plot.margin       = margin(3, 4, 2, 2, unit = "mm")
  )

# ---- 3. Exportação (Bioresource Technology / Elsevier) ----------------------
out_dir <- "figures/final"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
w_mm <- 165; h_mm <- 177   # 145 mm para 18 barras -> 177 mm para 22 (19/10/2026)

# TIFF LZW, 1000 dpi (mesmo formato da versão original do script).
ggsave(file.path(out_dir, "FigS2_DRAM_vs_BLAST.tiff"), p,
       device = ragg::agg_tiff, width = w_mm, height = h_mm, units = "mm",
       dpi = 1000, compression = "lzw", bg = "white")
# encoding WinAnsi (CP1252): inclui "\u2021" (19/10/2026), como na Figura 1.
ggsave(file.path(out_dir, "FigS2_DRAM_vs_BLAST.pdf"), p,
       device = grDevices::pdf, width = w_mm, height = h_mm, units = "mm",
       bg = "white", encoding = "WinAnsi")
ggsave(file.path(out_dir, "FigS2_DRAM_vs_BLAST.png"), p,
       device = ragg::agg_png, width = w_mm, height = h_mm, units = "mm",
       dpi = 300, bg = "white")

sessionInfo()
