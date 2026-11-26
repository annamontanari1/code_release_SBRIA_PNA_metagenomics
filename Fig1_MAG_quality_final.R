# =============================================================================
# Figura 1 — Qualidade dos 25 MAGs (completude × contaminação, CheckM v1.0.18)
# Projeto: SBRIA PN/A (anammox) — amostra NGS1140_6S
# Fonte dos dados: dados_brutos_figuras.md (Seção 8) + bacteria.csv/archea.csv
# (classificação taxonômica completa, GTDB-Tk) do log de reprodutibilidade;
# rRNA e tRNA de "DRAM annotation/" (genome_stats.tsv, rrnas.tsv, trnas.tsv;
# barrnap e tRNAscan-SE dentro do DRAM).
# Revisão 16/10/2026: classificação MIMAG completa (inclui rRNA e tRNA).
# Ambiente testado: R >= 4.3, ggplot2 >= 3.4, ragg >= 1.2
# Compatível com macOS ARM64 (Apple Silicon): todos os pacotes têm binário
# nativo no CRAN; ragg dispensa XQuartz/Cairo do sistema.
# =============================================================================

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(ggplot2)
})

# ---- 1. Dados (CheckM v1.0.18 p/ compl./contam.; GTDB-Tk p/ filo) -----------
# O filo vem de bacteria.csv/archea.csv (campo "Classification", nível p__),
# não do "Marker Lineage" do CheckM: este último mistura rank taxonômico
# (domain/phylum/class/order/family, dependendo de até onde o CheckM resolveu
# cada bin) e não é uma variável categórica válida para uma única legenda.
mags <- tribble(
  ~bin,      ~phylum,             ~genus_gtdb,          ~compl, ~contam,
  "bin.001", "Actinomycetota",    "RHKX01",              51.38,   5.17,
  "bin.002", "Acidobacteriota",   "f__Fen-336",          87.54,   5.68,
  "bin.003", "Nanoarchaeota",     "f__GW2011-AR1",       64.28,   0.00,
  "bin.004", "Nitrospirota",      "Nitrospira_A",        93.86,   8.93,
  "bin.005", "Planctomycetota",   "f__UBA1924",          67.23,   2.15,
  "bin.006", "Pseudomonadota",    "Contendobacter",      77.86,   5.11,
  "bin.007", "Nitrospirota",      "Nitrospira_F",        92.22,  10.91,
  "bin.008", "Pseudomonadota",    "VBCG01",              97.37,   2.98,
  "bin.009", "Planctomycetota",   "Brocadia",           100.00,   2.75,
  "bin.010", "Actinomycetota",    "UBA10799",            98.82,   0.00,
  "bin.011", "Pseudomonadota",    "Desulfobacillus",     94.85,   0.70,
  "bin.012", "Armatimonadota",    "f__Fimbriimonadaceae",93.98,  13.43,
  "bin.013", "Bacteroidota",      "Kapabacteria",        93.72,   2.33,
  "bin.014", "Pseudomonadota",    "Methylocystis",       68.96,   4.13,
  "bin.015", "Pseudomonadota",    "QUBU01",              56.28,  20.85,
  "bin.016", "Pseudomonadota",    "Sphaerotilus",        99.28,   0.18,
  "bin.017", "Pseudomonadota",    "Thiothrix",           99.00,   1.15,
  "bin.018", "Pseudomonadota",    "VBCG01",              82.38,   4.93,
  "bin.019", "Zixibacteria",      "FEB-12",              67.23,   0.00,
  "bin.020", "Bacteroidota",      "OLB6",                54.08,   0.91,
  "bin.021", "Chloroflexota",     "UBA7227",             93.64,   1.27,
  "bin.022", "Pseudomonadota",    "Methylotetracoccus",  95.54,   5.84,
  "bin.023", "Myxococcota_A",     "MYXO",                79.25,   1.81,
  "bin.024", "Acidobacteriota",   "JAJVIN01",            93.49,   3.42,
  "bin.025", "Chloroflexota",     "Caldilinea",          96.12,   1.09
)

# ---- 2. Classificação MIMAG completa (Bowers et al., 2017) -----------------
# HQ: compl >= 90, contam < 5, rRNA 5S, 16S e 23S presentes (parciais aceitos)
#     e >= 18 tipos de tRNA | MQ: compl >= 50 e contam < 10 | LQ: demais.
# rRNA: presença no genome_stats.tsv do DRAM (barrnap). tRNA: tipos distintos
# de aminoácido no trnas.tsv (tRNAscan-SE), excluindo Undet, Pseudo, SeC, Sup,
# fMet e Ile2 (Ile já contado).
gs <- read.delim("DRAM annotation/genome_stats.tsv", check.names = FALSE,
                 colClasses = "character")
gs$bin <- sub("\\.fasta.*$", "", gs$genome)
gs$all_rrna <- nzchar(trimws(gs$`5S rRNA`)) & nzchar(trimws(gs$`16S rRNA`)) &
               nzchar(trimws(gs$`23S rRNA`))
tr <- read.delim("DRAM annotation/trnas.tsv", check.names = FALSE,
                 colClasses = "character")
tr$bin  <- sub("\\.fasta.*$", "", tr$fasta)
tr$Type <- trimws(tr$Type)
tr <- tr[nzchar(tr$Type) & !tr$Type %in% c("Undet", "Pseudo", "SeC", "Sup",
                                            "fMet", "Ile2"), ]
trna <- aggregate(Type ~ bin, data = tr, FUN = function(x) length(unique(x)))
names(trna)[2] <- "trna_types"
mags <- mags |>
  left_join(gs[, c("bin", "all_rrna")], by = "bin") |>
  left_join(trna, by = "bin") |>
  mutate(trna_types = coalesce(trna_types, 0L))
stopifnot(!anyNA(mags$all_rrna))

mags <- mags |>
  mutate(
    mimag = case_when(
      compl >= 90 & contam < 5 & all_rrna & trna_types >= 18 ~ "High",
      compl >= 50 & contam < 10 ~ "Medium",
      TRUE                      ~ "Low"
    ),
    mimag = factor(mimag, levels = c("High", "Medium", "Low"))
  )

# Checagem de integridade: imprime a contagem e a lista efetiva de bins HQ.
stopifnot(nrow(mags) == 25, !anyDuplicated(mags$bin))
print(count(mags, mimag))
cat("High (MIMAG completo):", paste(mags$bin[mags$mimag == "High"], collapse = ", "), "\n")

# Conferência contra a Table S1 (results/tables_text/TableS1_bins.tsv):
# completude, contaminação e faixa devem bater bin a bin.
ts1 <- read.delim("results/tables_text/TableS1_bins.tsv", check.names = FALSE)
chk <- inner_join(mags, ts1, by = c("bin" = "Bin"))
stopifnot(nrow(chk) == 25,
          isTRUE(all.equal(chk$compl,  chk$`Completeness (%)`)),
          isTRUE(all.equal(chk$contam, chk$`Contamination (%)`)),
          all(as.character(chk$mimag) == chk$`Tier (MIMAG)`),
          all(chk$trna_types == chk$`tRNA types`))

# ---- 3. Filo (GTDB-Tk), único rank taxonômico em todos os 25 bins -----------
# Ordem por frequência decrescente (empates em ordem alfabética); Archaea
# (Nanoarchaeota, 1 bin) por último por ser de domínio distinto.
phylum_order <- c(
  "Pseudomonadota", "Acidobacteriota", "Actinomycetota", "Bacteroidota",
  "Chloroflexota", "Nitrospirota", "Planctomycetota",
  "Armatimonadota", "Myxococcota_A", "Zixibacteria", "Nanoarchaeota"
)
stopifnot(setequal(phylum_order, unique(mags$phylum)))

mags <- mags |>
  mutate(phylum = factor(phylum, levels = phylum_order))

# Paleta qualitativa (Okabe-Ito + complementos Paul Tol), colorblind-safe.
# Planctomycetota em preto: filo da Candidatus Brocadia (bin.009), a bactéria
# anammox de interesse central do reator.
pal <- c(
  "Pseudomonadota"   = "#0072B2",
  "Acidobacteriota"  = "#E69F00",
  "Actinomycetota"   = "#009E73",
  "Bacteroidota"     = "#CC79A7",
  "Chloroflexota"    = "#D55E00",
  "Nitrospirota"     = "#F0E442",
  "Planctomycetota"  = "#000000",
  "Armatimonadota"   = "#56B4E9",
  "Myxococcota_A"    = "#999999",
  "Zixibacteria"     = "#882255",
  "Nanoarchaeota"    = "#44AA99"
)

# Forma = qualidade MIMAG (variável independente do filo; cor e forma não
# codificam mais a mesma coisa).
shp <- c("High" = 21, "Medium" = 24, "Low" = 22)

# ---- 4. Rótulos dos bins-chave (posição fixa, plotmath) ---------------------
# Posições definidas manualmente (reprodutíveis, sem dependência de algoritmo
# de repulsão). Convenção: *Candidatus* em itálico e epíteto em romano;
# sufixos GTDB (_A, _F) em romano.
# Versão final (FigS1): rotulados apenas os bins de interesse do manuscrito
# (004, 007, 008, 009, 011, 016, 017). VBCG01 é nome provisório GTDB
# (alfanumérico), por isso em romano.
# 22/10/2026 (cópia *_v_antiga9.R): nenhum rótulo sobre as linhas de limiar (5%, 10%, 50%, 90%) e rótulos sem fundo, para as
# linhas não terem interrupções. Os rótulos ficam à direita da linha de 10%; a ordem vertical segue a dos pontos, para os traços não
# se cruzarem, e as posições de bin.011/bin.021 desviam os traços dos pontos de bin.013, bin.002 e bin.018. bin.016 tem traço em
# cotovelo (sobe e segue em y = 101,3, acima de bin.017 e bin.009).
key_labels <- tribble(
  ~bin,      ~lab,                                         ~lx,   ~ly,   ~hj, ~cotovelo,
  "bin.016", 'italic("Sphaerotilus")~"(bin.016)"',          10.4, 101.3,  0,  TRUE,
  "bin.009", 'italic("Candidatus")~"Brocadia"~"(bin.009)"', 10.4,  99.5,  0,  FALSE,
  "bin.017", 'italic("Thiothrix")~"(bin.017)"',             10.4,  97.7,  0,  FALSE,
  "bin.008", '"VBCG01 (bin.008)"',                          10.4,  95.9,  0,  FALSE,
  "bin.004", 'italic("Nitrospira")*"_A"~"(bin.004)"',       14.6,  92.4,  0,  FALSE,
  "bin.007", 'italic("Nitrospira")*"_F"~"(bin.007)"',       13.0,  88.0,  0,  FALSE,
  "bin.011", 'italic("Ca.")~"Desulfobacillus (bin.011)"',   10.6,  83.3,  0,  FALSE,
  "bin.021", '"UBA7227 (bin.021)"',                         10.6,  79.3,  0,  FALSE
)
lab_df <- inner_join(mags, key_labels, by = "bin") |>
  # o segmento termina ligeiramente antes do texto
  mutate(sx = lx - ifelse(hj == 0, 0.25, -0.25))
# traços: um segmento reto por rótulo, ou dois (vertical + horizontal) no cotovelo
seg_df <- bind_rows(
  lab_df |> filter(!cotovelo) |> transmute(x = contam, y = compl, xend = sx, yend = ly),
  lab_df |> filter(cotovelo)  |> transmute(x = contam, y = compl, xend = contam, yend = ly),
  lab_df |> filter(cotovelo)  |> transmute(x = contam, y = ly,    xend = sx,     yend = ly))

# ---- 5. Gráfico --------------------------------------------------------------
base_size <- 8   # pt no tamanho final impresso (Elsevier: >= 7 pt)
# Helvetica no macOS (fonte recomendada pela Elsevier); fallback métrico
# equivalente (Liberation Sans) em outros sistemas.
font_family <- Sys.getenv("FIG_FONT",
  if (Sys.info()[["sysname"]] == "Darwin") "Helvetica" else "Liberation Sans")

p <- ggplot(mags, aes(x = contam, y = compl)) +
  # Região HQ (compl >= 90, contam < 5)
  annotate("rect", xmin = -Inf, xmax = 5, ymin = 90, ymax = Inf,
           fill = "#1B9E77", alpha = 0.10) +
  # Limiares HQ (sólido) e MQ (tracejado)
  geom_hline(yintercept = 90, linewidth = 0.35, colour = "grey30") +
  geom_vline(xintercept = 5,  linewidth = 0.35, colour = "grey30") +
  geom_hline(yintercept = 50, linewidth = 0.35, colour = "grey50",
             linetype = "dashed") +
  geom_vline(xintercept = 10, linewidth = 0.35, colour = "grey50",
             linetype = "dashed") +
  # v_antiga8: rótulo "HQ thresholds" removido (a área verde e as linhas de
  # 90%/5% ficam; os limiares vão para a legenda da figura no Word).
  annotate("text", x = 9.7, y = 51.5, label = "Medium limit", hjust = 1,
           family = font_family,
           vjust = 0, size = (base_size - 1) / .pt, colour = "grey40") +
  geom_segment(data = seg_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.25, colour = "grey40") +
  # v_antiga6: rótulos dos bins com fundo branco e sem borda, para que as
  # linhas de referência (5%, 10%, 90%) passem por trás do texto. Desenhados
  # antes dos pontos, então nenhum fundo cobre um ponto.
  geom_text(data = lab_df,
            aes(x = lx, y = ly, label = lab, hjust = hj),
            parse = TRUE, size = (base_size - 1) / .pt,
            family = font_family) +
  geom_point(aes(fill = phylum, shape = mimag),
             size = 2.1, stroke = 0.3, colour = "black") +
  scale_fill_manual(values = pal, name = "Phylum (GTDB-Tk)", drop = FALSE) +
  scale_shape_manual(values = shp, name = "Draft quality",
                     drop = FALSE) +
  scale_x_continuous(name = "Contamination (%)",
                     limits = c(-0.3, 22), breaks = seq(0, 20, 5),
                     expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(name = "Completeness (%)",
                     limits = c(48, 102), breaks = seq(50, 100, 10),
                     expand = expansion(mult = c(0, 0))) +
  guides(fill  = guide_legend(ncol = 1, title.position = "top",
                              order = 1,
                              override.aes = list(shape = 21)),
         shape = guide_legend(ncol = 1, title.position = "top",
                              order = 2,
                              override.aes = list(fill = "grey60"))) +
  theme_classic(base_size = base_size, base_family = font_family) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    axis.line         = element_line(linewidth = 0.35),
    axis.ticks        = element_line(linewidth = 0.35),
    axis.text         = element_text(colour = "black"),
    legend.position   = "right",
    legend.box        = "vertical",
    legend.title      = element_text(size = base_size, face = "bold"),
    legend.text       = element_text(size = base_size - 1,
                                     margin = margin(l = 0.8, unit = "mm")),
    legend.key.size   = unit(3.2, "mm"),
    legend.key.width  = unit(2.6, "mm"),
    legend.spacing.y  = unit(0.5, "mm"),
    legend.margin     = margin(0, 0, 0, 0),
    legend.justification = "left",
    legend.box.just   = "left",
    plot.margin       = margin(3, 4, 2, 2, unit = "mm")
  )

# ---- 6. Exportação (Bioresource Technology / Elsevier) ----------------------
# Legenda à direita (1 coluna cada) em vez de embaixo: a figura passa a exigir
# largura de coluna dupla (~140 mm), não mais os 90 mm de coluna única usados
# quando a legenda ficava embaixo. Arte combinada (linha + texto) = mín. 500
# dpi; exporta-se em 1000 dpi (TIFF LZW) e vetorial (PDF) para a submissão.
out_dir <- "figures/final"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
w_mm <- 140; h_mm <- 95

ggsave(file.path(out_dir, "FigS1_MAG_quality.pdf"), p,
       device = grDevices::pdf, width = w_mm, height = h_mm, units = "mm",
       bg = "white")
ggsave(file.path(out_dir, "FigS1_MAG_quality.tiff"), p,
       device = ragg::agg_tiff, width = w_mm, height = h_mm, units = "mm",
       dpi = 1000, compression = "lzw", bg = "white")
# PNG apenas para pré-visualização e inserção no Word (não submeter)
ggsave(file.path(out_dir, "FigS1_MAG_quality_preview.png"), p,
       device = ragg::agg_png, width = w_mm, height = h_mm, units = "mm",
       dpi = 300, bg = "white")

sessionInfo()  # registrar no log de reprodutibilidade
