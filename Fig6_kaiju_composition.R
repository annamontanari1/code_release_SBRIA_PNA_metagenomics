# =============================================================================
# Figura 6 — Composição da comunidade por leitura (Kaiju, nível de genero)
# Projeto: SBRIA PN/A (anammox) — amostra NGS1140_6S
# Amostragem: dia 595, fim da Fase VI operacional, aeração DESLIGADA no
# momento da coleta (fase anoxica do ciclo intermitente).
# Fonte dos dados: Kaiju/kaiju_report/NGS1140_6S_trimmed_reads-genus.kaijuReport
# (2157 generos; percentuais de leitura classificada, nao de leitura total).
#
# Escopo: os oito organismos-chave da Table 1 (via gene-level MAG), reduzidos
# aos que tem correspondencia direta em nivel de genero no Kaiju. VBCG01
# (bin.008/bin.018) nao entra: e um nome de linhagem GTDB sem taxonomia de
# genero formal, nao aparece como tal no banco de referencia do Kaiju.
#
# DUAS RESSALVAS DE DADOS (mesma logica ja documentada no texto do artigo,
# secao "Unresolved findings and methodological limitations"):
#   1. Brocadia (bin.009) NAO aparece como genero proprio no relatorio Kaiju
#      (0 leituras). Usa-se Candidatus Kuenenia (1,838%) como proxy, pelo
#      mesmo argumento ja no texto: o mapeamento de 61,87% dessas leituras
#      contra bin.009 sugere que a maioria se origina de Brocadia e foi
#      atribuida a Kuenenia por vies de banco de referencia (K. stuttgartiensis
#      mais representada que Brocadia no banco do Kaiju). Fatia rotulada
#      "Candidatus Brocadia (via Ca. Kuenenia, vies de referencia)".
#   2. Desulfobacillus (bin.011) tambem NAO aparece como genero proprio
#      (0 leituras; a unica entrada proxima e Acididesulfobacillus, 0,008%,
#      genero nao relacionado). Consistente com a descricao de Okubo e Takami
#      (2021) do genoma como um "imperfect denitrifier" divergente o
#      suficiente para escapar da classificacao padrao por leitura. Fatia
#      mantida em 0% com rotulo e chamada de nota, nao omitida nem substituida
#      por proxy (ao contrario de Brocadia, nao ha um candidato de
#      classificacao errada documentado para este caso).
# =============================================================================

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(ggplot2)
})

# ---- 1. Dados (Kaiju, nivel de genero; ver caminho do arquivo no cabecalho) --
comp <- tribble(
  ~organismo,                                                  ~pct,       ~nota,
  "Candidatus Brocadia (via Ca. Kuenenia, vies de referencia)", 1.838476,  "proxy",
  "Nitrospira",                                                 13.448339, "",
  "Nitrosomonas",                                                0.783689, "",
  "Desulfobacillus",                                             0.000000, "nao_detectado",
  "Methylocystis",                                               2.370893, "",
  "Thiothrix",                                                   1.677262, "",
  "Methanothrix",                                                1.160403, ""
)

outros_pct <- 100 - sum(comp$pct)
comp <- comp |>
  add_row(organismo = "Others", pct = outros_pct, nota = "") |>
  mutate(
    organismo = factor(organismo, levels = rev(organismo)),
    label = sprintf("%.2f%%", pct)
  )

stopifnot(abs(sum(comp$pct) - 100) < 1e-6)

# ---- 2. Paleta (Okabe-Ito, colorblind-safe; consistente com Figura 1) -------
pal <- c(
  "Candidatus Brocadia (via Ca. Kuenenia, vies de referencia)" = "#000000",
  "Nitrospira"       = "#F0E442",
  "Nitrosomonas"     = "#0072B2",
  "Desulfobacillus"  = "#D55E00",
  "Methylocystis"    = "#009E73",
  "Thiothrix"        = "#CC79A7",
  "Methanothrix"     = "#56B4E9",
  "Others"           = "#CCCCCC"
)

# ---- 3. Grafico (donut; coord_polar sobre barra empilhada) -------------------
base_size <- 8
font_family <- Sys.getenv("FIG_FONT",
  if (Sys.info()[["sysname"]] == "Darwin") "Helvetica" else "Liberation Sans")

p <- ggplot(comp, aes(x = 2, y = pct, fill = organismo)) +
  geom_col(colour = "white", linewidth = 0.4, width = 1) +
  coord_polar(theta = "y") +
  xlim(0.4, 2.5) +
  scale_fill_manual(values = pal, name = NULL,
                     breaks = rev(levels(comp$organismo))) +
  labs(caption = paste0(
    "Dia 595 (fim da Fase VI); aeracao desligada no momento da coleta.\n",
    "Percentuais de leitura classificada (Kaiju, nivel de genero); genero\n",
    "\"Others\" agrega os demais 2150 generos classificados.\n",
    "Candidatus Brocadia: sem genero proprio no Kaiju (0 leituras);\n",
    "fatia usa Ca. Kuenenia como proxy por vies de banco de referencia\n",
    "(ver Unresolved findings, no texto). Desulfobacillus: tambem sem\n",
    "genero proprio no Kaiju (0 leituras); mantido como fatia nula, sem\n",
    "proxy documentado, consistente com o genoma \"imperfect denitrifier\"\n",
    "descrito por Okubo e Takami (2021)."
  )) +
  theme_void(base_size = base_size, base_family = font_family) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.position   = "right",
    legend.text       = element_text(size = base_size - 1),
    plot.caption      = element_text(hjust = 0, size = base_size - 2,
                                     colour = "grey30", lineheight = 1.15),
    plot.margin       = margin(3, 4, 18, 2, unit = "mm")
  )

# ---- 4. Exportacao (Bioresource Technology / Elsevier) -----------------------
out_dir <- "figures"
dir.create(out_dir, showWarnings = FALSE)
w_mm <- 140; h_mm <- 110

ggsave(file.path(out_dir, "Fig6_kaiju_composition.tiff"), p,
       device = ragg::agg_tiff, width = w_mm, height = h_mm, units = "mm",
       dpi = 1000, compression = "lzw", bg = "white")
ggsave(file.path(out_dir, "Fig6_kaiju_composition.pdf"), p,
       device = grDevices::pdf, width = w_mm, height = h_mm, units = "mm",
       bg = "white")
ggsave(file.path(out_dir, "Fig6_kaiju_composition_preview.png"), p,
       device = ragg::agg_png, width = w_mm, height = h_mm, units = "mm",
       dpi = 300, bg = "white")

sessionInfo()
