# =============================================================================
# Figure 1 (final) — Heatmap organismo x gene do ciclo do N, por nível de
# evidência, com barra lateral de profundidade média.
# Projeto: SBRIA PN/A (anammox) — amostra NGS1140_6S
# Atualização 19/10/2026 (cópia: *_v_antiga7.R): auditoria derivada da Table S2; marcador "\u2021" para
# amoA/amoB do bin.007 (NODE_13511, atribuição indireta); nota "*" de amoC ambíguo removida (amoC do bin.007 é
# comammox, 85,8% vs Ca. N. inopinata). Layout inalterado.
# Atualização 19/10/2026-b (cópia: *_v_antiga8.R): linha nova "Nitrospira (NODE_2432, unbinned)" (nxrA 85,8%;
# 43,0x); células de desnitrificação das unidades não binadas voltam a "not tested"; altura 185 -> 200 mm (uma linha
# a mais, mesmo tamanho de célula).
# Base: Fig5_heatmap_pathways.R (cópia de segurança: *_v_antiga.R). Mantidos
# tema (theme_minimal, base 8 pt, Helvetica), eixo x no topo a 30°, tiles
# brancos de 0,7, legenda à direita e dimensões 190 x 165 mm.
# Mudanças autorizadas pela autora (15/10/2026): cor contínua -> 4 níveis
# categóricos (+ "not tested"), paleta salmão/turquesa segura para daltônicos,
# painel lateral de profundidade. Caption interna removida (legenda no Word).
#
# Dados (nada digitado à mão):
#   - evidência por célula: results/tables_text/Fig1_evidence_audit.tsv
#     (gerado por results/tables_text/make_Fig1_audit.py)
#   - profundidade dos bins: results/tables_text/TableS1_bins.tsv
#   - profundidade de NODE_4969 e dos 9 contigs anammox:
#     results/contig_check/asm_coverage.tsv (samtools coverage, meandepth)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ---- 1. Dados -----------------------------------------------------------
aud <- read.delim("results/tables_text/Fig1_evidence_audit.tsv",
                  check.names = FALSE, colClasses = "character")

genes <- c("amoA", "amoB", "amoC", "nxrA", "hzsA", "hzsB", "hzsC", "hdh",
           "nirS", "narG/nxrA family", "napA", "norB", "norC", "nosZ")
units <- c("Nitrosomonas (NODE_4969)", "Nitrospira_A (bin.004)",
           "Nitrospira_F (bin.007)", "Nitrospira (NODE_2432, unbinned)", "Ca. Brocadia (bin.009)",
           "Anammox contigs (unassigned)", "bin.011", "bin.016", "bin.008",
           "bin.017")
# Rótulos (plotmath): gênero em itálico; bin.007 sinalizado pela contaminação.
unit_lab <- c(
  "Nitrosomonas (NODE_4969)"     = 'atop(italic("Nitrosomonas"), "(NODE_4969, unbinned)")',
  "Nitrospira_A (bin.004)"       = 'italic("Nitrospira")*"_A (bin.004)"',
  # três rótulos mais largos em duas linhas (atop), para caber em 190 mm a 10 pt
  "Nitrospira_F (bin.007)"       = 'atop(italic("Nitrospira")*"_F (bin.007;", "10.9% contamination)")',
  "Nitrospira (NODE_2432, unbinned)" = 'italic("Nitrospira")~"(NODE_2432, unbinned)"',
  "Ca. Brocadia (bin.009)"       = 'italic("Ca.")~"Brocadia (bin.009)"',
  # 21/10/2026 (cópia *_v_antiga10.R): rótulo alinhado à legenda do manuscrito
  "Anammox contigs (unassigned)" = 'atop("Unbinned contigs", "("*italic("Ca.")~"Brocadia-affiliated)")',
  "bin.011"                      = 'italic("Ca.")~"Desulfobacillus (bin.011)"',
  "bin.016"                      = 'italic("Sphaerotilus")~"(bin.016)"',
  "bin.008"                      = '"VBCG01 (bin.008)"',
  "bin.017"                      = 'italic("Thiothrix")~"(bin.017)"'
)
gene_lab <- setNames(ifelse(genes == "narG/nxrA family",
                            'italic("narG/nxrA")~"family"',
                            sprintf('italic("%s")', genes)), genes)

lev <- c("supported" = "BLAST: identity 40%+, coverage 70%+",
         "partial"   = "BLAST: coverage below 70%",
         "kofam"     = "KOfam call only",
         "not detected" = "Not detected",
         "not tested"   = "Not searched or not applicable")

cells <- expand_grid(unit = units, gene = genes) |>
  left_join(aud |> select(unit, gene, evidence, identity_pct, note),
            by = c("unit", "gene")) |>
  mutate(
    evidence = replace_na(evidence, "not detected"),  # nxrA fora de Nitrospira
    ev = factor(lev[evidence], levels = lev),
    # "*" = nota de afiliação ambígua registrada na auditoria (amoC de bin.007)
    ambig = grepl("ambiguous amo/pmo", replace_na(note, "")),
    # "†" = identidade contra a NxrA de Ca. K. stuttgartiensis (bin.009, 16/10/2026)
    kuen = grepl("kustd1700", replace_na(note, "")),
    # "\u2021" = atribuição indireta ao bin.007 (amoA/amoB de NODE_13511, contig não binado; 19/10/2026)
    indir = grepl("indirect attribution", replace_na(note, "")),
    lab = ifelse(evidence %in% c("supported", "partial"),
                 paste0(identity_pct, ifelse(ambig, "*", ""), ifelse(kuen, "\u2020", ""),
                        ifelse(indir, "\u2021", "")), ""),
    unit = factor(unit, levels = rev(units)),
    gene = factor(gene, levels = genes)
  )
stopifnot(nrow(cells) == length(units) * length(genes))

# Profundidade média
ts1 <- read.delim("results/tables_text/TableS1_bins.tsv", check.names = FALSE)
cov <- read.delim("results/contig_check/asm_coverage.tsv", check.names = FALSE) |>
  mutate(node = sub("_length.*", "", `#rname`))
anammox_nodes <- paste0("NODE_", c(5668, 5312, 7170, 17695, 29475, 47442,
                                    20914, 24842, 23806))
amx <- range(cov$meandepth[cov$node %in% anammox_nodes])
stopifnot(sum(cov$node %in% anammox_nodes) == 9)
bin_of <- c("Nitrospira_A (bin.004)" = "bin.004", "Nitrospira_F (bin.007)" = "bin.007",
            "Ca. Brocadia (bin.009)" = "bin.009", "bin.011" = "bin.011",
            "bin.016" = "bin.016", "bin.008" = "bin.008", "bin.017" = "bin.017")
depth <- tibble(unit = units) |>
  mutate(
    d = case_when(
      unit == "Nitrosomonas (NODE_4969)" ~ cov$meandepth[cov$node == "NODE_4969"],
      unit == "Nitrospira (NODE_2432, unbinned)" ~ cov$meandepth[cov$node == "NODE_2432"],
      unit %in% names(bin_of) ~
        ts1$`Mean depth (x)`[match(bin_of[unit], ts1$Bin)],
      TRUE ~ NA_real_),
    lo = ifelse(is.na(d), amx[1], NA), hi = ifelse(is.na(d), amx[2], NA),
    # 22/10/2026 (cópia *_v_antiga11.R): intervalo com travessão (en dash, U+2013)
    lab = ifelse(is.na(d), sprintf("%.1f\u2013%.1f\u00d7", floor(amx[1] * 10) / 10, amx[2]),
                 sprintf("%.1f\u00d7", d)),
    xlab = ifelse(is.na(d), amx[2], d),
    unit = factor(unit, levels = rev(units))
  )

# ---- 2. Gráfico -----------------------------------------------------------
base_size <- 8
# Tamanhos fora das células (ajuste 16/10/2026, v_antiga5): eixos 10 pt,
# painel de profundidade 9 pt, legenda 10/9 pt, notas 8 pt. Os números dentro
# das células: 7 pt (v_antiga6; antes base_size - 1.5 = 6.5 pt).
sz_cell <- 7; sz_axis <- 10; sz_depth <- 9; sz_leg_title <- 10; sz_leg_text <- 9; sz_note <- 8
font_family <- Sys.getenv("FIG_FONT",
  if (Sys.info()[["sysname"]] == "Darwin") "Helvetica" else "Liberation Sans")

# Turquesa (BLAST) x salmão (KOfam): contraste de matiz e de luminosidade,
# distinguíveis em deuteranopia/protanopia/tritanopia.
ev_pal <- setNames(c("#00897B", "#9ADBD2", "#F4A582", "#EDEDED", "white"), lev)

p_heat <- ggplot(cells, aes(x = gene, y = unit)) +
  # v_antiga6: separação branca 0.7 -> 0.3 (mais área útil para os números a 7 pt)
  geom_tile(aes(fill = ev), colour = "white", linewidth = 0.3) +
  geom_tile(data = filter(cells, evidence == "not tested"),
            fill = "white", colour = "grey75", linewidth = 0.3,
            width = 0.9, height = 0.9, linetype = "22") +
  geom_text(aes(label = lab), size = sz_cell / .pt, family = font_family,
            colour = ifelse(cells$evidence == "supported", "white", "black")) +
  scale_fill_manual(values = ev_pal, drop = FALSE, name = "Evidence") +
  scale_x_discrete(position = "top", labels = parse(text = gene_lab)) +
  scale_y_discrete(labels = parse(text = unit_lab[rev(units)])) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = base_size, base_family = font_family) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    panel.grid        = element_blank(),
    axis.ticks.x      = element_line(colour = "grey50", linewidth = 0.3),
    axis.ticks.length.x = unit(1.5, "mm"),
    axis.text.x       = element_text(angle = 45, hjust = 0, vjust = 0,
                                     colour = "black", size = sz_axis),
    axis.text.y       = element_text(colour = "black", size = sz_axis),
    legend.title      = element_text(size = sz_leg_title, face = "bold"),
    legend.text       = element_text(size = sz_leg_text),
    legend.key.height = unit(4.5, "mm"),
    legend.key.width  = unit(3.2, "mm"),
    legend.key        = element_rect(colour = "grey70", linewidth = 0.3),
    plot.caption      = element_text(hjust = 0, size = sz_note,
                                     colour = "grey20"),
    plot.caption.position = "plot"
  ) +
  guides(fill = guide_legend(ncol = 1))

p_depth <- ggplot(depth, aes(y = unit)) +
  geom_col(aes(x = d), fill = "grey55", width = 0.62, na.rm = TRUE) +
  geom_segment(aes(x = lo, xend = hi, yend = unit), colour = "grey35",
               linewidth = 0.9, na.rm = TRUE) +
  geom_text(aes(x = xlab, label = lab), hjust = -0.12,
            size = sz_depth / .pt, family = font_family) +
  scale_x_continuous(name = "Mean depth (\u00d7)", limits = c(0, 70),
                     breaks = c(0, 30, 60), expand = expansion(mult = c(0, 0))) +
  coord_cartesian(clip = "off") +
  labs(y = NULL) +
  theme_minimal(base_size = base_size, base_family = font_family) +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid       = element_blank(),
    axis.text.y      = element_blank(),
    axis.line.x      = element_line(colour = "grey50", linewidth = 0.3),
    axis.ticks.x     = element_line(colour = "grey50", linewidth = 0.3),
    axis.text.x      = element_text(colour = "black", size = sz_depth),
    axis.title.x     = element_text(size = sz_depth)
  )

# Notas de rodapé em plotmath (gênero e genes em itálico; "Ca." em itálico e o
# nome em romano), 8 pt, alinhadas à esquerda. Substituem o caption em texto.
notas <- c(
  if (any(cells$ambig))
    '"* "*italic("amoC")~"in bin.007: ambiguous"~italic("amo/pmo")~"affiliation."',
  if (any(cells$kuen)) c(
    '"\u2020 "*italic("nxrA")~"in bin.009: identity to the anammox nitrite oxidoreductase of"~italic("Ca.")~"Kuenenia stuttgartiensis (kustd1700);"',
    '"   other"~italic("nxrA")~"values are identities to"~italic("Nitrospira")~"NxrA."'),
  if (any(cells$indir)) c(
    '"\u2021 "*italic("amoA")*","~italic("amoB")~"in bin.007: unbinned contig NODE_13511 (comammox clade A), attributed to bin.007 by depth, GC,"',
    '"   tetranucleotide frequency, MaxBin2 co-binning and read pairs linking it to bin.007 contigs."'))
p_notes <- ggplot() +
  annotate("text", x = 0, y = -seq_along(notas), label = notas, parse = TRUE,
           hjust = 0, vjust = 0.5, size = sz_note / .pt, family = font_family,
           colour = "grey20") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-length(notas) - 0.5, -0.5), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_family = font_family)

# notas como elemento "full" (ocupa a largura toda da figura, não só a do painel)
p <- p_heat + p_depth + guide_area() + wrap_elements(full = p_notes) +
  plot_layout(design = "AB\nCC\nDD", widths = c(14, 2.2),
              heights = c(1, 0.28, 0.055 * length(notas)), guides = "collect") &
  theme(plot.margin = margin(3, 2, 2, 2, unit = "mm"),
        legend.position = "bottom", legend.justification = "left")
p[[2]] <- p[[2]] + theme(plot.margin = margin(3, 16, 2, 1, unit = "mm"))

# ---- 3. Exportação ------------------------------------------------------------
out_dir <- "figures/final"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
w_mm <- 190; h_mm <- 200   # 165 -> 185 mm (fontes maiores, v_antiga5); 185 -> 200 mm (linha NODE_2432, 19/10/2026)

# TIFF LZW, 1000 dpi (mesmo formato da versão original do script).
ggsave(file.path(out_dir, "Fig1_heatmap_SNAD.tiff"), p,
       device = ragg::agg_tiff, width = w_mm, height = h_mm, units = "mm",
       dpi = 1000, compression = "lzw", bg = "white")
# encoding WinAnsi (CP1252): inclui "†" (U+2020); o padrão ISOLatin1 do pdf()
# não tem esse caractere e o trocava por "." no PDF. Mesma fonte (Helvetica).
ggsave(file.path(out_dir, "Fig1_heatmap_SNAD.pdf"), p,
       device = grDevices::pdf, width = w_mm, height = h_mm, units = "mm",
       bg = "white", encoding = "WinAnsi")
ggsave(file.path(out_dir, "Fig1_heatmap_SNAD.png"), p,
       device = ragg::agg_png, width = w_mm, height = h_mm, units = "mm",
       dpi = 300, bg = "white")

sessionInfo()
