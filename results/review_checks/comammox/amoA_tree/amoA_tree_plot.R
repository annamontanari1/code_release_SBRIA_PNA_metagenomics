# =============================================================================
# Árvore de AmoA (comammox clados A/B, AOB beta/gama, PmoA) — revisão comammox
# Projeto: SBRIA PN/A (anammox) — amostra NGS1140_6S
# Fonte: IQ-TREE 3.1.3 (ModelFinder, BIC: Q.YEAST+F+G4; UFBoot 1000 + SH-aLRT
# 1000), alinhamento MAFFT L-INS-i aparado com trimAl v1.5.rev1 -automated1
# (329 -> 274 colunas). Estilo de Fig3_nirS_tree_final.R.
# Raiz: clado PmoA (Methylocystis, Methylococcus, Methylotetracoccus) +
# AmoA de Nitrosococcus (gama), usado apenas para orientar a leitura.
# Saída somente em results/review_checks/comammox/amoA_tree/ (figuras do artigo
# não são regeneradas).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ape)
  library(ggtree)
})

dir_tree <- path.expand("~/sbria_anammox_shotgun/results/review_checks/comammox/amoA_tree")
tree <- read.tree(file.path(dir_tree, "amoA_tree.treefile"))

# ---- 1. Enraizamento --------------------------------------------------------
outg <- grep("^pmoA\\||^AOB_gamma\\||bin.014_NODE_2726_PmoA", tree$tip.label, value = TRUE)
tree <- root(tree, outgroup = outg, resolve.root = TRUE, edgelabel = TRUE)
tree <- ladderize(tree, right = FALSE)
tree$node.label[tree$node.label %in% c("Root", "root")] <- ""

# ---- 2. Suporte dos nós de interesse (impresso no console) -------------------
sup <- function(tips) {
  n <- getMRCA(tree, tips)
  lab <- tree$node.label[n - Ntip(tree)]
  mono <- is.monophyletic(tree, tips)
  sprintf("%-55s monophyletic=%s support(SH-aLRT/UFBoot)=%s", paste(length(tips), "tips"), mono, lab)
}
tl <- tree$tip.label
grp <- list(
  "NODE_13511 + Ca. N. nitrosa"                  = c(grep("NODE_13511", tl, value = TRUE), grep("nitrosa", tl, value = TRUE)),
  "NODE_13511 + nitrosa + CMX1"                  = c(grep("NODE_13511|nitrosa|CMX1", tl, value = TRUE)),
  "NODE_13511 + nitrosa + CMX1 + inopinata"       = c(grep("NODE_13511|nitrosa|CMX1|inopinata", tl, value = TRUE)),
  "comammox clade A (all)"                        = c(grep("^comammox_A", tl, value = TRUE), grep("NODE_13511", tl, value = TRUE)),
  "comammox clade B"                              = grep("^comammox_B", tl, value = TRUE),
  "comammox A + B"                                = c(grep("^comammox_", tl, value = TRUE), grep("NODE_13511", tl, value = TRUE)),
  "NODE_4969 + N. oligotropha"                    = c(grep("NODE_4969", tl, value = TRUE), grep("oligotropha", tl, value = TRUE)),
  "AOB beta (+ NODE_4969)"                        = c(grep("^AOB_beta", tl, value = TRUE), grep("NODE_4969", tl, value = TRUE))
)
for (g in names(grp)) message(sprintf("%-42s %s", g, sup(grp[[g]])))

# ---- 3. Rótulos ---------------------------------------------------------------
lab_df <- tibble(label = tl) |>
  mutate(
    grupo = case_when(
      grepl("^STUDY", label) ~ "This study",
      grepl("^comammox_A", label) ~ "Comammox clade A",
      grepl("^comammox_B", label) ~ "Comammox clade B",
      grepl("^AOB_beta", label) ~ "AOB (Betaproteobacteria)",
      grepl("^AOB_gamma", label) ~ "AOB (Gammaproteobacteria)",
      TRUE ~ "PmoA (methanotrophs)"),
    nome = sub("^[^|]+\\|([^|]+)\\|.*$", "\\1", label),
    acc  = sub("^.*\\|", "", label),
    texto = case_when(
      grepl("NODE_13511", label) ~ "NODE_13511 AmoA (this study; unbinned, 26x)",
      grepl("NODE_4969", label)  ~ "NODE_4969 AmoA (this study; unbinned, 7.9x)",
      grepl("NODE_2726", label)  ~ "bin.014 PmoA (Methylocystis; this study)",
      TRUE ~ paste0(gsub("_", " ", nome), " (", acc, ")"))
  )

base_size <- 7
font_family <- Sys.getenv("FIG_FONT",
  if (Sys.info()[["sysname"]] == "Darwin") "Helvetica" else "Liberation Sans")
cores <- c("This study" = "black", "Comammox clade A" = "#1b6ca8", "Comammox clade B" = "#5aa0d0",
           "AOB (Betaproteobacteria)" = "#b2182b", "AOB (Gammaproteobacteria)" = "#e08a3c",
           "PmoA (methanotrophs)" = "grey45")

p <- ggtree(tree, linewidth = 0.35, colour = "black") %<+% lab_df +
  geom_tiplab(aes(label = texto, colour = grupo,
                  fontface = ifelse(grupo == "This study", "bold", "plain")),
              size = base_size / .pt, family = font_family, offset = 0.01, show.legend = FALSE) +
  geom_tippoint(aes(colour = grupo), size = 0, alpha = 0) +
  geom_text2(aes(subset = !isTip & label != "" & !is.na(label), label = label),
             hjust = 1.1, vjust = -0.45, size = 4.6 / .pt, family = font_family, colour = "grey30") +
  geom_treescale(x = 0, y = -0.5, width = 0.1, linesize = 0.35,
                 fontsize = (base_size - 1) / .pt, family = font_family, offset = -0.8) +
  scale_colour_manual(values = cores, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1, shape = 15), nrow = 2)) +
  hexpand(0.55, direction = 1) +
  theme_tree() +
  theme(text = element_text(size = base_size, family = font_family),
        legend.position = "bottom", legend.text = element_text(size = base_size),
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(3, 4, 2, 2, unit = "mm"))

ggsave(file.path(dir_tree, "amoA_tree_preview.png"), p, device = ragg::agg_png,
       width = 170, height = 150, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(dir_tree, "amoA_tree_preview.pdf"), p, device = grDevices::pdf,
       width = 170, height = 150, units = "mm", bg = "white")
sessionInfo()
