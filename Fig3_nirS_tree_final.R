# =============================================================================
# Figura 3 — Árvore filogenética do nirS (bin.009/Brocadia + 4 referências)
# Projeto: SBRIA PN/A (anammox) — amostra NGS1140_6S
# Fonte: IQ-TREE3 v3.1.3 (ModelFinder + bootstrap ultrarrápido, 1000 réplicas;
# modelo selecionado por BIC: WAG+G4), logs/reproducibility_log_NGS1140_6S.md,
# Seção final sobre validação filogenética do nirS em bin.009_refined.
#
# Contexto: achado inesperado de nirS (76% identidade, cobertura completa) em
# bin.009 (Ca. Brocadia), gênero cuja literatura publicada documenta ausência
# de nirS/nirK. Árvore construída para descartar contaminação/mau-binning,
# comparando contra nirS de outros gêneros anammox (Ca. Kuenenia, Ca.
# Scalindua) e de um desnitrificante clássico (Pseudomonas aeruginosa,
# Paracoccus denitrificans). Resultado: a sequência de Brocadia não se agrupa
# com o clado desnitrificante clássico (Pseudomonas+Paracoccus+Scalindua,
# 100% UFBoot nos dois nós internos) — consistente com gene retido de forma
# distinta na linhagem, não com contaminação cruzada.
#
# Preferência de fonte de dados: lê o .treefile local se existir; caso
# contrário usa a topologia Newick já registrada em dados_brutos_figuras.md
# como fallback (idêntica — verificado byte a byte contra o arquivo local).
#
# NOTA PARA A LEGENDA DO ARTIGO (não representada no gráfico): esta é uma
# árvore de máxima verossimilhança NÃO-ENRAIZADA (5 sequências, sem outgroup
# definido) — a trifurcação na base (Brocadia | Kuenenia | clado
# Pseudomonas+Paracoccus+Scalindua) é o resultado esperado do IQ-TREE nessas
# condições, não um erro de reconstrução. O layout em cladograma da esquerda
# para a direita é convenção visual comum mesmo para árvores não-enraizadas,
# mas não implica direção ancestral. Recomenda-se declarar explicitamente
# "unrooted maximum-likelihood tree" na legenda da Figura 3.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ape)
  library(ggtree)
})

# ---- 1. Árvore (arquivo local se existir; Newick do documento como fallback)
caminho_arvore <- path.expand(
  "~/sbria_anammox_shotgun/results/targeted_search/nirS_iqtree.treefile"
)

newick_fallback <- "(NOSSO_bin009_Brocadia:0.1152687642,
 ((sp|P24474.1|NIRS_PSEAE:0.3196143665,sp|Q51700.1|NIRS_PARDP:0.2262079926)100:0.3335266123,
  WP_203415591.1:0.3686542654)100:0.3167018626,
 WP_164995172.1:0.1567489551);"

# Versão final (FigS3): exige o .treefile existente em results/targeted_search/.
stopifnot(file.exists(caminho_arvore))
if (file.exists(caminho_arvore)) {
  tree <- read.tree(caminho_arvore)
  message("Arvore lida do arquivo local: ", caminho_arvore)
} else {
  tree <- read.tree(text = newick_fallback)
  message("Arquivo local nao encontrado; usando topologia Newick do documento do projeto.")
}

# ---- 2. Tradução dos códigos de acesso para nome de organismo ---------------
# "Ca." em itálico (convenção já usada na Figura 1 para Ca. Brocadia); demais
# binômios completos em itálico via plotmath (parse = TRUE nos geoms de rótulo).
rotulos <- tribble(
  ~label,                    ~organismo,
  "NOSSO_bin009_Brocadia",  'italic("Ca.")~"Brocadia bin.009 (NODE_696, this study)"',
  "sp|P24474.1|NIRS_PSEAE", 'italic("Pseudomonas aeruginosa")~"(reference)"',
  "sp|Q51700.1|NIRS_PARDP", 'italic("Paracoccus denitrificans")~"(reference)"',
  "WP_203415591.1",         'italic("Ca.")~"Scalindua japonica (reference)"',
  "WP_164995172.1",         'italic("Ca.")~"Kuenenia stuttgartiensis (reference)"'
)

# ---- 3. Gráfico --------------------------------------------------------------
base_size <- 8
font_family <- Sys.getenv("FIG_FONT",
  if (Sys.info()[["sysname"]] == "Darwin") "Helvetica" else "Liberation Sans")

tree_data <- ggtree(tree)$data |> left_join(rotulos, by = "label")
tip_focal <- filter(tree_data, isTip, label == "NOSSO_bin009_Brocadia")
tip_ref   <- filter(tree_data, isTip, label != "NOSSO_bin009_Brocadia")

p <- ggtree(tree, linewidth = 0.4, colour = "black") +
  geom_tiplab(data = tip_focal, aes(x = x, y = y, label = organismo),
              parse = TRUE, fontface = "bold", size = (base_size + 0.5) / .pt,
              family = font_family, offset = 0.012) +
  geom_tiplab(data = tip_ref, aes(x = x, y = y, label = organismo),
              parse = TRUE, size = base_size / .pt,
              family = font_family, offset = 0.012) +
  geom_text2(aes(subset = !isTip & label != "", label = label),
             hjust = 1.3, vjust = -0.6, size = 7 / .pt,  # v_antiga6: 6.5 -> 7 pt
             family = font_family, colour = "grey30") +
  geom_treescale(x = 0, y = 4.6, width = 0.1, linesize = 0.35,
                 fontsize = (base_size - 1) / .pt, family = font_family,
                 offset = -0.15) +
  # hexpand() reserva espaço à direita para o texto dos rótulos sem redefinir
  # os limites numéricos do eixo (evita marcações além da profundidade real
  # da árvore, ~0,97 substituições/sítio no ramo mais longo, Pseudomonas).
  hexpand(0.7, direction = 1) +
  scale_x_continuous(breaks = seq(0, 1.0, 0.2)) +
  theme_tree2() +
  theme(
    text              = element_text(size = base_size, family = font_family),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    # v_antiga5: eixo x graduado removido (linha, marcas, números e título);
    # a escala fica só na barra de 0.1 (geom_treescale acima).
    axis.line.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.text.x       = element_blank(),
    axis.title.x      = element_blank(),
    plot.margin       = margin(3, 4, 2, 2, unit = "mm")
  )

# ---- 4. Exportação (Bioresource Technology / Elsevier) ----------------------
out_dir <- "figures/final"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
w_mm <- 130; h_mm <- 85

# TIFF LZW, 1000 dpi (mesmo formato da versão original do script).
ggsave(file.path(out_dir, "FigS3_nirS_tree.tiff"), p,
       device = ragg::agg_tiff, width = w_mm, height = h_mm, units = "mm",
       dpi = 1000, compression = "lzw", bg = "white")
ggsave(file.path(out_dir, "FigS3_nirS_tree.pdf"), p,
       device = grDevices::pdf, width = w_mm, height = h_mm, units = "mm",
       bg = "white")
ggsave(file.path(out_dir, "FigS3_nirS_tree.png"), p,
       device = ragg::agg_png, width = w_mm, height = h_mm, units = "mm",
       dpi = 300, bg = "white")

sessionInfo()
