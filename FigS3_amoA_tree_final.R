# =============================================================================
# Figura S3 — Árvore filogenética de AmoA (comammox clados A e B, AOB, PmoA)
# Projeto: SBRIA PN/A (anammox) — amostra NGS1140_6S. Criada em 19/10/2026.
# Estilo de Fig3_nirS_tree_final.R (agora FigS4_nirS_tree_final.R): ggtree, base 8 pt, Helvetica, linhas pretas,
# sequências deste estudo em negrito.
#
# Dados (results/review_checks/comammox/amoA_tree/; ver README.md da pasta comammox):
#   - 27 AmoA/PmoA de referência do NCBI (panel.tsv) + AmoA de NODE_13511 e de NODE_4969 + PmoA do bin.014;
#   - MAFFT v7.526 --localpair --maxiterate 1000; trimAl v1.5.rev1 -automated1 (329 -> 274 colunas);
#   - IQ-TREE 3.1.3: -m MFP (BIC: Q.YEAST+F+G4) -B 1000 -alrt 1000 -seed 12345.
# Suporte nos nós: SH-aLRT (%) / UFBoot (%). Raiz: PmoA + AmoA de Nitrosococcus (Gammaproteobacteria) — só para
# orientar a leitura (enraizamento por grupo externo, não temporal).
# Barras de grupo calculadas por blocos contíguos de terminais (AOB beta e PmoA não são monofiléticos nesta árvore).
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ape)
  library(ggtree)
})

caminho_arvore <- path.expand("~/sbria_anammox_shotgun/results/review_checks/comammox/amoA_tree/amoA_tree.treefile")
stopifnot(file.exists(caminho_arvore))
tree <- read.tree(caminho_arvore)
outg <- grep("^pmoA\\||^AOB_gamma\\||bin.014_NODE_2726_PmoA", tree$tip.label, value = TRUE)
tree <- root(tree, outgroup = outg, resolve.root = TRUE, edgelabel = TRUE)
tree <- ladderize(tree, right = FALSE)
tree$node.label[tree$node.label %in% c("Root", "root")] <- ""

# ---- 1. Rótulos -------------------------------------------------------------
# Plotmath: gênero/espécie em itálico, "Ca." em itálico, códigos e acessos em romano.
rotulo <- function(lab) {
  if (grepl("NODE_13511", lab)) return('bold("NODE_13511 AmoA (this study; bin.007\u2021)")')
  if (grepl("NODE_4969", lab))  return('bold("NODE_4969 AmoA (this study; unbinned)")')
  # 24/10/2026 (cópia *_v_antiga5.R): bin.014 com o mesmo destaque (negrito + caixa cinza) dos outros terminais deste estudo
  if (grepl("NODE_2726", lab))  return('bold("bin.014 PmoA (")*bolditalic("Methylocystis")*bold("; this study)")')
  nome <- gsub("_", " ", sub("^[^|]+\\|([^|]+)\\|.*$", "\\1", lab)); acc <- sub("^.*\\|", "", lab)
  # WP_129883941.1: organismo do registro NCBI = Methylococcaceae (MULTISPECIES); "Methylotetracoccus" é o gênero GTDB de uma das fontes
  nome <- sub(" MULTISPECIES$", "", sub("^Methylotetracoccus Methylococcaceae", "Methylococcaceae", nome))
  if (startsWith(nome, "Ca. ")) return(sprintf('italic("Ca.")~"%s (%s)"', sub("^Ca\\. ", "", nome), acc))
  partes <- strsplit(nome, " ")[[1]]
  if (length(partes) >= 2 && partes[2] == "sp.")
    return(sprintf('italic("%s")~"%s (%s)"', partes[1], paste(partes[-1], collapse = " "), acc))
  if (grepl("^(Nitrospiraceae|Methylococcaceae$)", nome)) return(sprintf('"%s (%s)"', nome, acc))
  # binômio em itálico; linhagem/estirpe (ex.: C91, OBBP) em romano
  if (length(partes) > 2) return(sprintf('italic("%s")~"%s (%s)"', paste(partes[1:2], collapse = " "),
                                         paste(partes[-(1:2)], collapse = " "), acc))
  sprintf('italic("%s")~"(%s)"', nome, acc)
}
tl <- tree$tip.label
tip_df <- tibble(label = tl,
                 rot = vapply(tl, rotulo, character(1)),
                 foco = grepl("^STUDY\\|NODE_(13511|4969)|NODE_2726", tl),
                 grupo = case_when(
                   grepl("^comammox_A|NODE_13511", tl) ~ "Comammox clade A",
                   grepl("^comammox_B", tl)            ~ "Comammox clade B",
                   grepl("^AOB_beta|NODE_4969", tl)    ~ "AOB (Betaproteobacteria)",
                   grepl("^AOB_gamma", tl)             ~ "Nitrosococcus (AOB, Gammaproteobacteria)",
                   TRUE                                ~ "pmoA (methanotrophs)"))
invisible(parse(text = tip_df$rot))

base_size <- 8
font_family <- Sys.getenv("FIG_FONT",
  if (Sys.info()[["sysname"]] == "Darwin") "Helvetica" else "Liberation Sans")

g0 <- ggtree(tree)
td <- g0$data |> left_join(tip_df, by = "label")
xmax <- max(td$x[td$isTip])
# largura aproximada do rótulo em destaque (caracteres visíveis; negrito)
# trecho em bolditalic conta 0,8 caractere (letras inclinadas mais estreitas), para a caixa terminar como nos outros dois
nch_vis <- function(r) {
  it <- regmatches(r, gregexpr('bolditalic\\("[^"]*"\\)', r))[[1]]
  n_it <- sum(nchar(gsub('bolditalic\\("|"\\)', "", it)))
  nchar(gsub('(bold|bolditalic)\\("|"\\)|\\*', "", r)) - 0.2 * n_it
}
td <- td |> mutate(nch = ifelse(isTip & foco, vapply(rot, nch_vis, numeric(1)), NA))

# Barras de grupo: blocos contíguos de terminais (ordem vertical do gráfico)
barras <- td |> filter(isTip) |> arrange(y) |>
  mutate(bloco = cumsum(grupo != lag(grupo, default = first(grupo)))) |>
  group_by(grupo, bloco) |> summarise(y0 = min(y) - 0.3, y1 = max(y) + 0.3, ym = mean(y), .groups = "drop") |>
  mutate(rot = case_when(
    grupo == "Nitrosococcus (AOB, Gammaproteobacteria)" ~ 'italic("Nitrosococcus")~"(AOB)"',
    grupo == "pmoA (methanotrophs)" ~ '"PmoA"',
    grupo == "AOB (Betaproteobacteria)" ~ '"AOB (Beta)"',
    TRUE ~ sprintf('"%s"', grupo)))
x_bar <- xmax * 1.98
# Rótulos de suporte (SH-aLRT/UFBoot); nenhum rótulo é omitido, conteúdo inalterado.
# Posicionamento por checagem geométrica (19/10/2026-b): para cada nó testam-se, nesta ordem,
#   (1) acima do ramo que leva ao nó, no ponto médio; (2) à direita do nó, entre os dois ramos filhos;
#   (3) acima e à esquerda do nó; usa-se a primeira posição cuja caixa de texto não cruza ramos horizontais,
#   conectores verticais, rótulos de terminais ou rótulos de suporte já colocados (senão, a de menos conflitos).
# Altura da figura 165 -> 205 mm (19/10/2026-b): mais espaço entre linhas para os rótulos de suporte; fontes e largura
# inalteradas. Alturas convertidas de mm para "linhas" da árvore.
w_mm <- 180; h_mm <- 205
row_mm <- h_mm * 0.86 / (sum(td$isTip) + 2)
# 21/10/2026 (cópia *_v_antiga4.R): suporte 5,2 -> 5,7 pt, para >= 6 pt quando a figura é ampliada a 190 mm
# (figures/final_tiff); largura e altura usadas na checagem de colisão escalam com o tamanho.
sz_sup <- 5.7
cw <- xmax * 0.0145 * sz_sup / 5.2   # largura aproximada de um caractere (unidades do eixo x)
ch <- 2.0 * sz_sup / 5.2 / row_mm    # altura do texto de suporte (5,2 pt ~ 1,8 mm, com folga)
th <- 3.0 / row_mm / 2         # meia altura do texto dos terminais (8 pt ~ 2,8 mm)
nodes <- td |> mutate(px = td$x[match(parent, td$node)])
hseg <- nodes |> filter(node != parent) |> transmute(x0 = px, x1 = x, y)
vseg <- nodes |> filter(!isTip) |> rowwise() |>
  mutate(y0 = min(td$y[td$parent == node & td$node != node]), y1 = max(td$y[td$parent == node & td$node != node])) |>
  ungroup() |> transmute(x, y0, y1)
tipbox <- td |> filter(isTip) |> transmute(x0 = x, x1 = x + xmax * 0.9, y0 = y - th, y1 = y + th)
caixa <- function(lx, ly, hj, vj, w) c(lx - hj * w, lx + (1 - hj) * w, ly - vj * ch, ly + (1 - vj) * ch)
conflitos <- function(b, colocadas) {
  pad <- 0.004 * xmax
  n_h <- sum(hseg$y > b[3] & hseg$y < b[4] & pmax(hseg$x0, hseg$x1) > b[1] - pad & pmin(hseg$x0, hseg$x1) < b[2] + pad)
  n_v <- sum(vseg$x > b[1] - pad & vseg$x < b[2] + pad & vseg$y1 > b[3] & vseg$y0 < b[4])
  n_t <- sum(tipbox$x0 < b[2] & tipbox$x1 > b[1] & tipbox$y0 < b[4] & tipbox$y1 > b[3])
  n_l <- if (length(colocadas)) sum(vapply(colocadas, function(o) o[1] < b[2] && o[2] > b[1] && o[3] < b[4] && o[4] > b[3], logical(1))) else 0
  n_h + n_v + n_t + n_l
}
sup <- nodes |> filter(!isTip, !is.na(label), label != "") |> arrange(desc(x))
res <- vector("list", nrow(sup)); colocadas <- list()
for (i in seq_len(nrow(sup))) {
  r <- sup[i, ]; w <- nchar(r$label) * cw
  cand <- list(
    list(lx = (r$px + r$x) / 2, ly = r$y + 0.08, hj = 0.5, vj = 0),
    list(lx = r$x + xmax * 0.008, ly = r$y, hj = 0, vj = 0.5),
    list(lx = r$x - xmax * 0.006, ly = r$y + 0.08, hj = 1, vj = 0),
    list(lx = r$x - xmax * 0.006, ly = r$y - 0.08, hj = 1, vj = 1))
  nc <- vapply(cand, function(cd) conflitos(caixa(cd$lx, cd$ly, cd$hj, cd$vj, w), colocadas), numeric(1))
  k <- which(nc == 0)[1]; if (is.na(k)) k <- which.min(nc)
  cd <- cand[[k]]; colocadas[[length(colocadas) + 1]] <- caixa(cd$lx, cd$ly, cd$hj, cd$vj, w)
  res[[i]] <- tibble(label = r$label, lx = cd$lx, ly = cd$ly, hj = cd$hj, vj = cd$vj, pos = k, conflitos = nc[k])
}
sup <- bind_rows(res)
message("support labels: ", nrow(sup), "; positions used: ", paste(names(table(sup$pos)), table(sup$pos), sep = "=", collapse = " "),
        "; labels with residual conflicts: ", sum(sup$conflitos > 0))

p <- ggtree(tree, linewidth = 0.35, colour = "black") %<+% tip_df +
  geom_rect(data = filter(td, isTip, foco), inherit.aes = FALSE,
            aes(xmin = x + 0.01, xmax = x + 0.02 + xmax * 0.0240 * nch, ymin = y - 0.45, ymax = y + 0.45),
            fill = "grey88", colour = NA) +
  geom_tiplab(aes(label = rot), parse = TRUE, size = base_size / .pt, family = font_family, offset = 0.012) +
  geom_text(data = sup, inherit.aes = FALSE, aes(x = lx, y = ly, label = label, vjust = vj, hjust = hj),
            size = sz_sup / .pt, family = font_family, colour = "grey25") +
  geom_segment(data = barras, inherit.aes = FALSE, aes(x = x_bar, xend = x_bar, y = y0, yend = y1),
               linewidth = 0.6, colour = "black") +
  geom_text(data = barras, inherit.aes = FALSE, aes(x = x_bar + xmax * 0.02, y = ym, label = rot), parse = TRUE,
            hjust = 0, size = base_size / .pt, family = font_family) +
  geom_treescale(x = 0, y = -0.6, width = 0.1, linesize = 0.35, fontsize = (base_size - 1) / .pt,
                 family = font_family, offset = -0.9) +
  coord_cartesian(xlim = c(0, x_bar + xmax * 0.42), clip = "off") +
  theme_tree() +
  theme(text = element_text(size = base_size, family = font_family),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(3, 4, 4, 2, unit = "mm"))

# ---- 2. Exportação (TIFF 1000 dpi LZW + PDF; sem PNG desde 21/10/2026, como pedido pela autora) ----
out_dir <- "figures/final"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(out_dir, "FigS3_amoA_tree.tiff"), p,
       device = ragg::agg_tiff, width = w_mm, height = h_mm, units = "mm",
       dpi = 1000, compression = "lzw", bg = "white")
ggsave(file.path(out_dir, "FigS3_amoA_tree.pdf"), p, device = grDevices::pdf,
       width = w_mm, height = h_mm, units = "mm", bg = "white", encoding = "WinAnsi")
sessionInfo()
