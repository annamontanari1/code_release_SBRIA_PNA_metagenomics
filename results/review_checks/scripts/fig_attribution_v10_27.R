# =============================================================================
# Figura composta de atribuição — v10 (README seção 27). Em relação à v9: pink de bin.009/Brocadia #E05A8A (era #CC79A7) e barras de B
# com transparência (alfa 0,6). v1–v9 não são alteradas.
# Paleta v8 (mesmo hex = mesmo significado em todos os painéis): bin.009/Brocadia #CC79A7; bin.007 #1E9E3E; contigs hzs/hdh
# amarelo-canário #FFE81A com contorno preto de 0,5 pt; Kuenenia #C04A00; Jettenia #7A5230; Scalindua #666666.
# v6: hzs/hdh e contigs deste estudo em amarelo ouro (#E6AB02; texto #8C6A00); bin.007 em verde vivo (#1E9E3E); outros bins em azul
# royal (#2E5BDB). Escolhas testadas com simulação de daltonismo (Machado et al. 2009; ver README seção 23). v1–v5 não são alteradas.
# Saída: figures/composite_attribution/Figure_attribution_v10.{pdf,tiff} (190 mm; PDF Quartz com fontes embutidas; TIFF 600 dpi RGB LZW)
# Rodar a partir de ~/sbria_anammox_shotgun:  Rscript results/review_checks/scripts/fig_attribution_v10_27.R
# Dados:
#   A) 18_attribution/panelA_contigs_v2.tsv (attribution_18.py panelA_v2): contigs >= 5 kb (outros bins / não binados), todos os
#      contigs de bin.009 e bin.007, e 18 candidatos (os 17 da v1 + NODE_16231, nxrC/nxrB tipo Brocadia);
#   B) 18_attribution/candidates_links_v2.tsv (attribution_18.py v2): pares inter-contig (primário, MAPQ >= 20, mate em outro contig),
#      sem os pares com os contigs repetitivos NODE_17176 (148×) e NODE_19120 (472×); sem as linhas dos contigs que já estão no bin.009;
#   C) 19_hzsA_v2/tree_main/hzsA_v2_main.rooted_Scalindua.treefile (17 seq.; IQ-TREE 3.1.4, Q.PFAM+G4; SH-aLRT/UFBoot 1000),
#      enraizada em Scalindua; suporte mostrado só quando SH-aLRT >= 70 e UFBoot >= 70; ramo da raiz encurtado (marcado com //).
# Cores (Okabe-Ito): bin.009 e Brocadia #009E73; bin.007 #E69F00; Kuenenia #D55E00; outros bins #0072B2 (barras) e azul-acinzentado
# claro (fundo de A); contigs candidatos de hzs/hdh #CC79A7; NODE_2432 e NODE_4969 #56B4E9; Scalindua #0072B2; Jettenia cinza escuro.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork); library(ape); library(ggtree) })
RC <- "results/review_checks"; OUT <- "figures/composite_attribution"
base_size <- 8; small <- 7
font_family <- Sys.getenv("FIG_FONT", if (Sys.info()[["sysname"]] == "Darwin") "Helvetica" else "Liberation Sans")
GREEN <- "#009E73"; ORANGE <- "#E69F00"; VERM <- "#D55E00"; BLUE <- "#0072B2"; PINK <- "#E05A8A"; SKY <- "#56B4E9"
# v5: tons mais fortes para texto (contraste sobre branco: rosa 5,25:1; verde 5,48:1); #CC79A7 e #009E73 ficam para preenchimentos
PINK_TXT <- "#B0457F"; GREEN_TXT <- "#00785A"; BROWN <- "#8B5A2B"
GOLD <- "#E6AB02"; GOLD_TXT <- "#8C6A00"; GREEN007 <- "#1E9E3E"; ROYAL <- "#2E5BDB"   # v6
# v7, árvore (C): cores que não aparecem em A/B, a não ser por correspondência. NODE verde escuro (≠ verde vivo do bin.007),
# referências de Brocadia roxas (≠ rosa do bin.009), Scalindua cinza (≠ azul royal dos outros bins), Jettenia marrom mais claro
# (separação do verde escuro em daltonismo); pior ΔE entre as cores da árvore (normal/deutan/protan/tritan) = 19
DKGREEN <- "#14532D"; PURPLE <- "#7B3294"; BROWN2 <- "#996633"; SCAL <- "#666666"
# v8
CANARY <- "#FFE81A"; KUEN8 <- "#C04A00"; JETT8 <- "#8E44AD"; SCAL8 <- "#882255"   # v9: Jettenia roxo, Scalindua vinho (v8: #7A5230, #666666)
ST05 <- 0.5 / (96 / 25.4 / 2 * 72 / 96)   # stroke do ggplot que dá contorno de 0,5 pt (lwd = stroke * .stroke / 2, em 1/96 pol)
BLUEGREY <- "#AFC3D6"; LGREY <- "#D0D0D0"
tema <- theme_classic(base_size = base_size, base_family = font_family) +
  theme(plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA),
        axis.text = element_text(colour = "black", size = small), axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3), legend.text = element_text(size = small),
        legend.title = element_text(size = base_size, face = "bold"), plot.tag = element_text(size = 10, face = "bold"))

# ---- A ------------------------------------------------------------------------
a <- read.delim(file.path(RC, "18_attribution/panelA_contigs_v2.tsv"))
gc_m <- 42.53; gc_sd <- 1.95; d_mean9 <- 26.7      # GC: 96 contigs, sem ponderar; profundidade: média ponderada por comprimento (Table S1)
cat_of <- function(c) case_when(
  c %in% c("NODE_5312", "NODE_23806", "NODE_5668", "NODE_20914", "NODE_7170", "NODE_17695", "NODE_29475", "NODE_24575", "NODE_47442",
           "NODE_24842") ~ "hzs",
  c == "NODE_16231" ~ "nxrBC", c %in% c("NODE_4298", "NODE_1742", "NODE_696") ~ "in009",
  c %in% c("NODE_13511", "NODE_2578") ~ "amo007", c %in% c("NODE_2432", "NODE_4969") ~ "nitro_unb")
a_c <- a |> filter(class == "candidate") |> mutate(cat = cat_of(contig))
stopifnot(nrow(a_c) == 18, !anyNA(a_c$cat))
CAT_LAB <- c(hzs = "hzs/hdh contigs (unbinned)", nxrBC = "NODE_16231 (nxrB/nxrC)", in009 = "nxrA, nirS contigs in bin.009",
             amo007 = "amo contigs (bin.007)", nitro_unb = "NODE_2432 (nxrA), NODE_4969 (amo)")
CAT_FILL <- c(hzs = CANARY, nxrBC = "white", in009 = PINK, amo007 = GREEN007, nitro_unb = SKY)
CAT_STROKE <- c(hzs = ST05, nxrBC = 0.3, in009 = 0.3, amo007 = 0.3, nitro_unb = 0.3)   # v8: só o canário muda de contorno
CAT_SHAPE <- c(hzs = 21, nxrBC = 25, in009 = 22, amo007 = 23, nitro_unb = 24)
labs_a <- tribble(
  ~contig,      ~lab,                                 ~lx,   ~ly,   ~hj,
  "NODE_5312",  "NODE_5312",                          56.0,  92,    0,
  "NODE_23806", "NODE_23806",                         52.5,  122,   0,
  "NODE_16231", "NODE_16231 (nxrB/nxrC)",             45.0,  72,    1,   # v3: à esquerda, sem cruzar o traço de NODE_23806
  "NODE_5668",  "NODE_5668",                          44.0,  16.5,  1,
  "NODE_20914", "NODE_20914",                         44.0,  9.0,   1,
  "NODE_13511", "NODE_13511",                         50.5,  40,    1,
  "NODE_2578",  "NODE_2578",                          60.5,  15.0,  0,
  "NODE_2432",  "NODE_2432",                          63.0,  60,    0,
  "NODE_4969",  "NODE_4969",                          54.0,  6.0,   0,
  "NODE_4298",  "NODE_4298, NODE_1742,\nNODE_696",   31.0,  47,    0
) |> left_join(a_c |> select(contig, GC_pct, depth_x), by = "contig") |>
  # ponto final do traço: o próprio rótulo, exceto no rótulo de duas linhas (base do texto)
  mutate(sx = ifelse(contig == "NODE_4298", 36.5, lx), sy = ifelse(contig == "NODE_4298", 41.5, ly))
bg <- a |> filter(class == "background") |> mutate(bgk = ifelse(final_bin == "unbinned", "Unbinned contigs (>= 5 kb)", "Other bins (>= 5 kb)"))
pA <- ggplot() +
  annotate("rect", xmin = gc_m - gc_sd, xmax = gc_m + gc_sd, ymin = 3, ymax = 170, fill = PINK, alpha = 0.10) +
  geom_hline(yintercept = d_mean9, linetype = "22", linewidth = 0.3, colour = PINK) +
  geom_point(data = bg, aes(GC_pct, depth_x, colour = bgk), alpha = 0.22, size = 0.4, stroke = 0) +   # v7: menos destaque
  geom_point(data = filter(a, final_bin == "bin.007", class != "candidate"), aes(GC_pct, depth_x, colour = "bin.007"), alpha = 0.75, size = 0.55, stroke = 0) +
  geom_point(data = filter(a, final_bin == "bin.009", class != "candidate"), aes(GC_pct, depth_x, colour = "bin.009"), alpha = 0.9, size = 0.7, stroke = 0) +
  geom_segment(data = labs_a, aes(x = GC_pct, y = depth_x, xend = sx, yend = sy), linewidth = 0.25, colour = "grey30") +
  geom_point(data = a_c |> mutate(st = CAT_STROKE[cat]), aes(GC_pct, depth_x, shape = cat, fill = cat, stroke = st), colour = "black", size = 2.0) +
  geom_text(data = labs_a, aes(lx, ly, label = lab, hjust = hj), size = small / .pt, family = font_family, lineheight = 0.9) +
  annotate("text", x = 79.5, y = d_mean9 * 0.87, label = "bin.009 length-weighted mean 26.7×", hjust = 1, size = small / .pt,
           family = font_family, colour = PINK) +
  annotate("text", x = gc_m, y = 3.25, label = "bin.009 GC\n42.5 ± 2.0%", size = small / .pt, family = font_family,
           colour = PINK, lineheight = 0.9, vjust = 0) +
  scale_colour_manual(values = c("bin.009" = PINK, "bin.007" = GREEN007, "Other bins (>= 5 kb)" = ROYAL, "Unbinned contigs (>= 5 kb)" = LGREY),
                      breaks = c("bin.009", "bin.007", "Other bins (>= 5 kb)", "Unbinned contigs (>= 5 kb)"), name = "Contigs",
                      labels = c(expression("bin.009 ("*italic("Ca.")~"Brocadia)"), expression("bin.007 ("*italic("Nitrospira")*"_F)"),
                                 expression("Other bins (">="5 kb)"), expression("Unbinned (">="5 kb)"))) +
  scale_shape_manual(values = CAT_SHAPE, labels = CAT_LAB, breaks = names(CAT_LAB), name = "Candidate contigs") +
  scale_fill_manual(values = CAT_FILL, labels = CAT_LAB, breaks = names(CAT_LAB), name = "Candidate contigs") +
  scale_y_log10(name = "Mean depth (×)", breaks = c(5, 10, 20, 50, 100), limits = c(3, 170), expand = c(0, 0)) +
  scale_x_continuous(name = "GC (%)", limits = c(30, 80), breaks = seq(30, 80, 10)) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = 1.8, alpha = 1)),
         shape = guide_legend(order = 2, override.aes = list(size = 2, stroke = unname(CAT_STROKE[names(CAT_LAB)]))), fill = guide_legend(order = 2)) +
  labs(tag = "A") + tema + theme(legend.position = "right", legend.spacing.y = unit(2, "mm"), legend.key.size = unit(3.2, "mm"))

# ---- B ------------------------------------------------------------------------
genes <- read.delim(file.path(RC, "13_hzs_links/summary_links.tsv")) |> select(contig, genes)
b <- read.delim(file.path(RC, "18_attribution/candidates_links_v2.tsv")) |>
  filter(!contig %in% c("NODE_4298", "NODE_1742", "NODE_696")) |> left_join(genes, by = "contig") |>
  mutate(genes = case_when(contig == "NODE_16231" ~ "nxrB/nxrC", contig == "NODE_13511" ~ "amoAB", TRUE ~ genes),
         lab = sprintf('"%s ("*italic("%s")*")"', contig, genes),
         painel = factor(ifelse(target_bin == "bin.009", "Target:\nbin.009", "Target:\nbin.007"), c("Target:\nbin.009", "Target:\nbin.007")),
         ord = ifelse(contig == "NODE_16231", 2000, ifelse(target_bin == "bin.007", 3000, 1000)) - inter_q20)
stopifnot(nrow(b) == 12, b$q20_to_target[b$contig == "NODE_20914"] == 4, b$inter_q20[b$contig == "NODE_16231"] == 158)
b_long <- bind_rows(
  b |> transmute(contig, painel, inter_q20, cat = ifelse(target_bin == "bin.009", "bin.009 (target)", "bin.007 (target)"), n = q20_to_target),
  b |> transmute(contig, painel, inter_q20, cat = "Other bins", n = q20_to_other_bins),
  b |> transmute(contig, painel, inter_q20, cat = "Unbinned contigs", n = q20_to_unbinned)) |>
  mutate(frac = ifelse(inter_q20 > 0, n / inter_q20, 0),
         cat = factor(cat, c("Unbinned contigs", "Other bins", "bin.007 (target)", "bin.009 (target)")),
         contig = factor(contig, levels = b$contig[order(-b$ord)]))
lab_y <- setNames(b$lab, b$contig)
pB <- ggplot(b_long, aes(y = contig, x = frac, fill = cat)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.2, alpha = 0.6) +   # v10: barras com transparência
  # v3: ponto rosa ao lado dos contigs hzs/hdh; "no pairs" onde n = 0
  geom_point(data = filter(b, genes != "nxrB/nxrC", genes != "amoAB") |> mutate(contig = factor(contig, levels(b_long$contig))),
             inherit.aes = FALSE, aes(x = -0.035, y = contig), shape = 21, fill = CANARY, colour = "black", stroke = ST05, size = 1.6) +
  geom_text(data = filter(b, inter_q20 == 0) |> mutate(contig = factor(contig, levels(b_long$contig))), inherit.aes = FALSE,
            aes(x = 0.02, y = contig, label = "no pairs"), hjust = 0, size = small / .pt, family = font_family, colour = "grey30") +
  geom_text(data = distinct(b_long, contig, painel, inter_q20), inherit.aes = FALSE,
            aes(y = contig, x = 1.02, label = paste0("n = ", inter_q20)), hjust = 0, size = small / .pt, family = font_family) +
  facet_grid(painel ~ ., scales = "free_y", space = "free_y") +
  scale_y_discrete(labels = function(x) parse(text = lab_y[x])) +
  scale_x_continuous(name = expression("Inter-contig read pairs, MAPQ" >= "20 (%)"),
                     breaks = seq(0, 1, 0.25), labels = function(x) x * 100, expand = c(0, 0)) +
  scale_fill_manual(values = c("bin.009 (target)" = PINK, "bin.007 (target)" = GREEN007, "Other bins" = ROYAL, "Unbinned contigs" = LGREY),
                    breaks = c("bin.009 (target)", "bin.007 (target)", "Other bins", "Unbinned contigs"), name = NULL) +
  labs(y = NULL, tag = "B") + coord_cartesian(xlim = c(0, 1.32), expand = FALSE, clip = "off") +
  guides(fill = guide_legend(nrow = 2)) +
  tema + theme(legend.position = "bottom", strip.background = element_blank(), strip.text.y = element_text(size = small, angle = 0, hjust = 0),
               axis.text.y = element_text(size = small, margin = margin(r = 8)), axis.ticks.y = element_blank(),
               panel.spacing.y = unit(2, "mm"), legend.key.size = unit(3, "mm"),
               legend.margin = margin(0, 0, 0, 0))

# ---- C (v3: circular, mesma árvore e suportes da v2) ------------------------------
tr <- read.tree(file.path(RC, "19_hzsA_v2/tree_main/hzsA_v2_main.rooted_Scalindua.treefile"))
stopifnot(Ntip(tr) == 17)
rt <- Ntip(tr) + 1
e_long <- which(tr$edge[, 1] == rt & tr$edge.length == max(tr$edge.length[tr$edge[, 1] == rt]))
L_real <- tr$edge.length[e_long]; tr$edge.length[e_long] <- 0.25          # ramo da raiz encurtado (real 0,50, anotado); 0,15 afasta os nós internos do centro
# rosa das referências de Brocadia: #E8B4D0 pedido; contraste 1,77:1 sobre branco (pouco legível a 7 pt) -> #D88CB5 (2,52:1)
BROC_REF <- PINK; JETT <- BROWN                                  # v5: ramos; os rótulos usam os tons de texto (TXT_COL)
partial <- c(NODE_7170 = "570 aa", NODE_29475 = "408 aa", NODE_17695 = "115 aa")
abbr <- c(Brocadia = "B.", Jettenia = "J.", Kuenenia = "K.", Scalindua = "S.")
tip <- tibble(label = tr$tip.label) |>
  mutate(node_id = sub("^THISSTUDY_(NODE_\\d+)_\\d+$", "\\1", label),
         genus = ifelse(grepl("^THISSTUDY", label), "contig", sub("_.*$", "", label)),
         sp = sub("^[A-Za-z]+_([A-Za-z0-9]+)_.*$", "\\1", label),
         acc = sub("v(\\d+)$", ".\\1", sub("^.*_((WP|ODS|KHE|AEW|TLD|GAN|KKO|BBO)_?\\d+v\\d+)$", "\\1", label)),
         grp = recode(genus, Brocadia = "Brocadia"),
         rot = ifelse(genus == "contig",
                      sprintf('bold("%s%s")', node_id, ifelse(node_id %in% names(partial), " (partial)", "")),
                      sprintf('italic("Ca.")~"%s %s"', abbr[genus], ifelse(grepl("^sp", sp), sub("^sp", "sp. ", sp), sp))))
GRP_COL <- c(Brocadia = PINK, contig = "black", Jettenia = JETT8, Kuenenia = KUEN8, Scalindua = SCAL8)   # v8: contigs com ramo e rótulo pretos
TXT_COL <- GRP_COL                                            # v7: rótulos na mesma cor dos ramos (todas com contraste >= 3,9:1)
gC <- ggtree(tr, layout = "circular", linewidth = 0.4) %<+% tip
d <- gC$data
sup <- d |> filter(!isTip, !is.na(label), label != "") |>
  mutate(sh = suppressWarnings(as.numeric(sub("/.*", "", label))), uf = suppressWarnings(as.numeric(sub(".*/", "", label)))) |>
  filter(sh >= 70, uf >= 70, x > 0)   # x = 0: nó na raiz (mesma bipartição do grupo interno); mostrado uma vez só
brk <- d[d$node == tr$edge[e_long, 2], ]
xmax_tip <- max(d$x[d$isTip])
# barra de escala radial na direção de 3 h (ângulo 0 dos rótulos), fora dos rótulos; y obtido da relação ângulo ~ y das pontas
# posição de cada valor de suporte, nó a nó (identificado pelo conjunto de folhas): recuo radial dx (para dentro do nó) e lado
# do deslocamento angular (meio passo de ângulo). Ajustado à mão, comparando as opções, para não tocar ramos, pontas, rótulos
# nem o número do nó vizinho.
SUP_DY <- 0.5
pos_sup <- list(
  list(t = c("Scalindua_brodae_KHE92665v1", "Scalindua_sp15L_AEW50032v1"), dx = 0.095, lado = 1),   # v4: escolhido entre 6 testes
  list(t = c("Jettenia_asiatica_AEW50031v1", "Jettenia_ecosi_TLD39960v1"), dx = 0.02, lado = 1),
  list(t = c("Brocadia_fulgida_KKO20072v1", "THISSTUDY_NODE_7170_1", "Brocadia_pituitae_BBO18370v1"), dx = 0.05, lado = 1),
  list(t = c("Brocadia_fulgida_KKO20072v1", "THISSTUDY_NODE_7170_1"), dx = 0.042, lado = 1),
  list(t = c("Brocadia_sapporoensis_WP_080325085v1", "THISSTUDY_NODE_5312_3", "THISSTUDY_NODE_29475_2", "Brocadia_sinica_GAN34140v1",
             "Brocadia_spAMX2_WP_285442451v1"), dx = 0.05, lado = -1),
  list(t = c("Brocadia_sapporoensis_WP_080325085v1", "THISSTUDY_NODE_5312_3", "THISSTUDY_NODE_29475_2"), dx = 0.064, lado = -1),
  list(t = c("Brocadia_sapporoensis_WP_080325085v1", "THISSTUDY_NODE_5312_3"), dx = 0.1, lado = -1),   # dentro do círculo
  list(t = c("Brocadia_sinica_GAN34140v1", "Brocadia_spAMX2_WP_285442451v1"), dx = 0.035, lado = -1))
achar <- function(n) { tips <- extract.clade(tr, n)$tip.label
  for (q in pos_sup) if (setequal(tips, q$t)) return(c(q$dx, q$lado))
  c(0.045, 1) }                                                    # padrão (nó do grupo interno): dentro do arco da raiz
pp <- t(vapply(sup$node, achar, numeric(2))); sup$dx <- pp[, 1]; sup$lado <- pp[, 2]
ROT <- as.numeric(Sys.getenv("ROT_C", "160")); LEG_POS <- as.numeric(strsplit(Sys.getenv("LEG_C", "0,0"), ",")[[1]])                # rotação do layout (graus), escolhida por legibilidade
gC <- rotate_tree(gC + aes(colour = grp), ROT)
d <- gC$data
fit <- lm(angle ~ y, data = filter(d, isTip)); ang <- function(y) (coef(fit)[1] + coef(fit)[2] * y) %% 360
sb_a <- as.numeric(Sys.getenv("SB_ANG", "168")) + 360 * (-3:3)       # ângulo pedido e equivalentes (±360°)
y0 <- (sb_a - coef(fit)[1]) / coef(fit)[2]; y0 <- y0[y0 >= 0 & y0 < Ntip(tr) + 1][1]   # o que cai dentro da faixa de y da árvore                       # barra de escala na direção 0, fora dos rótulos
XL <- xmax_tip + 0.183                                              # margem radial (v3: + 0.2); com a legenda no canto, raio ~18% maior que na v3 (medido)
tl <- d |> filter(isTip) |> mutate(a = angle %% 360, flip = a > 90 & a < 270, a_txt = ifelse(flip, a + 180, a), hj = ifelse(flip, 1, 0))
SB_X <- XL - 0.065; SB_Y <- y0                                     # dentro do limite radial, além dos rótulos de Scalindua (9 h)
# suportes: texto na posição da v3 (pos_sup) e traço fino curto até o próprio nó
sup <- sup |> mutate(lx = x - dx, ly = y + SUP_DY * lado, sx = lx + 0.3 * (x - lx), sy = ly + 0.3 * (y - ly))
pC <- gC +
  # rótulos das pontas desenhados aqui (e não com geom_tiplab): depois de rotate_tree a inversão automática dos rótulos do lado
  # esquerdo deixa de funcionar; a orientação vem do ângulo de cada ponta já rotacionado
  geom_text(data = tl, aes(x = x + 0.035, y = y, label = rot, angle = a_txt, hjust = hj), colour = TXT_COL[tl$grp], parse = TRUE,   # v5: mais espaço até a ponta
            size = small / .pt, family = font_family, show.legend = FALSE) +
  geom_tippoint(aes(subset = genus == "contig"), shape = 21, fill = CANARY, colour = "black", stroke = ST05, size = 1.9, show.legend = FALSE) +
  geom_segment(data = sup, aes(x = sx, y = sy, xend = x, yend = y), inherit.aes = FALSE, linewidth = 0.2, colour = "grey45") +
  geom_text(data = sup, aes(x = lx, y = ly, label = uf), inherit.aes = FALSE, size = small / .pt, family = font_family, colour = "grey15") +
  annotate("text", x = brk$x * 0.45, y = brk$y, label = "//", size = small / .pt + 1, family = font_family, fontface = "bold", colour = "black") +
  annotate("segment", x = SB_X, xend = SB_X + 0.05, y = SB_Y, yend = SB_Y, linewidth = 0.4, colour = "black") +
  annotate("text", x = SB_X + 0.025, y = SB_Y + 0.35, label = "0.05", hjust = 0.5, size = small / .pt, family = font_family) +
  scale_colour_manual(values = GRP_COL, na.value = "black", breaks = c("Brocadia", "contig", "Jettenia", "Kuenenia", "Scalindua"), name = NULL,
                      labels = c(expression(italic("Ca.")~"Brocadia (references)"), "This study (contigs)", expression(italic("Ca.")~"Jettenia"),
                                 expression(italic("Ca.")~"Kuenenia"), expression(italic("Ca.")~"Scalindua (outgroup)"))) +
  scale_x_continuous(limits = c(0, XL), expand = c(0, 0)) +
  guides(colour = guide_legend(ncol = 1)) +
  labs(tag = "C", caption = sprintf(paste0("Node values: UFBoot (shown when SH-aLRT \u2265 70 and UFBoot \u2265 70).\n",
                                           "// root branch shortened (%.2f substitutions/site).\n",
                                           "Partial HzsA: NODE_7170 (570 aa), NODE_29475 (408 aa), NODE_17695 (115 aa)."), L_real)) +
  theme(legend.position = "inside", legend.position.inside = LEG_POS, legend.justification = LEG_POS,
        legend.text = element_text(size = small, family = font_family), legend.key.size = unit(3, "mm"), legend.background = element_blank(),
        plot.caption = element_text(size = small, family = font_family, hjust = 0), plot.caption.position = "plot",
        plot.tag = element_text(size = 10, face = "bold", family = font_family), plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(1, 1, 1, 1, "mm"))
# v8: legenda do C com chave própria. As camadas da árvore não entram na legenda; uma camada vazia (x = NA) gera as chaves:
# círculo canário contornado (0,5 pt) para os contigs deste estudo e linha colorida para os demais
key_v8 <- function(data, params, size) {
  if (identical(data$colour, "black"))
    grid::pointsGrob(0.5, 0.5, pch = 21, size = unit(2.2, "mm"), gp = grid::gpar(col = "black", fill = CANARY, lwd = ST05 * 96 / 25.4 / 2))
  else grid::segmentsGrob(0.05, 0.5, 0.95, 0.5, gp = grid::gpar(col = data$colour, lwd = 1.2 * .pt, lineend = "butt"))
}
for (i in seq_along(pC$layers)) pC$layers[[i]]$show.legend <- FALSE
pC <- pC + geom_point(data = tibble(grp = names(GRP_COL), x = NA_real_, y = NA_real_), aes(x = x, y = y, colour = grp),
                      inherit.aes = FALSE, na.rm = TRUE, key_glyph = key_v8)
if (nzchar(Sys.getenv("ONLY_C"))) {                                  # teste rápido: só o painel C, no tamanho que ele tem na figura
  ggsave(Sys.getenv("ONLY_C"), pC, device = ragg::agg_png, width = 104, height = 118, units = "mm", dpi = 200, bg = "white"); quit(save = "no")
}

# ---- Montagem -------------------------------------------------------------------
p <- wrap_elements(full = pA) / (wrap_elements(full = pB) + wrap_elements(full = pC) + plot_layout(widths = c(0.9, 1.1))) +
  plot_layout(heights = c(0.9, 1.25))
w_mm <- 190; h_mm <- 205
ggsave(file.path(OUT, "Figure_attribution_v10.tiff"), p, device = ragg::agg_tiff, width = w_mm, height = h_mm, units = "mm",
       dpi = 600, compression = "lzw", bg = "white")
quartz_pdf <- function(filename, width, height, ...) grDevices::quartz(type = "pdf", file = filename, width = width, height = height, bg = "white")
ggsave(file.path(OUT, "Figure_attribution_v10.pdf"), p, device = quartz_pdf, width = w_mm, height = h_mm, units = "mm")
sessionInfo()
