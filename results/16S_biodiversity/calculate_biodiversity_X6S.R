#!/usr/bin/env Rscript
# ============================================================================
# Cálculo de índices de biodiversidade - Amostra X6S (reator SBRIA anammox/PN-A)
# Projeto independente (não pertence à tese de doutorado)
#
# Input: rarefiedCounts.tsv (contagem de ASVs já normalizada por rarefação,
#        profundidade padronizada em 40.321 reads/amostra)
#        ASV_genus.nwt (árvore filogenética dos ASVs, para Faith's PD)
#
# Método: vegan (índices padrão de diversidade alfa) + picante (Faith's PD)
# Autora: Anna Flavia Pereira Montanari
# Data: 02/09/2026
# ============================================================================

library(tidyverse)
library(vegan)
library(picante)

# --- 1. Carregar dados ---
counts_raw <- read_tsv("rarefiedCounts.tsv", show_col_types = FALSE)

tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
first_col <- colnames(counts_raw)[1]
if (first_col == "...1" || first_col == "") {
  counts_raw <- counts_raw %>% rename(ASV_ID = 1)
} else {
  counts_raw <- counts_raw %>% rename(ASV_ID = all_of(first_col))
}
sample_cols <- setdiff(colnames(counts_raw), c("ASV_ID", tax_cols))

cat("Amostras disponíveis no arquivo:", length(sample_cols), "\n")
cat("Amostra de interesse: X6S\n\n")

stopifnot("X6S" %in% sample_cols)

# --- 2. Extrair vetor de abundância da amostra X6S ---
asv_table <- counts_raw %>%
  select(ASV_ID, X6S) %>%
  filter(!is.na(X6S), X6S > 0)

abund_vector <- asv_table$X6S
names(abund_vector) <- asv_table$ASV_ID

cat("Número de ASVs presentes em X6S:", length(abund_vector), "\n")
cat("Total de reads (rarefeito):", sum(abund_vector), "\n\n")

# --- 3. Índices de diversidade alfa via vegan ---
observed_richness <- specnumber(abund_vector)
shannon <- diversity(abund_vector, index = "shannon")
simpson <- diversity(abund_vector, index = "simpson")
invsimpson <- diversity(abund_vector, index = "invsimpson")
pielou_evenness <- shannon / log(observed_richness)

chao_est <- estimateR(abund_vector)
chao1 <- chao_est["S.chao1"]

# --- 4. Faith's Phylogenetic Diversity via picante ---
tree <- read.tree("ASV_genus.nwt")

# IMPORTANTE: drop = FALSE evita que o R derrube a matriz para vetor
# quando há apenas 1 amostra (1 linha), o que quebrava pd() do picante
comm_matrix <- matrix(abund_vector, nrow = 1, dimnames = list("X6S", names(abund_vector)))

asvs_in_tree <- intersect(colnames(comm_matrix), tree$tip.label)
cat("ASVs de X6S presentes na árvore filogenética:", length(asvs_in_tree),
    "de", ncol(comm_matrix), "\n\n")

if (length(asvs_in_tree) < ncol(comm_matrix)) {
  cat("AVISO: nem todos os ASVs da amostra estão na árvore.",
      "Faith's PD será calculado só com o subconjunto presente.\n\n")
}

comm_matched <- match.phylo.comm(tree, comm_matrix)

# Proteção extra: garantir que comm continua sendo matriz 2D antes de pd()
comm_for_pd <- comm_matched$comm
if (is.null(dim(comm_for_pd))) {
  comm_for_pd <- matrix(comm_for_pd, nrow = 1,
                         dimnames = list("X6S", names(comm_for_pd)))
}
comm_for_pd <- comm_for_pd[, colnames(comm_for_pd), drop = FALSE]

faith_pd <- pd(comm_for_pd, comm_matched$phy, include.root = TRUE)

# --- 5. Consolidar resultados ---
results <- tibble(
  Amostra = "X6S",
  Observed_Richness = observed_richness,
  Chao1 = as.numeric(chao1),
  Shannon = shannon,
  Simpson_GiniSimpson = simpson,
  InverseSimpson = invsimpson,
  Pielou_Evenness = pielou_evenness,
  Faith_PD = faith_pd$PD[1],
  Faith_SR = faith_pd$SR[1]
)

cat("=== RESULTADOS ===\n")
print(as.data.frame(results))

# --- 6. Salvar output ---
write_tsv(results, "biodiversity_indices_X6S_calculated.tsv")
cat("\nResultados salvos em: biodiversity_indices_X6S_calculated.tsv\n")

# --- 7. Registro de reprodutibilidade ---
cat("\n=== sessionInfo() ===\n")
print(sessionInfo())
