# =============================================================================
# Regenera as figuras do manuscrito somente em TIFF em figures/final_tiff/ (21/10/2026; README seção 16).
# Uso (a partir de ~/sbria_anammox_shotgun):
#   Rscript results/review_checks/scripts/render_final_tiff.R <script.R> <Figure_X.tiff> <largura_final_mm>
# Roda o script original sem alterá-lo, com duas substituições no ambiente global:
#   - ggsave(): só a chamada .tif/.tiff é atendida, e o arquivo vai para figures/final_tiff/<Figure_X.tiff>;
#     PDF/PNG são ignorados (nada é gravado em figures/final/);
#   - Sys.which(): devolve "" (o script da Fig 2 não roda o sips sobre figures/final/Fig2_mapa_SNAD.tiff).
# Tamanho: o layout do script (largura w_mm) é mantido e o dispositivo é reescalado para a largura final
# (ragg `scaling` = largura_final / w_mm): mesma figura, texto e linhas multiplicados pelo mesmo fator.
# Formato (Elsevier): TIFF LZW, RGB, 600 dpi, fundo branco opaco.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
SCRIPT <- args[1]; OUT <- file.path("figures/final_tiff", args[2]); FINAL_W <- as.numeric(args[3])
dir.create("figures/final_tiff", showWarnings = FALSE)
.n_tiff <- 0

ggsave <- function(filename, plot = ggplot2::last_plot(), device = NULL, path = NULL, scale = 1,
                   width = NA, height = NA, units = "mm", dpi = 300, ...) {
  if (!grepl("\\.tiff?$", filename, ignore.case = TRUE)) {
    message("[render_final_tiff] ignorado: ", filename); return(invisible(NULL))
  }
  stopifnot(units == "mm", is.finite(width), is.finite(height))
  s <- FINAL_W / width
  ragg::agg_tiff(OUT, width = FINAL_W, height = height * s, units = "mm", res = 600, scaling = s,
                 compression = "lzw", background = "white")
  print(plot)
  grDevices::dev.off()
  .n_tiff <<- .n_tiff + 1
  message(sprintf("[render_final_tiff] %s <- %s: layout %.1f x %.1f mm, escala %.4f, final %.1f x %.1f mm, 600 dpi",
                  OUT, filename, width, height, s, FINAL_W, height * s))
  invisible(OUT)
}
Sys.which <- function(names) setNames(rep("", length(names)), names)

source(SCRIPT, echo = FALSE)
stopifnot(.n_tiff == 1)
