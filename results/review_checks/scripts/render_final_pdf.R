# =============================================================================
# Versão em PDF das figuras do manuscrito em figures/final_pdf/ (README seção 22). Par do render_final_tiff.R.
# Uso (a partir de ~/sbria_anammox_shotgun):
#   Rscript results/review_checks/scripts/render_final_pdf.R <script.R> <Figure_X.pdf> <largura_final_mm>
# Roda o script original sem alterá-lo. O ggsave é interceptado: na chamada TIFF do script, o mesmo gráfico é desenhado em PDF
# (dispositivo Quartz do macOS, que embute as fontes) no tamanho de layout do script; as demais chamadas (PDF/PNG) são ignoradas.
# Depois o PDF é reescalado como vetor para a largura final (pypdf, scale_pdf.py), igual à escala usada nos TIFF.
# Sys.which() devolve "" (o script da Fig 2 não roda o sips sobre figures/final/).
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
SCRIPT <- args[1]; OUT <- file.path("figures/final_pdf", args[2]); FINAL_W <- as.numeric(args[3])
dir.create("figures/final_pdf", showWarnings = FALSE)
PY <- Sys.getenv("PDF_PY")                     # python com pypdf
SCALER <- "results/review_checks/scripts/scale_pdf.py"
.n_pdf <- 0

ggsave <- function(filename, plot = ggplot2::last_plot(), device = NULL, path = NULL, scale = 1,
                   width = NA, height = NA, units = "mm", dpi = 300, ...) {
  if (!grepl("\\.tiff?$", filename, ignore.case = TRUE)) {
    message("[render_final_pdf] ignorado: ", filename); return(invisible(NULL))
  }
  stopifnot(units == "mm", is.finite(width), is.finite(height))
  tmp <- tempfile(fileext = ".pdf")
  grDevices::quartz(type = "pdf", file = tmp, width = width / 25.4, height = height / 25.4, bg = "white")
  print(plot)
  grDevices::dev.off()
  st <- system2(PY, c(SCALER, shQuote(tmp), shQuote(OUT), FINAL_W), stdout = TRUE)
  message("[render_final_pdf] ", OUT, " <- ", filename, sprintf(": layout %.1f x %.1f mm; ", width, height), paste(st, collapse = " "))
  .n_pdf <<- .n_pdf + 1
  invisible(OUT)
}
Sys.which <- function(names) setNames(rep("", length(names)), names)

source(SCRIPT, echo = FALSE)
stopifnot(.n_pdf == 1)
