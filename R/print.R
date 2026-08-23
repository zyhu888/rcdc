print.robustCDC <- function(x, ...) {
  cat("RCDC fit\n")
  cat("Selected K:", x$K, "\n")
  cat("WAIC:", x$WAIC, "\n")
  if (!is.null(x$waic_table)) {
    cat("\nWAIC comparison:\n")
    print(x$waic_table)
  }
  invisible(x)
}
