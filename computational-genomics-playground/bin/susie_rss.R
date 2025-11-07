#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(susieR)
})

# ---------- helpers ----------
ts <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))

read_ld_square <- function(path, expect_n = NULL, verbose = TRUE) {
  if (verbose) cat(sprintf("[LD] reading: %s\n", path))
  dt <- data.table::fread(path, header = FALSE, fill = TRUE, na.strings = c("NA", "", "NaN"))
  if (verbose) cat(sprintf("[LD] raw dims from fread: %d x %d\n", nrow(dt), ncol(dt)))

  # drop trailing all-NA columns (sometimes appear from stray delimiters)
  all_na_cols <- which(vapply(dt, function(x) all(is.na(x)), logical(1)))
  if (length(all_na_cols)) {
    dt <- dt[, -all_na_cols, with = FALSE]
    if (verbose) cat(sprintf("[LD] dropped %d all-NA columns -> %d x %d\n",
                             length(all_na_cols), nrow(dt), ncol(dt)))
  }

  # if first col looks like an index/non-numeric, drop it
  if (ncol(dt) > 1) {
    c1 <- suppressWarnings(as.numeric(dt[[1]]))
    if (any(is.na(c1)) || identical(c1, seq_len(nrow(dt))) || all(c1 == nrow(dt))) {
      dt <- dt[, -1, with = FALSE]
      if (verbose) cat(sprintf("[LD] dropped first column (index/labels) -> %d x %d\n", nrow(dt), ncol(dt)))
    }
  }

  # coerce to numeric matrix
  mat <- as.matrix(as.data.frame(lapply(dt, function(x) as.numeric(x))))
  # drop any all-NA rows
  keep_rows <- which(rowSums(is.na(mat)) < ncol(mat))
  if (length(keep_rows) < nrow(mat)) {
    mat <- mat[keep_rows, , drop = FALSE]
    if (verbose) cat(sprintf("[LD] dropped %d all-NA rows -> %d x %d\n",
                             nrow(dt) - length(keep_rows), nrow(mat), ncol(mat)))
  }

  # ensure square (trim to top-left if needed, but warn)
  nr <- nrow(mat); nc <- ncol(mat)
  if (nr != nc) {
    k <- min(nr, nc)
    if (verbose) cat(sprintf("[LD] WARNING: non-square (%d x %d); trimming to %d x %d\n", nr, nc, k, k))
    mat <- mat[seq_len(k), seq_len(k), drop = FALSE]
  }

  if (!is.null(expect_n) && nrow(mat) != expect_n) {
    stop(sprintf("LD size (%d) does not match expected SNP count from Z/keep (%d).", nrow(mat), expect_n))
  }

  if (!all(is.finite(mat))) stop("LD contains non-finite values after cleaning.")
  diag(mat) <- 1
  mat
}

sanitize_ld_adaptive <- function(R) {
  # clamp, symmetrize, unit diag (dimension-safe)
  R[R >  1] <-  1
  R[R < -1] <- -1
  R <- (R + t(R)) * 0.5
  diag(R) <- 1

  Rs <- Matrix::forceSymmetric(R, uplo = "U")
  ok <- tryCatch({ chol(Rs); TRUE }, error = function(e) FALSE)
  if (ok) {
    cat("[LD] PD check: OK (lambda=0)\n")
    return(as.matrix(Rs))
  }

  lambda <- 1e-8
  repeat {
    ok <- tryCatch({ chol(Rs + Matrix::Diagonal(nrow(Rs), lambda)); TRUE }, error = function(e) FALSE)
    if (ok) break
    lambda <- lambda * 2
    if (lambda > 1e-1) stop("Unable to stabilize LD with ridge (lambda > 1e-1).")
  }
  cat(sprintf("[LD] PD check: used adaptive ridge lambda=%.2e\n", lambda))
  as.matrix(Rs + Matrix::Diagonal(nrow(Rs), lambda))
}

# ---------- args ----------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: susie_rss.R z.tsv[.gz] R.txt[.gz] kept_snps.txt out_prefix [L] [n]")
}
zfile   <- args[1]
rfile   <- args[2]
keepfile<- args[3]
outpref <- args[4]
L       <- ifelse(length(args) >= 5, as.integer(args[5]), 5)
n_arg   <- if (length(args) >= 6) as.integer(args[6]) else NA_integer_

# ---------- run ----------
ts("Starting")

# Z & keep
ts(paste("Reading Z:", zfile))
t1 <- Sys.time()
z <- data.table::fread(zfile)
t2 <- Sys.time(); ts(sprintf("Read Z in %.2f sec", as.numeric(t2 - t1, "secs")))

if (!("SNP" %in% names(z))) stop("Z file must have a 'SNP' column.")
if (!("Z" %in% names(z)))  stop("Z file must have a numeric 'Z' column.")

ts(paste("Reading keep list:", keepfile))
keep <- data.table::fread(keepfile, header = FALSE)$V1

ts("Aligning Z to 'keep' order (strict match)")
t1 <- Sys.time()
z_aligned <- merge(data.table(SNP = keep), z, by = "SNP", sort = FALSE)
if (nrow(z_aligned) != length(keep)) {
  missing <- setdiff(keep, z$SNP)
  stop(sprintf("Z alignment failed: %d SNPs from keep not found in Z. Example missing: %s",
               length(missing), paste(head(missing, 5), collapse = ",")))
}
if (!is.numeric(z_aligned$Z) || any(!is.finite(z_aligned$Z))) {
  stop("Z column has non-finite values after alignment.")
}
z <- z_aligned
t2 <- Sys.time(); ts(sprintf("Aligned in %.2f sec", as.numeric(t2 - t1, "secs")))

# infer n if not provided
n <- n_arg
if (is.na(n)) {
  if ("N" %in% names(z)) {
    n <- suppressWarnings(as.integer(stats::median(z$N, na.rm = TRUE)))
    ts(sprintf("Inferred n from z$N: n = %s", ifelse(is.na(n), "NA", n)))
  } else {
    ts("No n provided and no N column; proceeding with n=NA (SuSiE will warn).")
  }
}

# LD
ts(paste("Reading LD:", rfile))
t1 <- Sys.time()
R <- read_ld_square(rfile, expect_n = nrow(z), verbose = TRUE)
t2 <- Sys.time(); ts(sprintf("LD read & basic clean in %.2f sec", as.numeric(t2 - t1, "secs")))
cat(sprintf("[LD] after read: %d x %d (class=%s)\n", nrow(R), ncol(R), class(R)))

# sanitize to PD (adaptive ridge only if needed)
ts("Sanitizing LD (clamp/symmetrize + adaptive ridge if needed)")
t1 <- Sys.time()
R <- sanitize_ld_adaptive(R)
t2 <- Sys.time(); ts(sprintf("Sanitized LD in %.2f sec", as.numeric(t2 - t1, "secs")))

# SuSiE-RSS
zscores <- z$Z
ts(sprintf("Running SuSiE-RSS (L=%d, n=%s, estimate_residual_variance=FALSE)", L, ifelse(is.na(n),"NA",n)))
t1 <- Sys.time()
fit <- susie_rss(
  zscores, R,
  L = L,
  n = n,                              # recommended with reference LD
  coverage = 0.90,
  max_iter = 500,
  track_fit = FALSE,
  estimate_residual_variance = FALSE  # key for reference LD
)
t2 <- Sys.time(); ts(sprintf("Finished SuSiE in %.2f sec", as.numeric(t2 - t1, "secs")))

# Outputs
ts("Building PIP table")
pp <- data.table(SNP = z$SNP, PP = colSums(fit$alpha))
data.table::fwrite(pp, paste0(outpref, ".pip.tsv"), sep = "\t")


ts("Extracting credible sets")
pip <- colSums(fit$alpha)
fwrite(data.table(SNP = z$SNP, PP = pip), paste0(outpref, ".pip.tsv"), sep = "\t")

csobj <- susie_get_cs(fit, Xcorr = R, coverage = 0.95, min_abs_corr = 0.0)

if (length(csobj$cs) > 0) {
  # find component index for each CS, trying several field names for compatibility
  get_cs_comp <- function(j) {
    if (!is.null(csobj$cs_index))       return(csobj$cs_index[j])
    if (!is.null(csobj$index))          return(csobj$index[j])
    if (!is.null(fit$sets$cs_index))    return(fit$sets$cs_index[j])
    # fallback: assume CS j corresponds to component j (may be wrong, but better than crash)
    return(j)
  }

  # ensure we know p and L
  p <- nrow(as.matrix(fit$alpha))
  L <- ncol(as.matrix(fit$alpha))

  out_cs <- rbindlist(lapply(seq_along(csobj$cs), function(j) {
    idx <- csobj$cs[[j]]
    idx <- idx[idx >= 1 & idx <= p]  # guard against any stray indices
    comp <- get_cs_comp(j)
    comp <- max(1, min(L, comp))     # clamp just in case

    # extract PP for the component; handle L=1 (alpha can drop to a vector)
    alpha_mat <- as.matrix(fit$alpha) # safe matrix view
    pp_j <- alpha_mat[idx, comp, drop = TRUE]

    data.table(
      CS = j,
      Component = comp,
      Coverage = csobj$coverage[j],
      SNP = z$SNP[idx],
      PP = pp_j
    )
  }), fill = TRUE)

  fwrite(out_cs, paste0(outpref, ".credible_sets.tsv"), sep = "\t")
} else {
  # naive fallback by cumulative PP
  ord <- order(pip, decreasing = TRUE)
  cum <- cumsum(pip[ord]); k <- which(cum >= 0.95)[1]
  take <- ifelse(is.na(k), length(ord), k)
  naive <- data.table(
    CS = 1L, Coverage = 0.95,
    SNP = z$SNP[ord[seq_len(take)]],
    PP  = pip[ord[seq_len(take)]]
  )
  fwrite(naive, paste0(outpref, ".credible_sets.tsv"), sep = "\t")
}


# small model artifact (optional)
saveRDS(list(alpha = fit$alpha, sets = csobj), paste0(outpref, ".susie_fit.rds"))

ts("Done")
