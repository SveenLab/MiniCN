#' @importFrom magrittr %>%
#' @importFrom rlang .data

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(BSgenome)
  library(Biostrings)
  library(GenomicFeatures)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ggforce)
})

options(stringsAsFactors = FALSE)

# ----------------------------- Utilities -------------------------------------

# Read coverage file (requires target + total_coverage in header)
read_cov <- function(path){
  if(!file.exists(path)) stop(paste("File not found:", path))
  dt <- data.table::fread(path, sep = "\t", header = TRUE)
  orig <- names(dt)
  data.table::setnames(dt, tolower(names(dt)))
  # normalize Target column name to 'target'
  if ("interval_id" %in% names(dt)) data.table::setnames(dt, "interval_id", "target")
  if ("interval" %in% names(dt)) data.table::setnames(dt, "interval", "target")
  if (!("target" %in% names(dt)))
    stop("Coverage must have a 'Target' column (chr:start-end). Found columns: ", paste(orig, collapse=", "))
  # normalize coverage column to 'total_coverage'
  if (!("total_coverage" %in% names(dt))) {
    cand <- intersect(c("totalcount","totalcounts","count","counts","coverage","totalcoverage"), names(dt))
    if (length(cand)) data.table::setnames(dt, cand[1], "total_coverage")
  }
  if (!("total_coverage" %in% names(dt)))
    stop("Coverage must have a 'total_coverage' column (or a recognized alias).")
  dt
}

# Parse "chr:start-end" into data.frame(chr,start,end)
parse_target <- function(x){
  x <- trimws(x)
  m <- regexec("^([^:]+):([0-9]+)(?:-([0-9]+))?$", x)
  regm <- regmatches(x, m)

  bad <- which(lengths(regm) == 0L)
  if (length(bad)) {
    stop("Invalid Target format at rows: ", paste(bad, collapse = ", "),
         ". Expected 'chr:start' or 'chr:start-end'.")
  }

  chr <- vapply(regm, function(g) g[2], character(1))
  start <- as.integer(vapply(regm, function(g) g[3], character(1)))
  end <- as.integer(vapply(regm, function(g) if (length(g) >= 4 && nzchar(g[4])) g[4] else g[3], character(1)))

  data.frame(chr = chr, start = start, end = end, stringsAsFactors = FALSE)
}

# Vectorized GC fraction for GRanges
gc_for_ranges <- function(gr, genome){
  seqs <- BSgenome::getSeq(genome, gr)
  rowSums(Biostrings::letterFrequency(seqs, letters = c("G","C"), as.prob = TRUE))
}

# Read panel gene list (HGNC symbols only)
read_panel_genes <- function(path) {
  if (is.null(path) || !nzchar(path)) stop("Panel gene list is required (path not provided).")
  if (!file.exists(path)) stop(paste("Panel gene file not found:", path))
  # try as table
  dt <- tryCatch(suppressWarnings(data.table::fread(path, header = FALSE)), error = function(e) NULL)
  syms <- NULL
  if (!is.null(dt)) {
    nm <- tolower(names(dt))
    pick <- intersect(nm, c("gene","symbol","genes","symbols"))
    if (length(pick) >= 1) syms <- dt[[ pick[1] ]]
    else if (ncol(dt) == 1) syms <- dt[[1]]
  }
  if (is.null(syms)) {
    # fallback: one symbol per line
    syms <- scan(path, what = character(), quiet = TRUE, sep = "\n")
  }
  syms <- trimws(syms)
  if (!length(syms)) stop("Panel gene list parsed but empty.")

  # Return clean, unique SYMBOLS and their ENTREZ mapping
  syms <- unique(syms)
  # Enforce "symbols only": map SYMBOL -> ENTREZ; if any fail, error
  eg <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = syms,
                              keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
  bad <- is.na(eg)
  if (any(bad)) {
    stop("The following entries are not valid HGNC symbols (no mapping in org.Hs.eg.db): ",
         paste(unique(syms[bad]), collapse = ", "))
  }

  list(symbols = syms, entrez = unname(eg))
}

# Check if amplicon is within selected gene transcripts
annotate_to_panel <- function(gr, panel_entrez) {
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

  # all per-gene loci (can be multi-part across strands/seqlevels)
  gene_gl <- GenomicFeatures::genes(txdb, single.strand.genes.only = FALSE)
  keep <- names(gene_gl) %in% panel_entrez
  gene_gl <- gene_gl[keep]

  gene_gr <- unlist(gene_gl, use.names = FALSE)
  gene_id <- rep(names(gene_gl), lengths(gene_gl))

  out <- rep(NA_character_, length(gr))

  # Overlap to panel transcripts (max overlap per query)
  if (length(gene_gr)) {
    hits <- GenomicRanges::findOverlaps(gr, gene_gr, ignore.strand = TRUE)
    if (length(hits)) {
      ovw <- IRanges::pintersect(GenomicRanges::ranges(gr)[S4Vectors::queryHits(hits)],
                                 GenomicRanges::ranges(gene_gr)[S4Vectors::subjectHits(hits)])

      dfh <- data.frame(q = S4Vectors::queryHits(hits),
                       s = S4Vectors::subjectHits(hits),
                       ov = IRanges::width(ovw))

      dfh <- dfh[order(dfh$q, -dfh$ov), , drop = FALSE]
      best <- dfh[!duplicated(dfh$q), , drop = FALSE]
      out[best$q] <- gene_id[best$s]
    }
  }

  out # ENTREZ IDs or NA
}

# Build regions (metadata) from target + panel, and drop intervals not mapping to panel
build_regions_from_target_panel <- function(target_vec, genome, panel_symbols, panel_entrez) {
  coords <- parse_target(target_vec)
  gr <- GenomicRanges::GRanges(seqnames = coords$chr, ranges = IRanges::IRanges(coords$start, coords$end))
  # normalize naming style (add/remove 'chr' as needed)
  try({ GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC" }, silent = TRUE)

  # Annotate to panel (ENTREZ)
  eg <- annotate_to_panel(gr, panel_entrez)

  # Map ENTREZ -> SYMBOL for output
  sym_by_entrez <- suppressMessages(
    AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                          keys = unique(stats::na.omit(eg)),
                          keytype = "ENTREZID",
                          column = "SYMBOL",
                          multiVals = "first"))
  gene_sym <- ifelse(is.na(eg), NA_character_, unname(sym_by_entrez[eg]))

  # Compose DF
  df <- tibble::tibble(
    Target = target_vec,
    chr = as.character(GenomicRanges::seqnames(gr)),
    start = GenomicRanges::start(gr),
    end = GenomicRanges::end(gr),
    len = GenomicRanges::end(gr) - GenomicRanges::start(gr) + 1L,
    GC_content = gc_for_ranges(gr, genome),
    gene = gene_sym
  )

  # Drop NA genes (not in panel / too far)
  kept <- !is.na(df$gene)
  dropped_n <- sum(!kept)
  if (dropped_n) message("Excluded ", dropped_n, " targets not mapping to panel genes.")
  df <- df[kept, , drop = FALSE]

  if (!nrow(df)) stop("After panel filtering, no targets remain. Check panel list.")

  # Order targets within gene (for plotting)
  df <- df %>%
    dplyr::group_by(.data$gene) %>%
    dplyr::arrange(.data$chr, .data$start, .by_group = TRUE) %>%
    dplyr::mutate(exon_index = dplyr::row_number(),
                  amplicon = paste0("amplicon", .data$exon_index),
                  gene_exon = factor(paste0(.data$gene, "_", .data$amplicon),
                                     levels = unique(paste0(.data$gene, "_", .data$amplicon)))) %>%
    dplyr::ungroup()

  df
}

# longest contiguous run of TRUEs (for focal deletions)
rle_max_run <- function(x) {
  if (!length(x)) return(0L)
  r <- rle(as.logical(x))
  if (!any(r$values)) return(0L)
  max(r$lengths[r$values])
}

# read all normal coverage files in a folder and compute per-target median and
# build per-target baseline SD from a folder of normals using the same Targets
pon_baseline_from_folder <- function(folder, targets, pc = 0.5) {
  files <- list.files(folder, pattern = "_coverage\\.sample_interval_summary$", full.names = TRUE)
  if (length(files) < 2) stop("PON folder must contain >=2 normal coverage files: ", folder)

  # read and align coverage
  read_aligned <- function(f) {
    d <- read_cov(f)
    key <- if ("target" %in% names(d)) d$target else d$Target
    idx <- match(targets, key)
    if (anyNA(idx)) stop("Normal coverage missing some Targets: ", basename(f))
    d$total_coverage[idx]
  }

  mat <- sapply(files, read_aligned)  # rows=targets, cols=normals
  if (is.vector(mat)) mat <- matrix(mat, ncol = 1L)

  # PON median coverage per target (used for normalization)
  pon_median <- apply(mat, 1, stats::median, na.rm = TRUE)

  # per-normal log2 ratios vs PON median (centers mu around 0 per target)
  log2R_normals <- apply(mat, 2, function(x) {
    log2((x + pc) / (pon_median + pc))
  })

  # robust per-target SD with floor
  sd_norm <- apply(log2R_normals, 1, function(v)
    stats::mad(v, constant = 1.4826, na.rm = TRUE)
  )
  sd_floor <- 0.05
  sd_norm <- pmax(sd_norm, sd_floor)

  list(
    targets = targets,
    median  = pon_median,
    mu_norm = rep(0, length(targets)),
    sd_norm = sd_norm
  )
}

# When only a single matched normal is available, approximate sigma with a single
# within-sample SD from the tumor's neutral bins (10-90% middle)
single_sample_sigma <- function(log2R, targets) {
  qr <- stats::quantile(log2R, c(0.1, 0.9), na.rm = TRUE)
  neut <- which(is.finite(log2R) & log2R >= qr[1] & log2R <= qr[2])
  sd0  <- stats::mad(log2R[neut], constant = 1.4826, na.rm = TRUE)
  sd_floor <- 0.05
  data.frame(Target = targets,
             mu_norm = 0,
             sd_norm = rep.int(max(sd0, sd_floor), length(targets)),
             stringsAsFactors = FALSE)
}
# --------------------------- Main runner -------------------------------------

#' Run MiniCN copy-number caller
#'
#' @param sample_csv CSV with columns: Tumor Folder, Tumor File Name, Normal Folder, Normal File Name
#' @param outdir output directory
#' @param genome_pkg BSgenome package name, (defualt "BSgenome.Hsapiens.UCSC.hg38")
#' @param panel_path path to panel gene list (HGNC symbols, required)
#' @param pc pseudocount for normalization (default 0.5)
#' @param HIGH_AMP_T high-amplification copy number threshold (default log2(17/2) or ~3.09)
#' @param LOW_AMP_T low-amplification copy number threshold (default log2(7/2) or ~1.81)
#' @param GAIN_T gain copy number threshold (default log2(3/2) or ~0.58)
#' @param LOSS_T loss copy number threshold (default log2(1/2) or -1)
#' @param DEEP_DEL_T deletion copy number threshold (default -2.0)
#' @param MIN_AMPLICONS minimum number of amplicons supporting respective changes (default 2)
#' @param MIN_GENE_Z minimum z-score cutoff for gains and losses (default 3)
#' @param MAX_Q_VALUE maximum q-value cutoff for gains and losses (default 0.05)
#'
#' @return merged tsv files and per sample pdf files
#' @export
#'
#' @examples
#' # Package example files
#' sample_csv <- system.file("extdata", "sample_sheet.csv", package = "MiniCN")
#' panel_file <- system.file("extdata", "panel_genes.txt", package = "MiniCN")
#'
#' # Temporary output directory for examples
#' out <- file.path(tempdir(), "miniCN_example")
#' dir.create(out, showWarnings = FALSE, recursive = TRUE)
#'
#' # Skip if genome package is not available on this machine
#' if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
#'   run_caller(sample_csv = sample_csv, outdir = out, panel_path = panel_file)
#' }
run_caller <- function(sample_csv, outdir,
                       genome_pkg = "BSgenome.Hsapiens.UCSC.hg38",
                       panel_path = NULL,
                       pc = 0.5,
                       HIGH_AMP_T = log2(17/2),
                       LOW_AMP_T = log2(7/2),
                       GAIN_T = log2(3/2),
                       LOSS_T = log2(1/2),
                       DEEP_DEL_T = -2.0,
                       MIN_AMPLICONS  = 2L,
                       MIN_GENE_Z = 3L,
                       MAX_Q_VALUE = 0.05) {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # Load genome package and object
  suppressPackageStartupMessages(library(genome_pkg, character.only = TRUE))
  genome <- get(genome_pkg)

  # Panel genes (SYMBOLS only; error if any not valid symbols)
  panel <- suppressMessages(read_panel_genes(panel_path))
  panel_symbols <- panel$symbols
  panel_entrez  <- panel$entrez

  # Read sample sheet
  samples <- data.table::fread(sample_csv)
  csv_dir <- normalizePath(dirname(sample_csv), mustWork = TRUE)
  need_cols <- c("Tumor Folder", "Tumor File Name", "Normal Folder", "Normal File Name")
  if (!all(need_cols %in% names(samples))) {
    stop("Sample CSV must contain columns: ", paste(need_cols, collapse = ", "))
  }

  for (i in seq_len(nrow(samples))) {
    sample <- samples[i, ]
    sample_id <- sample$`Tumor File Name`

    tumor_path  <- file.path(csv_dir,
                             sample$`Tumor Folder`,
                             paste0(sample$`Tumor File Name`, "_coverage.sample_interval_summary"))
    normal_path <- file.path(csv_dir,
                             sample$`Normal Folder`,
                             paste0(sample$`Normal File Name`, "_coverage.sample_interval_summary"))

    tumor  <- read_cov(tumor_path)

    # Build regions from tumor Targets with strict panel filtering
    target_vec <- if ("target" %in% names(tumor)) tumor$target else tumor$Target
    regions <- build_regions_from_target_panel(target_vec, genome,
                                               panel_symbols = panel_symbols,
                                               panel_entrez  = panel_entrez)

    # Align tumor rows to kept Targets
    t_all <- if ("target" %in% names(tumor)) tumor$target else tumor$Target
    idxT <- match(regions$Target, t_all); if (anyNA(idxT)) stop("Tumor missing kept Targets.")
    t_keep <- tumor[idxT, ]
    tvec <- t_keep$total_coverage

    is_pon <- identical(sample$`Normal File Name`, "PON_median")

    if (is_pon) {
      # PON mode: build baseline once and use its median for ratios
      pon_folder <- file.path(csv_dir, sample$`Normal Folder`)
      pon <- pon_baseline_from_folder(pon_folder, regions$Target, pc = pc)
      nvec <- pon$median
    } else {
      # Matched-normal mode
      normal <- read_cov(normal_path)
      n_all <- if ("target" %in% names(normal)) normal$target else normal$Target
      idxN  <- match(regions$Target, n_all); if (anyNA(idxN)) stop("Normal missing kept Targets.")
      nvec <- normal$total_coverage[idxN]
    }

    # Record excluded targets (for audit)
    dropped <- setdiff(t_all, regions$Target)
    if (length(dropped)) {
      data.table::fwrite(data.table::data.table(Sample = sample_id, Target = dropped),
             file = file.path(outdir, "0_Excluded_Targets.tsv"),
             sep = "\t")
    }

    # Ratio-first normalization
    raw_ratio <- (tvec + pc) / (nvec + pc)
    qr <- stats::quantile(raw_ratio, c(0.1, 0.9), na.rm=TRUE)
    neut_idx <- which(is.finite(raw_ratio) & raw_ratio >= qr[1] & raw_ratio <= qr[2])
    sf <- stats::median(raw_ratio[neut_idx], na.rm=TRUE)
    log2_raw <- log2(raw_ratio / sf)

    df <- cbind(regions, log2_raw = log2_raw)

    # Robust LOESS GC correction
    mid <- abs(log2_raw - stats::median(log2_raw, na.rm=TRUE)) <= stats::mad(log2_raw, na.rm=TRUE) * 2.5
    fit_idx <- is.finite(df$GC_content) & is.finite(log2_raw) & mid

    lo <- stats::loess(log2_raw ~ GC_content,
                       data = df[fit_idx,],
                       family="symmetric",
                       control = stats::loess.control(surface="direct"),
                       degree=2,
                       span=0.65)

    trend <- stats::predict(lo, newdata = data.frame(GC_content = df$GC_content))
    trend[!is.finite(trend)] <- 0
    df$log2_gc <- log2_raw - trend

    df <- cbind(sample = sample_id, df)
    data.table::fwrite(df, file = file.path(outdir, "1_Amplicons_CN.tsv"), sep = "\t", append = TRUE)

    # QC summary
    qc <- data.frame(
      sample = sample_id,
      median_depth_tumor  = stats::median(tvec),
      median_depth_normal = stats::median(nvec),
      kept_targets = nrow(df),
      excluded_targets = length(dropped),
      MAD_raw = stats::mad(df$log2_raw, na.rm = TRUE),
      MAD_gc = stats::mad(df$log2_gc, na.rm = TRUE)
    )
    data.table::fwrite(qc, file = file.path(outdir, "2_Samples_QC.tsv"), sep = "\t", append = TRUE)

    # z-score baseline (per-target if PON, single-sigma if matched normal)
    kept_targets <- df$Target
    if (is_pon) {
      baseline <- data.frame(
        Target = pon$targets,
        mu_norm = pon$mu_norm,
        sd_norm = pon$sd_norm
      )
    } else {
      baseline <- single_sample_sigma(df$log2_gc, kept_targets)
    }

    df <- df %>%
      dplyr::left_join(baseline, by = "Target") %>%
      dplyr::mutate(sd_norm = pmax(.data$sd_norm, 0.05),
                    z_i = (.data$log2_gc - .data$mu_norm) / .data$sd_norm,
                    w_i = 1 / .data$sd_norm)

    # Gene-level aggregation (length-weighted)
    gene_stats <- df %>%
      dplyr::arrange(.data$gene, .data$start) %>%
      dplyr::group_by(.data$gene) %>%
      dplyr::summarize(
        sample = .data$sample[1],
        n_targets = dplyr::n(),
        mu = stats::weighted.mean(.data$log2_gc, w = .data$len),
        mad_gene = stats::mad(.data$log2_gc, constant = 1.4826),
        Z_gene = sum(.data$w_i * .data$z_i, na.rm = TRUE) /
          sqrt(sum((.data$w_i)^2, na.rm = TRUE)),

        # support counts & contiguity for amplifications & deletions
        n_sup_lowamp = sum(.data$log2_gc >= LOW_AMP_T),
        n_sup_highamp = sum(.data$log2_gc >= HIGH_AMP_T),
        max_run_del = rle_max_run(.data$log2_gc <= DEEP_DEL_T),

        # directional consistency for gains/losses
        n_sup_gain = sum(.data$log2_gc >= GAIN_T),
        n_sup_loss = sum(.data$log2_gc <= LOSS_T),

        .groups = "drop"
      ) %>%
      dplyr::mutate(
        p_two_sided = 2 * stats::pnorm(-abs(.data$Z_gene)),
        p_gain = 1 - stats::pnorm(.data$Z_gene),
        p_loss = stats::pnorm(.data$Z_gene),
        q_gain = stats::p.adjust(.data$p_gain, method = "BH"),
        q_loss = stats::p.adjust(.data$p_loss, method = "BH"),

        cn_ratio = 2^.data$mu,
        diff_copies = 2*(.data$cn_ratio - 1),

        call = dplyr::case_when(
          # Amplifications (require at least 2 amplicons over copy-based thresholds)
          .data$mu >= HIGH_AMP_T & .data$n_sup_highamp >= MIN_AMPLICONS ~ "High amplification",
          .data$mu >= LOW_AMP_T & .data$n_sup_lowamp >= MIN_AMPLICONS ~ "Low amplification",

          # Focal deletion (contiguous run of -2.0 across at least 2 adjacent amplicons)
          .data$max_run_del >= MIN_AMPLICONS ~ "Deletion",

          # Single-copy LOSS & GAIN
          (.data$Z_gene <= -MIN_GENE_Z & .data$q_loss < MAX_Q_VALUE &
             .data$n_sup_loss >= MIN_AMPLICONS & .data$mu <= LOSS_T) ~ "Loss",
          (.data$Z_gene >= MIN_GENE_Z & .data$q_gain < MAX_Q_VALUE &
             .data$n_sup_gain >= MIN_AMPLICONS & .data$mu >= GAIN_T) ~ "Gain",

          TRUE ~ "Neutral"
        )
      )

    # Save gene results
    data.table::fwrite(gene_stats, file = file.path(outdir, "3_Gene_CN_Final.tsv"), sep = "\t", append = TRUE)

    # Plot CN for all genes
    axisdf <- df %>%
      dplyr::arrange(.data$gene, .data$start) %>%
      dplyr::group_by(.data$gene) %>%
      dplyr::mutate(idx = dplyr::row_number()) %>%
      dplyr::summarize(center = .data$gene_exon[ceiling(max(.data$idx, na.rm = TRUE) / 2)],
                       .groups="drop")

    df$gene_order <- match(df$gene, axisdf$gene)
    df$alt <- (df$gene_order %% 3) + 1

    p1 <- ggplot2::ggplot(df,
                          ggplot2::aes(x=.data$gene_exon, y=.data$log2_gc, color=.data$gene)) +
      ggplot2::geom_point(ggplot2::aes(color = factor(.data$alt), shape = factor(.data$alt)),
                          alpha = 0.8, size = 2) +
      ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) +
      ggplot2::scale_shape_manual(values = c(16, 17, 18)) +
      ggplot2::scale_x_discrete(labels = axisdf$gene, breaks = axisdf$center) +
      ggplot2::labs(x=NULL, y="Log2 (GC-corrected CN Ratio)") +
      ggplot2::theme_bw(base_size=12) +
      ggplot2::theme(legend.position="none") +
      ggplot2::theme(legend.position = "none",
                     panel.border = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle=45, hjust=1, size=12),
                     axis.text.y = ggplot2::element_text(size=12)) +
      ggplot2::coord_cartesian(ylim = c(-6,6))

    # Plot CN for significant genes
    sign_CN <- gene_stats[gene_stats$call != "Neutral",]

    if (nrow(sign_CN) > 0) {
      sign_df <- df[df$gene %in% sign_CN$gene, , drop = FALSE]
      sign_df <- merge(sign_df, sign_CN[, c("gene","call")], by = "gene", all.x = TRUE)
      sign_df$copies_dev <- 2 * (2^sign_df$log2_gc - 1)
      sign_df$facet_lab <- paste0(sign_df$gene, " (", sign_df$call, ")")

      # Per-facet upper limit: max(2, max deviation for that gene)
      limits <- sign_df %>%
        dplyr::group_by(.data$facet_lab) %>%
        dplyr::summarise(L = max(2, max(.data$copies_dev, na.rm = TRUE)),
                         .groups="drop") %>%
        # pick any Target present in that facet to satisfy x mapping for geom_blank
        dplyr::left_join(
          sign_df %>%
            dplyr::group_by(.data$facet_lab) %>%
            dplyr::summarise(Target = dplyr::first(.data$Target), .groups="drop"),
          by = "facet_lab"
        )

      lim_low  <- limits; lim_low$y  <- -2
      lim_high <- limits; lim_high$y <- limits$L

      # colors
      call_levels <- c("Loss", "Gain", "Deletion", "Low amplification", "High amplification")
      sign_df$call <- factor(sign_df$call, levels = call_levels)

      call_colors <- c(
        "Loss" = "#377eb8",
        "Gain" = "#4daf4a",
        "Deletion" = "#984ea3",
        "Low amplification" = "#ff7f00",
        "High amplification"= "#e41a1c"
      )

      # page layout: 2 rows x 3 cols
      ncol_fac <- 3L
      nrow_fac <- 2L
      per_page <- ncol_fac * nrow_fac
      n_facets <- length(unique(sign_df$facet_lab))
      n_pages <- ceiling(n_facets / per_page)

      # Build base plot
      p2_base <- ggplot2::ggplot(sign_df,
                                 ggplot2::aes(x = .data$Target, y = .data$copies_dev,
                                     group = .data$gene, color = .data$call)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
        ggplot2::geom_line(linewidth = 0.3) +
        ggplot2::geom_point(size = 1.6) +
        ggplot2::scale_color_manual(values = call_colors, drop = FALSE) +
        ggplot2::geom_blank(data = lim_low,
                            ggplot2::aes(x = .data$Target, y = .data$y), inherit.aes = FALSE) +
        ggplot2::geom_blank(data = lim_high,
                            ggplot2::aes(x = .data$Target, y = .data$y), inherit.aes = FALSE) +
        ggforce::facet_wrap_paginate(~ facet_lab, ncol = ncol_fac, nrow = nrow_fac,
                                     page = 1, scales = "free") +
        ggplot2::theme_minimal(base_size = 12, base_family = "sans") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                                           vjust = 0.5, size = 9),
                       legend.position = "none",
                       panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::xlab("") +
        ggplot2::ylab("Copies deviation")

      grDevices::pdf(file.path(outdir, paste0(sample_id, "_CN_plots.pdf")), width = 11, height = 8)

      print(p1)

      for (pg in seq_len(n_pages)) {
        p2_page <- p2_base +
          ggforce::facet_wrap_paginate(~ facet_lab, ncol = ncol_fac, nrow = nrow_fac,
                                       page = pg, scales = "free")
        print(p2_page)
      }
      grDevices::dev.off()

    } else {
      # No significant genes -> only plot 1
      grDevices::pdf(file.path(outdir, paste0(sample_id, "_CN_plots.pdf")), width = 11, height = 8)
      print(p1)
      grDevices::dev.off()
    }
  }
}
