# ==========================================================================
# 01_prepare_data.R  --  BIO463 within-host evolution app
#
# Runs ONCE on your laptop (full R, all packages). Does every heavy or
# package-dependent step from Norovirus-Prolonged-Infection.Rmd and writes
# small plain-data files into app/data/.
#
# The Shiny app then needs ONLY ggplot2 + dplyr, so it can run in the
# browser via shinylive (no ape / ggtree / treedater / seqinr in WASM).
#
# Usage:  Rscript prep/01_prepare_data.R
# ==========================================================================

suppressPackageStartupMessages({
  library(ape)
  library(seqinr)
  library(lubridate)
  library(treedater)
})

set.seed(463)                      # reproducible Sydney subsample
OUT  <- "app/data"        # small .rds files the app reads (goes into the WASM bundle)
DATA <- "site/data"       # bigger files offered as plain downloads
dir.create(OUT,  recursive = TRUE, showWarnings = FALSE)
dir.create(DATA, recursive = TRUE, showWarnings = FALSE)

NORO <- path.expand("~/Library/CloudStorage/Dropbox/Norovirus_GII.4")
P33_FASTA <- file.path(NORO, "prolonged-infection-evolution-tutorial/P33_capsid.fasta")
POP_FASTA <- file.path(NORO, "Variants_Pop_Level/tempest_data/all",
                       "G24_allseqs_capsid_sstrainfull.selection.fasta")

stopifnot(file.exists(P33_FASTA), file.exists(POP_FASTA))

say <- function(...) cat(sprintf(...), "\n")

# --------------------------------------------------------------------------
# Helper: rectangular tree layout -> plain data frames
# Reproduces what ggtree draws, using only base R, so the app can plot it
# with ggplot2::geom_segment().
# --------------------------------------------------------------------------
tree_layout <- function(phy, mrsd = NULL) {
  n_tip <- length(phy$tip.label)
  tot   <- n_tip + phy$Nnode
  edge  <- phy$edge
  elen  <- if (is.null(phy$edge.length)) rep(1, nrow(edge)) else phy$edge.length
  root  <- n_tip + 1L

  kids <- split(edge[, 2], edge[, 1])          # parent -> children
  elen_of <- setNames(elen, edge[, 2])

  # x = distance from root (preorder)
  x <- numeric(tot)
  stack <- root
  while (length(stack)) {
    nd <- stack[1]; stack <- stack[-1]
    ch <- kids[[as.character(nd)]]
    if (!is.null(ch)) {
      x[ch] <- x[nd] + elen_of[as.character(ch)]
      stack <- c(ch, stack)
    }
  }

  # y = tip ladder order (depth-first), internal nodes = mean of children
  y <- numeric(tot)
  tip_counter <- 0
  assign_y <- function(nd) {
    ch <- kids[[as.character(nd)]]
    if (is.null(ch)) {
      tip_counter <<- tip_counter + 1
      y[nd] <<- tip_counter
      return(tip_counter)
    }
    ys <- vapply(ch, assign_y, numeric(1))
    y[nd] <<- mean(range(ys))
    y[nd]
  }
  assign_y(root)

  # convert x (substitutions/site OR years) to calendar time if asked
  time_x <- x
  if (!is.null(mrsd)) {
    max_x <- max(x)
    time_x <- x - max_x + decimal_date(as.Date(mrsd))
  }

  nodes <- data.frame(
    node    = seq_len(tot),
    x       = time_x,
    y       = y,
    is_tip  = seq_len(tot) <= n_tip,
    label   = c(phy$tip.label, rep(NA_character_, phy$Nnode)),
    stringsAsFactors = FALSE
  )

  # horizontal branch for every edge + vertical connector per parent
  hseg <- data.frame(
    x    = time_x[edge[, 1]],
    xend = time_x[edge[, 2]],
    y    = y[edge[, 2]],
    yend = y[edge[, 2]]
  )
  vseg <- do.call(rbind, lapply(names(kids), function(p) {
    p  <- as.integer(p)
    ys <- y[kids[[as.character(p)]]]
    data.frame(x = time_x[p], xend = time_x[p], y = min(ys), yend = max(ys))
  }))

  list(nodes = nodes, segs = rbind(hseg, vseg))
}

# treedater stores the strict-clock rate in $mean.rate (NOT $clockRate, which
# does not exist -- referencing it silently yields NULL).
td_rate <- function(td) {
  for (f in c("mean.rate", "adjusted.mean.rate", "meanRate", "clockRate")) {
    v <- td[[f]]
    if (!is.null(v) && length(v) == 1 && is.finite(v)) return(as.numeric(v))
  }
  stop("could not find clock rate in treedater object")
}

# Root-to-tip distances from a phylo
root_to_tip <- function(phy) {
  d <- node.depth.edgelength(phy)
  setNames(d[seq_along(phy$tip.label)], phy$tip.label)
}

# --------------------------------------------------------------------------
# 1. P33 within-host alignment
# --------------------------------------------------------------------------
say("Reading P33 alignment ...")
p33 <- read.dna(P33_FASTA, format = "fasta")

# order sequences by sampling date (earliest first)
p33_dates_chr <- vapply(strsplit(rownames(p33), "\\|"), `[`, character(1), 2)
p33 <- p33[order(as.Date(p33_dates_chr)), ]
p33_names     <- rownames(p33)
p33_dates_chr <- vapply(strsplit(p33_names, "\\|"), `[`, character(1), 2)
p33_acc       <- vapply(strsplit(p33_names, "\\|"), `[`, character(1), 1)
p33_dec       <- decimal_date(as.Date(p33_dates_chr))

# day-since-first label, which is what students actually reason about
day0     <- min(as.Date(p33_dates_chr))
p33_day  <- as.integer(as.Date(p33_dates_chr) - day0) + 1L
p33_lab  <- sprintf("Day %d (%s)", p33_day, p33_acc)

nt_mat <- toupper(as.character(p33))
rownames(nt_mat) <- p33_lab
n_codon <- ncol(nt_mat) %/% 3

say("  %d sequences x %d nt  (%d codons)", nrow(nt_mat), ncol(nt_mat), n_codon)

# translate every sequence
aa_mat <- matrix("", nrow(nt_mat), n_codon, dimnames = list(p33_lab, NULL))
for (i in seq_len(nrow(nt_mat))) {
  cod <- matrix(tolower(nt_mat[i, 1:(n_codon * 3)]), ncol = 3, byrow = TRUE)
  aa_mat[i, ] <- apply(cod, 1, seqinr::translate)
}

meta <- data.frame(
  name = p33_names, acc = p33_acc, date = p33_dates_chr,
  decimal_date = p33_dec, day = p33_day, label = p33_lab,
  stringsAsFactors = FALSE
)

var_nt <- which(apply(nt_mat, 2, function(z) length(unique(z)) > 1))
var_aa <- which(apply(aa_mat, 2, function(z) length(unique(z)) > 1))
say("  %d variable nucleotide sites, %d variable amino-acid sites",
    length(var_nt), length(var_aa))

file.copy(P33_FASTA, file.path(DATA, "P33_capsid.fasta"), overwrite = TRUE)

saveRDS(list(nt = nt_mat, aa = aa_mat, meta = meta,
             var_nt = var_nt, var_aa = var_aa),
        file.path(OUT, "p33_alignment.rds"))

# --------------------------------------------------------------------------
# 2. P33 dated tree + root-to-tip
# --------------------------------------------------------------------------
say("Building & dating P33 tree ...")
p33_nj <- nj(dist.dna(p33, model = "K80"))
p33_td <- dater(p33_nj, setNames(p33_dec, p33_names),
                s = ncol(p33), clock = "strict")

mrsd_p33 <- as.character(max(as.Date(p33_dates_chr)))
lay <- tree_layout(p33_td, mrsd = mrsd_p33)
# pretty tip labels
lay$nodes$label <- meta$label[match(lay$nodes$label, meta$name)]

p33_rtt_raw <- root_to_tip(p33_nj)
rtt_p33 <- data.frame(
  label = meta$label[match(names(p33_rtt_raw), meta$name)],
  day   = meta$day[match(names(p33_rtt_raw), meta$name)],
  decimal_date = meta$decimal_date[match(names(p33_rtt_raw), meta$name)],
  distance = as.numeric(p33_rtt_raw),
  stringsAsFactors = FALSE
)
fit_p33 <- lm(distance ~ decimal_date, rtt_p33)

p33_rate <- td_rate(p33_td)
say("  within-host clock rate = %.3e subs/site/year", p33_rate)
say("  root-to-tip slope      = %.3e   (R2 = %.3f)",
    coef(fit_p33)[2], summary(fit_p33)$r.squared)

saveRDS(list(layout = lay, rtt = rtt_p33,
             slope = unname(coef(fit_p33)[2]),
             intercept = unname(coef(fit_p33)[1]),
             r2 = summary(fit_p33)$r.squared,
             rate = p33_rate, mrsd = mrsd_p33,
             tmrca = min(lay$nodes$x)),
        file.path(OUT, "p33_tree.rds"))

# --------------------------------------------------------------------------
# 3. Pairwise dN/dS for every pair of the 9 sequences
# --------------------------------------------------------------------------
say("Computing pairwise dN/dS (%d pairs) ...", choose(nrow(nt_mat), 2))
p33_aln <- read.alignment(P33_FASTA, format = "fasta")
aln_idx <- match(tolower(p33_names), tolower(p33_aln$nam))

codon_changes <- list()
pair_rows <- list()

for (a in 1:(nrow(nt_mat) - 1)) for (b in (a + 1):nrow(nt_mat)) {
  s1 <- tolower(nt_mat[a, ]); s2 <- tolower(nt_mat[b, ])
  syn <- 0; non <- 0; rows <- list()
  for (k in seq_len(n_codon)) {
    pos <- ((k - 1) * 3 + 1):(k * 3)
    c1 <- s1[pos]; c2 <- s2[pos]
    if (!identical(c1, c2)) {
      a1 <- seqinr::translate(c1); a2 <- seqinr::translate(c2)
      if (a1 == a2) syn <- syn + 1 else non <- non + 1
      rows[[length(rows) + 1]] <- data.frame(
        codon = k, nt_pos = pos[1],
        codon_from = toupper(paste(c1, collapse = "")),
        codon_to   = toupper(paste(c2, collapse = "")),
        aa_from = a1, aa_to = a2,
        type = if (a1 == a2) "Synonymous" else "Non-synonymous",
        change = if (a1 == a2) sprintf("%s%d (silent)", a1, k)
                 else sprintf("%s%d%s", a1, k, a2),
        stringsAsFactors = FALSE)
    }
  }

  sub <- p33_aln
  sub$nb <- 2
  sub$nam <- p33_aln$nam[aln_idx[c(a, b)]]
  sub$seq <- p33_aln$seq[aln_idx[c(a, b)]]
  sub$com <- p33_aln$com[aln_idx[c(a, b)]]
  kk <- seqinr::kaks(sub)
  dN <- as.numeric(kk$ka); dS <- as.numeric(kk$ks)

  key <- paste(a, b, sep = "-")
  codon_changes[[key]] <- if (length(rows)) do.call(rbind, rows) else
    data.frame(codon = integer(), nt_pos = integer(), codon_from = character(),
               codon_to = character(), aa_from = character(), aa_to = character(),
               type = character(), change = character(), stringsAsFactors = FALSE)
  pair_rows[[key]] <- data.frame(
    i = a, j = b, key = key,
    label_i = p33_lab[a], label_j = p33_lab[b],
    days_apart = abs(p33_day[b] - p33_day[a]),
    syn = syn, nonsyn = non, total = syn + non,
    dN = dN, dS = dS, omega = dN / dS,
    stringsAsFactors = FALSE)
}

pairs_df <- do.call(rbind, pair_rows)
first_last <- pairs_df[pairs_df$i == 1 & pairs_df$j == nrow(nt_mat), ]
say("  day 1 vs day %d:  S=%d  N=%d  dN/dS=%.3f",
    max(p33_day), first_last$syn, first_last$nonsyn, first_last$omega)

saveRDS(list(pairs = pairs_df, changes = codon_changes),
        file.path(OUT, "p33_dnds.rds"))

# --------------------------------------------------------------------------
# 4. Sydney 2012 population-level subsample (n = 350)
# --------------------------------------------------------------------------
say("Reading population alignment ...")
pop <- read.dna(POP_FASTA, format = "fasta")
pop_var <- vapply(strsplit(rownames(pop), "\\|"), `[`, character(1), 3)
syd_all <- pop[pop_var == "Sydney_2012", ]
say("  %d Sydney 2012 sequences available", nrow(syd_all))

N_SUB <- 350
syd <- syd_all[sort(sample(nrow(syd_all), N_SUB)), ]
syd_names <- rownames(syd)
syd_dchr  <- vapply(strsplit(syd_names, "\\|"), `[`, character(1), 2)
syd_dec   <- decimal_date(as.Date(syd_dchr))

# save the exact subsample so the analysis is reproducible / inspectable
write.dna(syd, file.path(DATA, "Sydney_2012_capsid_n350.fasta"),
          format = "fasta", nbcol = -1, colsep = "")
say("  wrote subsample FASTA (%s to %s)", min(syd_dchr), max(syd_dchr))

say("Building & dating Sydney tree (this is the slow step) ...")
syd_nj <- nj(dist.dna(syd, model = "K80"))
syd_td <- dater(syd_nj, setNames(syd_dec, syd_names),
                s = ncol(syd), clock = "strict")

mrsd_syd <- as.character(max(as.Date(syd_dchr)))
syd_lay  <- tree_layout(syd_td, mrsd = mrsd_syd)
syd_lay$nodes$label <- NA_character_          # 350 tips: no labels

syd_rtt_raw <- root_to_tip(syd_nj)
rtt_syd <- data.frame(
  acc = vapply(strsplit(names(syd_rtt_raw), "\\|"), `[`, character(1), 1),
  decimal_date = decimal_date(as.Date(
    vapply(strsplit(names(syd_rtt_raw), "\\|"), `[`, character(1), 2))),
  distance = as.numeric(syd_rtt_raw),
  stringsAsFactors = FALSE
)
fit_syd  <- lm(distance ~ decimal_date, rtt_syd)
syd_rate <- td_rate(syd_td)
say("  population clock rate = %.3e subs/site/year", syd_rate)

saveRDS(list(layout = syd_lay, rtt = rtt_syd,
             slope = unname(coef(fit_syd)[2]),
             intercept = unname(coef(fit_syd)[1]),
             r2 = summary(fit_syd)$r.squared,
             rate = syd_rate, mrsd = mrsd_syd, n = N_SUB,
             n_available = nrow(syd_all),
             date_min = min(syd_dchr), date_max = max(syd_dchr),
             tmrca = min(syd_lay$nodes$x)),
        file.path(OUT, "sydney_tree.rds"))

# --------------------------------------------------------------------------
# 5. Headline numbers for the app
# --------------------------------------------------------------------------
saveRDS(list(
  p33_rate = p33_rate, syd_rate = syd_rate,
  fold = p33_rate / syd_rate,
  p33_n = nrow(nt_mat), p33_nt = ncol(nt_mat), p33_codons = n_codon,
  p33_span_days = max(p33_day),
  n_var_nt = length(var_nt), n_var_aa = length(var_aa),
  fl_syn = first_last$syn, fl_non = first_last$nonsyn,
  fl_omega = first_last$omega, fl_dN = first_last$dN, fl_dS = first_last$dS,
  built = as.character(Sys.Date())
), file.path(OUT, "summary.rds"))

say("")
say("DONE. Files in %s:", OUT)
print(data.frame(file = basename(list.files(OUT, full.names = TRUE)),
                 kb = round(file.size(list.files(OUT, full.names = TRUE)) / 1024, 1)))
