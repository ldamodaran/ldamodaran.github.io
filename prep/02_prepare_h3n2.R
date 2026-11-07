# ==========================================================================
# 02_prepare_h3n2.R  --  BIO463 molecular epidemiology of H3N2
#
# Precomputes everything in Molecular-Epi-H3N2.Rmd so the app can run in the
# browser via shinylive. Same approach as 01_prepare_data.R: the heavy and
# package-dependent work happens here, the app only reads data frames.
#
# Usage:  Rscript prep/02_prepare_h3n2.R
# ==========================================================================

suppressPackageStartupMessages({
  library(ape)
  library(lubridate)
  library(treedater)
})

OUT  <- "app-h3n2/data"
DATA <- "site/data"
dir.create(OUT,  recursive = TRUE, showWarnings = FALSE)
dir.create(DATA, recursive = TRUE, showWarnings = FALSE)

SRC <- path.expand(file.path(
  "~/Library/CloudStorage/Dropbox/scientific_reasources/teaching",
  "BIO463/BIO463-phylo-tutorial"))

FASTA  <- file.path(SRC, "Northern_hemisphere_sequences_2010_2026.fasta")
DATES  <- file.path(SRC, "north_date_df.csv")
CASES  <- file.path(SRC, "WHO-H3N2-global-cases.csv")
stopifnot(file.exists(FASTA), file.exists(DATES), file.exists(CASES))

say <- function(...) cat(sprintf(...), "\n")

# treedater keeps the strict-clock rate in $mean.rate ($clockRate does not
# exist and silently returns NULL).
td_rate <- function(td) {
  for (f in c("mean.rate", "adjusted.mean.rate", "meanRate", "clockRate")) {
    v <- td[[f]]
    if (!is.null(v) && length(v) == 1 && is.finite(v)) return(as.numeric(v))
  }
  stop("could not find clock rate in treedater object")
}

# Rectangular tree layout -> segments (same helper as script 01)
tree_layout <- function(phy, mrsd = NULL) {
  n_tip <- length(phy$tip.label)
  tot   <- n_tip + phy$Nnode
  edge  <- phy$edge
  elen  <- if (is.null(phy$edge.length)) rep(1, nrow(edge)) else phy$edge.length
  root  <- n_tip + 1L
  kids  <- split(edge[, 2], edge[, 1])
  elen_of <- setNames(elen, edge[, 2])

  x <- numeric(tot); stack <- root
  while (length(stack)) {
    nd <- stack[1]; stack <- stack[-1]
    ch <- kids[[as.character(nd)]]
    if (!is.null(ch)) {
      x[ch] <- x[nd] + elen_of[as.character(ch)]
      stack <- c(ch, stack)
    }
  }

  y <- numeric(tot); tip_counter <- 0
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

  time_x <- x
  if (!is.null(mrsd)) time_x <- x - max(x) + decimal_date(as.Date(mrsd))

  nodes <- data.frame(node = seq_len(tot), x = time_x, y = y,
                      is_tip = seq_len(tot) <= n_tip,
                      label = c(phy$tip.label, rep(NA_character_, phy$Nnode)),
                      stringsAsFactors = FALSE)
  hseg <- data.frame(x = time_x[edge[, 1]], xend = time_x[edge[, 2]],
                     y = y[edge[, 2]], yend = y[edge[, 2]])
  vseg <- do.call(rbind, lapply(names(kids), function(p) {
    p <- as.integer(p); ys <- y[kids[[as.character(p)]]]
    data.frame(x = time_x[p], xend = time_x[p], y = min(ys), yend = max(ys))
  }))
  list(nodes = nodes, segs = rbind(hseg, vseg))
}

# --------------------------------------------------------------------------
# 1. WHO sentinel case data
# --------------------------------------------------------------------------
say("Reading WHO case data ...")
cases <- read.csv(CASES)
cases$Week <- as.Date(cases$Week)
cases <- cases[order(cases$Week), ]

cases_out <- data.frame(
  week         = cases$Week,
  decimal_date = decimal_date(cases$Week),
  hemisphere   = cases$Hemisphere,
  tested       = cases$Specimen.tested,
  positive     = cases$Influenza.positive,
  negative     = cases$Influenza.negative,
  h3           = cases$A.H3,
  stringsAsFactors = FALSE
)
cases_out$pct_positive <- with(cases_out,
  ifelse(tested > 0, 100 * positive / tested, NA_real_))

say("  %d weekly records, %s to %s, hemispheres: %s",
    nrow(cases_out), min(cases_out$week), max(cases_out$week),
    paste(unique(cases_out$hemisphere), collapse = " / "))
saveRDS(cases_out, file.path(OUT, "h3n2_cases.rds"))

# --------------------------------------------------------------------------
# 2. Sequences, tree, dating
# --------------------------------------------------------------------------
say("Reading HA alignment ...")
aln <- read.dna(FASTA, format = "fasta")
dd  <- read.csv(DATES)

stopifnot(setequal(rownames(aln), dd$taxa))
dd <- dd[match(rownames(aln), dd$taxa), ]
dd$decimal_date <- decimal_date(as.Date(dd$date))
say("  %d sequences x %d nt (HA), %s to %s",
    nrow(aln), ncol(aln), min(dd$date), max(dd$date))

# pull a readable strain name out of the long FASTA header
#   "PX982508.1 |Influenza A virus (A/Rhode Island/.../2026(H3N2)) ...|date"
strain_of <- function(h) {
  m <- regmatches(h, regexpr("A/[^()]+?\\(H3N2\\)", h))
  s <- ifelse(lengths(m) > 0 | nzchar(m), m, NA_character_)
  s <- sub("\\(H3N2\\)$", "", s)
  acc <- sub("^([A-Z0-9.]+).*", "\\1", h)
  ifelse(is.na(s) | !nzchar(s), acc, trimws(s))
}
dd$strain <- strain_of(dd$taxa)
dd$acc    <- sub("^([A-Z0-9.]+).*", "\\1", dd$taxa)
dd$year   <- as.integer(substr(dd$date, 1, 4))

say("Building & dating the tree ...")
tre <- nj(dist.dna(aln, model = "K80"))
td  <- dater(tre, setNames(dd$decimal_date, dd$taxa),
             s = ncol(aln), clock = "strict")

mrsd <- as.character(max(as.Date(dd$date)))
lay  <- tree_layout(td, mrsd = mrsd)
# attach season/year to tips so the app can color by sampling year
lay$nodes$year <- NA_integer_
ti <- match(lay$nodes$label, dd$taxa)
lay$nodes$year[!is.na(ti)]   <- dd$year[ti[!is.na(ti)]]
lay$nodes$strain <- NA_character_
lay$nodes$strain[!is.na(ti)] <- dd$strain[ti[!is.na(ti)]]
lay$nodes$label  <- NULL          # drop the long headers

rate <- td_rate(td)
say("  clock rate = %.3e subs/site/year", rate)

# --------------------------------------------------------------------------
# 3. Root-to-tip regression
# --------------------------------------------------------------------------
rtt_raw <- node.depth.edgelength(tre)[seq_along(tre$tip.label)]
names(rtt_raw) <- tre$tip.label
k <- match(names(rtt_raw), dd$taxa)
rtt <- data.frame(
  acc = dd$acc[k], strain = dd$strain[k], year = dd$year[k],
  decimal_date = dd$decimal_date[k],
  distance = as.numeric(rtt_raw),
  stringsAsFactors = FALSE
)
fit <- lm(distance ~ decimal_date, rtt)
say("  root-to-tip slope = %.3e  (R2 = %.3f)",
    coef(fit)[2], summary(fit)$r.squared)

# --------------------------------------------------------------------------
# 4. Lineages through time
# --------------------------------------------------------------------------
say("Computing lineages through time ...")
co <- ltt.plot.coords(td, backward = TRUE)
ltt <- as.data.frame(co)
# co$time is negative, measured back from the most recent tip
ltt$decimal_date <- ltt$time + decimal_date(as.Date(mrsd))
ltt <- ltt[, c("decimal_date", "N")]
ltt$week <- as.Date(date_decimal(ltt$decimal_date))
say("  %d steps, %s to %s, up to %d lineages",
    nrow(ltt), min(ltt$week), max(ltt$week), max(ltt$N))

saveRDS(list(layout = lay, rtt = rtt, ltt = ltt,
             slope = unname(coef(fit)[2]),
             intercept = unname(coef(fit)[1]),
             r2 = summary(fit)$r.squared,
             rate = rate, mrsd = mrsd,
             tmrca = min(lay$nodes$x),
             n_tip = nrow(dd), n_nt = ncol(aln),
             date_min = min(dd$date), date_max = max(dd$date),
             years = sort(unique(dd$year))),
        file.path(OUT, "h3n2_tree.rds"))

# --------------------------------------------------------------------------
# 5. Headline numbers
# --------------------------------------------------------------------------
peak <- cases_out[which.max(cases_out$positive), ]
saveRDS(list(
  rate = rate, slope = unname(coef(fit)[2]), r2 = summary(fit)$r.squared,
  n_tip = nrow(dd), n_nt = ncol(aln),
  date_min = min(dd$date), date_max = max(dd$date),
  tmrca = min(lay$nodes$x),
  max_lineages = max(ltt$N),
  n_weeks = nrow(cases_out),
  cases_min = as.character(min(cases_out$week)),
  cases_max = as.character(max(cases_out$week)),
  peak_week = as.character(peak$week),
  peak_n = peak$positive, peak_hemi = peak$hemisphere,
  built = as.character(Sys.Date())
), file.path(OUT, "summary.rds"))

# copy the inputs for students
file.copy(FASTA, file.path(DATA, "H3N2_HA_northern_2010_2026.fasta"),
          overwrite = TRUE)
file.copy(CASES, file.path(DATA, "WHO-H3N2-global-cases.csv"),
          overwrite = TRUE)

say("")
say("DONE. Files in %s:", OUT)
f <- list.files(OUT, full.names = TRUE)
print(data.frame(file = basename(f), kb = round(file.size(f) / 1024, 1)))
