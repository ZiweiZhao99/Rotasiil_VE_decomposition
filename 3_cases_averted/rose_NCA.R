###############################################################################
# rose_NCA.R
#
# Number of cases averted (NCA) per 1000 children for the Rotasiil trial,
# computed separately for the COMPOSITE and PRIMARY time-varying effects,
# and their difference (the secondary effect on a case-count scale).
#
# Coefficients are NOT refitted here. They are the Stata estimates from the
# main analysis (stcox for the Andersen-Gill model, strmcure for the
# conditional frailty model) and are entered below as constants. Only the
# control-arm event counts and risk time are read from the data.
#
# Method (following the RTS,S implementation in Vivli_01_MasterFunctions.R):
#   For each monthly interval [lo, hi):
#     incidence_c = control events in interval / control risk time in interval
#     ve_int      = INTEGRAL of VE(t) over [lo, hi)      <- raw, NOT normalized
#     cases_averted = incidence_c * ve_int * 1000
#   NCA = sum over intervals.
#
#   Control incidence is computed per interval because rotavirus incidence
#   falls steeply with age; a single pooled rate would misallocate averted
#   cases across time. Risk time accumulates partial-interval overlap.
#
# IMPORTANT: ve_int is the RAW integral of VE(t) over the interval, in units
# of months. This differs from the normalized AUC reported in the manuscript
# (integral over 0-20 divided by 20, an average VE). The raw integral is
# required here because it multiplies a rate to give a count.
#
# INTERPRETATION: the composite NCA is the expected real-world reduction in
# episodes. The primary NCA is a counterfactual in which the natural-immunity
# pathway is blocked; it is a decomposition device, not a projected impact.
# Their difference expresses the secondary effect as cases per 1000 children.
#
# Coefficients are from the strmcure/stcox fit on ROSE_split_tt.csv
# (2_main_analysis_stata/ch2roseCF_tt.do).
#
# Usage:  ROSE_DATA_DIR=/path/to/data Rscript rose_NCA.R
# Requires: ROSE_ph.csv (from 1_data_preparation/prepare_data.R); package dplyr.
###############################################################################

suppressPackageStartupMessages({ library(dplyr) })

# ---------------- Settings ----------------
data_file <- file.path(Sys.getenv("ROSE_DATA_DIR", "../data"), "ROSE_ph.csv")
Tau       <- 20
breaks    <- seq(0, Tau, by = 1)
per_n     <- 1000            # report per 1000 children

# ---------------- Coefficients from Stata (main analysis) ----------------
# Composite effect: stcox with tt(group), g_c(t) = b0c + b1c * log(t)
ag_b0 <- -0.5488418
ag_b1 <-  0.1858125

# Primary effect: strmcure, g_p(t) = b0p + b1p * log(t)
cf_b0 <- -0.571524
cf_b1 <-  0.1776118

cat("Coefficients (from Stata, not refitted):\n")
cat(sprintf("  Composite (AG, stcox)    : b0 = %+.4f, b1 = %+.4f\n", ag_b0, ag_b1))
cat(sprintf("  Primary   (CF, strmcure) : b0 = %+.4f, b1 = %+.4f\n\n", cf_b0, cf_b1))

# VE(t) as integrable functions
ve_composite <- function(t) 1 - exp(ag_b0 + ag_b1 * log(t))
ve_primary   <- function(t) 1 - exp(cf_b0 + cf_b1 * log(t))

# ---------------- Data (control arm only is used) ----------------
d <- read.csv(data_file); names(d) <- tolower(names(d))
d <- d[, !grepl("^unnamed|^x$", names(d)), drop = FALSE]
stopifnot(all(c("id","group","start","end","event") %in% names(d)))
cat(sprintf("Loaded %d rows, %d persons, %d events (%d in control arm).\n\n",
            nrow(d), length(unique(d$id)), sum(d$event),
            sum(d$event[d$group == 0])))

# ---------------- Risk-time helper (partial-interval overlap) ----------------
risktime <- function(start, end, lo, hi) {
  dplyr::case_when(
    end   <= lo                         ~ 0,
    start <= lo & end <= hi             ~ end - lo,
    start <= lo & end >  hi             ~ hi - lo,
    start >  lo & end <= hi             ~ end - start,
    start >  lo & start < hi & end > hi ~ hi - start,
    TRUE                                ~ 0
  )
}

# ---------------- NCA ----------------
calculate_NCA <- function(data, ve_func, breaks, per_n = 1000) {
  control  <- dplyr::filter(data, group == 0)
  n_int    <- length(breaks) - 1
  events_c <- integer(n_int); risktime_c <- numeric(n_int)

  for (i in seq_len(n_int)) {
    lo <- breaks[i]; hi <- breaks[i + 1]
    events_c[i]   <- sum(control$event == 1 & control$end > lo & control$end <= hi,
                         na.rm = TRUE)
    risktime_c[i] <- sum(risktime(control$start, control$end, lo, hi), na.rm = TRUE)
  }

  # RAW integral of VE(t) over each interval; guard the log singularity at t = 0
  ve_int <- sapply(seq_len(n_int), function(i) {
    lo <- breaks[i]; hi <- breaks[i + 1]
    lo_adj <- if (lo == 0) 1e-6 else lo
    tryCatch(integrate(ve_func, lower = lo_adj, upper = hi)$value,
             error = function(e) NA_real_)
  })

  incidence_c   <- events_c / risktime_c
  cases_averted <- incidence_c * ve_int * per_n
  cases_averted[is.nan(cases_averted) | is.infinite(cases_averted)] <- 0
  cases_averted[risktime_c == 0] <- 0

  data.frame(month = seq_len(n_int), events_c = events_c,
             risktime_c = risktime_c, incidence_c = incidence_c,
             ve_int = ve_int, cases_averted = cases_averted)
}

nca_c <- calculate_NCA(d, ve_composite, breaks, per_n)
nca_p <- calculate_NCA(d, ve_primary,   breaks, per_n)

NCA_composite <- sum(nca_c$cases_averted, na.rm = TRUE)
NCA_primary   <- sum(nca_p$cases_averted, na.rm = TRUE)
NCA_secondary <- NCA_composite - NCA_primary

# ---------------- Report ----------------
out <- data.frame(
  month         = nca_c$month,
  events_c      = nca_c$events_c,
  risktime_c    = round(nca_c$risktime_c, 1),
  incidence_c   = round(nca_c$incidence_c, 5),
  ve_int_comp   = round(nca_c$ve_int, 4),
  ve_int_prim   = round(nca_p$ve_int, 4),
  nca_composite = round(nca_c$cases_averted, 2),
  nca_primary   = round(nca_p$cases_averted, 2)
)
out$nca_diff <- round(out$nca_composite - out$nca_primary, 2)

cat("Monthly detail (control arm; NCA per", per_n, "children):\n")
print(out, row.names = FALSE)

cat(sprintf("\n================ NCA over 0-%d months, per %d children ================\n",
            Tau, per_n))
cat(sprintf("  Composite effect : %.1f cases averted\n", NCA_composite))
cat(sprintf("  Primary effect   : %.1f cases averted (counterfactual)\n", NCA_primary))
cat(sprintf("  Difference       : %.1f  (secondary effect, cases per %d)\n",
            NCA_secondary, per_n))
cat(sprintf("\n  Control episodes observed: %d over %.0f child-months\n",
            sum(nca_c$events_c), sum(nca_c$risktime_c)))

write.csv(out, "rose_NCA_monthly.csv", row.names = FALSE)
cat("\nWrote rose_NCA_monthly.csv\n")
