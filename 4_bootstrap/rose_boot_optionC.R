###############################################################################
# rose_boot_optionC.R   (parallel, timed, estimability-screened)
#
# Bootstrap CI for the ROSE secondary effect (AG - CF) using Option C:
#   - Point estimates, including the frailty variance theta, come from the
#     STATA strmcure fit on the full data (run separately). theta is held
#     FIXED at that value in every resample.
#   - In each resample, AG and CF are refit with coxph + tt() (exact-time
#     evaluation), theta fixed. Persons are the resampling unit.
#
# Reports the AUC-scale secondary effect with:
#   a single 95% CI (no screen needed; see pooling note below).
#
# EVENT-ORDER POOLING (fixed two strata):
#   Episodes are pooled into two strata: first episodes and all recurrences.
#   This matches the full-data fit (recurrences already form one stratum).
#   Because the recurrence stratum is always retained, the conditional frailty
#   model never collapses onto Andersen-Gill in sparse resamples, so every
#   resample is estimable and no screening is applied.
#

#
# Parallel (PSOCK); crash-safe checkpointing; resumable; wall-clock timing
# printed to console AND logged.
#
# Usage:  ROSE_DATA_DIR=/path/to/data Rscript rose_boot_optionC.R
# Requires: ROSE_ph.csv (from 1_data_preparation/prepare_data.R),
#           theta_fixed from the STATA strmcure run (ch2roseCF_tt.do).
#           packages survival, parallel.
###############################################################################

suppressPackageStartupMessages({ library(survival); library(parallel) })

# ---------------- Settings ----------------
data_file   <- file.path(Sys.getenv("ROSE_DATA_DIR", "../data"), "ROSE_ph.csv")
out_dir     <- "./rose_boot_C_out"
B           <- 1000
n_cores     <- 8
Tau         <- 20
theta_fixed <- 0.8236            # <-- STATA strmcure theta (ROSE_split_tt.csv; ch2roseCF_tt.do)
base_seed   <- 20260101
checkpoint_every <- 100
log_file    <- file.path(out_dir, "rose_boot_runlog.txt")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
progress_file <- file.path(out_dir, "boot_C_progress.rds")

log_msg <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...))
  cat(line, "\n"); cat(line, "\n", file = log_file, append = TRUE)
}
cat(sprintf("=== rose_boot_optionC.R run log - started %s ===\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = log_file, append = FALSE)
t_start <- Sys.time()

# ---------------- Data ----------------
df0 <- read.csv(data_file); names(df0) <- tolower(names(df0))
if ("x" %in% names(df0)) df0$x <- NULL
stopifnot(all(c("id","group","start","end","event","event_count") %in% names(df0)))
uid0 <- unique(df0$id)
log_msg("Loaded %d rows, %d persons, %d events. theta_fixed=%.4f, B=%d, cores=%d",
        nrow(df0), length(uid0), sum(df0$event), theta_fixed, B, n_cores)

tt_fun  <- function(x, t, ...) x * log(t)
ve_fun  <- function(t, b0, b1) 1 - exp(b0 + b1 * log(t))
auc_fun <- function(b0, b1) 1 - exp(b0) * Tau^b1 / (b1 + 1)   # normalized AUC
make_pool <- function(dat) {
  factor(
    pmin(dat$event_count, 1),
    levels = c(0, 1),
    labels = c("First event", "Recurrent events")
  )
}

# ---------------- Point estimates on original data ----------------
df0$ec_pool <- make_pool(df0)
log_msg("Failure counts by pooled event stratum:")
print(table(df0$ec_pool[df0$event == 1]))
ag0 <- coxph(Surv(start,end,event) ~ group + tt(group) + cluster(id),
             data=df0, ties="breslow", tt=tt_fun, control=coxph.control(timefix=FALSE))
cf0 <- coxph(Surv(start,end,event) ~ group + tt(group) + strata(ec_pool) +
               frailty(id, theta=theta_fixed),
             data=df0, ties="breslow", tt=tt_fun, control=coxph.control(timefix=FALSE))
ag_b0<-coef(ag0)[1]; ag_b1<-coef(ag0)["tt(group)"]
cf_b0<-coef(cf0)["group"]; cf_b1<-coef(cf0)["tt(group)"]
sec_auc_pt <- auc_fun(ag_b0, ag_b1) - auc_fun(cf_b0, cf_b1)
log_msg("Point: AG b1=%.4f, CF b1=%.4f (theta %.3f). Secondary AUC=%.4f",
        ag_b1, cf_b1, theta_fixed, sec_auc_pt)

# ---------------- Worker ----------------
boot_one <- function(b, df0, Tau, theta_fixed, base_seed) {
  library(survival)
  set.seed(base_seed + b)
  tt_fun  <- function(x, t, ...) x * log(t)
  auc_fun <- function(b0, b1) 1 - exp(b0) * Tau^b1 / (b1 + 1)
  make_pool <- function(dat) {
    factor(
      pmin(dat$event_count, 1),
      levels = c(0, 1),
      labels = c("First event", "Recurrent events")
    )
  }
  uid <- unique(df0$id)
  rows_by_id <- split(seq_len(nrow(df0)), df0$id)
  samp <- sample(uid, length(uid), replace = TRUE)
  idx  <- unlist(rows_by_id[as.character(samp)])
  d <- df0[idx, , drop = FALSE]
  d$id <- rep(seq_along(samp), lengths(rows_by_id[as.character(samp)]))
  d <- d[order(d$id, d$start), ]
  d$event_count <- ave(d$event, d$id, FUN = function(e) cumsum(e) - e)
  d$ec_pool <- make_pool(d)
  n_second <- sum(d$event == 1 & d$event_count >= 1)
  # Number of strata that actually contain failures
  n_strata <- length(unique(d$ec_pool[d$event == 1]))
  
  ag <- tryCatch(coxph(Surv(start,end,event) ~ group + tt(group) + cluster(id),
                       data=d, ties="breslow", tt=tt_fun,
                       control=coxph.control(timefix=FALSE)), error=function(e) NULL)
  cf <- tryCatch(coxph(Surv(start,end,event) ~ group + tt(group) + strata(ec_pool) +
                         frailty(id, theta=theta_fixed),
                       data=d, ties="breslow", tt=tt_fun,
                       control=coxph.control(timefix=FALSE)), error=function(e) NULL)
  if (is.null(ag) || is.null(cf)) {
    return(data.frame(rep=b, ag_b1=NA, cf_b1=NA, coef_gap=NA, sec_auc=NA,
                      n_second=n_second, n_strata=n_strata, failed=TRUE))
  }
  a0<-coef(ag)[1]; a1<-coef(ag)["tt(group)"]; p0<-coef(cf)["group"]; p1<-coef(cf)["tt(group)"]
  data.frame(rep=b, ag_b1=a1, cf_b1=p1, coef_gap=a1-p1,
             sec_auc=auc_fun(a0,a1)-auc_fun(p0,p1),
             n_second=n_second, n_strata=n_strata, failed=FALSE)
}

# ---------------- Resume + run in checkpointed chunks ----------------
done <- integer(0); acc <- NULL
if (file.exists(progress_file)) { acc <- readRDS(progress_file); done <- acc$rep
log_msg("Resuming: %d draws already done.", length(done)) }
todo <- setdiff(1:B, done)

if (length(todo) > 0) {
  cl <- makeCluster(n_cores)
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
  chunks <- split(todo, ceiling(seq_along(todo)/checkpoint_every))
  for (ci in seq_along(chunks)) {
    t0 <- Sys.time()
    res <- parLapply(cl, chunks[[ci]], boot_one,
                     df0=df0, Tau=Tau, theta_fixed=theta_fixed, base_seed=base_seed)
    acc <- rbind(acc, do.call(rbind, res)); saveRDS(acc, progress_file)
    log_msg("Chunk %d/%d: %d draws in %.1f s. Total %d/%d.",
            ci, length(chunks), length(chunks[[ci]]),
            as.numeric(difftime(Sys.time(), t0, units="secs")), nrow(acc), B)
  }
  stopCluster(cl)
}

# ---------------- Summarize ----------------
# Fixed two-strata pooling: recurrences are always kept as a stratum, so the
# conditional frailty model never collapses onto Andersen-Gill and every
# resample is estimable. No screen is applied; a single CI is reported.
res  <- readRDS(progress_file)
n_fail <- sum(res$failed)

pc <- function(x) quantile(x[is.finite(x)], c(.025, .975))
auc_ci <- pc(res$sec_auc[!res$failed])

log_msg("=============================================================")
log_msg("ROSE secondary effect - Option C bootstrap (theta fixed=%.3f)", theta_fixed)
log_msg("Replicates: %d total, %d failed. (Fixed two-strata pooling; all estimable.)",
        B, n_fail)
log_msg("Recurrences per resample: median %.0f.", median(res$n_second[!res$failed]))
log_msg("Secondary effect (AUC): point=%.4f  95%% CI (%.4f, %.4f)",
        sec_auc_pt, auc_ci[1], auc_ci[2])

write.csv(res, file.path(out_dir, "boot_C_draws.csv"), row.names = FALSE)
saveRDS(list(theta_fixed=theta_fixed, sec_auc_pt=sec_auc_pt,
             auc_ci=auc_ci, n_total=B, n_failed=n_fail),
        file.path(out_dir, "boot_C_summary.rds"))
log_msg("Wrote boot_C_draws.csv and boot_C_summary.rds.")
log_msg("TOTAL wall-clock: %.1f min.", as.numeric(difftime(Sys.time(), t_start, units="mins")))