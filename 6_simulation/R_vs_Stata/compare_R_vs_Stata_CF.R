###############################################################################
# compare_R_vs_Stata_CF.R
#
# Same-data comparison of the conditional frailty (CF) estimator in
# R (coxph + frailty) versus Stata (strmcure).
#
# Reads the datasets saved by ../stata/run_save_datasets.do (same DGP and
# seeds as the Stata estimation runs), located in SIM_DATA_DIR:
#   out_tt_full_1/simdata_001.csv ... simdata_100.csv   (theta = 1)
#   out_tt_full_0.2/simdata_001.csv ... simdata_100.csv (theta = 0.2)
# fits FE, AG, CF in R on EACH dataset (R's own tt() splitting), and
# summarizes bias by model and theta. Because the datasets are identical
# to those Stata analyzed, the only thing that differs is the software.
#
# This revision:
#   * per-dataset PROGRESS with elapsed time (so you can see it running)
#   * per-dataset TIMING of each model fit (FE / AG / CF)
#   * CONVERGENCE capture for the R frailty fit — a fit that emits a
#     "failed to converge" warning is flagged, so non-converged datasets
#     can be reported and excluded from the R bias summary.
#
# Datasets are fitted in parallel (forked workers; no random numbers are
# drawn, so results do not depend on the number of cores).
#
# Also reports (Supplementary Material 2 and Results):
#   * share of datasets in which R's frailty variance collapsed to ~0
#   * proportion of children with >=1 episode who had >=2 episodes
#     (simulation_recurrence_proportion.csv)
#
# Usage:  SIM_DATA_DIR=/path/to/simulated/data Rscript compare_R_vs_Stata_CF.R
#
# Output:
#   R_sim_fits.csv            per-dataset R estimates + timing + convergence
#   R_vs_Stata_CF_bias.csv    tidy bias table, both packages
#   simulation_recurrence_proportion.csv
#
# Requires: survival, dplyr, readr, purrr, tidyr.
###############################################################################

suppressPackageStartupMessages({
  library(survival); library(dplyr); library(readr); library(purrr); library(tidyr)
})

# ---------------- Settings ----------------
Tau       <- 20
true_b0   <- -0.57
true_b1   <-  0.17
true_auc  <- 1 - exp(true_b0) * Tau^true_b1 / (true_b1 + 1)   # 0.196
tail_frac <- 0.05

sim_dir <- Sys.getenv("SIM_DATA_DIR", "../simdata")
n_cores <- as.integer(Sys.getenv("N_CORES", "4"))
dirs <- c("1" = file.path(sim_dir, "out_tt_full_1"), "0.2" = file.path(sim_dir, "out_tt_full_0.2"))
# Stata (strmcure) results on the same datasets, from ../stata
stata_csv <- c("1"   = "../stata/results_theta1/results_all_STATA.csv",
               "0.2" = "../stata/results_theta0.2/results_all_theta0.2.csv")
cat(sprintf("R %s; survival %s; %d cores\n", getRversion(), packageVersion("survival"), n_cores))
cf_method <- NULL   # "em" to mimic Stata; set to NULL for R's default (faster, more stable)

# ---------------- Helpers ----------------
tt_fun  <- function(x, t, ...) x * log(t)
auc_fun <- function(b0, b1) 1 - exp(b0) * Tau^b1 / (b1 + 1)
make_pool <- function(dat) {
  ft <- table(dat$event_count[dat$event == 1] + 1)
  if (length(ft) <= 1) return(factor(rep(1, nrow(dat))))
  cumf <- rev(cumsum(rev(ft))) / sum(ft)
  elig <- as.integer(names(cumf[cumf >= tail_frac]))
  J <- if (length(elig) > 0) max(elig) else max(as.integer(names(ft)))
  factor(pmin(dat$event_count + 1, J))
}

fit_one <- function(path) {
  d <- readr::read_csv(path, show_col_types = FALSE); names(d) <- tolower(names(d))
  d <- d %>% arrange(id, start)
  
  out <- tibble(file = basename(path),
                fe_b0=NA_real_, fe_b1=NA_real_, fe_auc=NA_real_,
                ag_b0=NA_real_, ag_b1=NA_real_, ag_auc=NA_real_,
                cf_b0=NA_real_, cf_b1=NA_real_, cf_auc=NA_real_, theta=NA_real_,
                cf_failed=1L, cf_converged=NA, cf_warn="",
                t_fe=NA_real_, t_ag=NA_real_, t_cf=NA_real_)
  
  # --- FE ---
  tf <- system.time({
    d_fe <- d %>% group_by(id) %>% mutate(ce = cumsum(event)) %>%
      filter(ce == 0 | (ce == 1 & event == 1)) %>% slice_tail(n = 1) %>% ungroup()
    fe <- tryCatch(coxph(Surv(start, end, event) ~ z + tt(z), data = d_fe,
                         tt = tt_fun, ties = "breslow",
                         control = coxph.control(timefix = FALSE, iter.max = 500, eps = 1e-6)), error = function(e) NULL)
  })["elapsed"]
  out$t_fe <- as.numeric(tf)
  if (!is.null(fe)) { out$fe_b0 <- coef(fe)[1]; out$fe_b1 <- coef(fe)[2]
  out$fe_auc <- auc_fun(out$fe_b0, out$fe_b1) }
  
  # --- AG ---
  ta <- system.time({
    ag <- tryCatch(coxph(Surv(start, end, event) ~ z + tt(z) + cluster(id), data = d,
                         tt = tt_fun, ties = "breslow",
                         control = coxph.control(timefix = FALSE, iter.max = 500, eps = 1e-6)), error = function(e) NULL)
  })["elapsed"]
  out$t_ag <- as.numeric(ta)
  if (!is.null(ag)) { out$ag_b0 <- coef(ag)[1]; out$ag_b1 <- coef(ag)[2]
  out$ag_auc <- auc_fun(out$ag_b0, out$ag_b1) }
  
  # --- CF with convergence capture ---
  d$ec_pool <- make_pool(d)
  conv <- TRUE; warn_txt <- ""
  frailty_term <- if (is.null(cf_method)) quote(frailty(id)) else quote(frailty(id, method = "em"))
  tc <- system.time({
    cf <- withCallingHandlers(
      tryCatch(
        coxph(reformulate(c("z", "tt(z)", "strata(ec_pool)", deparse(frailty_term)),
                          response = "Surv(start, end, event)"),
              data = d, tt = tt_fun, ties = "breslow",
              control = coxph.control(timefix = FALSE, iter.max = 500, eps = 1e-6)),
        error = function(e) { conv <<- FALSE; warn_txt <<- conditionMessage(e); NULL }),
      warning = function(w) {
        m <- conditionMessage(w)
        if (grepl("converge", m, ignore.case = TRUE)) { conv <<- FALSE; warn_txt <<- m }
        invokeRestart("muffleWarning")
      })
  })["elapsed"]
  out$t_cf <- as.numeric(tc)
  out$cf_converged <- conv
  out$cf_warn <- warn_txt
  if (!is.null(cf)) {
    out$cf_b0 <- coef(cf)["z"]; out$cf_b1 <- coef(cf)["tt(z)"]
    out$cf_auc <- auc_fun(out$cf_b0, out$cf_b1)
    out$theta <- tryCatch(cf$history[[1]]$theta,
                          error = function(e) tryCatch(cf$fvar, error = function(e2) NA_real_))
    out$cf_failed <- 0L
  }
  out
}

# ---------------- Run over all datasets, with progress + timing ----------------
# Set REUSE_FITS=1 to skip refitting and summarize the saved R_sim_fits.csv.
reuse_fits <- Sys.getenv("REUSE_FITS") == "1" && file.exists("R_sim_fits.csv")
t_all <- Sys.time()
if (reuse_fits) {
  r_results <- readr::read_csv("R_sim_fits.csv", show_col_types = FALSE)
  cat("REUSE_FITS=1: summarizing saved R_sim_fits.csv (no refitting).\n")
} else {
r_results <- imap_dfr(dirs, function(dir, th) {
  files <- sort(list.files(dir, pattern = "^simdata_\\d+\\.csv$", full.names = TRUE))
  if (length(files) == 0) { message("No simdata CSVs in ", dir); return(tibble()) }
  message("\n=== theta = ", th, " : ", length(files), " datasets in ", dir, " ===")
  t_dir <- Sys.time()
  res <- bind_rows(parallel::mclapply(seq_along(files), mc.cores = n_cores, function(i) {
    r <- fit_one(files[i])
    cat(sprintf("  [%3d/%3d] %s  fe=%.1fs ag=%.1fs cf=%.1fs  %s\n",
                i, length(files), basename(files[i]), r$t_fe, r$t_ag, r$t_cf,
                ifelse(isTRUE(r$cf_converged), "ok", "CF-NONCONV")))
    r
  }))
  message(sprintf("theta=%s done in %.1f min. CF converged %d/%d.",
                  th, as.numeric(difftime(Sys.time(), t_dir, units="mins")),
                  sum(res$cf_converged %in% TRUE), nrow(res)))
  res %>% mutate(theta_true = as.numeric(th))
})
message(sprintf("\nAll fits done in %.1f min total.",
                as.numeric(difftime(Sys.time(), t_all, units="mins"))))

readr::write_csv(r_results, "R_sim_fits.csv")
}

# ---------------- Convergence report ----------------
cat("\n================ R frailty convergence ================\n")
r_results %>% group_by(theta_true) %>%
  summarise(n = n(),
            converged = sum(cf_converged %in% TRUE),
            nonconverged = sum(!(cf_converged %in% TRUE)),
            .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\n================ Mean fit time per model (seconds) ================\n")
r_results %>% group_by(theta_true) %>%
  summarise(fe = mean(t_fe, na.rm=TRUE), ag = mean(t_ag, na.rm=TRUE),
            cf = mean(t_cf, na.rm=TRUE), .groups="drop") %>%
  mutate(across(c(fe,ag,cf), ~round(.x,2))) %>% as.data.frame() %>% print(row.names = FALSE)

# ---------------- R bias summary (converged CF only) ----------------
r_summary <- r_results %>%
  filter(cf_failed == 0, cf_converged %in% TRUE) %>%
  group_by(theta_true) %>%
  summarise(across(c(fe_b0,fe_b1,fe_auc, ag_b0,ag_b1,ag_auc,
                     cf_b0,cf_b1,cf_auc, theta), ~mean(.x, na.rm=TRUE)), .groups="drop") %>%
  rename(cf_theta = theta)

r_tidy <- r_summary %>%
  pivot_longer(-theta_true, names_to="col", values_to="mean") %>%
  separate(col, into=c("model","quantity"), sep="_", extra="merge") %>%
  mutate(model = recode(model, fe="First event", ag="Composite (AG)", cf="Primary (CF)"),
         truth = case_when(quantity=="b0"~true_b0, quantity=="b1"~true_b1,
                           quantity=="auc"~true_auc, quantity=="theta"~theta_true, TRUE~NA_real_),
         bias = mean - truth, package = "R")

# ---------------- Stata side (optional) ----------------
read_stata_side <- function(path, th) {
  if (!file.exists(path)) return(tibble())
  s <- readr::read_csv(path, show_col_types = FALSE)
  if ("cf_failed" %in% names(s)) s <- s %>% filter(cf_failed == 0)
  tibble(theta_true = as.numeric(th),
         `First event_b0`=mean(s$fe_b0hat,na.rm=TRUE), `First event_b1`=mean(s$fe_b1hat,na.rm=TRUE),
         `First event_auc`=mean(s$fe_auc,na.rm=TRUE),
         `Composite (AG)_b0`=mean(s$ag_b0hat,na.rm=TRUE), `Composite (AG)_b1`=mean(s$ag_b1hat,na.rm=TRUE),
         `Composite (AG)_auc`=mean(s$ag_auc,na.rm=TRUE),
         `Primary (CF)_b0`=mean(s$cf_b0hat,na.rm=TRUE), `Primary (CF)_b1`=mean(s$cf_b1hat,na.rm=TRUE),
         `Primary (CF)_auc`=mean(s$cf_auc,na.rm=TRUE), `Primary (CF)_theta`=mean(s$theta,na.rm=TRUE)
  ) %>% pivot_longer(-theta_true, names_to="col", values_to="mean") %>%
    separate(col, into=c("model","quantity"), sep="_", extra="merge") %>%
    mutate(truth = case_when(quantity=="b0"~true_b0, quantity=="b1"~true_b1,
                             quantity=="auc"~true_auc, quantity=="theta"~theta_true, TRUE~NA_real_),
           bias = mean - truth, package = "Stata")
}
stata_tidy <- imap_dfr(stata_csv, ~read_stata_side(.x, .y))

combined <- bind_rows(r_tidy, stata_tidy) %>% arrange(desc(theta_true), model, quantity, package)
readr::write_csv(combined, "R_vs_Stata_CF_bias.csv")

cat("\n================ CF (primary) AUC bias: R vs Stata ================\n")
combined %>% filter(model=="Primary (CF)", quantity=="auc") %>%
  mutate(mean=sprintf("%.4f",mean), bias=sprintf("%+.4f",bias)) %>%
  select(theta_true, package, mean, bias) %>% as.data.frame() %>% print(row.names=FALSE)

cat("\n================ CF theta estimate: R vs Stata ================\n")
combined %>% filter(model=="Primary (CF)", quantity=="theta") %>%
  mutate(mean=sprintf("%.3f",mean), bias=sprintf("%+.3f",bias)) %>%
  select(theta_true, package, mean, bias) %>% as.data.frame() %>% print(row.names=FALSE)

cat("\n================ R frailty variance collapsing to ~0 ================\n")
r_results %>% filter(cf_failed == 0) %>% group_by(theta_true) %>%
  summarise(n = n(), theta_lt_1e3 = sum(theta < 1e-3),
            pct_theta_near0 = round(100 * mean(theta < 1e-3), 1),
            max_theta = round(max(theta), 2), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)

cat("\nWrote R_sim_fits.csv (per-dataset) and R_vs_Stata_CF_bias.csv (summary).\n")

# ---------------- Recurrence proportion in simulated data ----------------
# Proportion of children with >=1 episode who had >=2 episodes (Results:
# 5.1% at theta = 0.2 and 8.1% at theta = 1; observed trial value 6.3%).

event_dispersion <- function(dir, theta_lab) {
  files <- list.files(dir, pattern = "^simdata_\\d+\\.csv$", full.names = TRUE)
  map_dfr(files, function(f) {
    d <- read.csv(f); names(d) <- tolower(names(d))
    # episodes per subject = number of event==1 rows per id
    per <- d %>% filter(event == 1) %>% count(id, name = "k")
    n_ge1 <- sum(per$k >= 1)
    tibble(
      file  = basename(f),
      n_ge1 = n_ge1,
      p_ge2 = sum(per$k >= 2) / n_ge1,
      p_ge3 = sum(per$k >= 3) / n_ge1
    )
  }) %>% mutate(theta_true = theta_lab)
}

disp <- bind_rows(
  event_dispersion(dirs[["1"]],   1),
  event_dispersion(dirs[["0.2"]], 0.2)
)

# average across the 100 replicates per theta
recur <- disp %>% group_by(theta_true) %>%
  summarise(
    n_datasets = n(),
    p_ge2_mean = mean(p_ge2), p_ge2_sd = sd(p_ge2),
    p_ge3_mean = mean(p_ge3), p_ge3_sd = sd(p_ge3),
    .groups = "drop"
  )
cat("\n================ Recurrence among children with >=1 episode ================\n")
print(as.data.frame(recur), row.names = FALSE)
readr::write_csv(recur, "simulation_recurrence_proportion.csv")
cat("Wrote simulation_recurrence_proportion.csv\n")
