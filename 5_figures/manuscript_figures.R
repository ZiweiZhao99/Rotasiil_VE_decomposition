###############################################################################
# manuscript_figures.R   (parallel, timed, fixed two-strata pooling)
#
# Reproduces the two main-text figures and the AUC-scale secondary-effect CI:
#   Figure 1 - Three VE(t) curves (First event/Cox, Composite/AG, Primary/CF)
#   Figure 2 - Secondary effect VE_s(t) = VE_c(t) - VE_p(t) with bootstrap band
#   Table 2  - AUC-scale secondary effect with 95% CI, on the normalized
#              (average-VE) scale, and pointwise VE at 1, 6, 12, 20 months.
#
# Models are fit from ROSE_ph.csv with coxph + tt() (exact-time evaluation of
# the log-t term). The frailty variance is taken from the STATA strmcure fit
# (theta_hat) and held fixed in the bootstrap.
#
# EVENT-ORDER POOLING (fixed two strata):
#   Episodes are pooled into exactly two event-order strata: first episodes
#   and all recurrences (2nd, 3rd, ...). This matches the full-data fit, in
#   which recurrences already form a single stratum. Fixing the strata means
#   the conditional frailty model never collapses onto the Andersen-Gill
#   model in low-recurrence resamples, so every resample is estimable and no
#   screening is needed. (The earlier adaptive 5%-threshold pooling could drop
#   to one stratum in sparse resamples; that issue does not arise here.)
#
# Point estimates used in the manuscript (tt splitting):
#   FE : b0 = -0.532, b1 = 0.156
#   AG : b0 = -0.549, b1 = 0.186
#   CF : b0 = -0.572, b1 = 0.178, theta = 0.824  (strmcure, STATA)
#
# Bootstrap: B = 1000 resamples of children. Replicate b uses
# set.seed(base_seed + b) with the same resampling code as
# 4_bootstrap/rose_boot_optionC.R, so the AUC-scale draws (and CI) are
# identical to that run; replicates run in parallel (forked workers), which
# does not change the results because each replicate sets its own seed.
#
# Writes: Figure1_VE_curves.(png|pdf), Figure2_secondary_effect.(png|pdf),
#         Table2_pointwise_VE.csv, Figure2_secondary_band.csv, run log.
#
# Usage:  ROSE_DATA_DIR=/path/to/data Rscript manuscript_figures.R
# Requires: ROSE_ph.csv (from 1_data_preparation/prepare_data.R);
#           packages survival, ggplot2, dplyr, scales, parallel.
###############################################################################

suppressPackageStartupMessages({
  library(survival); library(ggplot2); library(dplyr); library(scales); library(parallel)
})

# ---------------- Settings ----------------
data_file   <- file.path(Sys.getenv("ROSE_DATA_DIR", "../data"), "ROSE_ph.csv")
Tau         <- 20
theta_hat   <- 0.8236            # strmcure frailty variance (2_main_analysis_stata)
B           <- 1000
n_cores     <- as.integer(Sys.getenv("N_CORES", "4"))
base_seed   <- 20260101
out_prefix  <- "Figure"
log_file    <- "manuscript_figures_runlog.txt"

# ---------------- Logging ----------------
log_msg <- function(...) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...))
  cat(line, "\n"); cat(line, "\n", file = log_file, append = TRUE)
}
cat(sprintf("=== manuscript_figures.R run log - started %s ===\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = log_file, append = FALSE)
t_start <- Sys.time()
log_msg("Settings: B=%d, parallel bootstrap (%d cores), theta_hat=%.4f, Tau=%d, base_seed=%d",
        B, n_cores, theta_hat, Tau, base_seed)
log_msg("R %s; survival %s", getRversion(), packageVersion("survival"))

# ---------------- Data ----------------
d <- read.csv(data_file); names(d) <- tolower(names(d))
if ("x" %in% names(d)) d$x <- NULL
stopifnot(all(c("id","group","start","end","event","event_count") %in% names(d)))
log_msg("Loaded %s: %d rows, %d persons, %d events.",
        basename(data_file), nrow(d), length(unique(d$id)), sum(d$event))

# ---------------- Functions ----------------
tt_fun  <- function(x, t, ...) x * log(t)
ve_fun  <- function(t, b0, b1) 1 - exp(b0 + b1 * log(t))
auc_fun <- function(b0, b1) 1 - exp(b0) * Tau^b1 / (b1 + 1)   # normalized (average-VE) AUC
make_pool <- function(dat) {
  # Fixed two strata: 0 = first episode, 1 = any recurrence (2nd, 3rd, ...).
  # Never collapses to one stratum, so every resample is estimable.
  factor(pmin(dat$event_count, 1), levels = c(0, 1),
         labels = c("First event", "Recurrent events"))
}
d$ec_pool <- make_pool(d)

# ---------------- Full-data fits ----------------
t_fit0 <- Sys.time()
d_fe <- d %>% arrange(id, start) %>% group_by(id) %>%
  mutate(ce = cumsum(event)) %>%
  filter(ce == 0 | (ce == 1 & event == 1)) %>%
  slice_tail(n = 1) %>% ungroup() %>% mutate(start = 0)
fe <- coxph(Surv(start, end, event) ~ group + tt(group), data = d_fe,
            ties = "breslow", tt = tt_fun, control = coxph.control(timefix = FALSE, iter.max = 500, eps = 1e-6))
fe_b0 <- coef(fe)[1]; fe_b1 <- coef(fe)[2]
ag <- coxph(Surv(start, end, event) ~ group + tt(group) + cluster(id), data = d,
            ties = "breslow", tt = tt_fun, control = coxph.control(timefix = FALSE, iter.max = 500, eps = 1e-6))
ag_b0 <- coef(ag)[1]; ag_b1 <- coef(ag)[2]
cf <- coxph(Surv(start, end, event) ~ group + tt(group) + strata(ec_pool) +
              frailty(id, theta = theta_hat), data = d,
            ties = "breslow", tt = tt_fun, control = coxph.control(timefix = FALSE, iter.max = 500, eps = 1e-6))
cf_b0 <- coef(cf)["group"]; cf_b1 <- coef(cf)["tt(group)"]
log_msg("Full-data fits done in %.1f s.", as.numeric(difftime(Sys.time(), t_fit0, units="secs")))
log_msg("FE: b0=%.3f b1=%.3f AUC=%.3f", fe_b0, fe_b1, auc_fun(fe_b0, fe_b1))
log_msg("AG: b0=%.3f b1=%.3f AUC=%.3f", ag_b0, ag_b1, auc_fun(ag_b0, ag_b1))
log_msg("CF: b0=%.3f b1=%.3f AUC=%.3f (theta fixed %.3f)", cf_b0, cf_b1, auc_fun(cf_b0, cf_b1), theta_hat)
sec_auc_pt <- auc_fun(ag_b0, ag_b1) - auc_fun(cf_b0, cf_b1)
log_msg("Secondary AUC point (AG-CF, normalized): %.4f", sec_auc_pt)

# ---------------- Sequential bootstrap ----------------
t_grid <- seq(0, Tau, by = 0.2)

boot_one <- function(b, d, rows_by_id, uid, t_grid, theta_hat, Tau, base_seed,
                     tt_fun, ve_fun, auc_fun, make_pool) {
  set.seed(base_seed + b)
  samp <- sample(uid, length(uid), replace = TRUE)
  idx  <- unlist(rows_by_id[as.character(samp)])
  db   <- d[idx, , drop = FALSE]
  db$id <- rep(seq_along(samp), lengths(rows_by_id[as.character(samp)]))
  db <- db[order(db$id, db$start), ]
  db$event_count <- ave(db$event, db$id, FUN = function(e) cumsum(e) - e)
  db$ec_pool <- make_pool(db)
  n_second <- sum(db$event == 1 & db$event_count >= 1)
  n_strata <- nlevels(db$ec_pool)
  
  agb <- tryCatch(coxph(Surv(start,end,event) ~ group + tt(group) + cluster(id),
                        data=db, ties="breslow", tt=tt_fun,
                        control=coxph.control(timefix=FALSE)), error=function(e) NULL)
  cfb <- tryCatch(coxph(Surv(start,end,event) ~ group + tt(group) + strata(ec_pool) +
                          frailty(id, theta=theta_hat),
                        data=db, ties="breslow", tt=tt_fun,
                        control=coxph.control(timefix=FALSE)), error=function(e) NULL)
  if (is.null(agb) || is.null(cfb)) {
    return(list(curve=rep(NA_real_, length(t_grid)), auc=NA_real_,
                n_second=n_second, n_strata=n_strata, coef_gap=NA_real_, failed=TRUE))
  }
  a0 <- coef(agb)[1]; a1 <- coef(agb)[2]
  p0 <- coef(cfb)["group"]; p1 <- coef(cfb)["tt(group)"]
  list(curve   = ve_fun(t_grid, a0, a1) - ve_fun(t_grid, p0, p1),
       auc     = auc_fun(a0, a1) - auc_fun(p0, p1),
       n_second= n_second, n_strata = n_strata, coef_gap = a1 - p1, failed = FALSE)
}

uid <- unique(d$id)
rows_by_id <- split(seq_len(nrow(d)), d$id)

log_msg("Starting parallel bootstrap: B=%d on %d cores...", B, n_cores)
t_boot0 <- Sys.time()
progress_every <- max(1L, floor(B / 20L))   # report approximately every 5%

boot_list <- mclapply(seq_len(B), mc.cores = n_cores, function(b) {
  ans <- boot_one(
    b = b,
    d = d,
    rows_by_id = rows_by_id,
    uid = uid,
    t_grid = t_grid,
    theta_hat = theta_hat,
    Tau = Tau,
    base_seed = base_seed,
    tt_fun = tt_fun,
    ve_fun = ve_fun,
    auc_fun = auc_fun,
    make_pool = make_pool
  )
  
  if (b %% progress_every == 0L || b == B) {
    elapsed <- as.numeric(difftime(Sys.time(), t_boot0, units = "secs"))
    log_msg(
      "Bootstrap progress: %d/%d (%.1f%%); elapsed %.1f s.",
      b, B, 100 * b / B, elapsed
    )
  }
  
  ans
})

boot_secs <- as.numeric(difftime(Sys.time(), t_boot0, units = "secs"))
log_msg(
  "Bootstrap done in %.1f s (%.2f s per replicate, parallel wall-clock).",
  boot_secs, boot_secs / B
)

# ---------------- Assemble ----------------
boot_curves <- do.call(rbind, lapply(boot_list, `[[`, "curve"))
boot_auc    <- vapply(boot_list, `[[`, numeric(1), "auc")
n_second    <- vapply(boot_list, `[[`, numeric(1), "n_second")
failed      <- vapply(boot_list, `[[`, logical(1), "failed")

n_fail <- sum(failed)
log_msg("Replicates: %d total, %d failed. (Fixed two-strata pooling; all estimable.)",
        B, n_fail)
log_msg("Recurrences (2nd+ episodes) per resample: median %.0f.", median(n_second[!failed]))

# Pointwise band for Figure 2 (all valid resamples)
sec_lo <- apply(boot_curves[!failed, , drop=FALSE], 2, quantile, 0.025, na.rm = TRUE)
sec_hi <- apply(boot_curves[!failed, , drop=FALSE], 2, quantile, 0.975, na.rm = TRUE)
sec_pt <- ve_fun(t_grid, ag_b0, ag_b1) - ve_fun(t_grid, cf_b0, cf_b1)

# Pointwise values reported in the text and Table 2
t_rep <- c(1, 6, 12, 20)
ix <- match(t_rep, round(t_grid, 6))
tab2 <- data.frame(
  month = t_rep,
  VE_first_event = ve_fun(t_rep, fe_b0, fe_b1),
  VE_composite   = ve_fun(t_rep, ag_b0, ag_b1),
  VE_primary     = ve_fun(t_rep, cf_b0, cf_b1),
  VE_secondary   = sec_pt[ix], sec_lo = sec_lo[ix], sec_hi = sec_hi[ix])
write.csv(tab2, "Table2_pointwise_VE.csv", row.names = FALSE)
for (i in seq_len(nrow(tab2))) with(tab2[i, ], log_msg(
  "t=%2d mo: VE FE=%.3f AG=%.3f CF=%.3f | secondary=%.4f (95%% band %.4f, %.4f)",
  month, VE_first_event, VE_composite, VE_primary, VE_secondary, sec_lo, sec_hi))
log_msg("Pointwise band upper limit, max over (0, %d]: %.4f (band excludes zero if < 0)",
        Tau, max(sec_hi[t_grid > 0]))
write.csv(data.frame(t = t_grid, sec = sec_pt, lo = sec_lo, hi = sec_hi),
          "Figure2_secondary_band.csv", row.names = FALSE)

# AUC-scale CI (single interval; no screen needed)
auc_ci <- quantile(boot_auc[!failed], c(0.025, 0.975), na.rm = TRUE)
log_msg("Table 2 AUC secondary effect: point=%.4f  95%% CI (%.4f, %.4f)",
        sec_auc_pt, auc_ci[1], auc_ci[2])

# ---------------- FIGURE 1 ----------------
curve_df <- bind_rows(
  data.frame(t=t_grid, VE=ve_fun(t_grid, fe_b0, fe_b1), Model="First event (Cox)"),
  data.frame(t=t_grid, VE=ve_fun(t_grid, ag_b0, ag_b1), Model="Composite (AG)"),
  data.frame(t=t_grid, VE=ve_fun(t_grid, cf_b0, cf_b1), Model="Primary (CF)"))
curve_df$Model <- factor(curve_df$Model, levels=c("Composite (AG)","Primary (CF)","First event (Cox)"))
pal <- c("Composite (AG)"="#C0392B","Primary (CF)"="#1E8449","First event (Cox)"="#2E86C1")
fig1 <- ggplot(curve_df, aes(t, VE, colour=Model, linetype=Model)) +
  geom_hline(yintercept=0, colour="grey70", linewidth=0.3) +
  geom_line(linewidth=1.1) +
  scale_colour_manual(values=pal) +
  scale_linetype_manual(values=c("Composite (AG)"="solid","Primary (CF)"="longdash","First event (Cox)"="dotted")) +
  scale_x_continuous(breaks=seq(0,20,2), expand=c(0.01,0)) +
  scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), 
                     labels=scales::percent_format(accuracy=1), 
                     limits=c(-0.01,0.55)) +
  labs(x="Months since vaccination", y="Vaccine efficacy VE(t)", colour=NULL, linetype=NULL) +
  theme_classic(base_size=12) +
  theme(legend.position=c(0.78,0.85), legend.background=element_rect(fill="white", colour="grey80"))
ggsave(paste0(out_prefix,"1_VE_curves.png"), fig1, width=7, height=5, dpi=300)
ggsave(paste0(out_prefix,"1_VE_curves.pdf"), fig1, width=7, height=5)
log_msg("Wrote %s1_VE_curves.png / .pdf", out_prefix)

# ---------------- FIGURE 2 ----------------
sec_df <- data.frame(t=t_grid, sec=sec_pt, lo=sec_lo, hi=sec_hi)
fig2 <- ggplot(sec_df, aes(t, sec)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), fill="#C0392B", alpha=0.15) +
  geom_hline(yintercept=0, colour="grey60", linetype="dashed", linewidth=0.4) +
  geom_line(colour="#922B21", linewidth=1.1) +
  scale_x_continuous(breaks=seq(0,20,2), expand=c(0.01,0)) +
  scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
  labs(x="Months since vaccination",
       y=expression("Secondary effect  " * VE[s](t) == VE[c](t) - VE[p](t))) +
  theme_classic(base_size=12)
ggsave(paste0(out_prefix,"2_secondary_effect.png"), fig2, width=7, height=5, dpi=300)
ggsave(paste0(out_prefix,"2_secondary_effect.pdf"), fig2, width=7, height=5)
log_msg("Wrote %s2_secondary_effect.png / .pdf", out_prefix)

# ---------------- Summary ----------------
total_secs <- as.numeric(difftime(Sys.time(), t_start, units="secs"))
log_msg("SECONDARY EFFECT (AUC): point=%.4f  95%% CI (%.4f, %.4f)",
        sec_auc_pt, auc_ci[1], auc_ci[2])
log_msg("Failed resamples: %d of %d.", n_fail, B)
log_msg("TOTAL wall-clock: %.1f s (%.1f min).", total_secs, total_secs/60)