###############################################################################
# model_selection_timeform.R
#
# Compares candidate forms of the time-varying term g(t) in
#   log HR(t) = b0 + b1 * g(t),   VE(t) = 1 - exp(b0 + b1 * g(t))
# for the Andersen-Gill (composite) and conditional frailty (primary) models,
# using information criteria (Methods: "Other candidate forms for g(t)
# (t, sqrt(t) and 1/t) were compared by information criteria and log t
# fitted best"). A proportional-hazards model (constant VE) is included as
# a reference. The first-event Cox model is shown for completeness.
#
# Implementation: coxph with tt() (exact evaluation of g(t) at event times),
# Breslow ties. The CF model uses two event-order strata (first episode /
# recurrences) and a gamma frailty with variance fixed at the strmcure
# estimate (theta = 0.8236), as in the bootstrap. All time-varying forms have
# the same number of parameters, so AIC and BIC rank them identically to the
# log-likelihood; BIC uses the number of events.
#
# Usage:  ROSE_DATA_DIR=/path/to/data Rscript model_selection_timeform.R
# Output: model_selection_timeform.csv (and console log)
# Requires: ROSE_ph.csv (from 1_data_preparation/prepare_data.R); survival, dplyr.
###############################################################################

suppressPackageStartupMessages({ library(survival); library(dplyr) })

data_file <- file.path(Sys.getenv("ROSE_DATA_DIR", "../data"), "ROSE_ph.csv")
theta_hat <- 0.8236
cat(sprintf("R %s; survival %s\n", getRversion(), packageVersion("survival")))

d <- read.csv(data_file); names(d) <- tolower(names(d))
if ("x" %in% names(d)) d$x <- NULL
d$ec_pool <- factor(pmin(d$event_count, 1))
n_ev <- sum(d$event)
cat(sprintf("Loaded %d rows, %d children, %d events.\n\n",
            nrow(d), length(unique(d$id)), n_ev))

# first-event data: follow-up to the first episode (or censoring)
d_fe <- d %>% arrange(id, start) %>% group_by(id) %>%
  mutate(ce = cumsum(event)) %>%
  filter(ce == 0 | (ce == 1 & event == 1)) %>%
  slice_tail(n = 1) %>% ungroup() %>% mutate(start = 0)

forms <- list(
  "PH (constant VE)" = NULL,
  "log(t)"           = function(t) log(t),
  "t"                = function(t) t,
  "sqrt(t)"          = function(t) sqrt(t),
  "1/t"              = function(t) 1 / t
)
ctrl <- coxph.control(timefix = FALSE, iter.max = 500, eps = 1e-6)

fit_form <- function(model, g) {
  rhs <- switch(model,
                FE = "group",
                AG = "group + cluster(id)",
                CF = "group + strata(ec_pool) + frailty(id, theta = theta_hat)")
  dat <- if (model == "FE") d_fe else d
  if (!is.null(g)) rhs <- paste(rhs, "+ tt(group)")
  f <- as.formula(paste("Surv(start, end, event) ~", rhs))
  fit <- if (is.null(g)) {
    coxph(f, data = dat, ties = "breslow", control = ctrl)
  } else {
    coxph(f, data = dat, ties = "breslow", control = ctrl,
          tt = function(x, t, ...) x * g(t))
  }
  k  <- if (is.null(g)) 1 else 2           # fixed-effect parameters
  ll <- fit$loglik[2]
  cf <- coef(fit)
  data.frame(b0 = unname(cf["group"]),
             b1 = if (is.null(g)) NA_real_ else unname(cf["tt(group)"]),
             loglik = ll, k = k,
             AIC = -2 * ll + 2 * k,
             BIC = -2 * ll + k * log(n_ev))
}

res <- do.call(rbind, lapply(c("AG", "CF", "FE"), function(m) {
  do.call(rbind, lapply(names(forms), function(nm) {
    cbind(model = m, g_t = nm, fit_form(m, forms[[nm]]))
  }))
}))

res <- res %>% group_by(model) %>%
  mutate(dAIC = AIC - min(AIC), best = ifelse(dAIC == 0, "<- best", "")) %>%
  ungroup()

for (m in c("AG", "CF", "FE")) {
  lab <- c(AG = "Andersen-Gill (composite effect)",
           CF = sprintf("Conditional frailty (primary effect; theta fixed %.4f)", theta_hat),
           FE = "First-event Cox model")[m]
  cat("==== ", lab, " ====\n", sep = "")
  out <- res %>% filter(model == m) %>% select(-model) %>% as.data.frame()
  num <- sapply(out, is.numeric)
  out[num] <- lapply(out[num], function(x) round(x, 3))
  print(out, row.names = FALSE)
  cat("\n")
}

write.csv(res, "model_selection_timeform.csv", row.names = FALSE)
cat("Wrote model_selection_timeform.csv\n")
