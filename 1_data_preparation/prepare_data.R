###############################################################################
# prepare_data.R
#
# Builds the analysis datasets for the Rotasiil (ROSE) trial re-analysis from
# the two Epicentre extraction files, and tabulates Table 1.
#
# Inputs  (in ROSE_DATA_DIR; trial data are NOT distributed with this code):
#   Baseline.xlsx   one row per child (arm, dose dates, last visit)
#   Outcome.xlsx    one row per diarrhoea episode (dates, ELISA result)
#
# Outputs (in ROSE_OUT_DIR, default = ROSE_DATA_DIR):
#   ROSE_ph.csv        counting-process data, one row per at-risk interval
#                      (ID, group, end_dof, start, end, event, event_count);
#                      used by all R scripts (coxph + tt()).
#   ROSE_split_tt.csv  the same data split at every observed end time, with
#                      vac.time = group * log(interval midpoint); used by the
#                      Stata stcox / strmcure analysis (imported as `vactime`).
# Output (in the working directory; aggregate counts only):
#   table1_episode_counts.csv   Table 1 of the manuscript
#
# Definitions (Methods):
#   * Time origin: 28 days after dose 3; time scale in months (365.25/12 days).
#   * Episode: gastroenteritis with a rotavirus-positive stool (ELISA_RESULT == 2).
#   * After an episode a child is not at risk for 5 days; the next at-risk
#     interval starts at max(end of previous episode, previous onset + 5 days).
#   * Follow-up is truncated at 20 months.
#   * event_count = number of previous episodes (event-order stratum).
#
# Usage:  ROSE_DATA_DIR=/path/to/data Rscript prepare_data.R
# Requires: tidyverse, readxl, survival.
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(survival)
})

in_dir  <- Sys.getenv("ROSE_DATA_DIR", "../data")
out_dir <- Sys.getenv("ROSE_OUT_DIR", in_dir)
a <- 0    # start of follow-up (months)
b <- 20   # end of follow-up (months)

# ---------------- Read and join ----------------
baseline <- read_excel(file.path(in_dir, "Baseline.xlsx")) %>%
  mutate(d1 = as.Date(DATE_DOSE_1),
         d2 = as.Date(DATE_DOSE_2),
         d3 = as.Date(DATE_DOSE_3),
         end_dof = as.Date(DATE_LAST_VISIT),
         group = case_when(GROUP == "BRV-PV Vaccine" ~ 1, .default = 0),
         age_d1 = AGE_DOSE_1,
         age_d2 = AGE_DOSE_2,
         age_d3 = AGE_DOSE_3)

outcome <- read_excel(file.path(in_dir, "Outcome.xlsx")) %>%
  mutate(start.epi = as.Date(DATE_BEGIN_EPI),
         end.epi = as.Date(DATE_END_EPI)) %>%
  dplyr::select(-DATE_BEGIN_EPI, -DATE_END_EPI)

# rotavirus gastroenteritis episodes of any severity
outcome.rota <- filter(outcome, ELISA_RESULT == 2)

data.rota <- left_join(x = baseline, y = outcome.rota, by = "ID") %>%
  # months since start of follow-up (28 days after dose 3)
  mutate(end_dof   = (as.double(difftime(end_dof,   as.Date(d3), units = "days")) - 28) * 12 / 365.25,
         start.epi = (as.double(difftime(start.epi, as.Date(d3), units = "days")) - 28) * 12 / 365.25,
         end.epi   = (as.double(difftime(end.epi,   as.Date(d3), units = "days")) - 28) * 12 / 365.25,
         event = case_when(is.na(end.epi) == TRUE ~ 0, .default = 1),
         end   = case_when(event == 0 ~ end_dof, .default = start.epi)) %>%
  arrange(ID, end) %>%
  group_by(ID) %>%
  # next interval starts at the later of previous episode end or onset + 5 days
  mutate(start = pmax(lag(end.epi), lag(end) + 5 * 12 / 365.25, na.rm = TRUE),
         start = case_when(is.na(start) == TRUE ~ 0, .default = start)) %>%
  filter(start < end) %>%
  mutate(event_count = row_number() - 1) %>%
  ungroup() %>%
  dplyr::select(ID, group, age_d1:age_d3, d1, d2, d3, end_dof, start, end, event, event_count, end.epi)

# ---------------- Restrict to (a, b] months ----------------
events <- filter(data.rota, event == 1 & end > a & end <= b) %>%
  group_by(ID) %>%
  mutate(start = pmax(lag(end.epi), lag(end) + 5 * 12 / 365.25, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(start = case_when(is.na(start) == TRUE ~ a, .default = start),
         end_dof = pmin(b, end_dof)) %>%
  dplyr::select(ID, group, age_d1:age_d3, d1, d2, d3, end_dof, start, end, event, event_count)

# censoring interval after each child's last episode
event.censoring <- events %>%
  group_by(ID) %>%
  arrange(ID, event_count) %>%
  summarise(group = last(group), age_d1 = last(age_d1), age_d2 = last(age_d2),
            age_d3 = last(age_d3), d1 = last(d1), d2 = last(d2), d3 = last(d3),
            end_dof = last(end_dof), end = last(end), event_count = last(event_count)) %>%
  mutate(start = end + 5 * 12 / 365.25,
         end = end_dof,
         event = 0,
         event_count = event_count + 1) %>%
  dplyr::select(ID, group, age_d1:age_d3, d1, d2, d3, end_dof, start, end, event, event_count) %>%
  filter(start < end)

# children with no episode in (a, b]
other <- anti_join(data.rota, events, by = c("ID", "end", "event")) %>%
  filter(!ID %in% pull(events, ID)) %>%
  mutate(start = a, end = pmin(b, end_dof), event = 0) %>%
  group_by(ID) %>%
  arrange(ID, event_count) %>%
  summarise(group = last(group), age_d1 = last(age_d1), age_d2 = last(age_d2),
            age_d3 = last(age_d3), d1 = last(d1), d2 = last(d2), d3 = last(d3),
            end_dof = last(end_dof), start = last(start), end = last(end),
            event = last(event), event_count = last(event_count)) %>%
  ungroup() %>%
  mutate(event = 0) %>%
  dplyr::select(ID, group, age_d1:age_d3, d1, d2, d3, end_dof, start, end, event, event_count)

analysis_a2b <- rbind(events, event.censoring, other) %>%
  arrange(ID, start) %>%
  filter(start < end)

analysis_initial <- analysis_a2b %>%
  dplyr::select(ID, group, end_dof, start, end, event, event_count)

cat(sprintf("Analysis data: %d rows, %d children (%d vaccine, %d placebo), %d episodes.\n",
            nrow(analysis_initial), n_distinct(analysis_initial$ID),
            n_distinct(analysis_initial$ID[analysis_initial$group == 1]),
            n_distinct(analysis_initial$ID[analysis_initial$group == 0]),
            sum(analysis_initial$event)))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(analysis_initial, file.path(out_dir, "ROSE_ph.csv"))

# ---------------- Table 1: children by number of episodes ----------------
kids <- analysis_initial %>%
  group_by(ID, group) %>% summarise(k = sum(event), .groups = "drop")
tab_arm <- function(x) c(`0` = sum(x == 0), `1` = sum(x == 1), `2` = sum(x == 2),
                         `3` = sum(x == 3), `>=4` = sum(x >= 4),
                         `>=1` = sum(x >= 1), `>=2` = sum(x >= 2),
                         `Total episodes` = sum(x), N = length(x))
table1 <- data.frame(
  row       = names(tab_arm(kids$k)),
  vaccine   = tab_arm(kids$k[kids$group == 1]),
  placebo   = tab_arm(kids$k[kids$group == 0]),
  total     = tab_arm(kids$k), row.names = NULL)
cat("\nTable 1. Children by number of rotavirus gastroenteritis episodes:\n")
print(table1, row.names = FALSE)
cat(sprintf("Recurrence among children with >=1 episode: %d/%d = %.1f%%\n",
            sum(kids$k >= 2), sum(kids$k >= 1), 100 * mean(kids$k[kids$k >= 1] >= 2)))
cat(sprintf("Followed >=19 months: vaccine %.1f%%, placebo %.1f%%\n",
            100 * mean(analysis_initial$end_dof[analysis_initial$event_count == 0 & analysis_initial$group == 1] >= 19),
            100 * mean(analysis_initial$end_dof[analysis_initial$event_count == 0 & analysis_initial$group == 0] >= 19)))
write.csv(table1, "table1_episode_counts.csv", row.names = FALSE)

# ---------------- Split at every observed end time (for Stata) ----------------
# Splitting at all observed times lets Stata's stcox/strmcure evaluate the
# time-varying term group*log(t) on a fine grid, mimicking R's tt().
cut_points <- sort(analysis_initial$end, decreasing = FALSE)

analysis_split <-
  survSplit(Surv(start, end, event) ~ ., data = analysis_initial, cut = cut_points,
            end = "end", event = "event", start = "start") %>%
  mutate(time = (start + end) / 2,
         log_t = log(time),
         sqrt_t = sqrt(time),
         vac.time.t = group * time,
         vac.time.sqrt = group * sqrt(time),
         vac.time = group * log(time),
         interval = as.factor(findInterval(start, cut_points) + 1)) %>%
  filter(!(end - start < 1e-8 & event == 0)) %>%
  arrange(ID, start) %>%
  group_by(ID) %>%
  # move an event recorded on a zero-length interval to the preceding interval
  mutate(event_to_move = (end - start < 1e-8),
         event_above = lead(event_to_move, default = FALSE),
         event = case_when(event_above == TRUE ~ 1, .default = event)) %>%
  ungroup() %>%
  filter(!(event_to_move == TRUE)) %>%
  arrange(ID, start) %>%
  dplyr::select(-event_to_move, -event_above)

cat(sprintf("\nSplit data: %d rows, %d events.\n", nrow(analysis_split), sum(analysis_split$event)))
write.csv(analysis_split, file.path(out_dir, "ROSE_split_tt.csv"))
cat("Wrote ROSE_ph.csv and ROSE_split_tt.csv to ROSE_OUT_DIR.\n")
