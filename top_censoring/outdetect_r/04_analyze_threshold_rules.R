#!/usr/bin/env Rscript

# Recalculate the 5% welfare-share and 1% population-share rules using
# exact R outdetect results and the corresponding LIS results.

COMPARISON <- file.path("top_censoring", "outdetect_r", "outputs", "lis_vs_outdetect_comparison.csv")
LIS_POPULATION <- file.path("top_censoring", "outdetect_r", "outputs", "lis_population_rule_comparison.csv")
OUT_DIR <- file.path("top_censoring", "outdetect_r", "outputs")

comparison <- read.csv(COMPARISON, stringsAsFactors = FALSE, check.names = FALSE)
lis_population <- read.csv(LIS_POPULATION, stringsAsFactors = FALSE, check.names = FALSE)

keys <- c("Country", "Year", "Survey")
rule_data <- merge(
  comparison,
  lis_population[, c(keys, "LIS population share above ceiling (%)")],
  by = keys,
  all.x = TRUE
)

rule_data$`LIS: welfare >= 5%` <-
  as.numeric(rule_data[["Welfare share above ceiling (%)"]]) >= 5
rule_data$`outdetect: welfare >= 5%` <-
  as.numeric(rule_data[["Welfare share flagged (%)"]]) >= 5
rule_data$`LIS: population >= 1%` <-
  as.numeric(rule_data[["LIS population share above ceiling (%)"]]) >= 1
rule_data$`outdetect: population >= 1%` <-
  as.numeric(rule_data[["Population share flagged (%)"]]) >= 1

summary <- data.frame(
  method = c("LIS", "outdetect", "LIS", "outdetect"),
  rule = c("Welfare share >= 5%", "Welfare share >= 5%", "Population share >= 1%", "Population share >= 1%"),
  surveys_flagged = c(
    sum(rule_data$`LIS: welfare >= 5%`, na.rm = TRUE),
    sum(rule_data$`outdetect: welfare >= 5%`, na.rm = TRUE),
    sum(rule_data$`LIS: population >= 1%`, na.rm = TRUE),
    sum(rule_data$`outdetect: population >= 1%`, na.rm = TRUE)
  )
)

top10 <- rule_data[rule_data$`LIS: welfare >= 5%`, ]
top10 <- top10[order(-as.numeric(top10[["Welfare share above ceiling (%)"]])), ]
top10 <- top10[seq_len(min(10, nrow(top10))), c(
  keys,
  "Welfare share above ceiling (%)",
  "Welfare share flagged (%)",
  "LIS population share above ceiling (%)",
  "Population share flagged (%)",
  "LIS: welfare >= 5%",
  "outdetect: welfare >= 5%",
  "LIS: population >= 1%",
  "outdetect: population >= 1%"
)]

four_cases <- rule_data[
  rule_data$`LIS: welfare >= 5%` & rule_data$`outdetect: population >= 1%`,
  c(
    keys,
    "Welfare share above ceiling (%)",
    "Welfare share flagged (%)",
    "LIS population share above ceiling (%)",
    "Population share flagged (%)"
  )
]
four_cases <- four_cases[order(four_cases$Country, four_cases$Year), ]

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(summary, file.path(OUT_DIR, "threshold_rule_summary_exact_r.csv"), row.names = FALSE)
write.csv(top10, file.path(OUT_DIR, "threshold_rule_top10_exact_r.csv"), row.names = FALSE)
write.csv(four_cases, file.path(OUT_DIR, "threshold_rule_four_cases_exact_r.csv"), row.names = FALSE)

message("Wrote exact-R threshold rule analysis to ", OUT_DIR)
print(summary, row.names = FALSE)
message("\nFour cross-method cases: ")
print(four_cases, row.names = FALSE)