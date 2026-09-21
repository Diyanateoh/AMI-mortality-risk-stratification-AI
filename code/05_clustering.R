library(clustMixType)
library(dplyr)
library(pheatmap)

analysis_seed <- 123L
d <- readRDS("derived/analysis_data.rds")

vars_30d <- c(
  "age", "sex", "stemi", "bmi", "multivessel_disease",
  "residential_area", "onset_season"
)

vars_1y <- c(
  "age", "treatment_within_12h", "residential_area",
  "hypertension", "diabetes"
)

vars_5y <- c(
  "age", "treatment_within_12h", "diabetes", "sex", "dyslipidemia"
)

prepare_cluster_data <- function(data, variables) {
  x <- data[, variables, drop = FALSE]

  for (v in intersect(c("age", "bmi"), variables)) {
    x[[v]] <- as.numeric(scale(x[[v]]))
  }

  x
}

fit_kproto <- function(x, seed = 123L) {
  k_values <- 2:6
  elbow <- numeric(length(k_values))

  for (i in seq_along(k_values)) {
    set.seed(seed)
    elbow[i] <- kproto(x, k = k_values[i], verbose = FALSE)$tot.withinss
  }

  set.seed(seed)
  fit <- kproto(x, k = 3, verbose = FALSE)

  list(
    fit = fit,
    elbow = data.frame(k = k_values, total_withinss = elbow)
  )
}

profile_matrix <- function(data, variables, cluster, mortality) {
  mm <- model.matrix(~ . - 1, data = data[, variables, drop = FALSE])
  profile <- aggregate(mm, by = list(cluster = cluster), FUN = mean)

  mat <- t(as.matrix(profile[, -1, drop = FALSE]))
  colnames(mat) <- paste0("Cluster ", profile$cluster)

  mortality_rate <- as.numeric(tapply(mortality, cluster, mean))
  mat <- rbind(mat, observed_mortality = mortality_rate)

  scaled <- t(scale(t(mat)))
  scaled[is.na(scaled)] <- 0
  scaled
}

x_30d <- prepare_cluster_data(d$short, vars_30d)
x_1y <- prepare_cluster_data(d$long, vars_1y)
x_5y <- prepare_cluster_data(d$long, vars_5y)

cluster_30d <- fit_kproto(x_30d, analysis_seed)
cluster_1y <- fit_kproto(x_1y, analysis_seed)
cluster_5y <- fit_kproto(x_5y, analysis_seed)

d$short$cluster_30d <- factor(cluster_30d$fit$cluster)
d$long$cluster_1y <- factor(cluster_1y$fit$cluster)
d$long$cluster_5y <- factor(cluster_5y$fit$cluster)

summary_30d <- d$short %>%
  group_by(cluster_30d) %>%
  summarise(n = n(), mortality_rate = mean(mortality_30d) * 100, .groups = "drop")

summary_1y <- d$long %>%
  group_by(cluster_1y) %>%
  summarise(n = n(), mortality_rate = mean(event_1y) * 100, .groups = "drop")

summary_5y <- d$long %>%
  group_by(cluster_5y) %>%
  summarise(n = n(), mortality_rate = mean(event_5y) * 100, .groups = "drop")

heat_30d <- profile_matrix(
  d$short, vars_30d, d$short$cluster_30d, d$short$mortality_30d
)

heat_1y <- profile_matrix(
  d$long, vars_1y, d$long$cluster_1y, d$long$event_1y
)

heat_5y <- profile_matrix(
  d$long, vars_5y, d$long$cluster_5y, d$long$event_5y
)

dir.create("outputs", showWarnings = FALSE)

write.csv(cluster_30d$elbow, "outputs/elbow_30d.csv", row.names = FALSE)
write.csv(cluster_1y$elbow, "outputs/elbow_1y.csv", row.names = FALSE)
write.csv(cluster_5y$elbow, "outputs/elbow_5y.csv", row.names = FALSE)

write.csv(summary_30d, "outputs/cluster_summary_30d.csv", row.names = FALSE)
write.csv(summary_1y, "outputs/cluster_summary_1y.csv", row.names = FALSE)
write.csv(summary_5y, "outputs/cluster_summary_5y.csv", row.names = FALSE)

palette <- colorRampPalette(c("green", "yellow", "red"))(100)

pheatmap(heat_30d, cluster_rows = FALSE, cluster_cols = FALSE, color = palette)
pheatmap(heat_1y, cluster_rows = FALSE, cluster_cols = FALSE, color = palette)
pheatmap(heat_5y, cluster_rows = FALSE, cluster_cols = FALSE, color = palette)

print(summary_30d)
print(summary_1y)
print(summary_5y)
