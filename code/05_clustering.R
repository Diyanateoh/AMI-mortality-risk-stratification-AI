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

fit_kproto <- function(x, selected_k = 3L, seed = 123L) {
  k_values <- 2:6
  elbow <- numeric(length(k_values))

  for (i in seq_along(k_values)) {
    set.seed(seed)
    elbow[i] <- kproto(
      x,
      k = k_values[i],
      iter.max = 100,
      nstart = 10,
      verbose = FALSE
    )$tot.withinss
  }

  set.seed(seed)
  fit <- kproto(
    x,
    k = selected_k,
    iter.max = 100,
    nstart = 10,
    verbose = FALSE
  )

  list(
    fit = fit,
    selected_k = selected_k,
    elbow = data.frame(
      k = k_values,
      total_withinss = elbow
    )
  )
}

cluster_profile <- function(data, variables, cluster, mortality) {
  out <- data.frame(
    cluster = levels(factor(cluster)),
    n = as.integer(table(factor(cluster))),
    mortality_rate = 100 * as.numeric(tapply(mortality, cluster, mean))
  )

  for (v in variables) {
    if (is.numeric(data[[v]])) {
      out[[v]] <- as.numeric(tapply(data[[v]], cluster, mean))
    } else {
      for (level in levels(data[[v]])) {
        name <- paste(v, level, sep = "__")
        out[[name]] <- 100 * as.numeric(
          tapply(data[[v]] == level, cluster, mean)
        )
      }
    }
  }

  out
}

cluster_tests <- function(data, variables, cluster) {
  p <- sapply(variables, function(v) {
    if (is.numeric(data[[v]])) {
      summary(aov(data[[v]] ~ factor(cluster)))[[1]][["Pr(>F)"]][1]
    } else {
      suppressWarnings(
        chisq.test(table(data[[v]], cluster), correct = FALSE)$p.value
      )
    }
  })

  data.frame(variable = variables, p_value = as.numeric(p))
}

mortality_test <- function(cluster, mortality) {
  suppressWarnings(
    chisq.test(table(cluster, mortality), correct = FALSE)$p.value
  )
}

heatmap_matrix <- function(profile) {
  mat <- t(as.matrix(profile[, setdiff(names(profile), c("cluster", "n")), drop = FALSE]))
  colnames(mat) <- paste0("Cluster ", profile$cluster)

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

profile_30d <- cluster_profile(
  d$short,
  vars_30d,
  d$short$cluster_30d,
  d$short$mortality_30d
)

profile_1y <- cluster_profile(
  d$long,
  vars_1y,
  d$long$cluster_1y,
  d$long$event_1y
)

profile_5y <- cluster_profile(
  d$long,
  vars_5y,
  d$long$cluster_5y,
  d$long$event_5y
)

tests_30d <- bind_rows(
  cluster_tests(
    d$short,
    vars_30d,
    d$short$cluster_30d
  ),
  data.frame(
    variable = "mortality_30d",
    p_value = mortality_test(d$short$cluster_30d, d$short$mortality_30d)
  )
)

tests_1y <- bind_rows(
  cluster_tests(
    d$long,
    vars_1y,
    d$long$cluster_1y
  ),
  data.frame(
    variable = "mortality_1y",
    p_value = mortality_test(d$long$cluster_1y, d$long$event_1y)
  )
)

tests_5y <- bind_rows(
  cluster_tests(
    d$long,
    vars_5y,
    d$long$cluster_5y
  ),
  data.frame(
    variable = "mortality_5y",
    p_value = mortality_test(d$long$cluster_5y, d$long$event_5y)
  )
)

heat_30d <- heatmap_matrix(profile_30d)
heat_1y <- heatmap_matrix(profile_1y)
heat_5y <- heatmap_matrix(profile_5y)

dir.create("outputs", showWarnings = FALSE)

write.csv(cluster_30d$elbow, "outputs/elbow_30d.csv", row.names = FALSE)
write.csv(cluster_1y$elbow, "outputs/elbow_1y.csv", row.names = FALSE)
write.csv(cluster_5y$elbow, "outputs/elbow_5y.csv", row.names = FALSE)

pdf("outputs/elbow_kprototypes.pdf", width = 9, height = 3)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
plot(cluster_30d$elbow$k, cluster_30d$elbow$total_withinss,
     type = "b", pch = 19, xlab = "k", ylab = "Total within-cluster variation",
     main = "30-day")
plot(cluster_1y$elbow$k, cluster_1y$elbow$total_withinss,
     type = "b", pch = 19, xlab = "k", ylab = "Total within-cluster variation",
     main = "1-year")
plot(cluster_5y$elbow$k, cluster_5y$elbow$total_withinss,
     type = "b", pch = 19, xlab = "k", ylab = "Total within-cluster variation",
     main = "5-year")
dev.off()

write.csv(profile_30d, "outputs/cluster_profile_30d.csv", row.names = FALSE)
write.csv(profile_1y, "outputs/cluster_profile_1y.csv", row.names = FALSE)
write.csv(profile_5y, "outputs/cluster_profile_5y.csv", row.names = FALSE)

write.csv(tests_30d, "outputs/cluster_tests_30d.csv", row.names = FALSE)
write.csv(tests_1y, "outputs/cluster_tests_1y.csv", row.names = FALSE)
write.csv(tests_5y, "outputs/cluster_tests_5y.csv", row.names = FALSE)

palette <- colorRampPalette(c("green", "yellow", "red"))(100)

pdf("outputs/cluster_heatmaps_fig5.pdf", width = 7, height = 10)
pheatmap(
  heat_30d,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = palette,
  main = "30-day"
)
pheatmap(
  heat_1y,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = palette,
  main = "1-year"
)
pheatmap(
  heat_5y,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = palette,
  main = "5-year"
)
dev.off()

print(profile_30d)
print(profile_1y)
print(profile_5y)
