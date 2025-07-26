# ------------------------------------------------------------------------------
# Module: 5_clustering_kproto.R
# Purpose: Unsupervised clustering of AMI patients using AIC-selected predictors
# Method: k-prototypes algorithm for mixed-type data
# ------------------------------------------------------------------------------

# Description:
# This module implements unsupervised clustering based on variables selected
# via AIC in the supervised modeling pipeline. The goal is to identify distinct
# patient subgroups among AMI cohorts using both continuous and categorical features

library(clustMixType)   # for k-prototypes algorithm
library(cluster)        # for silhouette
library(factoextra)     # for visualization
library(dplyr)

# ---------------------------- Load and Prepare Data ----------------------------
# Assuming preprocessed and cleaned data from previous modules
# Replace below with actual balanced dataset if needed

data_cluster <- train_data_balanced %>% 
  select(Age_C, Time, Residential_area, DM, HTN)

# Ensure correct data types
cat_vars <- c("Time", "Residential_area", "DM", "HTN")
data_cluster[cat_vars] <- lapply(data_cluster[cat_vars], as.factor)

# ------------------------------ Clustering -------------------------------------
# Determine optimal number of clusters using silhouette
sil_width <- numeric()
for (k in 2:6) {
  set.seed(123)
  kpres <- kproto(data_cluster, k)
  clust <- kpres$cluster
  sil <- silhouette(clust, dist(data_cluster, method = "euclidean"))
  sil_width[k] <- mean(sil[, 3])
}

optimal_k <- which.max(sil_width)
cat("Optimal number of clusters:", optimal_k, "\n")

# Final clustering using optimal k
set.seed(123)
kproto_fit <- kproto(data_cluster, optimal_k)
data_cluster$cluster <- as.factor(kproto_fit$cluster)

# ------------------------------ Output -----------------------------------------
# Cluster sizes
print(table(data_cluster$cluster))

# Visualize clustering result (PCA-based for mixed data)
fviz_cluster(list(data = data_cluster, cluster = kproto_fit$cluster),
             ellipse.type = "convex",
             geom = "point",
             main = "K-prototypes Clustering of AMI Patients")

# Optional: Save cluster assignments for further analysis
# write.csv(data_cluster, "clustered_patients.csv", row.names = FALSE)
