suppressPackageStartupMessages({
  # Bioconductor packages
  library(structToolbox)
  library(pmp)
  library(ropls)
  library(mixOmics)
  
  # CRAN libraries
  library(ggpubr)
  library(patchwork)
  library(openxlsx)
  library(struct)
  library(rstatix)
  library(reshape2)
  library(caret)
  library(kableExtra)
  library(viridis)
  library(tidyr)
  library(tidyverse)
  library(dplyr)
  library(rsample)
  library(janitor)
})

## ---------------------------------------------------------------
# prepare data
data <- read.csv("data.csv", row.names = 1, check.names = TRUE)
sample_meta <- read.csv("sample_meta.csv", row.names = 1, check.names = TRUE)
var_meta <- read.csv("var_meta.csv", row.names = 1, check.names = TRUE)
var_meta <- clean_names(var_meta)

data[data == 0] <- NA

# update sample meta
common <- intersect(rownames(var_meta), rownames(data))
data <- data[common, ]
var_meta <- var_meta[common, ]

data <- as.data.frame(t(data)) %>% clean_names()
data <- as.data.frame(t(data))
var_meta <- as.data.frame(t(var_meta)) %>% clean_names()
var_meta <- as.data.frame(t(var_meta))

# create DatasetExperiment
lc_DE <- DatasetExperiment(
  data = as.data.frame(t(data)),
  sample_meta = sample_meta,
  variable_meta = as.data.frame(var_meta, row.names = colnames(as.data.frame(t(
    data
  )))),
  description = 'Mandible study',
  name = "LC-MS"
)

lc_DE

# convert to factors
lc_DE$sample_meta$Month = factor(lc_DE$sample_meta$Month)
lc_DE$sample_meta$Pre_Post = factor(lc_DE$sample_meta$Pre_Post)
lc_DE$sample_meta$Depth = factor(lc_DE$sample_meta$Depth)
lc_DE$sample_meta$Run_order = factor(lc_DE$sample_meta$Run_order)
lc_DE$sample_meta$Type = factor(lc_DE$sample_meta$Type)

lc_DE

## ---------------------------------------------------------------
# matrix processing
# blank filter
blk_filter <- blank_filter(
  fold_change = 20,
  blank_label = "Blank",
  qc_label = "QC",
  factor_name = 'Type'
)

lc_blk <- model_apply(blk_filter, lc_DE)

lc_blk <- predicted(lc_blk)
lc_blk

nc = ncol(lc_DE) - ncol(lc_blk)
cat(paste0('Number of features removed: ', nc))

## ---------------------------------------------------------------
# missing feature and % value
perc_features <-
  mv_feature_filter(threshold = 90,
                    method = 'across',
                    factor_name = 'Type')

lc_blk_perc <- model_apply(perc_features, lc_blk)

lc_blk_perc <- predicted(lc_blk_perc)
lc_blk_perc

nc = ncol(lc_blk) - ncol(lc_blk_perc)
cat(paste0('Number of features removed: ', nc))

## ---------------------------------------------------------------
# rsd qc by feature
qc_features <-
  rsd_filter(rsd_threshold = 10,
             qc_label = 'QC',
             factor_name = 'Type')

lc_blk_perc_qc <- model_apply(qc_features , lc_blk_perc)

lc_blk_perc_qc <- predicted(lc_blk_perc_qc)
lc_blk_perc_qc

nc = ncol(lc_blk_perc) - ncol(lc_blk_perc_qc)
cat(paste0('Number of features removed: ', nc))

## ---------------------------------------------------------------
# prepare the model sequence
process <- knn_impute(neighbours = 5, sample_max = 100) +
  pqn_norm(qc_label = 'QC', factor_name = 'Type') +
  glog_transform(qc_label = 'QC', factor_name = 'Type')

lc_pr <- model_apply(process, lc_blk_perc_qc)
# get the transformed, scaled and imputed matrix
lc_pr <- predicted(lc_pr)

## ---------------------------------------------------------------
# PCA
PCA <- mean_centre() + PCA(number_components = 3)
lc_pr_PCA <- model_apply(PCA, lc_pr)

# scores plot - the QC cluster correctly around the zero cohordinates
score <- pca_scores_plot(factor_name = 'Type', points_to_label = "outliers")

chart_plot(score, lc_pr_PCA[2]) + theme_bw()

## ---------------------------------------------------------------
# Subset to only keep sample and remove QC for further processing
TT = filter_smeta(mode = 'include',
                  factor_name = 'Type',
                  levels = c('sample'))
# +
#   filter_smeta(mode = 'exclude',
#                factor_name = 'Month',
#                levels = c('0'))
# apply model
TT = model_apply(TT, lc_pr)

# keep the data filtered by group for later
lc_pr_filtered = predicted(TT)
lc_pr_filtered

## ---------------------------------------------------------------
# PCA
PCA <- mean_centre() + PCA(number_components = 3)
lc_filtered_PCA <- model_apply(PCA, lc_pr_filtered)

# scores plot - the QC cluster correctly around the zero coordinates
score <- pca_scores_plot(
  factor_name = 'Month',
  points_to_label = "all",
  label_factor = "PMI",
  ellipse = 'none'
)

pca_lc <- chart_plot(score, lc_filtered_PCA[2]) + theme_bw() + 
  viridis::scale_color_viridis(option = "D", discrete = TRUE, name = "PMI")

## ---------------------------------------------------------------
# Univariate statistics (using rstatix package)
df <- lc_pr_filtered$data
meta <- lc_pr_filtered$sample_meta

# add variable of interest
df_melted <- melt(df)
df_melted$PMI <- lc_pr_filtered$sample_meta$Month
df_melted$Depth <- lc_pr_filtered$sample_meta$Depth

df_melted %>%
  group_by(variable) %>%
  kruskal_test(value ~ PMI) %>%
  add_significance() %>%
  write.xlsx("kruskal_test_PMI_lc.xlsx")

df_melted %>%
  group_by(variable) %>%
  dunn_test(value ~ PMI,
            p.adjust.method = "fdr") %>%
  add_significance() %>%
  write.xlsx("dunn_test_PMI_lc.xlsx")

df_melted %>%
  group_by(variable) %>%
  kruskal_test(value ~ Depth) %>%
  add_significance() %>%
  write.xlsx("kruskal_test_Depth_lc.xlsx")

df_melted %>%
  group_by(variable) %>%
  dunn_test(value ~ Depth,
            p.adjust.method = "fdr") %>%
  add_significance() %>%
  write.xlsx("dunn_test_Depth_lc.xlsx")

## ---------------------------------------------------------------
# heatmap
df <- lc_pr_filtered$data
df$month <- lc_pr_filtered$sample_meta$Month
df$depth <- lc_pr_filtered$sample_meta$Depth

hm <- df[1:147]
hm_t <- t(hm) 
annotation_col = data.frame(df[,148:149]) 
colnames(annotation_col) <- c("Month","Depth (cm)")
row.names(annotation_col) <- colnames(hm_t)

heat_lc <- pheatmap::pheatmap(
  hm_t,
  annotation_col = annotation_col,
  cellheight = 7,
  border_color = "grey60",
  # scale = "both",
  annotation_legend = TRUE,
  main = "",
  show_rownames = FALSE,
  angle_col =  "45"
) 

ggsave('heat_lc.pdf', width = 22.5, height = 16, plot = heat_lc)

## ---------------------------------------------------------------
# Machine learning - prepare data
df <- lc_pr_filtered$data
df$PMI <- as.numeric(lc_pr_filtered$sample_meta$PMI)

set.seed(123)
split <- df  %>%
  initial_split(prop = 0.70, PMI)

set.seed(123)
train <- training(split)
test <- testing(split)

## ---------------------------------------------------------------
# PLS model
set.seed(123)
#train model
pls_model <- train(
  PMI ~ .,
  data = train,
  method = "pls",
  tuneLength = 20,
  trControl = trainControl("repeatedcv", number = 5, repeats = 2),
  preProc = c("center", "scale")
)

pls_model

# Make predictions
pls_model_predictions_test <- pls_model %>% predict(test)
# error in test set
Metrics::rmse(test$PMI, pls_model_predictions_test)
Metrics::mae(test$PMI, pls_model_predictions_test)

data_gr_pls <- train %>%
  mutate(set = "train") %>%
  bind_rows(test %>% mutate(set = "test"))

data_gr_pls$fit <- predict(pls_model, data_gr_pls)

plsr_lc <- ggplot(data = data_gr_pls, mapping = aes(y = fit, x = PMI)) +
  geom_point(aes(colour = set, shape = set), alpha = 0.7, size = 4) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    aes(colour = set),
    alpha = 0.9,
    linetype = 2.5,
  ) +
  theme_bw() +
  labs(
    title = "PLS model",
    subtitle = "MAE = 0.49 ," ~ R ^ 2 ~ "= 0.82",
    color = "Sample set",
    shape = "Sample set"
  ) +
  xlab("PMI (month)") + ylab("Estimated PMI (month)") + 
  scale_colour_manual(values = c("#443a83", "#73d056")) 

# normalised variable importance
pls_vip <- vip::vi(pls_model) 
var_meta <- lc_pr_filtered$variable_meta

# extract variable metadata
r <- rownames(var_meta) %in% pls_vip$Variable
var_meta[r,] -> vip_bone_mandible_lcms_plsr_var_meta
# add importance
vip_bone_mandible_lcms_plsr_var_meta$Importance <- pls_vip$Importance

# save VIP
write.xlsx(vip_bone_mandible_lcms_plsr_var_meta,"vip_bone_mandible_lcms_plsr_var_meta.xlsx")


# box plots ---------------------------------------------------------------
df %>%
  dplyr::select(PMI,'x3_46_509_2715m_z', 'x7_74_301_2160m_z' ) %>%
  gather(Measure, Value,-PMI) %>%
  ggplot(aes(x=factor(PMI), y=Value,color=factor(PMI), fill=factor(PMI))) + 
  geom_boxplot(width=0.4,lwd=.7)   + facet_wrap(~Measure, scales = "free_y")+ 
  geom_jitter(aes(color = factor(PMI)), height = 0, width = .2, alpha= .9)+
  theme_bw()+ 
  viridis::scale_color_viridis(option = "D", discrete = TRUE)+
  viridis::scale_fill_viridis(option = "D", discrete = TRUE, alpha = 0.4)+
  xlab('PMI (month)') + ylab('Normalised intensity') + theme(legend.position = 'none')

ggsave('boxplot_lc.pdf', width = 7, height = 3)

# combine plots -----------------------------------------------------------
(pca_gc | pca_lc) + 
  plot_annotation(tag_levels = 'A')

ggsave('pca_tot.pdf', width = 10, height = 4)

(plsr_gc | plsr_lc) + 
  plot_annotation(tag_levels = 'A')

ggsave('plsr_tot.pdf', width = 10, height = 4)

