suppressPackageStartupMessages({
  # Bioconductor packages
  library(structToolbox)
  library(erah)
  library(mixOmics)
  
  # CRAN libraries
  library(ggpubr)
  library(ggprism)
  library(patchwork)
  library(rstatix)
  library(caret)
  library(tibble)
  library(dplyr)
  library(reshape2)
  library(rsample)
  library(kableExtra)
  library(openxlsx)
  library(RColorBrewer)
})

# load files --------------------------------------------------------------
# object created to reduce run time 
load("/Users/andreabonicelli/Documents/GitHub/metabolomics-mandible/data/mandible_al.rda")
load("/Users/andreabonicelli/Documents/GitHub/metabolomics-mandible/data/mandible_dec.rda")
load("/Users/andreabonicelli/Documents/GitHub/metabolomics-mandible/data/mandible_exp.rda")
load("/Users/andreabonicelli/Documents/GitHub/metabolomics-mandible/data/mandible_id.rda")
load("/Users/andreabonicelli/Documents/GitHub/metabolomics-mandible/data/golmdatabase.rda")

# plan(future::multisession, workers = 7)

instrumental <- read.csv('/Users/andreabonicelli/Documents/GitHub/metabolomics-mandible/data/inst.csv')
phenotype <- read.csv('/Users/andreabonicelli/Documents/GitHub/metabolomics-mandible/data/pheno.csv')

# # create object
# ex <- newExp(instrumental = instrumental,
#              phenotype = phenotype,
#              info = 'GC Mandible Experiment')
# 
# save(ex, file = "mandible_exp.rda")
# 
# plotChr(ex, N.sample = 5, type = "BIC", xlim=c(15,16))
# 
# # perform deconvolution
# ex.dec.par <- setDecPar(min.peak.width = 2,
#                         avoid.processing.mz = c(35:69,73:75,147:149),
#                         noise.threshold = 500)
# 
# ex_dec <- deconvolveComp(ex, ex.dec.par)
# # save deconvolution
# save(ex_dec, file = "mandible_dec.rda")
# 
# # perform alignment
# ex.al.par <- setAlPar(min.spectra.cor = 0.1, max.time.dist = 3)
# ex_al <- alignComp(ex_dec, alParameters = ex.al.par)
# # save alignment
# save(ex_al, file = "mandible_al.rda")
# 
# # load metabolite library 
# load("golmdatabase.rda")

mslib <- golm.database
# perform identification add "id.database = mslib" for GOLM
ex_id <- identifyComp(ex_al, n.putative = 1, id.database = mslib)
# save processed data
# save(ex_id, file = "mandible_id.rda")

# extract data
id.list <- idList(ex_id)
align.list <- alignList(ex_id)
data.list <- dataList(ex_id)

# plot identified features
plotSpectra(ex_id, 53, draw.color = "red") 

# prepare identification matrix -------------------------------------------
id.list$AlignID <- sub("^", "ft", id.list$AlignID )

# prepare intensity matrix
align.list.proc <- align.list[,-2:-4]
align.list.proc$AlignID <- sub("^", "ft", align.list.proc$AlignID )
align.list.proc[align.list.proc == 0] <- NA
align.list.proc <- align.list.proc %>% remove_rownames %>% 
  column_to_rownames(var="AlignID")

# create DatasetExperiment ------------------------------------------------
gc_DE <- DatasetExperiment(
  data = as.data.frame(t(align.list.proc)),
  sample_meta = phenotype,
  variable_meta = as.data.frame(id.list,
                                row.names = colnames(align.list.proc)),
  description = 'Mandible study',
  name = "GC-MS"
)

# convert to factors
gc_DE$sample_meta$class = factor(gc_DE$sample_meta$class)
gc_DE$sample_meta$type = factor(gc_DE$sample_meta$type)
gc_DE$sample_meta$depth = factor(gc_DE$sample_meta$depth)
gc_DE$sample_meta$run_order = as.numeric(gc_DE$sample_meta$run_order)
gc_DE

# matrix processing -------------------------------------------------------
# blank filter
blk_filter <- blank_filter(
  fold_change = 25,
  blank_label = "Blank",
  qc_label = "QC",
  factor_name = 'type'
)

gc_blk <- model_apply(blk_filter, gc_DE)

blank_plot = blank_filter_hist()
A <- chart_plot(blank_plot, gc_blk) + theme_bw(14)

gc_blk <- predicted(gc_blk)
gc_blk

nc = ncol(gc_DE) - ncol(gc_blk)
cat(paste0('Number of features removed: ', nc))

# missing feature and % value - remove eveery feature with more than 10% missing 
# data
perc_features <-
  mv_feature_filter(threshold = 80,
                    method = 'across',
                    factor_name = 'type')

gc_blk_perc <- model_apply(perc_features, gc_blk)

mv_plot = mv_feature_filter_hist()
B <- chart_plot(mv_plot, gc_blk_perc) + theme_bw(14)

gc_blk_perc <- predicted(gc_blk_perc)
gc_blk_perc

nc = ncol(gc_blk) - ncol(gc_blk_perc)
cat(paste0('Number of features removed: ', nc))

# removing samples based on missingness
perc_sample <-
  mv_sample_filter(mv_threshold = 50)

gc_blk_perc_samp <- model_apply(perc_sample, gc_blk_perc)

mv_plot_samp = mv_sample_filter_hist()
C <- chart_plot(mv_plot_samp, gc_blk_perc_samp) + theme_bw(14)

gc_blk_perc_samp <- predicted(gc_blk_perc_samp)
gc_blk_perc_samp

nc = ncol(gc_blk) - ncol(gc_blk_perc_samp)
cat(paste0('Number of samples removed: ', nc))

# rsd QC by feature 
qc_features <-
  rsd_filter(rsd_threshold = 30,
             qc_label = 'QC',
             factor_name = 'type')

gc_blk_perc_qc <- model_apply(qc_features , gc_blk_perc_samp)

rsd_plot = rsd_filter_hist()
D <- chart_plot(rsd_plot, gc_blk_perc_qc) + theme_bw(14)

gc_blk_perc_qc <- predicted(gc_blk_perc_qc)
gc_blk_perc_qc

nc = ncol(gc_blk_perc) - ncol(gc_blk_perc_qc)
cat(paste0('Number of features removed: ', nc))

((A+B) / (C+D))

# perform drift correction  -----------------------------------------------
M = # batch correction
  sb_corr(
    order_col='run_order',
    batch_col='Batch', 
    qc_col='type', 
    qc_label='QC'
  )

M = model_apply(M,gc_blk_perc_qc)

C = tic_chart(factor_name='type',run_order='run_order')

A = chart_plot(C,gc_blk_perc_qc) + theme_bw(14) + ggtitle('Original')
B = chart_plot(C,predicted(M)) + theme_bw(14) + ggtitle('Adjusted')

(A | B) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect") &
  theme(legend.position = 'bottom')

M = predicted(M)



# data normalisation and imputation ---------------------------------------
# prepare the model sequence
process <- knn_impute(neighbours = 5, sample_max = 50) +
  pqn_norm(qc_label = 'QC', factor_name = 'type') +
  glog_transform(qc_label='QC', factor_name='type')

gc_pr <- model_apply(process, M)

# get the transformed, scaled and imputed matrix
gc_pr <- predicted(gc_pr)
gc_pr

nrow(gc_pr$variable_meta)

box_sample <- DatasetExperiment_boxplot(
  factor_name = 'type',
  number = 20,
  by_sample = TRUE,
  per_class = TRUE
)

chart_plot(box_sample, gc_pr) + theme_bw(14)

# PCA - check data quality ------------------------------------------------
PCA <- autoscale() + PCA(number_components = 3)
gc_pr_PCA <- model_apply(PCA, gc_pr)

# scores plot - the QC cluster correctly around the zero coordinates
score <- pca_scores_plot(
  factor_name = 'type',
  points_to_label = 'outliers'
)
 
# prepare and print the plot
chart_plot(score, gc_pr_PCA[2]) + 
  theme_bw(16)

# subset object -----------------------------------------------------------
# Subset to only keep sample and remove QC for further processing and
# we remove those few outliers shown in the score plot above
TT = filter_smeta(mode = 'include',
                  factor_name = 'type',
                  levels = c('sample')) 
# apply model
TT = model_apply(TT, gc_pr)

# keep the data filtered by group for later
gc_pr_filtered = predicted(TT)
gc_pr_filtered

gc_pr_filtered$sample_meta$depth = factor(gc_pr_filtered$sample_meta$depth)

# remove duplicated features with lower match score -----------------------
M = filter_by_name(mode='exclude',dimension='variable',
                   names=c("ft160","ft622","ft254"))

# apply model
gc_pr_filtered = model_apply(M, gc_pr_filtered)

# keep the data filtered by group for later
gc_pr_filtered = predicted(gc_pr_filtered)
gc_pr_filtered

gc_pr_filtered_PCA <- model_apply(PCA, gc_pr_filtered)

# scores plot - the QC cluster correctly around the zero coordinates
score <- pca_scores_plot(
  factor_name = 'month',
  points_to_label = "all",
  label_factor = "pmi",
  ellipse = 'none'
)
# prepare and print the plot
pca_gc <- chart_plot(score, gc_pr_filtered_PCA[2]) + 
  theme_bw() + viridis::scale_color_viridis(option = "A",
                                              discrete = TRUE,
                                              name = "PMI") 

# Univariate statistics (using rstatix package) ---------------------------
df <- gc_pr_filtered$data
meta <- gc_pr_filtered$sample_meta

# add variable of interest
df_melted <- melt(df)
df_melted$PMI <- gc_pr_filtered$sample_meta$pmi
df_melted$Depth <- gc_pr_filtered$sample_meta$depth

df_melted %>%
  group_by(variable) %>%
  kruskal_test(value ~ PMI) %>%
  add_significance() %>%
  write.xlsx("/Users/andreabonicelli/Documents/GitHub/metabolomics-mandible/tables/kruskal_test_PMI_gcms.xlsx")

df_melted %>%
  group_by(variable) %>%
  dunn_test(value ~ PMI,
            p.adjust.method = "fdr") %>%
  add_significance() %>%
  write.xlsx("/Users/andreabonicelli/Documents/GitHub/metabolomics-mandible/tables/dunn_test_PMI_gcms.xlsx")

df_melted %>%
  group_by(variable) %>%
  kruskal_test(value ~ Depth) %>%
  add_significance() %>%
  write.xlsx("/Users/andreabonicelli/Documents/GitHub/metabolomics-mandible/tables/kruskal_test_Depth_gcms.xlsx")

df_melted %>%
  group_by(variable) %>%
  dunn_test(value ~ Depth,
            p.adjust.method = "fdr") %>%
  add_significance() %>%
  write.xlsx("/Users/andreabonicelli/Documents/GitHub/metabolomics-mandible/tables/dunn_test_Depth_gcms.xlsx")

# heatmap -----------------------------------------------------------------
df <- gc_pr_filtered$data
df$month <- gc_pr_filtered$sample_meta$month
df$depth <- gc_pr_filtered$sample_meta$depth

hm <- df[1:48]
hm_t <- t(hm) 
annotation_col = data.frame(df[,49:50]) 
colnames(annotation_col) <- c("PMI","Depth (cm)")
row.names(annotation_col) <- colnames(hm_t)

heat_gc = pheatmap::pheatmap(
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

ggsave('/Users/andreabonicelli/Documents/GitHub/metabolomics-mandible/figures/heat_gcms.pdf', width = 12, height = 8, plot= heat_gc)

# data split --------------------------------------------------------------
df <- gc_pr_filtered$data
df$PMI <- as.numeric(gc_pr_filtered$sample_meta$pmi)

set.seed(123)
split <- df  %>%
  initial_split(prop = 0.75, PMI)

set.seed(123)
train <- training(split)
test <- testing(split)

# PLS regression model ----------------------------------------------------
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

# Define the residual plot function
residual_plot <- function(model, data, response) {
  
  # Check if model is trained with caret
  if (!inherits(model, "train")) {
    stop("The model should be a 'train' object from the caret package.")
  }
  
  # Predict values
  predicted_values <- predict(model, data)
  
  # Extract actual values
  actual_values <- data[[response]]
  
  # Compute residuals
  residuals <- actual_values - predicted_values
  
  # Create a data frame for plotting
  residual_data <- data.frame(Predicted = predicted_values, Residuals = residuals)
  
  # Plot residuals
  ggplot(residual_data, aes(x = Predicted, y = Residuals)) +
    geom_point(alpha = 0.5, color = "blue") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
    labs(title = "Residual Plot",
         x = "Predicted Values",
         y = "Residuals") + theme_bw(12)
}

# Generate the residual plot
resid_gc <- residual_plot(pls_model, test, "PMI") + ggtitle('GC-MS model')

# Make predictions
pls_model_predictions_test <- pls_model %>% predict(test)
# error in test set
Metrics::rmse(test$PMI, pls_model_predictions_test)
Metrics::mae(test$PMI, pls_model_predictions_test)

# plot predicted vs actual values
data_gr_pls <- train %>%
  mutate(set = "train") %>%
  bind_rows(test %>% mutate(set = "test"))

data_gr_pls$fit <- predict(pls_model, data_gr_pls)

plsr_gc <- ggplot(data = data_gr_pls, mapping = aes(y = fit, x = PMI)) +
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
    subtitle = "MAE = 0.82," ~ R ^ 2 ~ "= 0.78",
    color = "Sample set",
    shape = "Sample set"
  ) +
  xlab("PMI (month)") + ylab("Estimated PMI (month)") + 
  scale_colour_manual(values = c("#56146DFF", "#F4685CFF")) 

plsr_gc

vip::vip(
  pls_model,
  num_features = 10,
  geom = "point",
  aesthetics = list(color = "#F4685CFF", size = 3)
) +
  theme_bw(14)

# normalised variable importance
pls_vip <- vip::vi(pls_model) 
var_meta <- gc_pr_filtered$variable_meta

# extract variable metadata
r <- rownames(var_meta) %in% pls_vip$Variable
var_meta[r,] -> vip_bone_mandible_gcms_plsr_var_meta
# add importance
vip_bone_mandible_gcms_plsr_var_meta$Importance <- pls_vip$Importance

# save VIP
write.xlsx(vip_bone_mandible_gcms_plsr_var_meta,
           "/Users/andreabonicelli/Documents/GitHub/metabolomics-mandible/tables/vip_bone_mandible_plsr_gcms.xlsx")

# box plots ---------------------------------------------------------------
df %>%
  dplyr::select(PMI, ft159, ft403, ft624, ft625 ) %>%
  gather(Measure, Value,-PMI) %>%
  ggplot(aes(x=factor(PMI), y=Value,color=factor(PMI), fill=factor(PMI))) + 
  geom_boxplot(width=0.4,lwd=.7)   +facet_wrap(~Measure, scales = "free_y")+ 
  geom_jitter(aes(color = factor(PMI)), height = 0, width = .2, alpha= .9)+
  theme_bw()+ 
  viridis::scale_color_viridis(option = "A", discrete = TRUE)+
  viridis::scale_fill_viridis(option = "A", discrete = TRUE, alpha = 0.4)+
  xlab('PMI (month)') + ylab('Normalised intensity') +theme(legend.position = 'none')

ggsave('/Users/andreabonicelli/Documents/GitHub/metabolomics-mandible/figures/boxplot_gcms.pdf', width = 7, height = 3)

