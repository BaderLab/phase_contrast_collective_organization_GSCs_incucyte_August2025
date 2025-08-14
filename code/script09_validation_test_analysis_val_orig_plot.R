### PACKAGES TO LOAD
###------------------------------------------------------------------#### 
library(ggplot2)  
library(colorspace)
library(tidyr)
library(dplyr)
library(ggthemes)
library(ggpubr)
library(ggrepel)
library(effectsize)
library(ggthemes)
library(scales)
library(forcats)

### FILE/FOLDER PATHS
###------------------------------------------------------------------#### 
path_to_repo <- '/Users/mystique27m/Documents/Professional/research/PostdoctoralResearch_2020/Projects/PatternRecognitionInGSCs_UsingCV/'
script_path <- paste0(path_to_repo,'scripts_final/code/')
out <- paste0(path_to_repo, 'out/')
save_dir_path <- paste0(out, 'script09_output_files/')
dir.create(save_dir_path)

figures_filepath <- paste0(out, '/figures_final/')

#tables_path <-  paste0(out,'/tables/')




### source data 
###------------------------------------------------------------------#### 
#source(paste0(script_path,"/script04b_.R"))
source(paste0(script_path,"source_scripts/sourceData_General_variables_and_functions.R"))
GSC.gsva_orig <- read.csv(paste0(path_to_repo,'/datasets_final/manuscript_analysis_data/bulk_RNA/GSC.gsva.csv'), row.names=1)
GSC.gsva_new <- read.csv(paste0(path_to_repo,'/incucyte_validation/out/GSC.gsva_including_validation_cohort.csv'), row.names=1)
load(paste0('/Users/mystique27m/Documents/Professional/research/PostdoctoralResearch_2020/Projects/PatternRecognitionInGSCs_UsingCV/out/script08_output_files/feature_correlations.Rdata'))



### predict using top loading 
## the common feature that can score for neurodevelopmental InfoMeas1
## the common feature that can score for injury is DifferenceVariance

#new_df_list <- readRDS(paste0(out, 'normalized_validation.Rds'))
sample_names <- unique(unlist(lapply(pca_list, function(x) x$Sample)))
colors_to_map <- c(cols_vivid, col_bold)
names(colors_to_map) <- sample_names


#glist <- list()

#for (i in 1:length(pca_list)){  
#  glist[[i]] <- pca_list[[i]] %>% 
#    ggplot(aes(x=reorder(Sample, InfoMeas1_00), y=InfoMeas1_00, fill=Sample)) + 
#    geom_boxplot() + 
#    geom_point(size=0.1) +
#    theme_minimal() + 
#    scale_fill_manual(values=colors_to_map) +
#    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=21),
#          axis.text.y = element_text(size=21),
#          axis.title.x = element_blank(),
#          axis.title.y = element_text(size=24), 
#          legend.position='none') + 
#    labs(x='Sample', y='Informational Measure 1') #+
#    #ylim(0,1)
#}


###------------------------------------------------------------------####
### FIGURE S2A
###------------------------------------------------------------------####
#dir.create(paste0(figures_filepath, '/Sup_Fig_S2_repeat/'))

#pdf(file = paste0(figures_filepath,'Sup_Fig_S2_repeat/S2_C1_infoMeas1_orig.pdf'), width=45, height=6)
#ggarrange(plotlist=glist, ncol=9, nrow=1)
#dev.off()


#glist <- list()

#for (i in 1:length(pca_list)){  
#  glist[[i]] <- pca_list[[i]] %>% 
#    ggplot(aes(x=reorder(Sample, PC2), y=InfoMeas1_00, fill=Sample)) + 
#    geom_boxplot() + 
#    geom_point(size=0.1) +
#    theme_minimal() + 
#    scale_fill_manual(values=colors_to_map) +
#    theme(axis.text.x = element_text(angle = 90, hjust = 1, size=21),
#          axis.text.y = element_text(size=21),
#          axis.title.x = element_blank(),
#          axis.title.y = element_text(size=24), 
#          legend.position='none') + 
#    labs(x='Sample', y='Informational Measure 1') #+
#  #ylim(0,1)
#}


###------------------------------------------------------------------####
### FIGURE S2A
###------------------------------------------------------------------####
#dir.create(paste0(figures_filepath, '/Sup_Fig_S2_repeat/'))

#pdf(file = paste0(figures_filepath,'Sup_Fig_S2_repeat/S2_C1_infoMeas1_pc2_orig.pdf'), width=48, height=6)
#ggarrange(plotlist=glist, ncol=9, nrow=1)
#dev.off()


################################################################################
### Regression model by each confluency group
################################################################################

library(tidyverse)
library(caret)
library(glmnet)
library(rsample)
library(ggplot2)
library(ggbiplot)

# Set seed for reproducibility
set.seed(42)

# Transpose GSVA and extract PC1 axis
gsc.gsva_t <- t(GSC.gsva_orig)
gsva_pca <- prcomp(gsc.gsva_t, center = TRUE, scale. = TRUE)

# Compute proportion of variance explained
pve <- (gsva_pca$sdev)^2 / sum(gsva_pca$sdev^2)
pve_df <- data.frame(
  PC = paste0("PC", 1:length(pve)),
  Variance = pve * 100  # in percentage
)

#g <- ggbiplot(gsva_pca,
#         labels.size = 9,
#         var.factor = 1,
#         var.scale = 1,
#         varname.size = 5.5,
#         varname.adjust = 1,
#         ellipse=T,
#         ellipse.linewidth = 0.1,
#         circle=T,
#         )+
#  theme_minimal(base_size = 9) +
#  theme(legend.direction = 'horizontal', legend.position = 'top')


#pdf(file = paste0(figures_filepath,'Sup_Fig_S2_repeat/S2_C1_gsva_pc1_2.pdf'), width=18, height=18)
#g
#dev.off()


# Ensure correct numeric ordering of PCs
pve_df$PC <- factor(pve_df$PC, levels = paste0("PC", 1:length(pve)))


gsva_axis_score_pc1 <- scale(gsva_pca$x[, 1])[, 1]

# Flip GSVA PC1 and loadings so that neurodevelopmental = high PC1
gsva_axis_score_pc1 <- -gsva_axis_score_pc1
gsva_pca$rotation[, 1] <- -gsva_pca$rotation[, 1]

#library(dplyr)
#library(ggplot2)
#####################################
#### Extract loadings for PC1 and PC2
#####################################

### PC1
# 1. Extract the PC1 loadings into a data.frame
loadings_pc1 <- gsva_pca$rotation[, 1]
loadings_df_pc1 <- data.frame(
  GeneSet = names(loadings_pc1),
  Loading = as.numeric(loadings_pc1),
  row.names = NULL
)

# 2. Select top10 and bottom10 by Loading
top10_pc1    <- loadings_df_pc1 %>% top_n(10,    wt = Loading)
bottom10_pc1 <- loadings_df_pc1 %>% top_n(-10,   wt = Loading)
plot_df_pc1  <- bind_rows(top10_pc1, bottom10_pc1) %>%
  # keep them in order so that the barplot orders correctly
  arrange(desc(Loading))

# 3. Make a barplot (horizontal) for PC1
ggplot(plot_df_pc1, aes(
  x = reorder(GeneSet, Loading), 
  y = Loading,
  fill = Loading > 0
)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Top 10 and Bottom 10 Loadings for PC1",
    x = "Gene Set",
    y = "Loading"
  ) +
  ylab('PC1 Loadings') +
  xlab('Gene signatures') +
  theme(legend.position = "none",
        title = element_text(size=18, face='bold', hjust=.5),
        plot.title = element_text(hjust=1.5),
        axis.text = element_text(size=21))

main_figure_path <- paste0(figures_filepath, 'Figure_04/')
dir.create(main_figure_path)

ggsave(filename='S4_A_training_gene_signature_pc1_gsva_top_loadings.png', path = main_figure_path, width=12, height=12)

### PC2
# 1. Extract the PC2 loadings into a data.frame
loadings_pc2 <- gsva_pca$rotation[, 2]
loadings_df_pc2 <- data.frame(
  GeneSet = names(loadings_pc2),
  Loading = as.numeric(loadings_pc2),
  row.names = NULL
)

# 2. Select top 10 and bottom 10 by Loading
top10_pc2    <- loadings_df_pc2 %>% top_n(10,    wt = Loading)
bottom10_pc2 <- loadings_df_pc2 %>% top_n(-10,   wt = Loading)
plot_df_pc2  <- bind_rows(top10_pc2, bottom10_pc2) %>%
  # keep them in order so that the barplot orders correctly
  arrange(desc(Loading))

# 3. Make a barplot (horizontal) for PC 1
ggplot(plot_df_pc2, aes(
  x = reorder(GeneSet, Loading), 
  y = Loading,
  fill = Loading > 0
)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Top 10 and Bottom 10 Loadings for PC2",
    x = "Gene Set",
    y = "Loading",
  ) +
  ylab('PC1 Loadings') +
  xlab('Gene signatures') +
  theme(legend.position = "none",
        title = element_text(size=18, face='bold', hjust=.5),
        plot.title = element_text(hjust=1.5),
        axis.text = element_text(size=21))
#ggsave(filename='Sup_Fig_S2_repeat/S2_training_gene_signature_pc2_gsva_top_loadings.png', path = figures_filepath, width=12, height=12)



# Create plotting df
df_train_pc1 <- data.frame(
  Sample = rownames(gsc.gsva_t),
  PC1 = gsva_axis_score_pc1
)


ggplot(df_train_pc1, aes(x = reorder(Sample, PC1), y = PC1)) +
  geom_col(fill = "gray50") +
  theme_minimal(base_size = 16) +
  labs(
    title = "Training Samples Ordered by GSVA gradient (PC1)",
    x = "",
    y = "GSVA PC1 (Scaled)"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12))

#ggsave(filename='Sup_Fig_S2_repeat/S2_training_pc1_gsva_gradient.png', path = figures_filepath, width = 9, height=6)


# Extract loadings
loadings_train_pc1 <- gsva_pca$rotation[, 1]

loadings_train_df <- data.frame(
  GeneSet = names(loadings_train_pc1),
  Loading = loadings_train_pc1
)

# Top and bottom drivers
#top_gs_train <- loadings_train_df %>% top_n(10, Loading)
#bottom_gs_train <- loadings_train_df %>% top_n(-10, Loading)
#top_bottom_train <- rbind(top_gs_train, bottom_gs_train)

# Plot
#ggplot(top_bottom_train, aes(x = reorder(GeneSet, Loading), y = Loading, fill = Loading > 0)) +
#  geom_col() +
#  coord_flip() +
#  theme_minimal(base_size = 12) +
#  labs(
#    title = "Top GSVA pathways driving PC1 (training gradient)",
#    x = "Gene Set",
#    y = "PC1 Loading"
#  ) +
#  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato")) +
#  theme(legend.position = "none",
#        axis.text.y = element_text(size=12))

#ggsave(filename='Sup_Fig_S2_repeat/S2_training_pc2_gsva_top_loadings.png', path = figures_filepath, width=9, height=9)


# Summarize image features
features <- colnames(pca_list[[1]])[36:64]
image_feature_list <- lapply(pca_list, function(df) {
  df <- df[, c("Sample", features)]
  df %>%
    group_by(Sample) %>%
    summarise(across(all_of(features), mean, na.rm = TRUE)) %>%
    column_to_rownames("Sample")
})


# Count unique samples per confluency group BEFORE summarization
raw_sample_counts <- sapply(pca_list, function(df) {
  length(unique(df$Sample))
})

names(raw_sample_counts) <- paste0("Confluency ", seq_along(pca_list))
raw_sample_counts


# Function to train + evaluate a regularized model
train_model_with_alpha <- function(df, alpha_val) {
  x <- as.matrix(df[, -1])
  y <- df$index
  
  cv_model <- cv.glmnet(x, y, alpha = alpha_val)
  best_lambda <- cv_model$lambda.min
  model <- glmnet(x, y, alpha = alpha_val, lambda = best_lambda)
  
  preds <- predict(model, newx = x)
  
  list(
    model = model,
    cor = cor(preds, y),
    rmse = RMSE(preds, y),
    r2 = R2(preds, y),
    preds = data.frame(index = y, predicted = as.numeric(preds), Sample = rownames(df))
  )
}

# Store performance
performance_df <- data.frame()
prediction_df <- data.frame()

for (alpha_val in c(0, 0.5, 1)) {
  for (i in seq_along(image_feature_list)) {
    features_df <- image_feature_list[[i]]
    common_samples <- intersect(rownames(features_df), names(gsva_axis_score_pc1))
    df <- data.frame(index = gsva_axis_score_pc1[common_samples], features_df[common_samples, ])
    df <- df[complete.cases(df), ]
    
    out <- train_model_with_alpha(df, alpha_val)
    
    performance_df <- rbind(performance_df, data.frame(
      model_type = paste0("alpha_", alpha_val),
      confluency_group = paste0("Confluency ", i),
      correlation = out$cor,
      rmse = out$rmse,
      r2 = out$r2
    ))
    
    pred_df <- out$preds
    pred_df$model_type <- paste0("alpha_", alpha_val)
    pred_df$con_grp <- paste0("Confluency ", i)
    prediction_df <- rbind(prediction_df, pred_df)
  }
}

# Update factor labels
prediction_df$model_type <- factor(prediction_df$model_type,
                                   levels = c("alpha_0", "alpha_0.5", "alpha_1"),
                                   labels = c("Ridge (α=0)", "Elastic Net (α=0.5)", "Lasso (α=1)"))
performance_df$model_type <- factor(performance_df$model_type,
                                    levels = c("alpha_0", "alpha_0.5", "alpha_1"),
                                    labels = c("Ridge (α=0)", "Elastic Net (α=0.5)", "Lasso (α=1)"))

# --- Bar Plots ---
#ggplot(performance_df, aes(x = confluency_group, y = correlation, fill = model_type)) +
#  geom_bar(stat = "identity", position = position_dodge()) +
#  labs(title = "Performance by Model and Confluency", y = "Correlation", fill = "Model") +
#  theme_minimal(base_size = 13) +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=18),
#        axis.text.y = element_text(size=15),
#        axis.title.x = element_blank(),
#        title = element_text(size=21, face = 'bold'))
#ggsave(filename='Sup_Fig_S2_repeat/S2_training_pc2_gsva_image_pc2_performance_barplot.png', path = figures_filepath, width=, height=6)



# --- Predicted vs Actual ---
#ggplot(prediction_df, aes(x = index, y = predicted, color = model_type)) +
#  geom_point(alpha = 0.7, size = 2) +
#  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black", size = 0.4) +
#  facet_wrap(~ model_type + con_grp, scales = "free_y", nrow=3, ) +
#  labs(
#    x = "Groundtruth Neurodev/Inj. Res gradient",
#    y = "Predicted Neurodev/Inj. Res gradient",
#    title = "Predicted vs Groundtruth by Model and Confluency",
#    color = "Model",
#    size = 21,
#    face = 'bold'
#  ) +
#  theme_minimal(base_size = 15) +
#  theme(strip.text = element_text(size = 12, face = "bold"),
#        axis.text.x = element_text(angle = 45, hjust = 1, size=9),
#        axis.text.y = element_text(size=9),
##        axis.title.x = element_blank(),
#        legend.text = element_text(size=15),
#        legend.key.spacing.y = unit(0.45, "cm"),
#        legend.position = "top",
#        title = element_text(size = 18, face = 'bold')) +
#  scale_color_manual(values = c("Ridge (α=0)" = "#F8766D",
#                                "Elastic Net (α=0.5)" = "#00BA38",
#                                "Lasso (α=1)" = "#619CFF"))

#ggsave(filename='Sup_Fig_S2_repeat/S2_predicted_gradient_training_gradient_gsva.png', path = figures_filepath, width=18, height=9)


### plot R2 and RMSE values 

ggplot(performance_df, aes(x=confluency_group, y=s0, color=model_type, group=model_type)) + 
  labs(
    x = "Confluency groups",
    y = "R2 value",
    title = "Variance explained (R2)",
    color = "Model",
    size = 21,
    face = 'bold'
  ) +
  geom_line(linewidth=3) +
  theme_minimal(base_size = 15) +
  theme(strip.text = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=21, face='bold'),
        axis.text.y = element_text(size=12, face='bold'),
        #        axis.title.x = element_blank(),
        legend.text = element_text(size=27),
        legend.key.spacing.y = unit(0.45, "cm"),
        legend.position = "top",
        title = element_text(size = 21, face = 'bold')) +
  scale_color_manual(values = c("Ridge (α=0)" = "#F8766D",
                                "Elastic Net (α=0.5)" = "#00BA38",
                                "Lasso (α=1)" = "#619CFF"))

#ggsave(filename='Sup_Fig_S2_repeat/S2_r2_variance_explained.png', path = figures_filepath, width=15, height=12)


performance_df_model_ridge <- subset(performance_df, model_type=='Ridge (α=0)')

ggplot(performance_df_model_ridge, aes(x=confluency_group, y=s0, color=model_type, group=1)) + 
  labs(
    x = "Confluency groups",
    y = "R2 value",
    title = "Variance explained (R2)",
    color = "Model",
    size = 27,
    face = 'bold'
  ) +
  geom_line(linewidth=3) +
  theme_minimal(base_size = 15) +
  theme(strip.text = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=21, face='bold'),
        axis.text.y = element_text(size=12, face='bold'),
        #        axis.title.x = element_blank(),
        legend.text = element_text(size=27),
        legend.key.spacing.y = unit(0.45, "cm"),
        legend.position = "top",
        title = element_text(size = 21, face = 'bold')) 

#ggsave(filename='Sup_Fig_S2_repeat/S2_r2_variance_explained_ridge.png', path = figures_filepath, width=15, height=12)



##########################
### END
########################



#########################
### COMBINED - combine image features from confluency group 1-9 
####


# -- Load Packages
# --- Load Required Libraries ---
library(tidyverse)
library(caret)
library(glmnet)
library(rsample)
library(ggplot2)

# --- Set Seed for Reproducibility ---
set.seed(42)

# --- Compute GSVA PC1 Axis ---
#gsc.gsva_t <- t(GSC.gsva_orig)
#gsva_pca <- prcomp(gsc.gsva_t, center = TRUE, scale. = TRUE)
#gsva_axis_score_pc1 <- scale(gsva_pca$x[, 1])[, 1]  # Use PC1
# Flip GSVA PC1 and loadings so that neurodevelopmental = high PC1
#gsva_axis_score_pc1 <- -gsva_axis_score_pc1
#gsva_pca$rotation[, 1] <- -gsva_pca$rotation[, 1]


# --- Extract Image Feature Names ---
features <- colnames(pca_list[[1]])[36:64]

# --- Combine and Aggregate All Samples Across All Confluency Groups ---
combined_features_df <- bind_rows(pca_list) %>%
  select(Sample, all_of(features)) %>%
  group_by(Sample) %>%
  summarise(across(all_of(features), mean, na.rm = TRUE), .groups = 'drop') %>%
  column_to_rownames("Sample")

# --- Find Common Samples with GSVA ---
common_samples <- intersect(rownames(combined_features_df), names(gsva_axis_score_pc1))

# --- Create Data Frame for Modeling ---
df_combined <- data.frame(
  index = gsva_axis_score_pc1[common_samples],
  combined_features_df[common_samples, ]
)
df_combined <- df_combined[complete.cases(df_combined), ]
df_combined <- df_combined %>% rename(PC1=index)


# --- Model Training Function ---
train_model_with_alpha <- function(df, alpha_val) {
  x <- as.matrix(df[, -1])
  y <- df$PC1
  
  cv_model <- cv.glmnet(x, y, alpha = alpha_val)
  best_lambda <- cv_model$lambda.min
  model <- glmnet(x, y, alpha = alpha_val, lambda = best_lambda)
  
  preds <- predict(model, newx = x)
  
  list(
    model = model,
    cor = cor(preds, y),
    rmse = RMSE(preds, y),
    r2 = R2(preds, y),
    preds = data.frame(index = y, predicted = as.numeric(preds), Sample = rownames(df))
  )
  }

library(tidyverse)
library(tidyverse)

# --- Step 0: Rename 'index' to 'PC1' ---
#df_combined <- df_combined %>%
#  rename(PC1 = index)

# --- Step 1: Z-normalize the response only ---
#df_combined$PC1 <- scale(df_combined$PC1)

# --- Step 2: Train Ridge, Elastic Net, and Lasso ---
model_types <- c(0, 0.5, 1)
performance_combined <- data.frame()
prediction_combined <- data.frame()
lasso_model_final <- NULL

for (alpha_val in model_types) {
  out <- train_model_with_alpha(df_combined, alpha_val)
  
  if (alpha_val == 1) {
    lasso_model_final <- out$model
  }
  
  performance_combined <- rbind(performance_combined, data.frame(
    model_type = paste0("alpha_", alpha_val),
    correlation = out$cor,
    rmse = out$rmse,
    r2 = out$r2
  ))
  
  pred_df <- out$preds
  pred_df$model_type <- paste0("alpha_", alpha_val)
  prediction_combined <- rbind(prediction_combined, pred_df)
}

# --- Step 3: Train Linear Regression separately ---
model_lm <- lm(PC1 ~ ., data = df_combined)
preds_lm <- predict(model_lm, newdata = df_combined)
cor_val <- cor(preds_lm, df_combined$PC1)
rmse_val <- sqrt(mean((preds_lm - df_combined$PC1)^2))
r2_val <- summary(model_lm)$r.squared



pred_df_lm <- data.frame(
  index = df_combined$PC1,
  predicted = preds_lm,
  Sample = rownames(df_combined),
  model_type = "alpha_lm"
)

names(performance_combined)
names(data.frame(
  model_type = "alpha_lm",
  correlation = cor_val,
  rmse = rmse_val,
  r2 = r2_val
))
names(prediction_combined)
names(pred_df_lm)

performance_combined <- performance_combined %>% rename(r2=s0)



# --- Step 4: Combine all results ---
performance_combined <- rbind(performance_combined, data.frame(
  model_type = "alpha_lm",
  correlation = cor_val,
  rmse = rmse_val,
  r2 = r2_val

))

#prediction_combined <- rbind(prediction_combined, pred_df_lm)

# --- Step 5: Relabel model names ---
performance_combined$model_type <- factor(performance_combined$model_type,
                                          levels = c("alpha_0", "alpha_0.5", "alpha_1"),
                                          labels = c("Ridge (α=0)", "Elastic Net (α=0.5)", "Lasso (α=1)"))

prediction_combined$model_type <- factor(prediction_combined$model_type,
                                         levels = c("alpha_0", "alpha_0.5", "alpha_1"),
                                         labels = c("Ridge (α=0)", "Elastic Net (α=0.5)", "Lasso (α=1)"))

# --- Step 6: Plot ---
ggplot(prediction_combined, aes(x = index, y = predicted, color = model_type)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed",
              aes(group = model_type, color = model_type)) +
  labs(
    title = "Predicted vs Actual",
    x = "Groundtruth",
    y = "Predicted",
    color = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = 'bold', size = 21),
    axis.text.x = element_text(size = 21),
    axis.text.y = element_text(size=21),
    legend.text = element_text(size = 18),
    legend.key.spacing.y = unit(0.45, "cm"),
    title = element_text(face = 'bold', size = 30)
  ) +
  scale_color_manual(values = c(
    "Ridge (α=0)" = "#F8766D",
    "Elastic Net (α=0.5)" = "#00BA38",
    "Lasso (α=1)" = "#619CFF",
    "Linear Regression" = "#C77CFF"
  ))


pred_ridge <- subset(prediction_combined, subset=prediction_combined$model_type=='Ridge (α=0)')
pred_lasso <- subset(prediction_combined, subset=prediction_combined$model_type=='Lasso (α=1)')
pred_elastic_net <- subset(prediction_combined, subset=prediction_combined$model_type=='Elastic Net (α=0.5)')

ridge_pval <- cor.test(pred_ridge$predicted, pred_ridge$index)$p.value
print(paste0('Ridge p_value ',ridge_pval))

lasso_pval <- cor.test(pred_lasso$predicted, pred_lasso$index)$p.value
print(paste0('Lasso p_value ',round(lasso_pval, digits=5)))

en_pval <- cor.test(pred_elastic_net$predicted, pred_elastic_net$index)$p.value
print(paste0('Elastic net p_value ',round(en_pval, digits=5)))



# --- Step 7: Save Plot ---
#ggsave(
#  filename = 'S2_predicted_all_confluency_combined_gradient_training_gradient_gsva.png',
#  path = figures_filepath,
#  width = 15,
#  height = 9
#)
ggsave(filename='Fig_4B_predicted_all_confluency_combined_gradient_training_gradient_gsva.png', path = main_figure_path, width=9, height=9)



# 1) PCA on the new-only GSVA matrix to get true PC1 for new samples
#gsc.gsva_new_t <- t(GSC.gsva_new)
#pca_new        <- prcomp(gsc.gsva_new_t, center = TRUE, scale. = TRUE)

# extract & flip PC1
#pc1_new_scores <- scale(pca_new$x[,1])[,1]

#gsva_axis_score_pc1 <- scale(gsva_pca$x[, 1])[, 1]

# Flip GSVA PC1 and loadings so that neurodevelopmental = high PC1
#pc1_new_scores <- -pc1_new_scores
#pca_new$rotation[, 1] <- -pca_new$rotation[, 1]

# 1) PCA on the new-only GSVA matrix to get true PC1 for new samples
gsc.gsva_new_t <- t(GSC.gsva_new)
pca_new        <- prcomp(gsc.gsva_new_t, center = TRUE, scale. = TRUE)

# extract & flip PC1
pc1_new_scores <- scale(pca_new$x[,1])[,1]

#gsva_axis_score_pc1 <- scale(gsva_pca$x[, 1])[, 1]

# Flip GSVA PC1 and loadings so that neurodevelopmental = high PC1
pc1_new_scores <- -pc1_new_scores
pca_new$rotation[, 1] <- -pca_new$rotation[, 1]


# make lookup data.frame
pc1_new_df <- data.frame(
  Sample       = rownames(pca_new$x),
  Computed_PC1 = pc1_new_scores,
  stringsAsFactors = FALSE
)


load(paste0('/Users/mystique27m/Documents/Professional/research/PostdoctoralResearch_2020/Projects/PatternRecognitionInGSCs_UsingCV/out/script08_output_files/feature_correlations.Rdata'))
combined_features_df_new <- bind_rows(pca_list) %>%
  select(Sample, all_of(features)) %>%
  group_by(Sample) %>%
  summarise(across(all_of(features), mean, na.rm = TRUE), .groups = 'drop') %>%
  column_to_rownames("Sample")


# --- Predict on New Validation Samples ---
new_samples <- c("G620", "G637", "G683", "G411")
combined_features_df_new <- combined_features_df_new[new_samples,]

x_new <- as.matrix(combined_features_df_new)
#lasso_model_final <- train_model_with_alpha(df_combined, 1)$model
predicted_lasso <- predict(lasso_model_final, newx = x_new)

predicted_lasso <- predict(lasso_model_final, newx = x_new)

pc1_new_only <- pc1_new_scores[which(names(pc1_new_scores) %in% rownames(combined_features_df_new))]

predicted_lasso

df_new <- data.frame(
  samples=rownames(combined_features_df_new),
  Computed_PC1=pc1_new_only,
  Predicted_PC1=predicted_lasso
)

df_new <- df_new %>% rename('Predicted_PC1'='s0')


r2_df <- df_new %>%
  summarize(
    R2    = cor(Computed_PC1, Predicted_PC1, use = "complete.obs")^2,
    x_pos = min(Computed_PC1, na.rm = TRUE),
    y_pos = max(Predicted_PC1, na.rm = TRUE)
  ) %>%
  mutate(label = paste0("R² = ", round(R2, 3)))





# 2) Prepare training set (all features → original PC1)
# 1) Find the samples in common
#common_train <- intersect(rownames(combined_features_df),
#                          names(gsva_axis_score_pc1))

# 2) Build your X and y
#x_train <- as.matrix(combined_features_df[common_train, ])

# since gsva_axis_score_pc1 is a named vector, subset it by name:
#y_train <- gsva_axis_score_pc1[common_train]

# 3) (Optional) drop any all-NA columns in X
#x_train <- x_train[, colSums(is.na(x_train)) == 0]

# 3) Prepare new-feature matrix
#x_new <- as.matrix(combined_features_df_new[, colnames(x_train)])
#validation_results <- do.call(rbind, all_results)



#library(dplyr)
#library(ggplot2)

# 1) Compute per‐Model R² and a common x/y position
#r2_df <- validation_results %>%
#  group_by(Model) %>%
#  summarize(
#    R2    = cor(Computed_PC1, Predicted_PC1, use = "complete.obs")^2,
#    x_pos = min(Computed_PC1, na.rm = TRUE),
#    y_pos = max(Predicted_PC1, na.rm = TRUE)
#  ) %>%
#  mutate(label = paste0("R² = ", round(R2, 3)))

# 2) Plot and add those labels
#ggplot(validation_results, aes(x = Computed_PC1, y = Predicted_PC1)) +
#  geom_point(size=4.5) +
#  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
#  facet_wrap(~ Model) +
#  geom_text(
#    data = r2_df,
#    aes(x = x_pos, y = .9, label = label),
#    inherit.aes = T,
#    hjust = 0,
#    vjust = 1,
#    size = 9
#  ) +
#  theme_minimal(base_size = 14) +
#  theme(#legend.position = "none",
#    axis.title.y = element_text(face='bold', size=21),
#    #axis.title.x = element_blank(),
#    axis.text.x = element_text(size=21),
#    legend.text = element_text(size=18),
#    legend.key.spacing.y = unit(0.45, "cm"),
#    title = element_text(face = 'bold',
#                         size = 21)) +
#  labs(
#    title = "Predicted vs Actual GSVA PC1 (New Samples)",
#    x     = "GSVA PC1",
#    y     = "Predicted GSVA PC1"
#  )

#ggsave(filename='Sup_Fig_S2_repeat/S2_predicted_new_samples_allfeatures_correlation_lineplot.png', path = figures_filepath, width=15, height=6)

#validation_Results_lasso = validation_results[9:12,]
#r2_df_lasso <- subset(r2_df, Model=='Lasso')


ggplot(df_new, aes(x = Computed_PC1 , y = Predicted_PC1)) +
  geom_point(size=4.5) +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
#  facet_wrap(~ Model) #+
  geom_text(
    data = r2_df,
    aes(x = x_pos, y = 2.75, label = label),
    inherit.aes = T,
    hjust = 0,
    vjust = 1,
    size = 9
  ) +
  theme_minimal(base_size = 14) +
  theme(#legend.position = "none",
    axis.title = element_text(face='bold', size=21),
    #axis.title.x = element_blank(),
    axis.text = element_text(size=18),
    legend.text = element_text(size=18),
    legend.key.spacing.y = unit(0.45, "cm"),
    title = element_text(face = 'bold',
                         size = 15)) +
  labs(
    title = "Predicted vs Actual GSVA PC1 (New Samples)",
    x     = "GSVA PC1",
    y     = "Predicted GSVA PC1"
  )

lasso_model_final$beta

cor_values_lasso_model <- cor.test(df_new$Computed_PC1, df_new$Predicted_PC1)
cor(df_new$Computed_PC1, df_new$Predicted_PC1, use='complete.obs')

pval <- cor_values_lasso_model$p.value
r2 <- cor_values_lasso_model$estimate

supplementary_fig_path <- paste0(figures_filepath, 'Sup_Fig_S7/')
dir.create(supplementary_fig_path)

ggsave(filename='S7_predicted_new_samples_allfeatures_correlation_lineplot_lasso.png', path = supplementary_fig_path, width=9, height=9)



#ggplot(r2_df, aes(x=Model, y=R2)) + 
#  geom_bar(stat='identity') + 
#  ylim(0,.25) + 
#  theme_minimal() +
#  theme(axis.text = element_text(size=15),
#        axis.title = element_text(size=21))

#ggsave(filename='Sup_Fig_S2_repeat/S2_predicted_new_samples_allfeatures_correlation_model_performance.png', path = figures_filepath, width=6, height=6)


# make lookup data.frame
#pc1_new_df <- data.frame(
#  Sample       = rownames(pca_new$x),
#  Computed_PC1 = pc1_new_scores,
#  stringsAsFactors = FALSE
#)


# 2) Prepare training set (all features → original PC1)
# 1) Find the samples in common
#common_train <- intersect(rownames(combined_features_df),
#                          names(gsva_axis_score_pc1))

# 2) Build your X and y
#x_train <- as.matrix(combined_features_df[common_train, ])

# since gsva_axis_score_pc1 is a named vector, subset it by name:
#y_train <- gsva_axis_score_pc1[common_train]

# 3) (Optional) drop any all-NA columns in X
#x_train <- x_train[, colSums(is.na(x_train)) == 0]

# 3) Prepare new-feature matrix
#x_new <- as.matrix(combined_features_df_new[, colnames(x_train)])
#validation_results <- do.call(rbind, all_results)



#library(dplyr)
#library(ggplot2)

# 1) Compute per‐Model R² and a common x/y position
#r2_df <- validation_results %>%
#  group_by(Model) %>%
#  summarize(
#    R2    = cor(Computed_PC1, Predicted_PC1, use = "complete.obs")^2,
#    x_pos = min(Computed_PC1, na.rm = TRUE),
#    y_pos = max(Predicted_PC1, na.rm = TRUE)
#  ) %>%
#  mutate(label = paste0("R² = ", round(R2, 3)))

# 2) Plot and add those labels
#ggplot(validation_results, aes(x = Computed_PC1, y = Predicted_PC1)) +
#  geom_point(size=4.5) +
#  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
#  facet_wrap(~ Model) +
#  geom_text(
#    data = r2_df,
#    aes(x = x_pos, y = .9, label = label),
#    inherit.aes = T,
#    hjust = 0,
#    vjust = 1,
#    size = 9
#  ) +
#  theme_minimal(base_size = 14) +
#  theme(#legend.position = "none",
#    axis.title.y = element_text(face='bold', size=21),
#    #axis.title.x = element_blank(),
#    axis.text.x = element_text(size=21),
#    legend.text = element_text(size=18),
#    legend.key.spacing.y = unit(0.45, "cm"),
#    title = element_text(face = 'bold',
#                         size = 21)) +
#  labs(
#    title = "Predicted vs Actual GSVA PC1 (New Samples)",
#    x     = "GSVA PC1",
#    y     = "Predicted GSVA PC1"
#  )

#ggsave(filename='Sup_Fig_S2_repeat/S2_predicted_new_samples_allfeatures_correlation_lineplot.png', path = figures_filepath, width=15, height=6)

#validation_Results_lasso = validation_results[9:12,]
#r2_df_lasso <- subset(r2_df, Model=='Lasso')

#ggplot(validation_Results_lasso, aes(x = Computed_PC1, y = Predicted_PC1)) +
#  geom_point(size=4.5) +
#  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
#  #facet_wrap(~ Model) +
##  geom_text(
#    data = r2_df_lasso,
#    aes(x = x_pos, y = 0.75, label = label),
#    inherit.aes = T,
#    hjust = 0,
#    vjust = 1,
#    size = 9
#  ) +
#  theme_minimal(base_size = 14) +
#  theme(#legend.position = "none",
#    axis.title = element_text(face='bold', size=21),
#    #axis.title.x = element_blank(),
#    axis.text = element_text(size=18),
#    legend.text = element_text(size=18),
#    legend.key.spacing.y = unit(0.45, "cm"),
#    title = element_text(face = 'bold',
#                         size = 15)) +
#  labs(
#    title = "Predicted vs Actual GSVA PC1 (New Samples)",
#    x     = "GSVA PC1",
#    y     = "Predicted GSVA PC1"
#  )

#ggsave(filename='Sup_Fig_S2_repeat/S2_predicted_new_samples_allfeatures_correlation_lineplot_lasso.png', path = figures_filepath, width=9, height=9)


#ggplot(r2_df, aes(x=Model, y=R2)) + 
#  geom_bar(stat='identity') + 
#  ylim(0,.25) + 
#  theme_minimal() +
#  theme(axis.text = element_text(size=15),
#        axis.title = element_text(size=21))

#ggsave(filename='Sup_Fig_S2_repeat/S2_predicted_new_samples_allfeatures_correlation_model_performance.png', path = figures_filepath, width=6, height=6)





#prediction_combined_ridge <- subset(prediction_combined, model_type=='Ridge (α=0)')

# --- Predicted vs Actual Plot ---
#ggplot(prediction_combined_ridge, aes(x = index, y = predicted, color = model_type)) +
#  geom_point(size = 4.5, alpha = 0.9) +
#  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
#  labs(
#    title = "Predicted vs Actual",
#    x = "Groundtruth Neurodev/Inj. Res gradient",
#    y = "Predicted Neurodev/Inj. Res gradient",
#    color = "Model"
#  ) +
#  theme_minimal(base_size = 14) +
#  theme(#legend.position = "none",
#    axis.title.y = element_text(face='bold', size=21),
#    #axis.title.x = element_blank(),
#    axis.text.x = element_text(size=21),
#    legend.text = element_text(size=18),
#    legend.key.spacing.y = unit(0.45, "cm"),
#    title = element_text(face = 'bold',
#                         size = 30)) 

#ggsave(filename='Sup_Fig_S2_repeat/S2_predicted_all_confluency_combined_gradient_training_gradient_gsva_ridge.png', path = figures_filepath, width=15, height=9)



### IMPORT NEW GSVA OF NEW SAMPLES TO BE PREDICTED
# Transpose GSVA matrix
#gsva_new_t <- t(GSC_gsva_new)

GSC.gsva_new <- read.csv(paste0(path_to_repo,'/incucyte_validation/out/GSC.gsva_including_validation_cohort.csv'), row.names=1)

# --- Transpose GSVA matrix and subset new samples ---
gsva_new_t <- t(GSC.gsva_new)
new_samples <- c("G620", "G637", "G683", "G411")
gsva_new_subset <- gsva_new_t[new_samples, ]

# --- Create GSVA DataFrame (top 3 pathways per sample) ---
#top_gsva_long <- as.data.frame(gsva_new_subset) %>%
#  rownames_to_column("Sample") %>%
#  pivot_longer(-Sample, names_to = "Pathway", values_to = "GSVA_Score") %>%
#  group_by(Sample) %>%
#  top_n(10, wt = GSVA_Score)

#################
### Apply above model to new samples to predict
### Apply lasso regression
#################

### NEW data
load(paste0('/Users/mystique27m/Documents/Professional/research/PostdoctoralResearch_2020/Projects/PatternRecognitionInGSCs_UsingCV/incucyte_validation/out/feature_correlations.Rdata'))
combined_features_df_new <- bind_rows(pca_list) %>%
  select(Sample, all_of(features)) %>%
  group_by(Sample) %>%
  summarise(across(all_of(features), mean, na.rm = TRUE), .groups = 'drop') %>%
  column_to_rownames("Sample")


# --- Predict on New Validation Samples ---
new_samples <- c("G620", "G637", "G683", "G411")
combined_features_df_new <- combined_features_df_new[new_samples,]

x_new <- as.matrix(combined_features_df_new)

#for (alpha_val in c(0, 0.5, 1.0)){
#  print(alpha_val)
#  model_allfeat <- train_model_with_alpha(df_combined, alpha_val=alpha_val)
#  model_ridge <- model_allfeat$model
#  pred_df <- model_allfeat$preds
  
#  #predicted_new_scores <- predict(model_2feat, newdata = combined_features_df_new[new_samples,])
#  predicted_new_all_features <- predict(model_ridge, newx=x_new)
  
#  # --- Create DataFrame with Predictions ---
#  validation_results <- data.frame(
#    Sample = rownames(combined_features_df_new),
#    Predicted_PC1 = as.numeric(predicted_new_all_features)
#  )
  
#  #print(validation_results)
#  # --- Plot Prediction Results ---
#  ggplot(validation_results, aes(x = Sample, y = Predicted_PC1)) +
#    geom_bar(stat = "identity", fill = "#619CFF") +
#    theme_minimal(base_size = 14) +
#    labs(title = "Predicted GSVA PC1 for New Samples (2 Feature Model)",
#         x = "Sample", y = "Predicted GSVA PC1") +
#    theme(axis.text.x = element_text(angle = 45, hjust = 1))
#  ggsave(filename=paste0('Sup_Fig_S2_repeat/S2_predicted_new_samples_all_features_training_gradient_gsva_model_', alpha_val,'.png'), path = figures_filepath, width=12, height=6)
  
  
  
  
#  # --- Add Predicted PC1 Scores ---
#  top_gsva_long_new <- top_gsva_long %>%
#    left_join(validation_results, by = "Sample") %>%
#    mutate(Source = "Top GSVA Pathway") %>%
#    rename(Value = GSVA_Score)
  
#  ggplot(top_gsva_long_new, 
#         aes(x=reorder(Sample, Predicted_PC1), 
#             y=Value, 
#             color=Pathway)) + 
#    geom_point() + 
#    geom_label_repel(label=top_gsva_long$Pathway, size=9, label.size=1.25, face='bold') +
#    theme_minimal() +
#    theme(axis.text.x = element_text(size=27, face='bold'),
#          axis.text.y = element_text(size=21, face='bold'),
#          axis.title.x = element_blank(),
#          axis.title.y = element_text(size=42, face='bold'),
#          legend.position = 'none') + 
#    ylab('Gene set enrichment score')
#  ggsave(filename=paste0('Sup_Fig_S2_repeat/S2_predicted_new_samples_confluency_combined_lasso_gsva_gradient_top_gene_signature_modules', alpha_val,'.png'), path = figures_filepath, width=27, height=21)



#  # --- Add Predicted PC1 as its own row per sample ---
#  pred_df <- validation_results %>%
#    mutate(Pathway = "Predicted_PC1",
#           Source = "Predicted PC1",
#           Value = Predicted_PC1)
#  
#  # --- Combine both datasets ---
#  combined_plot_df <- bind_rows(top_gsva_long, pred_df)
#  
#  # --- Plot ---
#  ggplot(combined_plot_df, aes(x = reorder(Pathway, Value), y = Value, fill = Source)) +
#    geom_bar(stat = "identity", position = "dodge") +
#    coord_flip() +
#    facet_wrap(~Sample, scales = "free") +
#    theme_minimal(base_size = 13) +
#    scale_fill_manual(values = c("Top GSVA Pathway" = "#F8766D", "Predicted PC1" = "#619CFF")) +
#    labs(title = "Predicted PC1 vs Top GSVA Pathways (New Validation Samples)",
#         y = "Score", x = "Feature", fill = "Source")
  
  
#  ggsave(filename='Sup_Fig_S2_repeat/S2_predicted_new_samples_all_features_gsva_gradient_top.png', path = figures_filepath, width=12, height=6)
  
  
#  gsva_axis_score_pc1 <- as.data.frame(gsva_axis_score_pc1)
#  gsva_axis_score_pc1$Samples <- rownames(gsva_axis_score_pc1)

#  predicted_PC1 <- as.data.frame(predicted_new_all_features)
#  predicted_PC1$Samples <- rownames(predicted_PC1)
#  colnames(predicted_PC1) <- c('gsva_axis_score_pc1', 'Samples')

#  combined_prediction_trained <- rbind(predicted_PC1, gsva_axis_score_pc1)
  
#  ggplot(combined_prediction_trained, aes(x=reorder(Samples, gsva_axis_score_pc1), 
#                                  y=gsva_axis_score_pc1, 
#                                  group=1)) + 
#    geom_point(size=1) + 
#    geom_line(data=gsva_axis_score_pc1, aes(x=Samples, y=gsva_axis_score_pc1)) +
#    geom_point(data=predicted_PC1 , 
#               aes(x=Samples, 
#                   y=gsva_axis_score_pc1, 
#                   color=Samples), 
#               size=6, 
#               ) +
#    theme_minimal() +
#    theme(axis.text.x = element_text(size=21, hjust=1, angle=45),
#          axis.text.y = element_text(size=15),
#          axis.title.x = element_blank(),
#          axis.title.y = element_text(size=27, face='bold'),
#          legend.key.spacing.y = unit(0.5, 'cm'),
#          legend.text = element_text(size=24),
#          legend.title = element_blank(),
#          line = element_line(linewidth = 0.25)) + 
#    ylab('Neurodev/Inj. Res gradient')
#  ggsave(filename=paste0('Sup_Fig_S2_repeat/S2_predicted_new_samples_confluency_combined_lasso_gsva_gradient_lineplot_', alpha_val,'.png'), path = figures_filepath, width=12, height=9)

#}

###########################

##################
### Validation using only top 2 features from image-ased PC2
### Do not perform regularization
##################

# --- Subset original training data with only the two key features ---

# Ensure sample matching with GSVA PC1 scores
# 1) Define the two features you actually want to use
#features_to_use <- c("InfoMeas1_00", "DifferenceVariance_00")
features_to_use <- c("InfoMeas1_00", "Granularity_13")

# 2) Subset your combined feature matrix
df_small <- combined_features_df[, features_to_use, drop = FALSE]

# 3) Find the samples that exist in both df_small and your PC1 vector
common_samples_small <- intersect(rownames(df_small), names(gsva_axis_score_pc1))
gsva_axis_score_pc1

gsva_axis_score_pc1_z <- scale(gsva_axis_score_pc1)[, 1]  # scale returns a matrix

# 3. Build training data
df_train <- data.frame(
  Sample   = common_samples_small,
  GSVA_PC1 = gsva_axis_score_pc1_z[common_samples_small],
  df_small[common_samples_small, , drop = FALSE]
)


# Inspect
str(df_train)

# 0) Make sure df_train’s rownames are the sample IDs:
#    (Assumes you built df_train exactly as in the previous step)
rownames(df_train) <- df_train$Sample
df_train$Sample   <- NULL

# 1) Extract predictors and response
x_train <- as.matrix(df_train[, features_to_use, drop = FALSE])
y_train <- df_train$GSVA_PC1



# Fit model
#model_2feat <- lm(
#  GSVA_PC1 ~ InfoMeas1_00 + DifferenceVariance_00,
#  data = df_train
#)

model_2feat <- lm(
  GSVA_PC1 ~ InfoMeas1_00 + Granularity_13,
  data = df_train
)


# Make predictions
pred <- predict(model_2feat)

# Calculate R²
r2 <- summary(model_2feat)$r.squared

# Create dataframe
trained_df <- data.frame(
  Sample    = rownames(df_train),
  Actual    = df_train$GSVA_PC1,
  Predicted = pred
)

# Plot
library(ggplot2)

ggplot(trained_df, aes(x = Actual, y = Predicted)) +
  geom_point(color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  xlim(-2.5, 2.5) +
  ylim(-2, 2) +
  annotate("text", x = -1.8, y = 1.8, label = paste("R² =", round(r2, 3)), hjust = 0, size = 5) +
  theme_minimal() +
  #labs(
  #  title = "Predicted vs Actual GSVA_PC1",
  #  x = "Actual GSVA_PC1",
  #  y = "Predicted GSVA_PC1"
  #)  + 
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 21),
    axis.text.y = element_text(size=21),
    legend.text = element_text(size = 18),
    legend.key.spacing.y = unit(0.45, "cm"),
    title = element_text(face = 'bold', size = 30)
  )

ggsave(filename='Sup_Fig_S2_repeat/S2_predicted_new_samples_top2loadings_correlation_lineplot_training.png', path = figures_filepath, width=7, height=7)


t <- cor.test(trained_df$Actual, trained_df$Predicted)
t$p.value

#### Train on new samples using 2 features
# --- Predict on New Validation Samples ---
#new_samples <- c("G620", "G637", "G683", "G411")
#combined_features_df_new <- combined_features_df_new[new_samples,]

#x_new <- as.matrix(combined_features_df_new[, features_to_use])
#predicted_new_scores <- predict(model_2feat, newdata = combined_features_df_new[, features_to_use])

# --- Create DataFrame with Predictions ---
#validation_results <- data.frame(
#  Sample = rownames(combined_features_df_new),
#  Predicted_PC1 = as.numeric(predicted_new_scores)
#)

# --- Plot Prediction Results ---
#ggplot(validation_results, aes(x = Sample, y = Predicted_PC1)) +
#  geom_bar(stat = "identity", fill = "#619CFF") +
#  theme_minimal(base_size = 14) +
#  labs(title = "Predicted GSVA PC1 for New Samples (2 Feature Model)",
#       x = "Sample", y = "Predicted GSVA PC1") +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#ggsave(filename='Sup_Fig_S2_repeat/S2_predicted_new_sampled_top2loadings_training_gradient_gsva.png', path = figures_filepath, width=12, height=6)


### IMPORT NEW GSVA OF NEW SAMPLES TO BE PREDICTED
# Transpose GSVA matrix
#gsva_new_t <- t(GSC_gsva_new)

#GSC.gsva_new <- read.csv(paste0(path_to_repo,'/incucyte_validation/out/GSC.gsva_including_validation_cohort.csv'), row.names=1)

# --- Transpose GSVA matrix and subset new samples ---
#gsva_new_t <- t(GSC.gsva_new)
#new_samples <- c("G620", "G637", "G683", "G411")
#gsva_new_subset <- gsva_new_t[new_samples, ]

# --- Create GSVA DataFrame (top 3 pathways per sample) ---
#top_gsva_long <- as.data.frame(gsva_new_subset) %>%
#  rownames_to_column("Sample") %>%
#  pivot_longer(-Sample, names_to = "Pathway", values_to = "GSVA_Score") %>%
#  group_by(Sample) %>%
#  top_n(10, wt = GSVA_Score)

# --- Add Predicted PC1 Scores ---
#top_gsva_long <- top_gsva_long %>%
#  left_join(validation_results, by = "Sample") %>%
#  mutate(Source = "Top GSVA Pathway") %>%
#  rename(Value = GSVA_Score)

# --- Add Predicted PC1 as its own row per sample ---
#pred_df <- validation_results %>%
#  mutate(Pathway = "Predicted_PC1",
#         Source = "Predicted PC1",
#         Value = Predicted_PC1)

# --- Combine both datasets ---
#combined_plot_df <- bind_rows(top_gsva_long, pred_df)

# --- Plot ---
#ggplot(combined_plot_df, aes(x = reorder(Pathway, Value), y = Value, fill = Source)) +
#  geom_bar(stat = "identity", position = "dodge") +
#  coord_flip() +
#  facet_wrap(~Sample, scales = "free") +
#  theme_minimal(base_size = 13) +
#  scale_fill_manual(values = c("Top GSVA Pathway" = "#F8766D", "Predicted PC1" = "#619CFF")) +
#  labs(title = "Predicted PC1 vs Top GSVA Pathways (New Validation Samples)",
#       y = "Score", x = "Feature", fill = "Source")


#ggsave(filename='Sup_Fig_S2_repeat/S2_predicted_new_samples_top2loadings_gsva_gradient_top_loadings.png', path = figures_filepath, width=12, height=6)

###########################

#gsva_axis_score_pc1 <- as.data.frame(gsva_axis_score_pc1)
#gsva_axis_score_pc1$Samples <- rownames(gsva_axis_score_pc1)
#gsva_axis_score_pc1

#predicted_PC1 <- as.data.frame(predicted_new_scores)
#predicted_PC1$Samples <- rownames(predicted_PC1)
#colnames(predicted_PC1) <- c('gsva_axis_score_pc1', 'Samples')
#predicted_PC1

#combined_prediction_trained <- rbind(predicted_PC1, gsva_axis_score_pc1)

#ggplot(combined_prediction_trained, aes(x=reorder(Samples, gsva_axis_score_pc1), 
#                                        y=gsva_axis_score_pc1, 
#                                        group=1)) + 
#  geom_point(size=1) + 
#  geom_line(data=gsva_axis_score_pc1, aes(x=Samples, y=gsva_axis_score_pc1)) +
#  geom_point(data=predicted_PC1 , 
#             aes(x=Samples, 
#                 y=gsva_axis_score_pc1, 
#                 color=Samples), 
#             size=6, 
#  ) +
#  theme_minimal() +
#  theme(axis.text.x = element_text(size=24, angle=45),
#        axis.title.x = element_blank(),
#        axis.text.y = element_text(size=21),
#        axis.title.y = element_text(size = 30, face = 'bold'),
#        line = element_line(linewidth = 0.25),
#        legend.text = element_text(size = 21),
#        legend.title = element_blank(),
#        legend.key.spacing.y = unit(0.5, 'cm')) + 
#  ylab('Neurodev/Inj. Res gradient')
##ggsave(filename='Sup_Fig_S2_repeat/S2_predicted_new_samples_confluency_combined_top2loadings_gsva_gradient_lineplot.png', path = figures_filepath, width=12, height=9)


### compute actual PC1 score and compare

# Transpose GSVA and extract PC1 axis
#gsc.gsva_combined_t <- t(GSC.gsva_new)
#gsva_pca_combined <- prcomp(gsc.gsva_combined_t, center = TRUE, scale. = TRUE)

# Compute proportion of variance explained
#pve <- (gsva_pca$sdev)^2 / sum(gsva_pca$sdev^2)
#pve_df <- data.frame(
#  PC = paste0("PC", 1:length(pve)),
#  Variance = pve * 100  # in percentage
#)

#gsva_axis_score_pc1_combined <- scale(gsva_pca_combined$x[, 1])[, 1]

# Flip GSVA PC1 and loadings so that neurodevelopmental = high PC1
#gsva_axis_score_pc1 <- -gsva_axis_score_pc1_combined
#gsva_pca$rotation[, 1] <- -gsva_pca$rotation[, 1]


# Flip GSVA PC1 and loadings so that neurodevelopmental = high PC1
#gsva_axis_score_pc1_combined <- -gsva_axis_score_pc1_combined
#gsva_pca_combined$rotation[, 1] <- -gsva_pca_combined$rotation[, 1]

#df <- data.frame(
#  Sample = names(gsva_axis_score_pc1_combined),
#  Computed_PC1 = as.numeric(gsva_axis_score_pc1_combined)
#)

#df_validation <- subset(df, Sample %in% validation_results$Sample)
#validation_results$Computed_PC1 <- df_validation$Computed_PC1



#ggplot(validation_results, aes(x=Computed_PC1, y=Predicted_PC1)) + geom_point()


#model <- lm(Predicted_PC1 ~ Computed_PC1, data = validation_results)
#r2 <- summary(model)$r.squared

#ggplot(validation_results, aes(x = Computed_PC1, y = Predicted_PC1)) +
#  geom_point(size=9) +
#  geom_smooth(method = "lm", se = FALSE, color = "red") +
#  annotate("text", x = min(validation_results$Computed_PC1),
#           y = 0.90,
#           label = paste0("R² = ", round(r2, 3)), size = 12, hjust = 0) + 
#  theme_minimal(base_size = 14) +
#  theme(#legend.position = "none",
#    axis.title = element_text(face='bold', size=21),
#    #axis.title.x = element_blank(),
#    axis.text = element_text(size=21),
#    legend.text = element_text(size=18),
#    legend.key.spacing.y = unit(0.45, "cm"),
#    title = element_text(face = 'bold',
#                         size = 21)) +
#  labs(
#    title = "Predicted vs Actual GSVA PC1 (New Samples) using top 2 image features",
#    x     = "GSVA PC1",
#    y     = "Predicted GSVA PC1"
#  )


#ggsave(filename='Sup_Fig_S2_repeat/S2_predicted_new_samples_top2loadings_correlation_lineplot.png', path = figures_filepath, width=12, height=9)


#######
### USING ALL FEATURES
#library(glmnet)
#library(ggplot2)



## 4) Loop over models, predict, and merge in true PC1
#alpha_vals    <- c(Ridge = 0, ElasticNet = 0.5, Lasso = 1)
#all_results   <- list()

#for (model_name in names(alpha_vals)) {
#  cvm   <- cv.glmnet(x_train, y_train, alpha = alpha_vals[model_name])
#  preds <- predict(cvm, newx = x_new, s = "lambda.min")
  
#  df <- data.frame(
#    Sample        = rownames(x_new),
#    Predicted_PC1 = as.numeric(preds),
#    Model         = model_name,
#    stringsAsFactors = FALSE
#  )
  
  # merge in the true PC1 from the new-only PCA
#  df <- merge(df, pc1_new_df, by = "Sample", all.x = TRUE)
#  #df$Correlation <- with(df, cor(Predicted_PC1, Computed_PC1, use = "complete.obs"))
#  test_result <- cor.test(df$Predicted_PC1, df$Computed_PC1, use = "complete.obs")
#  df$Correlation <- test_result$estimate
#  df$P_value     <- test_result$p.value
  
#  all_results[[model_name]] <- df
#}












# --- Plot ---
# filter out any rows without a Computed_PC1
#to_plot <- subset(validation_results, !is.na(Computed_PC1))

#ggplot(to_plot, aes(x = Computed_PC1, y = Predicted_PC1)) +
#  geom_point(alpha = 0.6) +
#  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
#  facet_wrap(~ Model) +
#  theme_minimal() +
#  labs(
#    title = "Predicted vs Actual GSVA PC1 (New Samples)",
#    x = "GSVA PC1 (Ground Truth)",
#    y = "Predicted GSVA PC1"
#  )

#library(ggplot2)

#ggplot(validation_results, aes(x = Computed_PC1, y = Predicted_PC1)) +
#  geom_point(alpha = 0.6) +
#  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
#  facet_wrap(~ Model) +
#  theme_minimal() +
#  labs(title = "Predicted vs Actual GSVA PC1 (New Samples)",
#       x = "GSVA PC1 (Ground Truth)",
#       y = "Predicted GSVA PC1")


# 1) Compute R² and label positions per model
#r2_df <- validation_results %>%
#  group_by(Model) %>%
#  summarize(
#    R2    = cor(Predicted_PC1, Computed_PC1, use = "complete.obs")^2,
#    x_pos = min(Computed_PC1, na.rm = TRUE),
#    y_pos = max(Predicted_PC1, na.rm = TRUE)
#  )

# 2) Plot with facets and per‐facet R²
#ggplot(validation_results, aes(x = Computed_PC1, y = Predicted_PC1)) +
#  geom_point(alpha = 0.6) +
#  geom_smooth(method = "lm", se = FALSE, color = "red") +
#  facet_wrap(~ Model) +
#  geom_text(
#    data = r2_df,
#    aes(x = x_pos, y = 0.95, label = paste0("R² = ", round(R2, 3))),
#    hjust = 0, vjust = 1, size = 5,
#    inherit.aes = FALSE
#  ) +
#  theme_minimal(base_size = 14) +
#  labs(
#    title = "Predicted vs Actual GSVA PC1 (New Samples)",
#    x     = "GSVA PC1 (Ground Truth)",
#    y     = "Predicted GSVA PC1"
#  )


###


#ggplot(top_gsva_long, 
#       aes(x=reorder(Sample, Predicted_PC1), 
#           y=Value, 
#           color=Pathway)) + 
#  geom_point() + 
#  geom_label_repel(label=top_gsva_long$Pathway, size=9, label.size=1, face='bold') +
#  theme_minimal() +
#  theme(axis.text.x = element_text(size=27, face='bold'),
#        axis.text.y = element_text(size=21, face='bold'),
#        axis.title.x = element_blank(),
#        axis.title.y = element_text(size=42, face='bold'),
#        legend.position = 'none') + 
#  ylab('Gene set enrichment score')

#ggsave(filename='Sup_Fig_S2_repeat/S2_predicted_new_samples_confluency_combined_top2loadings_gsva_gradient_top_gene_signature_modules.png', path = figures_filepath, width=24, height=24)


###########################




remove(list=ls())

