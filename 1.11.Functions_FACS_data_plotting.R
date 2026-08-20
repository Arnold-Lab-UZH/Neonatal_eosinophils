### Function for plotting of FACS data along time for each tissue 
FACS_plotting_heatmap_line_graph_and_stats_along_time_points <- function(
    input_path,
    measurement,
    output_path
){
  ### load data
  df1 <- read.csv(paste0(input_path,"BM_",measurement,".csv"), sep = ";")
  df2 <- read.csv(paste0(input_path,"BLD_",measurement,".csv"), sep = ";")
  df3 <- read.csv(paste0(input_path,"SPL_",measurement,".csv"), sep = ";")
  df4 <- read.csv(paste0(input_path,"CO_",measurement,".csv"), sep = ";")
  df5 <- read.csv(paste0(input_path,"SI_",measurement,".csv"), sep = ";")
  
  # remove the first column which has the mouse number information 
  df1 <- df1[,2:length(colnames(df1))]
  df2 <- df2[,2:length(colnames(df2))]
  df3 <- df3[,2:length(colnames(df3))]
  df4 <- df4[,2:length(colnames(df4))]
  df5 <- df5[,2:length(colnames(df5))]
  
  # change the , to dots 
  df1[] <- lapply(df1, function(x) {if (is.character(x)) {as.numeric(gsub(",", ".", x))} else {x}})
  df2[] <- lapply(df2, function(x) {if (is.character(x)) {as.numeric(gsub(",", ".", x))} else {x}})
  df3[] <- lapply(df3, function(x) {if (is.character(x)) {as.numeric(gsub(",", ".", x))} else {x}})
  df4[] <- lapply(df4, function(x) {if (is.character(x)) {as.numeric(gsub(",", ".", x))} else {x}})
  df5[] <- lapply(df5, function(x) {if (is.character(x)) {as.numeric(gsub(",", ".", x))} else {x}})
  
  # add a tissue column
  df1$Tissue <- "BM"
  df2$Tissue <- "BLD"
  df3$Tissue <- "SPL"
  df4$Tissue <- "CO"
  df5$Tissue <- "SI"
  
  ## combine all to one dataframe 
  dfs <- list(df1, df2, df3, df4, df5)
  
  clean_df <- function(df, name){
    df_long <- df %>%
      pivot_longer(
        cols = -Tissue,
        names_to = "Time_point",
        values_to = "Value"
      )
    return(df_long)
  }
  
  all_data <- bind_rows(
    clean_df(df1, "BM"),
    clean_df(df2, "BLD"),
    clean_df(df3, "SPL"),
    clean_df(df4, "CO"),
    clean_df(df5, "SI")
  )
  
  ### statistical test 
  all_data_df <- as.data.frame(all_data)
  all_data_df <- na.omit(all_data_df)
  
  P0 <- all_data_df[all_data_df$Time_point %in% "P0",]
  anova <- aov(Value ~ Tissue, data = P0)
  summary(anova)
  results <- TukeyHSD(anova)
  results_P0 <- as.data.frame(results$Tissue)
  results_P0$Time_point <- "P0"
  results_P0$Comparison <- rownames(results_P0)
  rownames(results_P0) <- NULL
  
  P4 <- all_data_df[all_data_df$Time_point %in% "P4",]
  anova <- aov(Value ~ Tissue, data = P4)
  summary(anova)
  results <- TukeyHSD(anova)
  results_P4 <- as.data.frame(results$Tissue)
  results_P4$Time_point <- "P4"
  results_P4$Comparison <- rownames(results_P4)
  rownames(results_P4) <- NULL
  
  P7 <- all_data_df[all_data_df$Time_point %in% "P7",]
  anova <- aov(Value ~ Tissue, data = P7)
  summary(anova)
  results <- TukeyHSD(anova)
  results_P7 <- as.data.frame(results$Tissue)
  results_P7$Time_point <- "P7"
  results_P7$Comparison <- rownames(results_P7)
  rownames(results_P7) <- NULL
  
  P14 <- all_data_df[all_data_df$Time_point %in% "P14",]
  anova <- aov(Value ~ Tissue, data = P14)
  summary(anova)
  results <- TukeyHSD(anova)
  results_P14 <- as.data.frame(results$Tissue)
  results_P14$Time_point <- "P14"
  results_P14$Comparison <- rownames(results_P14)
  rownames(results_P14) <- NULL
  
  P21 <- all_data_df[all_data_df$Time_point %in% "P21",]
  anova <- aov(Value ~ Tissue, data = P21)
  summary(anova)
  results <- TukeyHSD(anova)
  results_P21 <- as.data.frame(results$Tissue)
  results_P21$Time_point <- "P21"
  results_P21$Comparison <- rownames(results_P21)
  rownames(results_P21) <- NULL
  
  Adult <- all_data_df[all_data_df$Time_point %in% "Adult.",]
  anova <- aov(Value ~ Tissue, data = Adult)
  summary(anova)
  results <- TukeyHSD(anova)
  results_Adult <- as.data.frame(results$Tissue)
  results_Adult$Time_point <- "Adult"
  results_Adult$Comparison <- rownames(results_Adult)
  rownames(results_Adult) <- NULL
  
  results_all <- bind_rows(results_P0, results_P4, results_P7,results_P14,results_P21,results_Adult)
  results_all <- results_all[,colnames(results_all) %in% c("p adj","Time_point","Comparison")]
  write.csv(results_all,paste0(output_path,measurement,"_one_way_anova_tukeyHSD_along_time_points.csv"))
  
  # calculate the median for plotting 
  df <- all_data %>%
    group_by(Tissue, Time_point) %>%
    summarise(
      Median = median(Value, na.rm = TRUE),
      .groups = "drop"
    )
  df[is.na(df)] <- 0
  df <- as.data.frame(df)
  
  # set the order of Time points 
  df$Time_point <- factor(
    df$Time_point,
    levels = c("P0", "P4", "P7", "P14", "P21", "Adult.")
  )
  
  # Plot line graph 
  p <- ggplot(df, aes(x = Time_point, y = Median,
                      group = Tissue, color = Tissue)) +
    geom_point(size = 3) +
    geom_smooth(se = FALSE, method = "loess") +
    scale_color_manual(values = c(
      "BM"  = "blue",
      "BLD" = "red",
      "SPL" = "purple",
      "CO"  = "green",
      "SI"  = "black"
    )) +   
    labs(y = paste0("Median ", measurement)) +
    labs(x = "Time point") +
    theme_classic()
  ggsave(paste0(output_path, measurement,"_line_graph_along_time_points.pdf"), width = 8, height = 6, plot = p)
  
  ### plot a heatmap 
  #df$Tissue <- factor(df$Tissue, levels = c("SI","CO","SPL","BLD","BM"))
  #p <- ggplot(df, aes(x = Time_point, y = Tissue, fill = Median)) +
  #  geom_tile() +
  #  scale_fill_viridis_c(option = "C") +
  #  theme_classic() +
  #  labs(
  #    x = "Time Point",
  #    y = "Tissue",
  #    fill = paste0("Median ", measurement)
  #  )
  #ggsave(paste0(output_path,measurement,"_heatmap_along_time_points.pdf"), width = 8, height = 6, plot = p)
}

FACS_plotting_heatmap_line_graph_and_stats_along_pseudotime <- function(
    input_path,
    measurement,
    output_path
){
  ### load data
  df1 <- read.csv(paste0(input_path,"BM_",measurement,".csv"), sep = ";")
  df2 <- read.csv(paste0(input_path,"BLD_",measurement,".csv"), sep = ";")
  df3 <- read.csv(paste0(input_path,"SPL_",measurement,".csv"), sep = ";")
  df4 <- read.csv(paste0(input_path,"CO_",measurement,".csv"), sep = ";")
  df5 <- read.csv(paste0(input_path,"SI_",measurement,".csv"), sep = ";")
  
  # remove the first column which has the mouse number information 
  df1 <- df1[,2:length(colnames(df1))]
  df2 <- df2[,2:length(colnames(df2))]
  df3 <- df3[,2:length(colnames(df3))]
  df4 <- df4[,2:length(colnames(df4))]
  df5 <- df5[,2:length(colnames(df5))]
  
  # change the , to dots 
  df1[] <- lapply(df1, function(x) {if (is.character(x)) {as.numeric(gsub(",", ".", x))} else {x}})
  df2[] <- lapply(df2, function(x) {if (is.character(x)) {as.numeric(gsub(",", ".", x))} else {x}})
  df3[] <- lapply(df3, function(x) {if (is.character(x)) {as.numeric(gsub(",", ".", x))} else {x}})
  df4[] <- lapply(df4, function(x) {if (is.character(x)) {as.numeric(gsub(",", ".", x))} else {x}})
  df5[] <- lapply(df5, function(x) {if (is.character(x)) {as.numeric(gsub(",", ".", x))} else {x}})
  
  # add a tissue column
  df1$Tissue <- "BM"
  df2$Tissue <- "BLD"
  df3$Tissue <- "SPL"
  df4$Tissue <- "CO"
  df5$Tissue <- "SI"
  
  ## combine all to one dataframe 
  dfs <- list(df1, df2, df3, df4, df5)
  
  clean_df <- function(df, name){
    df_long <- df %>%
      pivot_longer(
        cols = -Tissue,
        names_to = "Time_point",
        values_to = "Value"
      )
    return(df_long)
  }
  
  all_data <- bind_rows(
    clean_df(df1, "BM"),
    clean_df(df2, "BLD"),
    clean_df(df3, "SPL"),
    clean_df(df4, "CO"),
    clean_df(df5, "SI")
  )
  
  ### statistical test 
  all_data_df <- as.data.frame(all_data)
  all_data_df <- na.omit(all_data_df)
  
  BM <- all_data_df[all_data_df$Tissue %in% "BM",]
  anova <- aov(Value ~ Time_point, data = BM)
  summary(anova)
  results <- TukeyHSD(anova)
  results_BM <- as.data.frame(results$Time_point)
  results_BM$Tissue <- "BM"
  results_BM$Comparison <- rownames(results_BM)
  rownames(results_BM) <- NULL
  
  BLD <- all_data_df[all_data_df$Tissue %in% "BLD",]
  anova <- aov(Value ~ Time_point, data = BLD)
  summary(anova)
  results <- TukeyHSD(anova)
  results_BLD <- as.data.frame(results$Time_point)
  results_BLD$Tissue <- "BLD"
  results_BLD$Comparison <- rownames(results_BLD)
  rownames(results_BLD) <- NULL
  
  SPL <- all_data_df[all_data_df$Tissue %in% "SPL",]
  anova <- aov(Value ~ Time_point, data = SPL)
  summary(anova)
  results <- TukeyHSD(anova)
  results_SPL <- as.data.frame(results$Time_point)
  results_SPL$Tissue <- "SPL"
  results_SPL$Comparison <- rownames(results_SPL)
  rownames(results_SPL) <- NULL
  
  CO <- all_data_df[all_data_df$Tissue %in% "CO",]
  anova <- aov(Value ~ Time_point, data = CO)
  summary(anova)
  results <- TukeyHSD(anova)
  results_CO <- as.data.frame(results$Time_point)
  results_CO$Tissue <- "CO"
  results_CO$Comparison <- rownames(results_CO)
  rownames(results_CO) <- NULL
  
  SI <- all_data_df[all_data_df$Tissue %in% "SI",]
  anova <- aov(Value ~ Time_point, data = SI)
  summary(anova)
  results <- TukeyHSD(anova)
  results_SI <- as.data.frame(results$Time_point)
  results_SI$Tissue <- "SI"
  results_SI$Comparison <- rownames(results_SI)
  rownames(results_SI) <- NULL
  
  results_all <- bind_rows(results_BM, results_BLD, results_SPL,results_CO,results_SI)
  results_all <- results_all[,colnames(results_all) %in% c("p adj","Tissue","Comparison")]
  write.csv(results_all,paste0(output_path,measurement,"_one_way_anova_tukeyHSD_along_pseudotime.csv"))
  
  # calculate the median for plotting 
  df <- all_data %>%
    group_by(Tissue, Time_point) %>%
    summarise(
      Median = median(Value, na.rm = TRUE),
      .groups = "drop"
    )
  df[is.na(df)] <- 0
  df <- as.data.frame(df)
  
  # set the order of Time points 
  df$Time_point <- factor(
    df$Time_point,
    levels = c("P0", "P4", "P7", "P14", "P21", "Adult.")
  )
  
  df$Tissue <- factor(
    df$Tissue,
    levels = c("BM", "BLD", "SPL", "CO", "SI")
  )
  
  # Plot line graph 
  p <- ggplot(df, aes(x = Tissue, y = Median,
                      group = Time_point, color = Time_point)) +
    geom_point(size = 3) +
    geom_smooth(se = FALSE, method = "loess") +
    scale_color_manual(values = c(
      "P0"  = "red",
      "P4" = "blue",
      "P7" = "green",
      "P14"  = "grey",
      "P21"  = "orange",
      "Adult."  = "black"
    )) +   
    labs(y = paste0("Median ", measurement)) +
    labs(x = "Pseudotime") +
    theme_classic()
  ggsave(paste0(output_path, measurement,"_line_graph_along_pseudotime.pdf"), width = 8, height = 6, plot = p)
  
  ### plot a heatmap 
  #df$Time_point <- factor(
  #  df$Time_point,
  #  levels = c("Adult", "P21", "P14", "P7", "P4", "P0")
  #)
  #p <- ggplot(df, aes(x = Tissue, y = Time_point, fill = Median)) +
  #  geom_tile() +
  #  scale_fill_viridis_c(option = "C") +
  #  theme_classic() +
  #  labs(
  #    x = "Time Point",
  #    y = "Tissue",
  #    fill = paste0("Median ", measurement)
  #  )
  #ggsave(paste0(output_path,measurement,"_heatmap_along_pseudotime.pdf"), width = 8, height = 6, plot = p)
}

### Function to plot measurements of interest for one tissue along time in a heatmap 
FACS_plotting_heatmap_tissue_of_interest_non_scaled <- function(
    input_path,
    tissue_oi, 
    measurements,
    output_path,
    id
){
  ### load data
  df_list  <- list()
  for (i in measurements) { 
    df <- read.csv(paste0(input_path,tissue_oi,i,".csv"), sep = ";")
    df <- df[,2:length(colnames(df))]
    df[] <- lapply(df, function(x) {if (is.character(x)) {as.numeric(gsub(",", ".", x))} else {x}})
    df$Measurement <- i 
    df_list[[i]] <- df
  }
  
  clean_df <- function(df, name){
    df_long <- df %>%
      pivot_longer(
        cols = -Measurement,
        names_to = "Time_point",
        values_to = "Value"
      )
    return(df_long)
  }
  
  a <- lapply(df_list, clean_df)
  
  all_data <- bind_rows(a)
  
  # calculate the median for plotting 
  df <- all_data %>%
    group_by(Measurement, Time_point) %>%
    summarise(
      Median = median(Value, na.rm = TRUE),
      .groups = "drop"
    )
  df[is.na(df)] <- 0
  df <- as.data.frame(df)
  
  # set the order of Time points 
  df$Time_point <- factor(
    df$Time_point,
    levels = c("P0", "P4", "P7", "P14", "P21", "Adult.")
  )
  
  ### plot a heatmap 
  p <- ggplot(df, aes(x = Time_point, y = Measurement, fill = Median)) +
    geom_tile() +
    scale_fill_viridis_c(option = "C") +
    theme_classic() +
    labs(
      x = "Time Point",
      y = "Measurement",
      fill = paste0("Median ", id)
    )
  ggsave(paste0(output_path,tissue_oi,id, "_heatmap_along_time_points.pdf"), width = 8, height = 6, plot = p)
}

FACS_plotting_heatmap_tissue_of_interest_scaled <- function(
    input_path,
    tissue_oi, 
    measurements,
    output_path,
    id
){
  ### load data
  df_list  <- list()
  for (i in measurements) { 
    df <- read.csv(paste0(input_path,tissue_oi,i,".csv"), sep = ";")
    df <- df[,2:length(colnames(df))]
    df[] <- lapply(df, function(x) {if (is.character(x)) {as.numeric(gsub(",", ".", x))} else {x}})
    df$Measurement <- i 
    df_list[[i]] <- df
  }
  
  clean_df <- function(df, name){
    df_long <- df %>%
      pivot_longer(
        cols = -Measurement,
        names_to = "Time_point",
        values_to = "Value"
      )
    return(df_long)
  }
  
  a <- lapply(df_list, clean_df)
  
  all_data <- bind_rows(a)
  
  # calculate the median for plotting 
  df <- all_data %>%
    group_by(Measurement, Time_point) %>%
    summarise(
      Median = median(Value, na.rm = TRUE),
      .groups = "drop"
    )
  df[is.na(df)] <- 0
  df <- as.data.frame(df)
  
  # set the order of Time points 
  df$Time_point <- factor(
    df$Time_point,
    levels = c("P0", "P4", "P7", "P14", "P21", "Adult.")
  )
  
  ### scale 0 to 1 
  df_scaled <- df %>%
    group_by(Measurement) %>%
    mutate(Median_scaled = (Median - min(Median, na.rm = TRUE)) /
             (max(Median, na.rm = TRUE) - min(Median, na.rm = TRUE))) %>%
    ungroup()
  
  ### plot a heatmap 
  p <- ggplot(df_scaled, aes(x = Time_point, y = Measurement, fill = Median_scaled)) +
    geom_tile() +
    scale_fill_viridis_c(option = "C") +
    theme_classic() +
    labs(
      x = "Time Point",
      y = "Measurement",
      fill = paste0("Median scaled ", id)
    )
  ggsave(paste0(output_path,tissue_oi,id, "_heatmap_along_time_points_scaled.pdf"), width = 8, height = 6, plot = p)
}
