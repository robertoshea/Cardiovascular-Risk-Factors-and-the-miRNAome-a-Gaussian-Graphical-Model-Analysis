#Final 23_02_21

#load libraries
if(T){
  library(GEOquery)
  library(igraph)
  library(miRBaseConverter)
  library(ggnetwork)
  library(ggplot2)
  library(ggpubr)
  library(SILGGM)
  library(huge)
  library(caret)
}

#set working dir
if(T){
  setwd( "D:/Projects/R/Vascular Risk Factors Lasso")
}

#set output folders
if(T){
  figures_folder <- "vrf_figures"
  tables_folder <- "vrf_tables"
}

#download data from Gene Expression Omnibus
if(T){
  
  GSE117064 <- getGEO("GSE117064",GSEMatrix=TRUE)
  save(GSE117064, file = "GSE117064.Rdata")
  non_cvd_patients <- !GSE117064[[1]]$`group:ch1`%in%c("3A", "3B") # get subset of patients in control cohort
  X <- t(exprs(GSE117064[[1]]))[non_cvd_patients,] #get expression data
  colnames(X) <- gsub(",.*", "", colnames(X))
  colnames(X) <- miRNA_AccessionToName(colnames(X),
                                       targetVersion = "v21")$TargetName #convert names to miR format
}

#preprocess data
if(T){
  
  #read in mirtarbase data
  mirtarbase_df <- read.csv("hsa_MTI.csv")
  mirtarbase_df <- mirtarbase_df[-grep("Weak", mirtarbase_df$Support.Type),]
  
  all_mirna <- intersect(colnames(X),mirtarbase_df$miRNA)
  
  X <- X[,all_mirna]
  mirtarbase_df <- mirtarbase_df[mirtarbase_df$miRNA %in% all_mirna,]
  
  #phenotypic variables
  age <- as.numeric(GSE117064[[1]]$`age:ch1`[non_cvd_patients])
  bmi <- as.numeric(GSE117064[[1]]$`bmi:ch1`[non_cvd_patients])
  sbp <- as.numeric(GSE117064[[1]]$`systolic bp:ch1`[non_cvd_patients])
  dbp <- as.numeric(GSE117064[[1]]$`diastolic bp:ch1`[non_cvd_patients])
  map <- ((2*dbp)+sbp)/3
  hba1c <-  as.numeric(GSE117064[[1]]$`hb-a1c:ch1`[non_cvd_patients])
  male <- (GSE117064[[1]]$`Sex:ch1`=="Male")[non_cvd_patients]*1
  smoking <- as.numeric(GSE117064[[1]]$`smoking:ch1`[non_cvd_patients])
  
  X <- do.call(cbind, list(Age=age,
                           BMI=bmi,
                           MAP=map,
                           HbA1c=hba1c,
                           Male=male,
                           Smoking=smoking,
                           X))
  
  binary_vars <- colnames(X)%in%c("Male")
  X[,!binary_vars] <- scale(apply(X[,!binary_vars],2,rank))
  
}

#risk factor prevalence
if(F){
  
  age_hist <- gghistogram(age,
                          bins=10,
                          fill = "darkgoldenrod3", xlab="Age",ylab="N")
  bmi_hist <- gghistogram(bmi, bins=10, fill="deepskyblue3",
                          xlab=expression(Body ~ Mass ~ Index ~ (kg/m^2)),ylab="N")
  map_hist <- gghistogram(map, bins=10, fill="darkorchid1",
                          xlab="Mean Arterial Pressure (mmHg)",ylab="N")
  hba1c_hist <- gghistogram(hba1c, bins=20, fill="chocolate1",
                            xlab="HbA1c (%)",ylab="N")
  smoking_hist <- gghistogram(smoking, bins=5, fill="lightcoral",
                              xlab="Smoking Rate (Cigarettes/Day)",ylab="N")
  
  
  sex_tbl <- as.data.frame(table(ifelse(male==1, "Male", "Female")))
  sex_barplot <- ggbarplot(data=sex_tbl,
                           x="Var1",
                           y="Freq",
                           fill="darkseagreen",
                           xlab="Sex",
                           ylab="N"
  )
  
  all_pheno_plots <- egg::ggarrange(plots = list(
    age_hist,
    sex_barplot,
    bmi_hist,
    map_hist,
    hba1c_hist,
    smoking_hist
  ),nrow=2, ncol=3)
  
  ggsave(plot=all_pheno_plots,
         filename = paste0(figures_folder,"/all_pheno_plots.pdf"),
         dpi=1500,
         units="mm",
         width = 250,
         height = 125)
  
  ggsave(plot=all_pheno_plots,
         filename = paste0(figures_folder,"/all_pheno_plots.jpg"),
         dpi=1500,
         units="mm",
         width = 250,
         height = 125)
  
  iqr_func <- function(vec, digits=1){
    quartiles <- round(quantile(vec, c(0.025,0.975)),digits=digits)
    paste0("(Interquartile Range=", quartiles[1], "-", quartiles[2], ")")
  }
  med_func <- function(vec){
    round(median(vec),1)
  }
  
}

#infer graphical model
if(T){
  
  S <- cor(X)
  alpha_level <- 0.05
  
  h <- SILGGM(X, method="D-S_NW_SL",
              global = TRUE,
              alpha=alpha_level)
  am <- h$global_decision[[1]]
  dimnames(am)<- list(colnames(X), colnames(X))
  g <- graph_from_adjacency_matrix(am, mode="undirected")
  V(g)$name <- gsub("hsa-", "", V(g)$name)
  pheno_vars <- 1:6
  
  #adjust p values with false discovery rate control
  h$p_precision_adj <- h$p_precision
  h$p_precision_adj[upper.tri(h$p_precision_adj)] <-
    p.adjust(h$p_precision_adj[upper.tri(h$p_precision_adj)], method="fdr")
  h$p_precision_adj[lower.tri(h$p_precision_adj)] <- 0
  h$p_precision_adj <- h$p_precision_adj + t(h$p_precision_adj)
  
  #record correlation types
  cor_type_mat <- S
  cor_type_mat[]<- "Spearman"
  cor_type_mat[!binary_vars,binary_vars]<-
    cor_type_mat[binary_vars,!binary_vars] <- "Point-Biserial"
}

#extract ego graph of phenotypic variables
if(T){
  
  #exclude miRNA which don't interact directly with the phenotypic variables
  ego_rf <- ego(g, order=1, nodes = pheno_vars)
  names(ego_rf)<- colnames(X)[pheno_vars]
  
  included_nodes <- unique(unlist(lapply(ego_rf, as.vector)))
  g_rf <- induced.subgraph(g, vids=included_nodes)
  g_df <- as_data_frame(g_rf)
  e_idx <- cbind(match(g_df$from,  V(g)$name),
                 match(g_df$to,   V(g)$name)
  )
  
  g_df$partial_cor <- h$partialCor[e_idx]
  g_df$precision <- h$precision[e_idx]
  g_df$precision_lo <- h$CI_low_precision[e_idx]
  g_df$precision_hi <- h$CI_high_precision[e_idx]
  g_df$precision_z <- h$z_score_precision[e_idx]
  g_df$p <- h$p_precision_adj[e_idx]
  g_df$Interaction <- ifelse(g_df$partial_cor>0, "Upregulation", "Inhibition")
  g_df$"Bivariate Corr." <- S[e_idx]
  g_df$"Corr. Type" <- cor_type_mat[e_idx]
  
  g_rf <- graph_from_data_frame(g_df[,c("from", "to", "Interaction")], directed=F)
}

#plot ego graph
if(T){
  
  g_rf_ggnetwork_df <- ggnetwork(g_rf)
  myColors <- c("magenta","springgreen")
  names(myColors) <- na.omit(unique(g_rf_ggnetwork_df$Interaction))
  colScale <- scale_colour_manual(name = "Interaction",values = myColors)

  g_rf_plot <- ggplot(g_rf_ggnetwork_df, aes(x,y,xend=xend, yend=yend))+
    geom_edges(aes(colour=Interaction), size=1)+
    geom_nodes(alpha=1, size = 2.5, colour="slategray")+
    geom_nodetext_repel(aes(label=name), size=3.5,
                        data =function(x) { x[ !x$name%in%V(g)$name[pheno_vars], ]})+
    geom_nodelabel(aes(label=name),
                   data =function(x) { x[ x$name%in%V(g)$name[pheno_vars], ]},
                   fontface = "bold",)+
    colScale+
    theme_void()+
    theme(plot.margin=unit(c(1,1,1,1),"cm"))
  
  ggsave(plot=g_rf_plot,
         filename = paste0(figures_folder,"/g_rf_plot.jpg"),
         dpi=1500,
         units="in",
         width = 8.5,
         height = 7.5)
  
  ggsave(plot=g_rf_plot,
         filename = paste0(figures_folder,"/g_rf_plot.pdf"),
         dpi=1500,
         units="in",
         width = 8.5,
         height = 7.5)
  
}

#bootstrap analysis
if(T){
  
  signed_am <- am*sign(h$partialCor)
  
  n_bs <- 500
  bs_am_list <- list()
  bs_pcor_list <- list()
  for(bs_i in 1:n_bs){
    message(bs_i)
    X_bs <- X[sample(nrow(X), replace=T),]
    h_bs <- SILGGM(X_bs, method="D-S_NW_SL",
                   global = TRUE,
                   alpha=alpha_level)
    bs_pcor_list[[bs_i]] <- h_bs$partialCor
    bs_am_list[[bs_i]]<- h_bs$global_decision[[1]]*sign(h_bs$partialCor)
  }
  bs_am_mat <- sapply(bs_am_list, function(i){
    i[e_idx]
  })
  bs_pcor_mat <- sapply(bs_pcor_list, function(i){
    i[e_idx]
  })
  g_df$Stability <- rowMeans(sign(bs_am_mat)==signed_am[e_idx])
  g_df$partial_cor_lo <- apply(bs_pcor_mat, 1, function(i){
    quantile(i, 0.025)
  })
  g_df$partial_cor_hi <- apply(bs_pcor_mat, 1, function(i){
    quantile(i, 0.975)
  })
  
}

#generating result outputs
if(T){
  
  rounding_func <- function(i){
    format(round(i, digits=2), nsmall=2)
  }
  p_rounding_func <- function(i){
    
    sapply(i, function(i_i){
      if(i_i<0.001){
        return("<.001")
      }else if(i_i<0.01){
        i_i <- format(round(i_i, 3), nsmall=3)
        i_i <- substr(i_i, start = 2, stop = nchar(i_i))
      }else{
        i_i <- rounding_func(i_i)
        i_i <- substr(i_i, start = 2, stop = nchar(i_i))
      }
      return(i_i)
    })
    
  }
  
  #produce result tables
  g_df2 <- g_df
  g_df2$"Partial Corr." <-  paste0(
    rounding_func(g_df2$partial_cor), " [",
    rounding_func(g_df2$partial_cor_lo), ",",
    rounding_func(g_df2$partial_cor_hi), "]"
  )
  g_df2$"P-value (adjusted)" <- p_rounding_func(g_df2$p)
  colnames(g_df2)[1:2]<- c("Variable A", "Variable B")
  g_df2$"Bivariate Corr." <- rounding_func( g_df2$"Bivariate Corr.")
  g_df2 <- g_df2[,c("Variable A","Variable B","Corr. Type",
                    "Bivariate Corr.", "Partial Corr.",
                    "P-value (adjusted)", "Stability")]
  g_df2 <- g_df2[order(g_df2[,1], g_df2[,2]),]
  
  g_df_pheno <- g_df2[g_df2$`Variable A`%in%colnames(X)[pheno_vars]&
                        g_df2$`Variable B`%in%colnames(X)[pheno_vars],
  ]
  g_df_pheno_mir <-  g_df2[xor(g_df2$`Variable A`%in%colnames(X)[pheno_vars],
                               g_df2$`Variable B`%in%colnames(X)[pheno_vars]),
  ]
  
  write.table(g_df_pheno, paste0(tables_folder,
                                 "/g_df_pheno.tsv"),
              sep="\t",
              row.names = F)
  
  write.table(g_df_pheno_mir, paste0(tables_folder,
                                     "/g_df_pheno_mir.tsv"),
              sep="\t",
              row.names = F)
  
  #full network table for supplementary data
  g_all_df <- as_data_frame(g)
  e_all_idx <- cbind(match(g_all_df$from,  V(g)$name),
                     match(g_all_df$to,   V(g)$name)
  )
  
  bs_all_am_mat <- sapply(bs_am_list, function(i){
    i[e_all_idx]
  })
  bs_all_pcor_mat <- sapply(bs_pcor_list, function(i){
    i[e_all_idx]
  })
  g_all_df$Stability <- rowMeans(sign(bs_all_am_mat)==signed_am[e_all_idx])
  g_all_df$partial_cor_lo <- apply(bs_all_pcor_mat, 1, function(i){
    quantile(i, 0.025)
  })
  g_all_df$partial_cor_hi <- apply(bs_all_pcor_mat, 1, function(i){
    quantile(i, 0.975)
  })

  g_all_df$"Partial Corr." <- paste0(
    rounding_func(h$partialCor[e_all_idx]), " [",
    rounding_func(g_all_df$partial_cor_lo), ",",
    rounding_func(g_all_df$partial_cor_hi), "]"
  )
    
  g_all_df$"P-value (adjusted)" <- p_rounding_func(h$p_precision_adj[e_all_idx])
  g_all_df$"Bivariate Corr." <- rounding_func(S[e_all_idx])
  g_all_df$"Corr. Type" <- cor_type_mat[e_all_idx]
  g_all_df$Precision <- paste0(
    rounding_func(h$precision[e_all_idx]), " [",
    rounding_func(h$CI_low_precision[e_all_idx]), ",",
    rounding_func(h$CI_high_precision[e_all_idx]), "]"
  )
  g_all_df$Stability <- h$stability[e_all_idx]
  colnames(g_all_df)[1:2]<- c("Variable A", "Variable B")
  
  g_all_df <- g_all_df[,c("Variable A","Variable B","Corr. Type",
                          "Bivariate Corr.", "Partial Corr.",
                          "P-value (adjusted)", "Stability")]
  g_all_df <- g_all_df[order(g_all_df[,1], g_all_df[,2]),]
  
  write.table(g_all_df, paste0(tables_folder,
                               "/g_all_df.tsv"),
              sep="\t",
              row.names = F)
}

#plot bootstrap confidence intervals
if(T){
  
  bs_df = data.frame(
    partial_cor = g_df$partial_cor,
    partial_cor_lo = g_df$partial_cor_lo,
    partial_cor_hi = g_df$partial_cor_hi
  )
  bs_df$Interaction <- ifelse(bs_df$partial_cor>0, "Upregulation", "Inhibition")
  bs_df$vertices <- tolower(paste(g_df$from, g_df$to))
  bs_df$Rank <- rank(bs_df$partial_cor)
  
  colScale_fill <- scale_fill_manual(name = "Interaction",values = myColors)
  bs_plot <- ggplot(bs_df, aes(x=Rank, y=partial_cor))+
    geom_line(size=1)+
    geom_ribbon(aes(ymin=partial_cor_lo,
                    ymax=partial_cor_hi,
                    fill=Interaction),
                alpha=0.3,  
                size=1)+
    facet_wrap(facets=Interaction~., scales="free_x")+
    colScale_fill+
    ylab("Edge Partial Correlation (95% CI)")+
    xlab("Edge Rank")+
    geom_hline(yintercept = 0, colour = "darkgrey")+
    theme(legend.position="None")
  
  
  network_and_bs_plot <- ggarrange(plotlist=list(bs_plot, g_rf_plot), widths=c(1,3))
  ggsave(plot=network_and_bs_plot,
         filename = paste0(figures_folder,"/network_and_bs_plot.jpg"),
         dpi=1500,
         units="in",
         width = 8.5,
         height = 6)
}

#save analysis
if(F){
  save.image("mirna_rf_network_23_02_21.RData")
}

#reload analysis
if(F){
  load("mirna_rf_network_23_02_21.RData")
}
