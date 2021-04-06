#Cardiometabolic risk factors and the miRNAome
#23_03_21

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

#set output folders
if(T){
  figures_folder <- "vrf_figures"
  tables_folder <- "vrf_tables"
}

#download data from Gene Expression Omnibus
if(T){
  GSE117064 <- getGEO("GSE117064",GSEMatrix=TRUE)
  non_cvd_patients <- !GSE117064[[1]]$`group:ch1`%in%c("3A", "3B") # get subset of patients in control cohort
  X <- t(exprs(GSE117064[[1]]))[non_cvd_patients,] #get expression data
  
  colnames(X) <- gsub(",.*", "", colnames(X))
  colnames(X) <- miRNA_AccessionToName(colnames(X),
                                       targetVersion = "v21")$TargetName #convert names to miR format
  
  accession_codes <- GSE117064[[1]]$geo_accession[non_cvd_patients]
  write.csv(accession_codes, "accession_codes.csv")
  
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
  binary_vars <- colnames(X)%in%c("Male")
  
  X_raw <- do.call(cbind, list(Age=age,
                           BMI=bmi,
                           MAP=map,
                           HbA1c=hba1c,
                           Male=male,
                           Smoking=smoking,
                           X))
  colnames(X_raw) <- gsub("hsa-", "", colnames(X_raw))
  
  #transform data with nonparanormal truncated ecdf
  X <- huge.npn(X_raw, npn.func = "truncation")
  
}

#cardiometabolic factor distributions
if(T){
  
  #obesity prevalence test
  obesity_test <- prop.test(x=sum(bmi>30),
                            n=nrow(X),
                            p=0.39)
  
  #generate histograms
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
         filename = paste0(figures_folder,"/all_pheno_plots.jpg"),
         dpi=1500,
         units="mm",
         width = 250,
         height = 125)
  

}

#infer graphical model
if(T){
  
  #measure bivariate correlation matrix
  S <- cor(X)
  
  #de-sparsified nodewise scaled lasso
  alpha_level <- 0.05
  h <- SILGGM(X, method="D-S_NW_SL",
              global = TRUE,
              alpha=alpha_level)
  am <- h$global_decision[[1]]
  dimnames(am)<- list(colnames(X), colnames(X))
  
  #extract graph
  g <- graph_from_adjacency_matrix(am, mode="undirected")
  pheno_vars <- 1:6
  
  #measure edge betweenness centrality
  E(g)$centrality <- rank(edge_betweenness(g))/ecount(g)
  h$centrality <- matrix(0, ncol(X), ncol(X),
                         dimnames = list(colnames(X), colnames(X)))
  g_df <- as_data_frame(g)
  h$centrality[cbind(g_df[,1], g_df[,2])] <- g_df[,3]
  h$centrality <- h$centrality + t(h$centrality)

  #adjust p values with false discovery rate control
  h$p_precision_adj <- h$p_precision
  h$p_precision_adj[upper.tri(h$p_precision_adj)] <-
    p.adjust(h$p_precision_adj[upper.tri(h$p_precision_adj)], method="fdr")
  h$p_precision_adj[lower.tri(h$p_precision_adj)] <- 0
  h$p_precision_adj <- h$p_precision_adj + t(h$p_precision_adj)
  
  #record correlation types
  binary_vars <- colnames(X)%in% "Male"
  cor_type_mat <- S
  cor_type_mat[]<- "Pearson"
  cor_type_mat[!binary_vars,binary_vars]<-
  cor_type_mat[binary_vars,!binary_vars] <- "Point-Biserial"
}

#bootstrap analysis
if(T){
  
  signed_am <- am*sign(h$partialCor)
  
  n_bs <- 500
  bs_am_list <- list()
  bs_pcor_list <- list()
  for(bs_i in 1:n_bs){
    message(bs_i)
    
    #resample X with replacement
    X_bs <- X_raw[sample(nrow(X_raw), replace=T),]
    
    #nonparanormal tranformaiton
    X_bs <- huge.npn(X_bs, npn.func = "truncation")
    
    #de-sparsified nodewise scaled lasso
    h_bs <- SILGGM(X_bs, method="D-S_NW_SL",
                   global = TRUE,
                   alpha=alpha_level)
    
    #record bootstrapped graph estimates
    bs_pcor_list[[bs_i]] <- h_bs$partialCor
    bs_am_list[[bs_i]]<- h_bs$global_decision[[1]]*sign(h_bs$partialCor)
  }
  
}

#utility functions
if(T){
  
  #interquartile range formatting
  iqr_func <- function(vec, digits=1){
    quartiles <- round(quantile(vec, c(0.025,0.975)),digits=digits)
    paste0("(Interquartile Range=", quartiles[1], "-", quartiles[2], ")")
  }
  
  #median formatting
  med_func <- function(vec){
    round(median(vec),1)
  }
  
  #rounding
  rounding_func <- function(i){
    format(round(i, digits=2), nsmall=2)
  }
  
  #format p values
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
  
  #function to make graph data frame
  make_g_df <- function(g_i){
    
    g_df <- as_data_frame(g_i)
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
    g_df$centrality <- h$centrality[e_idx]
    
    #bootstrapped results
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
    
    return(g_df)
    
  }
  
  #function to format graph data frame
  format_g_df <- function(g_df){
    
    g_df2 <- g_df
    g_df2$"Partial Corr." <-  paste0(
      rounding_func(g_df2$partial_cor), " [",
      rounding_func(g_df2$partial_cor_lo), ",",
      rounding_func(g_df2$partial_cor_hi), "]"
    )
    g_df2$"P-value (adjusted)" <- p_rounding_func(g_df2$p)
    colnames(g_df2)[1:2]<- c("Variable A", "Variable B")
    g_df2$"Bivariate Corr." <- rounding_func( g_df2$"Bivariate Corr.")
    g_df2$Centrality <- round(g_df2$centrality,2)
    g_df2 <- g_df2[,c("Variable A","Variable B","Corr. Type",
                      "Bivariate Corr.", "Partial Corr.",
                      "P-value (adjusted)", "Stability", "Centrality")]
    g_df2 <- g_df2[order(g_df2[,1], g_df2[,2]),]
    
    return(g_df2)
  }
  
}

#extract ego graph of phenotypic variables
if(T){
  
  #examine miRNA which interact directly with cardiometabolic factors
  ego_rf <- ego(g, order=1, nodes = pheno_vars)
  names(ego_rf)<- colnames(X)[pheno_vars]
  
  included_nodes <- unique(unlist(lapply(ego_rf, as.vector)))
  g_rf <- induced.subgraph(g, vids=included_nodes)
  g_df <- make_g_df(g_rf)

}

#plot ego graph
if(T){
  
  g_rf <- graph_from_data_frame(g_df[,c("from", "to", "Interaction")], directed=F)
  g_rf_ggnetwork_df <- ggnetwork(g_rf)
  myColors <- c("magenta", "springgreen")
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
    theme(plot.margin=unit(c(1,1,1,1),"cm"),
          legend.position = "left")
  
  ggsave(g_rf_plot,
         filename= paste0(figures_folder,"/g_rf_plot.jpg"),
         dpi=1500,
         units="in",
         width=7,
         height=7
         )
  
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
    geom_line(size=1, alpha = 0.7)+
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
  
  network_and_bs_plot <- ggarrange(plotlist=list(g_rf_plot, bs_plot), widths=c(3,1))
  ggsave(plot=network_and_bs_plot,
         filename = paste0(figures_folder,"/network_and_bs_plot.jpg"),
         dpi=1500,
         units="in",
         width = 9,
         height = 7.5)
}

#generating result outputs
if(T){
  
  #make induced subgraph of all vertices interacting with phenotypic variables
  g_df <- make_g_df(g_rf)
  g_df_formatted <- format_g_df(g_df)
  
  #interactions between cardiometabolic factors
  g_df_pheno <- g_df_formatted[g_df_formatted$`Variable A`%in%colnames(X)[pheno_vars]&
                                 g_df_formatted$`Variable B`%in%colnames(X)[pheno_vars],
  ]
  
  write.table(g_df_pheno, paste0(tables_folder,
                                 "/g_df_pheno.tsv"),
              sep="\t",
              row.names = F)
  
  #interactions between cardiometabolic factors and miRNA
  g_df_pheno_mir <-  g_df_formatted[xor(g_df_formatted$`Variable A`%in%colnames(X)[pheno_vars],
                                        g_df_formatted$`Variable B`%in%colnames(X)[pheno_vars]),
  ]
  write.table(g_df_pheno_mir, paste0(tables_folder,
                                     "/g_df_pheno_mir.tsv"),
              sep="\t",
              row.names = F)
  
  #full network table for supplementary data
  g_all_df <- make_g_df(g)
  g_all_df <- format_g_df(g_all_df)
  
  write.table(g_all_df, paste0(tables_folder,
                               "/g_all_df.tsv"),
              sep="\t",
              row.names = F)
  
 
}



