#Cardiovascular Risk Factors and the miRNAome: a Gaussian Graphical Model Analysis

#load libraries
if(T){
  library(GEOquery)
  library(igraph)
  library(miRBaseConverter)
  library(ggnetwork)
  library(ggplot2)
  library(ggpubr)
  library(GGMncv)
  library(polycor)
  library(huge)
}

#set output folders
if(T){
  figures_folder <- "vrf_figures"
  tables_folder <- "vrf_tables"
}

#download data from Gene Expression Omnibus
if(T){
  GSE117064 <- getGEO("GSE117064",GSEMatrix=TRUE)
  non_cvd_patients <- !GSE117064[[1]]$`group:ch1`%in%c("3A", "3B")
  X <- t(exprs(GSE117064[[1]]))[non_cvd_patients,]
  
  colnames(X) <- gsub(",.*", "", colnames(X))
  colnames(X) <- miRNA_AccessionToName(colnames(X),
                                       targetVersion = "v21")$TargetName
  
  #read in mirtarbase data
  mirtarbase_df <- read.csv("hsa_MTI.csv")
  mirtarbase_df <- mirtarbase_df[-grep("Weak", mirtarbase_df$Support.Type),]
  mirtarbase_df <- mirtarbase_df[mirtarbase_df$miRNA%in% colnames(X),]
  
  #extract most important mirna
  top_100_mirna <- names(sort(table(mirtarbase_df$miRNA), decreasing = T)[1:100])
  X <- X[,top_100_mirna]
  mirtarbase_df <- mirtarbase_df[mirtarbase_df$miRNA%in%top_100_mirna,]
  
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
  X[,!binary_vars] <- huge.npn(as.matrix(X[,!binary_vars]))
  
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
  
}

#reporting functions
if(T){
  #interquartile range
  iqr_func <- function(vec, digits=1){
    quartiles <- round(quantile(vec, c(0.025,0.975)),digits=digits)
    paste0("(Interquartile Range=", quartiles[1], "-", quartiles[2], ")")
  }
  #median
  med_func <- function(vec){
    round(median(vec),1)
  }
}

#infer graphical model
if(T){
  
  #infer pearson correlations
  S <- cor(X)
  
  #infer biserial correlations between binary variables and continuous variables
  for(i in which(!binary_vars)){
    for(j in which(binary_vars)){
      S[i,j]<- S[j,i]<- polyserial(X[,i], X[,j], ML=TRUE)
    }
  }
  
  #record correlation types
  cor_type_mat <- S
  cor_type_mat[]<- "Pearson"
  cor_type_mat[!binary_vars,binary_vars]<-
    cor_type_mat[binary_vars,!binary_vars] <- "Biserial"
  
  
  #infer gaussian graphical model with the graphical lasso
  ggm_mod <- GGMncv(x=S,
                    n=nrow(X),
                    penalty="lasso",
                    select=F)
  
  #estimate desparsified partial correlation matrix with Jankova's method
  desparsified_mod <- desparsify(ggm_mod)
  
  #perform false discovery rate controlled inference on graphical model
  inferred_mod <- inference(ggm_mod,
                            alpha=0.05)
  
  #generate adjacency matrix from significant interactions
  am <- inferred_mod$adj
  rownames(am)<- colnames(am)<- colnames(X)
  pheno_vars <- 1:6
  g <- graph_from_adjacency_matrix(am,
                                   mode="undirected")
  V(g)$name <- gsub("hsa-", "", V(g)$name)
  
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
  
  g_df$partial_cor <- desparsified_mod$P[e_idx]
  g_df$Interaction <- ifelse(g_df$partial_cor>0, "Upregulation", "Inhibition")
  g_df$"Bivariate Corr." <- S[e_idx]
  g_df$"Corr. Type" <- cor_type_mat[e_idx]
  
  g_rf <- graph_from_data_frame(g_df[,c("from", "to", "Interaction")], directed=F)
}

#plot ego graph
if(T){
  
  g_rf_ggnetwork_df <- ggnetwork(g_rf)
  myColors <- c("coral","springgreen")
  names(myColors) <- levels(g_rf_ggnetwork_df$Interaction)
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
        i_i <- substr(i_i, start = 1, stop = nchar(i_i))
      }else{
        i_i <- rounding_func(i_i)
        i_i <- substr(i_i, start = 1, stop = nchar(i_i))
      }
      return(i_i)
    })
    
  }
  
  #produce result tables
  g_df2 <- g_df
  g_df2$"Partial Corr." <- rounding_func(g_df2$partial_cor)
  p_adjusted <- inferred_mod$corrected
  g_df2$"P-value (adjusted)" <- p_rounding_func(p_adjusted[e_idx])
  colnames(g_df2)[1:2]<- c("Variable A", "Variable B")
  g_df2$"Bivariate Corr." <- rounding_func( g_df2$"Bivariate Corr.")
  g_df2 <- g_df2[,c("Variable A","Variable B","Corr. Type",
                    "Bivariate Corr.", "Partial Corr.",
                    "P-value (adjusted)", "Interaction")]
  g_df2 <- g_df2[order(g_df2[,1], g_df2[,2]),]
  
  write.table(g_df2, paste0(tables_folder,
                            "/result_table.tsv"),
              sep="\t",
              row.names = F)
  
  #full network table for supplementary data
  g_all_df <- as_data_frame(g)
  e_idx <- cbind(match(g_all_df$from,  V(g)$name),
                 match(g_all_df$to,   V(g)$name)
  )
  
  g_all_df$partial_cor <- desparsified_mod$P[e_idx]
  g_all_df$Interaction <- ifelse(g_all_df$partial_cor>0, "Upregulation", "Inhibition")
  g_all_df$"Partial Corr." <- rounding_func(g_all_df$partial_cor)
  p_adjusted <- inferred_mod$corrected
  g_all_df$"P-value (adjusted)" <- p_rounding_func(p_adjusted[e_idx])
  g_all_df$"Bivariate Corr." <- rounding_func(S[e_idx])
  g_all_df$"Corr. Type" <- cor_type_mat[e_idx]
  
  colnames(g_all_df)[1:2]<- c("Variable A", "Variable B")
  
  
  g_all_df <- g_all_df[,c("Variable A", "Variable B","Corr. Type",
                          "Bivariate Corr.","Partial Corr.",
                          "P-value (adjusted)", "Interaction")]
  g_all_df <- g_all_df[order(g_all_df[,1]),]
  
  write.table(g_all_df, paste0(tables_folder,
                               "/result_table_all.tsv"),
              sep="\t",
              row.names = F)
}



#save.image("mirna_rf_network.RData")
#load("mirna_rf_network.RData")
