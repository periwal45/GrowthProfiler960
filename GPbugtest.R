setwd("/Users/periwal/ShikiFactory/WP3/GrowthProfiler/")
library(dplyr)
library(tibble)
library(reshape2)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(ggsci) #fancy color palettes
library(Cairo) #high resolution images
library(pROC) #auc function in growthcurver
library(growthcurver)
library(pROC) #auc function in growthcurver
library(growthcurver)
library(MASS) #fitdistr
library(gtools) #mixedsort
library(qvalue) #estimates FDR from a list of input p-values
library(colr) #regex on columns


#read all files with .tab extension
file.names <- list.files(path = '/Users/periwal/ShikiFactory/WP3/GrowthProfiler/', recursive = TRUE, pattern = "SF100_all.*\\.tab$") #recursive reads through all subfolders and files
file.names

annot<-data.frame(read.table(file = "OD/SF100/all_bugs_growthtest/plate_layout", sep = '\t', header = TRUE))
annot<-annot %>% melt(measure.vars = c("rep1","rep2"), value.name = "col_name")
View(annot)

#create an empty data frame to store all raw points and fitted points
fitted_values<-data.frame(stringsAsFactors = FALSE)

#create an empty data frame to store all computed model parameters
gc_fit_params<-data.frame(stringsAsFactors = FALSE)


#read plates
for(i in 1:length(file.names)){
    
  data<-data.frame(read.table(file.names[i], header = TRUE, sep = '\t', stringsAsFactors = FALSE))
  head(data)
  colnames(data)[1]<-"time"
  data[,1]<-data[,1]/60


  #fit model to each well of a plate and save the fitted model
  for(col_name in names(data)){
  
  col_name
  
  if(col_name != "time"){
    
    current_well<-data[, c("time",col_name)]
    current_well
    min_OD<-min(current_well[, col_name])
    min_OD
    
    #do background correction using min OD of each well
    current_well[,col_name]<-current_well[,col_name] - min_OD
    current_well[,col_name]
    
    #each time create a new variable for each well
    gc_fit<-SummarizeGrowth(data_t = current_well[,"time"], data_n = current_well[,col_name])
    #saveRDS(gc_fit, file = paste0("/Users/vinitaperiwal/GrowthCurver/models/",ID,"_",rep,"_",plate,"_",drug,"_",col_name,".rds"))
    
    gc_fit
    
    #create a data frame of raw values and fitted values
    mod_t<-data.frame(matrix(unlist(gc_fit$data$t)))
    mod_t
    mod_N<-data.frame(matrix(unlist(gc_fit$data$N)))
    mod_N
    
    if(gc_fit$vals$k != 0 & gc_fit$vals$n0 != 0 & gc_fit$vals$r != 0){
      
      mod<-cbind(mod_t, mod_N, gc_fit$model$m$fitted()) #m is the model object with all fitted and residual values
      
    }else{
      
      mod<-cbind(mod_t, mod_N, gc_fit$vals$r)
      
    }
    
    colnames(mod)<-c("time","OD","fitted")
    mod
    
    #add annotation to each well
    col_name

    annot_mod<-cbind(file.names[i],col_name,mod)
    annot_mod
    
    merged<-merge(annot_mod,annot,by="col_name")
    head(merged)
    merged<-merged[,c(6,7,1,3,4,5,8,2)]
    
    fitted_values<-rbind(fitted_values, merged)
    
    fitted_param<-as.data.frame(cbind(merged$bug,merged$strain,col_name,gc_fit$vals$k,gc_fit$vals$k_se,gc_fit$vals$k_p,gc_fit$vals$n0,gc_fit$vals$n0_se,gc_fit$vals$n0_p,gc_fit$vals$r,gc_fit$vals$r_se,gc_fit$vals$r_p,gc_fit$vals$sigma,gc_fit$vals$df,gc_fit$vals$t_mid,gc_fit$vals$t_gen,gc_fit$vals$auc_l,gc_fit$vals$auc_e,gc_fit$vals$note))
    colnames(fitted_param)<-c("bug","bug_id","well","k","k_se","k_p","n0","n0_se","n0_p","r","r_se","r_p","sigma","df","t_mid","t_gen","auc_l","auc_e","note")
    fitted_param
    
    gc_fit_params<-rbind(gc_fit_params, fitted_param)
    tail(gc_fit_params)
    
    
  }
  }
}

View(fitted_values)
nrow(fitted_values) #6912

mean_fit<-fitted_values %>% dplyr::group_by(bug,strain,time) %>% summarise(avOD = mean(OD), avfit = mean(fitted))
View(mean_fit)

View(gc_fit_params)
nrow(gc_fit_params) #should be equal to number of wells: 96 for 1 plate

# set global theme for all plots

th<-theme(plot.title = element_text(size = 12, face = "bold"),axis.title=element_text(size=12,color = "black"),
          axis.text.x = element_text(size=10, color = "black"),axis.text.y = element_text(size=10, color = "black"))

#CairoSVG(file="/Users/periwal/ShikiFactory/WP3/GrowthProfiler/sucrose.svg", width = 5, height = 4, bg = "white")
pdf(file = "/Users/periwal/ShikiFactory/WP3/GrowthProfiler/OD/SF100/all_bugs_growthtest/GP_bugstestgrowth.pdf", width = 18, height = 10)
fitted_values %>% dplyr::group_by(bug,strain,time) %>% summarise(OD = mean(OD), fitted = mean(fitted)) %>% 
  dplyr::group_by(bug,strain) %>% 
ggplot(aes(x=time, y=fitted)) + geom_line(aes(color=strain)) + 
  geom_point(aes(y=OD,color=strain), size=0.05) + 
  scale_x_continuous(name = "Time (Hours)", limits = c(0,24)) + scale_y_continuous(name="bg corrected OD",) +
  th + theme_bw(base_rect_size = 0.1) + theme(legend.position="right") + facet_wrap("bug", scales = "free", ncol = 8)
dev.off()
