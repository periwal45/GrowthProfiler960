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

#read plate
data<-data.frame(read.table(file.names, header = TRUE, sep = '\t', stringsAsFactors = FALSE))
head(data)
colnames(data)[1]<-"time"
data[,1]<-data[,1]/60

#create an empty data frame to store all raw points and fitted points
fitted_values<-data.frame(stringsAsFactors = FALSE)

#create an empty data frame to store all computed model parameters
gc_fit_params<-data.frame(stringsAsFactors = FALSE)

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
    
    #add annotation to each row
    if(str_detect(col_name, "A")){
      
      bug = "ED1A"
      
      if(col_name %in% c("A1","A2","A3","A4","A5","A6")){
        
        comp = "bug-control"
      }else{
        
        comp = "sucrose"
      }

    }else if(str_detect(col_name, "B")){
      
      bug = "IAI1"
      
      if(col_name %in% c("B1","B2","B3","B4","B5","B6")){
        
        comp = "bug-control"
      }else{
        
        comp = "sucrose"
      }
      
    
    }else if(str_detect(col_name, "C")){
      
      bug = "lactis17"
      
      if(col_name %in% c("C1","C2","C3","C4","C5","C6")){
        
        comp = "bug-control"
      }else{
        
        comp = "sucrose"
      }
      
    
    }else if(str_detect(col_name, "D")){
      
      bug = "lactis261"
      
      if(col_name %in% c("D1","D2","D3","D4","D5","D6")){
        
        comp = "bug-control"
      }else{
        
        comp = "sucrose"
      }
      
    
    }else if(str_detect(col_name, "E")){
      
      bug = "rothia"
      
      if(col_name %in% c("E1","E2","E3","E4","E5","E6")){
        
        comp = "bug-control"
      }else{
        
        comp = "sucrose"
      }
      
    
    }else if(str_detect(col_name, "F")){
      
      bug = "luteus"
      
      if(col_name %in% c("F1","F2","F3","F4","F5","F6")){
        
        comp = "bug-control"
      }else{
        
        comp = "sucrose"
      }
      
    
    }else if(str_detect(col_name, "G")){
      
      bug = "infantis"
      
      if(col_name %in% c("G1","G2","G3","G4","G5","G6")){
        
        comp = "bug-control"
      }else{
        
        comp = "sucrose"
      }
      
    
    }else if(str_detect(col_name, "H")){
      
      bug = "simplex"
      
      if(col_name %in% c("H1","H2","H3","H4","H5","H6")){
        
        comp = "bug-control"
      }else{
        
        comp = "sucrose"
      }
      
    
    }
    
    annot_mod<-cbind(bug,comp,col_name,mod)
    annot_mod
    
    fitted_values<-rbind(fitted_values, annot_mod)
    
    fitted_param<-as.data.frame(cbind(bug,comp,col_name,gc_fit$vals$k,gc_fit$vals$k_se,gc_fit$vals$k_p,gc_fit$vals$n0,gc_fit$vals$n0_se,gc_fit$vals$n0_p,gc_fit$vals$r,gc_fit$vals$r_se,gc_fit$vals$r_p,gc_fit$vals$sigma,gc_fit$vals$df,gc_fit$vals$t_mid,gc_fit$vals$t_gen,gc_fit$vals$auc_l,gc_fit$vals$auc_e,gc_fit$vals$note))
    colnames(fitted_param)<-c("bug_name","comp","well","k","k_se","k_p","n0","n0_se","n0_p","r","r_se","r_p","sigma","df","t_mid","t_gen","auc_l","auc_e","note")
    fitted_param
    
    gc_fit_params<-rbind(gc_fit_params, fitted_param)
    tail(gc_fit_params)
    
  }
}

View(fitted_values)
nrow(fitted_values) #4512

mean_fit<-fitted_values %>% dplyr::group_by(bug,comp,time) %>% summarise(avOD = mean(OD), avfit = mean(fitted))
View(mean_fit)

View(gc_fit_params)
nrow(gc_fit_params) #should be equal to number of wells: 96 for 1 plate

# set global theme for all plots

th<-theme(plot.title = element_text(size = 12, face = "bold"),axis.title=element_text(size=12,color = "black"),
          axis.text.x = element_text(size=10, color = "black"),axis.text.y = element_text(size=10, color = "black"))

CairoSVG(file="/Users/periwal/ShikiFactory/WP3/GrowthProfiler/sucrose.svg", width = 5, height = 4, bg = "white")
fitted_values %>% dplyr::group_by(bug,comp,time) %>% summarise(OD = mean(OD), fitted = mean(fitted)) %>% 
  dplyr::group_by(bug,comp) %>%
ggplot(aes(x=time, y=fitted)) + geom_line(aes(color=comp)) + 
  geom_point(aes(y=OD,color=comp), size=0.05) + scale_color_jama() +
  scale_x_continuous(name = "Time (Hours)", limits = c(0,24)) + scale_y_continuous(name="bg corrected OD",) +
  th + theme_bw(base_rect_size = 0.1) + theme(legend.position="bottom") + facet_wrap("bug", scales = "free", ncol = 4)
dev.off()
