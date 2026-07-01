#install.packages('DescTools')
#install.packages('hash')
library(DescTools)
library(hash)
stopifnot(packageVersion('DescTools')>='0.99.60')
stopifnot(packageVersion('hash')>='2.2.6.4')

# Set to script directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source('plot_funcs.R')

optim_folder <- '../simulations_SANN'

wd <- paste0('./',optim_folder,'/figures_supp/to_etal_supp/')

##########################################################################
##################### Power Curves #######################################
##########################################################################

# Test statistic names
#'mean','ratio','tstat','cv','KST','MWWT','DT'

plot_plots <- c('mean','ratio','tstat','cv','MWWT','DT')
for (i in plot_plots) {
  test_stat <- i
  
  png(paste0(wd,'power_',test_stat,'.png'), width = 8000*0.3, height = 6000*0.3,res = 300)
  
  layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  
  # Define the heights and widths for layouts
  heights <- c(1, 1)  # Equal heights for rows
  widths <- c(1, 1)   # Equal widths for columns
  
  # Create a layout with zero inner spacing
  layout(layout_matrix, widths = widths, heights = heights)
  
  # Define outer margins, which affect space around the entire plot area
  par(oma = c(1, 1, 1, 1))  # No additional outer margins
  par(pty = 's')
  par(mfrow = c(2, 2))
  par(font.axis = 2, font.lab = 2)
  
  #par(mar = c(2.5, 4, 4, 2))  # Top space to reduce proximity to next plot
  par(mar = c(2.5, 3, 4, 3))  # Top space to reduce proximity to next plot
  
  #####################
  #####################
  
  dir <- paste0('./',optim_folder,'/to_etal_simulations/PowerAnalysis/RMC_750_3_25/kscale_highshape/')
  load(list.files(dir,full.names=TRUE)[grepl('pwrarrall',basename(list.files(dir,full.names=TRUE)),ignore.case=TRUE)])
  
  #par(mar = c(5, 5, 5, 8)) # 5,5,5,7
  
  power_to_plot <- power_select(test_stat)[[1]]
  test_type <- power_select(test_stat)[[2]]
  CI <- BinomCI(pwr_calc_array[,test_type,1], pwr_calc_array[,test_type,2])
  CI_low <- CI[,2]
  CI_high <- CI[,3]
  
  eff_sizes = klist
  
  ai = 1; bi = length(eff_sizes)#-1
  plot(eff_sizes[ai:bi],(power_to_plot[ai:bi]),yaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
       xlab='',#Effect Size',
       ylab='Power',main='750_3_25',cex.main=1.5,cex.lab=1.3)
  axis(side = 2,las = 2) 
  points(eff_sizes[ai:bi],(CI_low[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  points(eff_sizes[ai:bi],(CI_high[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  abline(h = 0.8,lty='dashed')
  
  #par(mar = c(2.5, 2, 4, 4))  # Top space to reduce proximity to next plot
  
  #####################
  #####################
  
  dir <- paste0('./',optim_folder,'/to_etal_simulations/PowerAnalysis/RMC_750_11_10/kscale_highshape/')
  load(list.files(dir,full.names=TRUE)[grepl('pwrarrall',basename(list.files(dir,full.names=TRUE)),ignore.case=TRUE)])
  
  # par(mar = c(5, 5, 5, 7)) # 5,5,5,7
  
  power_to_plot <- power_select(test_stat)[[1]]
  test_type <- power_select(test_stat)[[2]]
  CI <- BinomCI(pwr_calc_array[,test_type,1], pwr_calc_array[,test_type,2])
  CI_low <- CI[,2]
  CI_high <- CI[,3]
  
  eff_sizes = klist
  
  ai = 1; bi = length(eff_sizes)#-1
  plot(eff_sizes[ai:bi],(power_to_plot[ai:bi]),yaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
       xlab='',#Effect Size',
       ylab='',#Power',
       main='750_11_10',cex.main=1.5,cex.lab=1.3)
  axis(side = 2,las = 2) 
  points(eff_sizes[ai:bi],(CI_low[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  points(eff_sizes[ai:bi],(CI_high[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  abline(h = 0.8,lty='dashed')
  
  #par(mar = c(5, 4, 2, 2))
  par(mar = c(4, 3, 2.5, 3)) 
  
  #####################
  #####################
  
  dir <- paste0('./',optim_folder,'/to_etal_simulations/PowerAnalysis/RMC_995_3_25/kscale_highshape/')
  load(list.files(dir,full.names=TRUE)[grepl('pwrarrall',basename(list.files(dir,full.names=TRUE)),ignore.case=TRUE)])
  
  # par(mar = c(5, 5, 5, 7)) # 5,5,5,7
  
  power_to_plot <- power_select(test_stat)[[1]]
  test_type <- power_select(test_stat)[[2]]
  CI <- BinomCI(pwr_calc_array[,test_type,1], pwr_calc_array[,test_type,2])
  CI_low <- CI[,2]
  CI_high <- CI[,3]
  
  eff_sizes = klist
  
  ai = 1; bi = length(eff_sizes)#-1
  plot(eff_sizes[ai:bi],(power_to_plot[ai:bi]),yaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
       xlab='Effect Size',
       ylab='Power',main='995_3_25',cex.main=1.5,cex.lab=1.3)
  axis(side = 2,las = 2) 
  points(eff_sizes[ai:bi],(CI_low[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  points(eff_sizes[ai:bi],(CI_high[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  abline(h = 0.8,lty='dashed')
  
  #par(mar = c(5, 2, 2, 4))
  
  #####################
  #####################
  
  dir <- paste0('./',optim_folder,'/to_etal_simulations/PowerAnalysis/RMC_995_11_10/kscale_highshape/')
  load(list.files(dir,full.names=TRUE)[grepl('pwrarrall',basename(list.files(dir,full.names=TRUE)),ignore.case=TRUE)])
  
  # par(mar = c(5, 5, 5, 7)) # 5,5,5,7
  
  power_to_plot <- power_select(test_stat)[[1]]
  test_type <- power_select(test_stat)[[2]]
  CI <- BinomCI(pwr_calc_array[,test_type,1], pwr_calc_array[,test_type,2])
  CI_low <- CI[,2]
  CI_high <- CI[,3]
  
  eff_sizes = klist
  
  ai = 1; bi = length(eff_sizes)#-1
  plot(eff_sizes[ai:bi],(power_to_plot[ai:bi]),yaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
       xlab='Effect Size',
       ylab='',#Power',
       main='995_11_10',cex.main=1.5,cex.lab=1.3)
  axis(side = 2,las = 2) 
  points(eff_sizes[ai:bi],(CI_low[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  points(eff_sizes[ai:bi],(CI_high[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  abline(h = 0.8,lty='dashed')
  
  dev.off()
  
}

#####################
#####################

png(paste0(wd,'low_powers.png'), width = 8000*0.3, height = 6000*0.3,res = 300)

layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)

# Define the heights and widths for layouts
heights <- c(1, 1)  # Equal heights for rows
widths <- c(1, 1)   # Equal widths for columns

# Create a layout with zero inner spacing
layout(layout_matrix, widths = widths, heights = heights)

# Define outer margins, which affect space around the entire plot area
par(oma = c(1, 1, 1, 1))  # No additional outer margins
par(pty = 's')
par(mfrow = c(2, 2))
par(font.axis = 2, font.lab = 2)

#par(mar = c(2.5, 4, 4, 2))  # Top space to reduce proximity to next plot
par(mar = c(2.5, 3, 4, 3))  #Top space to reduce proximity to next plot

#####################
#####################

dir <- paste0('./',optim_folder,'/to_etal_simulations/PowerAnalysis/RMC_750_3_25/kscale_lowshape')
load(list.files(dir,full.names=TRUE)[grepl('pwrarrall',basename(list.files(dir,full.names=TRUE)),ignore.case=TRUE)])

eff_sizes = klist

ai = 1; bi = length(eff_sizes)#-1
plot(eff_sizes[ai:bi],(pwr_arr_mean[ai:bi]),yaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
     xlab='',
     ylab='Power',main='750_3_25',cex.main=1.5,cex.lab=1.3)
axis(side = 2,las = 2) 
points(eff_sizes[ai:bi],(pwr_arr_KST[ai:bi]),ylim=c(0,1),type='o',pch=20,col='orange')
points(eff_sizes[ai:bi],(pwr_arr_MWWT[ai:bi]),ylim=c(0,1),type='o',pch=20,col='brown')
points(eff_sizes[ai:bi],(pwr_arr_tstat[ai:bi]),ylim=c(0,1),type='o',pch=20,col='blue')
points(eff_sizes[ai:bi],(pwr_arr_ratio[ai:bi]),ylim=c(0,1),type='o',pch=20,col='green')
points(eff_sizes[ai:bi],(pwr_arr_cv[ai:bi]),ylim=c(0,1),type='o',pch=20,col='red')
points(eff_sizes[ai:bi],(pwr_arr_DT[ai:bi]),ylim=c(0,1),type='o',pch=20,col='purple')
abline(h = 0.8,lty='dashed')
#legend('topright', inset=c(-0.3,0), xpd=TRUE,legend=c('Mean-Difference','KS Test','MW Test','t-test','Mean-Ratio','CV','Dip Test'),col=c('black','orange','brown','blue','green','red','purple'),lty=1:1, cex=0.7)


#par(mar = c(2.5, 2, 4, 4))  # Top space to reduce proximity to next plot

#####################
#####################

dir <- paste0('./',optim_folder,'/to_etal_simulations/PowerAnalysis/RMC_750_11_10/kscale_lowshape')
load(list.files(dir,full.names=TRUE)[grepl('pwrarrall',basename(list.files(dir,full.names=TRUE)),ignore.case=TRUE)])

eff_sizes = klist

ai = 1; bi = length(eff_sizes)#-1
plot(eff_sizes[ai:bi],(pwr_arr_mean[ai:bi]),yaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
     xlab='',
     ylab='Power',main='750_11_10',cex.main=1.5,cex.lab=1.3)
axis(side = 2,las = 2) 
points(eff_sizes[ai:bi],(pwr_arr_KST[ai:bi]),ylim=c(0,1),type='o',pch=20,col='orange')
points(eff_sizes[ai:bi],(pwr_arr_MWWT[ai:bi]),ylim=c(0,1),type='o',pch=20,col='brown')
points(eff_sizes[ai:bi],(pwr_arr_tstat[ai:bi]),ylim=c(0,1),type='o',pch=20,col='blue')
points(eff_sizes[ai:bi],(pwr_arr_ratio[ai:bi]),ylim=c(0,1),type='o',pch=20,col='green')
points(eff_sizes[ai:bi],(pwr_arr_cv[ai:bi]),ylim=c(0,1),type='o',pch=20,col='red')
points(eff_sizes[ai:bi],(pwr_arr_DT[ai:bi]),ylim=c(0,1),type='o',pch=20,col='purple')
abline(h = 0.8,lty='dashed')
#legend('topright', inset=c(-0.3,0), xpd=TRUE,legend=c('Mean-Difference','KS Test','MW Test','t-test','Mean-Ratio','CV','Dip Test'),col=c('black','orange','brown','blue','green','red','purple'),lty=1:1, cex=0.7)

#par(mar = c(5, 4, 2, 2))
par(mar = c(4, 3, 2.5, 3)) 

#####################
#####################

dir <- paste0('./',optim_folder,'/to_etal_simulations/PowerAnalysis/RMC_995_3_25/kscale_lowshape')
load(list.files(dir,full.names=TRUE)[grepl('pwrarrall',basename(list.files(dir,full.names=TRUE)),ignore.case=TRUE)])

eff_sizes = klist

ai = 1; bi = length(eff_sizes)#-1
plot(eff_sizes[ai:bi],(pwr_arr_mean[ai:bi]),yaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
     xlab='Effect Size',
     ylab='Power',main='995_3_25',cex.main=1.5,cex.lab=1.3)
axis(side = 2,las = 2) 
points(eff_sizes[ai:bi],(pwr_arr_KST[ai:bi]),ylim=c(0,1),type='o',pch=20,col='orange')
points(eff_sizes[ai:bi],(pwr_arr_MWWT[ai:bi]),ylim=c(0,1),type='o',pch=20,col='brown')
points(eff_sizes[ai:bi],(pwr_arr_tstat[ai:bi]),ylim=c(0,1),type='o',pch=20,col='blue')
points(eff_sizes[ai:bi],(pwr_arr_ratio[ai:bi]),ylim=c(0,1),type='o',pch=20,col='green')
points(eff_sizes[ai:bi],(pwr_arr_cv[ai:bi]),ylim=c(0,1),type='o',pch=20,col='red')
points(eff_sizes[ai:bi],(pwr_arr_DT[ai:bi]),ylim=c(0,1),type='o',pch=20,col='purple')
abline(h = 0.8,lty='dashed')
#legend('topright', inset=c(-0.3,0), xpd=TRUE,legend=c('Mean-Difference','KS Test','MW Test','t-test','Mean-Ratio','CV','Dip Test'),col=c('black','orange','brown','blue','green','red','purple'),lty=1:1, cex=0.7)

#par(mar = c(5, 2, 2, 4))

#####################
#####################

dir <- paste0('./',optim_folder,'/to_etal_simulations/PowerAnalysis/RMC_995_11_10/kscale_lowshape')
load(list.files(dir,full.names=TRUE)[grepl('pwrarrall',basename(list.files(dir,full.names=TRUE)),ignore.case=TRUE)])

eff_sizes = klist

ai = 1; bi = length(eff_sizes)#-1
plot(eff_sizes[ai:bi],(pwr_arr_mean[ai:bi]),yaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
     xlab='Effect Size',
     ylab='Power',main='995_11_10',cex.main=1.5,cex.lab=1.3)
axis(side = 2,las = 2) 
points(eff_sizes[ai:bi],(pwr_arr_KST[ai:bi]),ylim=c(0,1),type='o',pch=20,col='orange')
points(eff_sizes[ai:bi],(pwr_arr_MWWT[ai:bi]),ylim=c(0,1),type='o',pch=20,col='brown')
points(eff_sizes[ai:bi],(pwr_arr_tstat[ai:bi]),ylim=c(0,1),type='o',pch=20,col='blue')
points(eff_sizes[ai:bi],(pwr_arr_ratio[ai:bi]),ylim=c(0,1),type='o',pch=20,col='green')
points(eff_sizes[ai:bi],(pwr_arr_cv[ai:bi]),ylim=c(0,1),type='o',pch=20,col='red')
points(eff_sizes[ai:bi],(pwr_arr_DT[ai:bi]),ylim=c(0,1),type='o',pch=20,col='purple')
abline(h = 0.8,lty='dashed')

dev.off()

#####################
#####################

png(paste0(wd,'low_powers_legend.png'), width = 8000*0.3, height = 6000*0.3,res = 300)

layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)

# Define the heights and widths for layouts
heights <- c(1, 1)  # Equal heights for rows
widths <- c(1, 1)   # Equal widths for columns

# Create a layout with zero inner spacing
layout(layout_matrix, widths = widths, heights = heights)

# Define outer margins, which affect space around the entire plot area
par(oma = c(1, 1, 1, 1))  # No additional outer margins

par(mfrow = c(2, 2))
par(font.axis = 2, font.lab = 2)

par(mar = c(2.5, 4, 4, 2))  # Top space to reduce proximity to next plot

plot(eff_sizes[ai:bi],(pwr_arr_mean[ai:bi]),bty='n',yaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='white',
     xlab='',
     ylab='Power',main='750_3_25',
     cex.main=1.5,cex.lab=1.3)
legend('topright', inset=c(-0.1,0), xpd=TRUE,legend=c('Mean-Difference','KS Test','MW Test','t-test','Mean-Ratio','CV','Dip Test'),col=c('black','orange','brown','blue','green','red','purple'),lty=1:1, cex=0.7)

dev.off()
