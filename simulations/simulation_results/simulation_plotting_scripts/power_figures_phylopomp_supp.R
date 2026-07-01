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

wd <- paste0('./',optim_folder,'/figures_supp/phylopomp_supp/')

##########################################################################
##################### Power Curves #######################################
##########################################################################

# Test statistic names
#'mean','ratio','tstat','cv','KST','MWWT','DT'

plot_plots <- c('ratio','tstat','cv','MWWT','DT')
for (i in plot_plots) {
  test_stat <- i
  
  png(paste0(wd,'power_',test_stat,'.png'), width = 7, height = 3.5,units='in',res = 400)
  
  # Define the heights and widths for layouts
  heights <- c(1, 1)  # Equal heights for rows
  widths <- c(1, 1)   # Equal widths for columns
  
  # layout matrix
  layout_matrix <- matrix(c(1, 2), nrow = 1, byrow = TRUE)
  layout(layout_matrix, widths = widths, heights = heights)
  
  par(mfrow = c(1, 2))
  par(pty = 's')
  par(oma = c(1, 1, 1, 1)) # Outer margins
  par(font.axis = 2, font.lab = 2, cex.axis = 0.6, cex.lab = 0.8, cex.main = 1, cex = 1) # axes
  par(mar = c(3, 2.5, 3, 2.5)) # figure margins
  par(mgp = c(1.5,0.4,0)) # axis to ticks
  par(tcl = -0.2) # tick label size
  
  #####################
  #####################
  
  dir <- paste0('./',optim_folder,'/phylopomp_simulations/PowerAnalysis/kscale_highshape/')
  load(list.files(dir,full.names=TRUE)[grepl('pwrarrall',basename(list.files(dir,full.names=TRUE)),ignore.case=TRUE)])
  
  # par(mar = c(5, 5, 5, 7)) # 5,5,5,7
  
  power_to_plot <- power_select(test_stat)[[1]]
  test_type <- power_select(test_stat)[[2]]
  CI <- BinomCI(pwr_calc_array[,test_type,1], pwr_calc_array[,test_type,2])
  CI_low <- CI[,2]
  CI_high <- CI[,3]
  
  eff_sizes = klist
  
  ai = 1; bi = length(eff_sizes)#-1
  plot(eff_sizes[ai:bi],(power_to_plot[ai:bi]),yaxt='n',xaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
       xlab='Effect Size',ylab='Power',main='Full sampling')
  axis(side = 2,las = 2) 
  ticks_to_label <- c(1, 2,5,10,20,40)
  axis(1, at = ticks_to_label, labels = ticks_to_label,font=2)
  points(eff_sizes[ai:bi],(CI_low[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  points(eff_sizes[ai:bi],(CI_high[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  abline(h = 0.8,lty='dashed')
  
  #####################
  #####################
  
  dir <- paste0('./',optim_folder,'/phylopomp_simulations/PowerAnalysis/kscale_highshape_sampled_third/')
  load(list.files(dir,full.names=TRUE)[grepl('pwrarrall',basename(list.files(dir,full.names=TRUE)),ignore.case=TRUE)])
  
  # par(mar = c(5, 5, 5, 7)) # 5,5,5,7
  
  power_to_plot <- power_select(test_stat)[[1]]
  test_type <- power_select(test_stat)[[2]]
  CI <- BinomCI(pwr_calc_array[,test_type,1], pwr_calc_array[,test_type,2])
  CI_low <- CI[,2]
  CI_high <- CI[,3]
  
  eff_sizes = klist
  
  ai = 1; bi = length(eff_sizes)#-1
  plot(eff_sizes[ai:bi],(power_to_plot[ai:bi]),yaxt='n',xaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
       xlab='Effect Size',ylab='',main='Het. sampling (33%)' )
  axis(side = 2,las = 2)  
  ticks_to_label <- c(1, 2,5,10,20,40)
  axis(1, at = ticks_to_label, labels = ticks_to_label,font=2)
  points(eff_sizes[ai:bi],(CI_low[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  points(eff_sizes[ai:bi],(CI_high[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  abline(h = 0.8,lty='dashed')
  
  dev.off()
  
}

#####################
#####################

png(paste0(wd,'Figure6_lowshape.png'), width = 7, height = 14,units='in',res = 400)

# # Define the heights and widths for layouts
# heights <- c(1, 1)  # Equal heights for rows
# widths <- c(1, 1)   # Equal widths for columns

# layout matrix
#layout_matrix <- matrix(c(1, 2, 3, 4, 5), nrow = 5, byrow = TRUE)
#layout(layout_matrix, widths = widths, heights = heights)

par(mfrow = c(5, 1))
par(pty = 's')
par(oma = c(2, 2, 1, 1)) # Outer margins
par(font.axis = 2, font.lab = 2, cex.lab = 1.5, cex.main = 1.6, cex.axis=1.2) # axes
par(mar = c(3.5, 3, 1.5, 1)) # figure margins
par(mgp = c(2.5,0.4,0)) # axis to ticks
par(tcl = -0.2) # tick label size

dir <- paste0('./',optim_folder,'/phylopomp_simulations/PowerAnalysis/kscale_lowshape/')
load(list.files(dir,full.names=TRUE)[grepl('pwrarrall',basename(list.files(dir,full.names=TRUE)),ignore.case=TRUE)])

# Test statistics to plot
#'mean','ratio','cv','DT','KST','KST','MWWT','tstat'

plot_plots <- c('mean','ratio','tstat','KST','MWWT') #,'DT','cv')
plot_title <- hash(keys=plot_plots,values=c('Mean-Difference', 'Mean-Ratio', 't-statistic', 'KS Test', 'MW Test'))
for (i in plot_plots) {
  test_stat <- i
  
  # par(mar = c(5, 5, 5, 7)) # 5,5,5,7
  
  power_to_plot <- power_select(test_stat)[[1]]
  test_type <- power_select(test_stat)[[2]]
  CI <- BinomCI(pwr_calc_array[,test_type,1], pwr_calc_array[,test_type,2])
  CI_low <- CI[,2]
  CI_high <- CI[,3]
  
  eff_sizes = klist
  
  ai = 1; bi = length(eff_sizes)#-1
  if (i == tail(plot_plots,1)) {
  plot(eff_sizes[ai:bi],(power_to_plot[ai:bi]),yaxt='n',xaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
       xlab='Effect Size',ylab='Power',main=plot_title[[i]])
  } else {
    plot(eff_sizes[ai:bi],(power_to_plot[ai:bi]),yaxt='n',xaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
         xlab='',ylab='Power',main=plot_title[[i]])
  }
  axis(side = 2,las = 2) 
  ticks_to_label <- c(1, 2,5,10,20,40)
  axis(1, at = ticks_to_label, labels = ticks_to_label,font=2)
  points(eff_sizes[ai:bi],(CI_low[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  points(eff_sizes[ai:bi],(CI_high[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  abline(h = 0.8,lty='dashed')
  
}

dev.off()

# #####################
# #####################

png(paste0(wd,'Figure6_scaled.png'), width = 7, height = 14,units='in',res = 400)

# # Define the heights and widths for layouts
# heights <- c(1, 1)  # Equal heights for rows
# widths <- c(1, 1)   # Equal widths for columns

# layout matrix
#layout_matrix <- matrix(c(1, 2, 3, 4, 5), nrow = 5, byrow = TRUE)
#layout(layout_matrix, widths = widths, heights = heights)

par(mfrow = c(5, 1))
par(pty = 's')
par(oma = c(2, 2, 1, 1)) # Outer margins
par(font.axis = 2, font.lab = 2, cex.lab = 1.5, cex.main = 1.6, cex.axis=1.2) # axes
par(mar = c(3.5, 3, 1.5, 1)) # figure margins
par(mgp = c(2.5,0.4,0)) # axis to ticks
par(tcl = -0.2) # tick label size

dir <- paste0('./',optim_folder,'/phylopomp_simulations/PowerAnalysis/kscale_highshape_scaled/')
load(list.files(dir,full.names=TRUE)[grepl('pwrarrall',basename(list.files(dir,full.names=TRUE)),ignore.case=TRUE)])

# Test statistics to plot
#'mean','ratio','cv','DT','KST','KST','MWWT','tstat'

plot_plots <- c('mean','ratio','tstat','KST','MWWT') #,'DT','cv')
plot_title <- hash(keys=plot_plots,values=c('Mean-Difference', 'Mean-Ratio', 't-statistic', 'KS Test', 'MW Test'))
for (i in plot_plots) {
  test_stat <- i
  
  # par(mar = c(5, 5, 5, 7)) # 5,5,5,7
  
  power_to_plot <- power_select(test_stat)[[1]]
  test_type <- power_select(test_stat)[[2]]
  CI <- BinomCI(pwr_calc_array[,test_type,1], pwr_calc_array[,test_type,2])
  CI_low <- CI[,2]
  CI_high <- CI[,3]
  
  eff_sizes = klist
  
  ai = 1; bi = length(eff_sizes)#-1
  if (i == tail(plot_plots,1)) {
    plot(eff_sizes[ai:bi],(power_to_plot[ai:bi]),yaxt='n',xaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
         xlab='Effect Size',ylab='Power',main=plot_title[[i]])
  } else {
    plot(eff_sizes[ai:bi],(power_to_plot[ai:bi]),yaxt='n',xaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
         xlab='',ylab='Power',main=plot_title[[i]])
  }
  axis(side = 2,las = 2) 
  ticks_to_label <- c(1, 2,5,10,20,40)
  axis(1, at = ticks_to_label, labels = ticks_to_label,font=2)
  points(eff_sizes[ai:bi],(CI_low[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  points(eff_sizes[ai:bi],(CI_high[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  abline(h = 0.8,lty='dashed')
  
}

dev.off()

#####################
#####################

png(paste0(wd,'Figure6_67samp.png'), width = 7, height = 14,units='in',res = 400)

# # Define the heights and widths for layouts
# heights <- c(1, 1)  # Equal heights for rows
# widths <- c(1, 1)   # Equal widths for columns

# layout matrix
#layout_matrix <- matrix(c(1, 2, 3, 4, 5), nrow = 5, byrow = TRUE)
#layout(layout_matrix, widths = widths, heights = heights)

par(mfrow = c(5, 1))
par(pty = 's')
par(oma = c(2, 2, 1, 1)) # Outer margins
par(font.axis = 2, font.lab = 2, cex.lab = 1.5, cex.main = 1.6, cex.axis=1.2) # axes
par(mar = c(3.5, 3, 1.5, 1)) # figure margins
par(mgp = c(2.5,0.4,0)) # axis to ticks
par(tcl = -0.2) # tick label size

dir <- paste0('./',optim_folder,'/phylopomp_simulations/PowerAnalysis/kscale_highshape_sampled_twothird/')
load(list.files(dir,full.names=TRUE)[grepl('pwrarrall',basename(list.files(dir,full.names=TRUE)),ignore.case=TRUE)])

# Test statistics to plot
#'mean','ratio','cv','DT','KST','KST','MWWT','tstat'

plot_plots <- c('mean','ratio','tstat','KST','MWWT') #,'DT','cv')
plot_title <- hash(keys=plot_plots,values=c('Mean-Difference', 'Mean-Ratio', 't-statistic', 'KS Test', 'MW Test'))
for (i in plot_plots) {
  test_stat <- i
  
  # par(mar = c(5, 5, 5, 7)) # 5,5,5,7
  
  power_to_plot <- power_select(test_stat)[[1]]
  test_type <- power_select(test_stat)[[2]]
  CI <- BinomCI(pwr_calc_array[,test_type,1], pwr_calc_array[,test_type,2])
  CI_low <- CI[,2]
  CI_high <- CI[,3]
  
  eff_sizes = klist
  
  ai = 1; bi = length(eff_sizes)#-1
  if (i == tail(plot_plots,1)) {
    plot(eff_sizes[ai:bi],(power_to_plot[ai:bi]),yaxt='n',xaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
         xlab='Effect Size',ylab='Power',main=plot_title[[i]])
  } else {
    plot(eff_sizes[ai:bi],(power_to_plot[ai:bi]),yaxt='n',xaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
         xlab='',ylab='Power',main=plot_title[[i]])
  }
  axis(side = 2,las = 2) 
  ticks_to_label <- c(1, 2,5,10,20,40)
  axis(1, at = ticks_to_label, labels = ticks_to_label,font=2)
  points(eff_sizes[ai:bi],(CI_low[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  points(eff_sizes[ai:bi],(CI_high[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  abline(h = 0.8,lty='dashed')
  
}

dev.off()

#####################
#####################

png(paste0(wd,'Figure6_cv.png'), width = 7, height = 14,units='in',res = 400)

# # Define the heights and widths for layouts
# heights <- c(1, 1)  # Equal heights for rows
# widths <- c(1, 1)   # Equal widths for columns

# layout matrix
#layout_matrix <- matrix(c(1, 2, 3, 4, 5), nrow = 5, byrow = TRUE)
#layout(layout_matrix, widths = widths, heights = heights)

par(mfrow = c(5, 1))
par(pty = 's')
par(oma = c(2, 2, 1, 1)) # Outer margins
par(font.axis = 2, font.lab = 2, cex.lab = 1.5, cex.main = 1.6, cex.axis=1.2) # axes
par(mar = c(3.5, 3, 2, 1)) # figure margins
par(mgp = c(2.5,0.4,0)) # axis to ticks
par(tcl = -0.2) # tick label size

dir <- c()
dir[1] <- paste0('./',optim_folder,'/phylopomp_simulations/PowerAnalysis/kscale_highshape/')
dir[2] <- paste0('./',optim_folder,'/phylopomp_simulations/PowerAnalysis/kscale_highshape_sampled_third/')
dir[3] <- paste0('./',optim_folder,'/phylopomp_simulations/PowerAnalysis/kscale_highshape_sampled_twothird/')
dir[4] <- paste0('./',optim_folder,'/phylopomp_simulations/PowerAnalysis/kscale_highshape_scaled/')
dir[5] <- paste0('./',optim_folder,'/phylopomp_simulations/PowerAnalysis/kscale_lowshape/')
plot_title <- c('High shape', 'Het. sampling (33%)', 'Het. sampling (67%)','Scaled','Low shape')

for (i in c(1:5)) {
  test_stat <- 'cv'
  
  load(list.files(dir[i],full.names=TRUE)[grepl('pwrarrall',basename(list.files(dir[i],full.names=TRUE)),ignore.case=TRUE)])
  
  # par(mar = c(5, 5, 5, 7)) # 5,5,5,7
  
  power_to_plot <- power_select(test_stat)[[1]]
  test_type <- power_select(test_stat)[[2]]
  CI <- BinomCI(pwr_calc_array[,test_type,1], pwr_calc_array[,test_type,2])
  CI_low <- CI[,2]
  CI_high <- CI[,3]
  
  eff_sizes = klist
  
  ai = 1; bi = length(eff_sizes)#-1
  if (i == 5) {
    plot(eff_sizes[ai:bi],(power_to_plot[ai:bi]),yaxt='n',xaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
         xlab='Effect Size',ylab='Power',main=plot_title[[i]])
  } else {
    plot(eff_sizes[ai:bi],(power_to_plot[ai:bi]),yaxt='n',xaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
         xlab='',ylab='Power',main=plot_title[[i]])
  }
  axis(side = 2,las = 2) 
  ticks_to_label <- c(1, 2,5,10,20,40)
  axis(1, at = ticks_to_label, labels = ticks_to_label,font=2)
  points(eff_sizes[ai:bi],(CI_low[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  points(eff_sizes[ai:bi],(CI_high[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  abline(h = 0.8,lty='dashed')
  
}

dev.off()

#####################
#####################

png(paste0(wd,'Figure6_dip.png'), width = 7, height = 14,units='in',res = 400)

# # Define the heights and widths for layouts
# heights <- c(1, 1)  # Equal heights for rows
# widths <- c(1, 1)   # Equal widths for columns

# layout matrix
#layout_matrix <- matrix(c(1, 2, 3, 4, 5), nrow = 5, byrow = TRUE)
#layout(layout_matrix, widths = widths, heights = heights)

par(mfrow = c(5, 1))
par(pty = 's')
par(oma = c(2, 2, 1, 1)) # Outer margins
par(font.axis = 2, font.lab = 2, cex.lab = 1.5, cex.main = 1.6, cex.axis=1.2) # axes
par(mar = c(3.5, 3, 2, 1)) # figure margins
par(mgp = c(2.5,0.4,0)) # axis to ticks
par(tcl = -0.2) # tick label size

dir <- c()
dir[1] <- paste0('./',optim_folder,'/phylopomp_simulations/PowerAnalysis/kscale_highshape/')
dir[2] <- paste0('./',optim_folder,'/phylopomp_simulations/PowerAnalysis/kscale_highshape_sampled_third/')
dir[3] <- paste0('./',optim_folder,'/phylopomp_simulations/PowerAnalysis/kscale_highshape_sampled_twothird/')
dir[4] <- paste0('./',optim_folder,'/phylopomp_simulations/PowerAnalysis/kscale_highshape_scaled/')
dir[5] <- paste0('./',optim_folder,'/phylopomp_simulations/PowerAnalysis/kscale_lowshape/')
plot_title <- c('High shape', 'Het. sampling (33%)', 'Het. sampling (67%)','Scaled','Low shape')

for (i in c(1:5)) {
  test_stat <- 'DT'
  
  load(list.files(dir[i],full.names=TRUE)[grepl('pwrarrall',basename(list.files(dir[i],full.names=TRUE)),ignore.case=TRUE)])
  
  # par(mar = c(5, 5, 5, 7)) # 5,5,5,7
  
  power_to_plot <- power_select(test_stat)[[1]]
  test_type <- power_select(test_stat)[[2]]
  CI <- BinomCI(pwr_calc_array[,test_type,1], pwr_calc_array[,test_type,2])
  CI_low <- CI[,2]
  CI_high <- CI[,3]
  
  eff_sizes = klist
  
  ai = 1; bi = length(eff_sizes)#-1
  if (i == 5) {
    plot(eff_sizes[ai:bi],(power_to_plot[ai:bi]),yaxt='n',xaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
         xlab='Effect Size',ylab='Power',main=plot_title[[i]])
  } else {
    plot(eff_sizes[ai:bi],(power_to_plot[ai:bi]),yaxt='n',xaxt='n',ylim=c(0,1),type='o',pch=20,log='x',col='black',
         xlab='',ylab='Power',main=plot_title[[i]])
  }
  axis(side = 2,las = 2) 
  ticks_to_label <- c(1, 2,5,10,20,40)
  axis(1, at = ticks_to_label, labels = ticks_to_label,font=2)
  points(eff_sizes[ai:bi],(CI_low[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  points(eff_sizes[ai:bi],(CI_high[ai:bi]),ylim=c(0,1),type='l',pch=20,col='black',lty=2)
  abline(h = 0.8,lty='dashed')
  
}

dev.off()
