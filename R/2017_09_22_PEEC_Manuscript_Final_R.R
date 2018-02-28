rm(list=ls())

######################################################################################
######################################################################################

# Load required packages
library(plyr) # 'split-apply-combine'
library(lubridate) # Deals with dates
library(ggplot2) # plotting
library(reshape2)
library(gridExtra)
library(cowplot)
library(scales)

# Use a function I wrote
source("R/lm_stats.R") # start from working directory
#source("../R/lm_stats.R") # for Drew
#source("R/process_LM_PEEC_compiled.R")
####Could not figure out how to process the data using the process_LM_PEEC_compiled.R file######## 

# Set graphical theme
theme_set(theme_bw() + theme(text=element_text(size=9)))

#load data
thirteen <- read.csv ("data/MSXUT_v2_2013_updated.csv")
fifteen <- read.csv ("data/MSXUT_2015_compiled.csv")
tnmine <- read.csv ("data/2015_TN_compiled_fixed.csv")
sixteen <- read.csv ("data/2016_12_12_PEEC_2016_data.csv")

all_data_list <- list(thirteen=thirteen, fifteen=fifteen, tnmine=tnmine, sixteen=sixteen)

raw_df <- ldply(all_data_list, identity)
#raw_df <- rename(raw_df, c(".id" = "location"))

#add substrate names and set unlabled reps to A
raw_df$substrate <- factor(raw_df$substrate, labels=c("Arg-AP", "Gly-AP", "Leu-AP", "Pyr-AP", "Trypsin"))
raw_df$rep[raw_df$rep==""] <- "A"

times <- paste("2015_06_01", raw_df$time)
raw_df$Rtime <- ymd_hm(times) #makes time into time object
raw_df <- ddply(raw_df, c("site"), mutate, elapsed=as.numeric(Rtime-min(Rtime))/3600)
attr(raw_df$elapsed, "units") <- "hours"

#generate plot of raw data
p_raw <- ggplot(raw_df, aes(elapsed, y=RFU, shape=treatment, colour=site, fill=treatment)) +
  geom_point() +
  geom_smooth(method="lm", se=TRUE) +
  facet_grid(when ~ fluorimeter ~ substrate, scales = "free") +
  theme(text=element_text(size=9), axis.text.x=element_text(angle=-70, hjust=0))
# print(p_raw)

# Calculate the slopes
uncalibrated_slopes <- ddply(raw_df, c(".id", "fluorimeter", "substrate", "when", "treatment", "site", "volume", "rep"), 
                             function(x) lm_stats(xvar="elapsed", yvar="RFU", d=x))

# Export the slopes to a .csv file
# write.csv(uncalibrated_slopes, "data/PEEC_TN_final_data_sheets/2016_12_12_PEEC_and_TN_uncalibrated_slopes.csv")

# Plot the slopes
p_slopes <- ggplot(uncalibrated_slopes, aes(x=substrate, y=slope, colour=site, shape=treatment)) +
  geom_pointrange(aes(ymin=slope-slope.se, ymax=slope+slope.se), position=position_jitter(width=0.1)) +
  facet_grid(when ~ fluorimeter, scales = "free_y") +
  theme_bw()+
  theme(text=element_text(size=9), axis.text.x=element_text(angle=-70, hjust=0)) +
  ggtitle("Uncalibrated Slopes")
# print(p_slopes)

#load calibration data
thirteen_calib <- read.csv("data/2013_all_calib_data.csv")
fifteen_calib <- read.csv("data/2016_12_13_PEEC_2015_calib_nM.csv")
tnmine_calib <- read.csv("data/2015_TN_compiled_calib_fixed.csv")
sixteen_calib <-read.csv("data/2016_12_12_2016_peec_calib.csv")

#bind calibration data
all_calib_list <- list(thirteen_calib=thirteen_calib, fifteen_calib=fifteen_calib, 
                       tnmine_calib=tnmine_calib, sixteen_calib=sixteen_calib)
all_sites_df_calib <- do.call(rbind.fill, all_calib_list)

#look at calibration data
p_calib <- ggplot(all_sites_df_calib, aes(x=conc.AMC.nM, y=RFU, colour= site)) +
  geom_point() + #scatterplot
  geom_smooth(method="lm", se=TRUE) +
  facet_wrap(when~fluorimeter , scales="free") +
  theme(text=element_text(size=9), axis.text.x=element_text(angle=90, hjust=0)) +
  theme_bw()+
  ggtitle("Total Calibration for AMC")
#print(p_calib)
#ggsave("plots/PEEC and TN Final/for_pub_2016_12_14_peec_and_tn_calibration_curve.tiff", height=2.5, width=7.5, units="in", dpi=300, compression='lzw')

# Calculate a slope of fluorescence vs concentration for each unique fluorimeter and site and year
calib_slopes <- ddply(all_sites_df_calib, c("site", "fluorimeter", "when", "fluorophore"), function(x) lm_stats(xvar="conc.AMC.nM", yvar="RFU", d=x))

# Merge the calibration slopes into the data of uncalibrated slopes of RFU vs time
calib_slopes_short <- plyr::rename(calib_slopes[ , c("site", "fluorimeter", "when", "slope")], c("slope" = "calib.slope"))
#write.csv(calib_slopes_short, "data/2016_12_19_calib_slopes_short_2.csv")

# Merge calibration slopes into main data frame of uncalibrate slopes 
slopes_merged <- merge(uncalibrated_slopes, calib_slopes_short, by=c("site", "fluorimeter", "when"), all=TRUE)

# Calculate the calibrated slopes
slopes_merged$v0 <- slopes_merged$slope / slopes_merged$calib.slope / 1000 # CHANGING FROM NM TO uM

# Set units correctly
attr(slopes_merged$v0, "units") <- "umol per liter per hour"

# Pull out slopes for the saturation curves
slopes_sat_curve <- subset(slopes_merged, .id="thirteen")

# Eliminate saturatino curve slopes from the main slopes data frame
slopes <- subset(slopes_merged, .id!="thirteen" | (.id == "thirteen" & volume == 40))

# Calculate average and standard deviation of slopes
slopes_summ <- ddply(slopes_merged, c("site", "fluorimeter", "when", "substrate", "treatment"), summarise, 
                     v0.mean = mean(v0, na.rm=TRUE),
                     v0.sd = sd(v0, na.rm=TRUE))

#Cast
live_v0 <- slopes_summ[slopes_summ$treatment == "live", ]
live_v0$substrate <- revalue(live_v0$substrate, c("Arg-AP"="ArgAP", "Leu-AP"="LeuAP", "Trypsin"="Trypsin", "Gly-AP"="GlyAP", "Pyr-AP"="PyrAP"))

# remove NA values
live_v0 <- live_v0[!is.na(live_v0$v0.mean), ]


# Calculate median, geometric mean, and IQR; make fig S1
library(EnvStats)
gmean.v0s <- live_v0$v0.mean[live_v0$v0.mean > 0]
g.mean <- geoMean(gmean.v0s); print(g.mean)
v0.IQR <- quantile(gmean.v0s, c(0.25, 0.5, 0.75)); print(v0.IQR)
p_dist <- ggplot(live_v0, aes(x=v0.mean)) + 
  geom_density() + 
  geom_rug(sides="t") +
  scale_x_log10() + 
  annotation_logticks(sides="b") +
  xlab(expression(paste(log[10], " ", v[0], ", ", mu, "M ", hr^{-1})))
print(p_dist)
# ggsave("plots/S1_v0_distribution.png", height=3, width=4, units="in", dpi=300)


v0_c <- dcast(live_v0, site+fluorimeter+when~substrate, value.var= 'v0.mean')
v0_c_sd <- dcast(live_v0, site+fluorimeter+when~substrate, value.var="v0.sd")
colnames(v0_c_sd)[4:ncol(v0_c_sd)] <- paste0(colnames(v0_c_sd)[4:ncol(v0_c_sd)], ".sd") # add .sd to the st dev data frame names, to prevent confusion in the merge

#merge data into one "wide" frame
v0_wide <- merge(v0_c, v0_c_sd, by=c("site", "fluorimeter", "when"))


# make a boxplot of live v0 data
p_boxes <- ggplot(live_v0, aes(x=substrate, y=log10(v0.mean))) + # trying log v0 for kicks
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(alpha=0.3) +
  xlab("Substrate") +
  ylab(expression(paste("log" [10], v[0], ", ", mu, "M ", hr^{-1})))  # see demo(plotmath)
print(p_boxes)
# ggsave("plots/v0_boxplot.tiff", height=2.5, width=7.08/2, units="in", dpi=300, compression="lzw")

# ANOVA with post-hoc testing of rates by substrate
# BUT WE NEED PAIRED ANOVA
# # Try it with un-transformed rates
# bad_aov <- aov(v0.mean ~ substrate, data=live_v0)
# # plot(bad_aov) # pretty bad!
# 
# live_v0$log10.v0 <- log10(live_v0$v0.mean) 
# live_v0_aov <- subset(live_v0, !is.na(log10.v0))
# v0_mod <- aov(log10.v0 ~ substrate, data=live_v0) # p<0.001, n=192
# # plot(v0_mod) # plots look fine I guess
# summary(v0_mod)
# v0_Tukey <- TukeyHSD(v0_mod)
# print(v0_Tukey) ## NOTES ON TUKEY:
# 
# # Pyr < ARG
# # Pyr < GLY
# # Pyr < Leu
# # Pyr < TRP

# #######hydrolysis rates
fs <- 8

ee_scatterplot <- function(y.var, y.var.name) {
  # Create strings for errorbars
  y.err.min <- paste0(y.var, "-", y.var, ".sd")
  y.err.max <- paste0(y.var, "+", y.var, ".sd")
  
  # Create string for axis label
  y.ax.label <- expression(paste("log" [10], " ", eval(y.var.name), ", ", mu, "M L ", hr^{-1}), sep="")
  
  # Make the plot
  p <- ggplot(v0_wide, aes_string(x="LeuAP", y=y.var)) +
    geom_point() + 
    geom_errorbar(aes_string(ymin=y.err.min, ymax=y.err.max)) +
    geom_errorbarh(aes(xmin=LeuAP-LeuAP.sd, xmax=LeuAP+LeuAP.sd)) +
    geom_smooth(colour="black", method="lm", se=TRUE) +
    scale_x_log10() + 
    scale_y_log10() +
    xlab(expression(paste("log" [10], " LeuAP", ", ", mu, "M L ", hr^{-1}))) +
    ylab(y.ax.label) +
    theme(axis.text.x=element_text(angle=-70, hjust=0), text=element_text(size=fs), axis.text=element_text(size=fs))
  p
}

p_arg_leu <- ee_scatterplot("ArgAP", "ArgAP")
print(p_arg_leu)

#for arg vs leu
hyd_rate_comp_arg_leu <- ggplot(v0_wide, aes(x=LeuAP, y=ArgAP)) +
  #geom_pointrange(aes(ymin=ArgAP-ArgAP.sd, ymax=ArgAP+ArgAP.sd)) +
  geom_point() +
  geom_errorbar(aes(ymin=ArgAP-ArgAP.sd, ymax=ArgAP+ArgAP.sd)) +
  geom_errorbarh(aes(xmin=LeuAP-LeuAP.sd, xmax=LeuAP+LeuAP.sd)) +
  geom_smooth(colour="black", method="lm", se=TRUE) +
  scale_x_log10() + 
  scale_y_log10() +
  ylab(expression(paste("log" [10], " ArgAP", ", ", mu, "M L ", hr^{-1}))) + # see demo(plotmath)
  xlab(expression(paste("log" [10], " LeuAP", ", ", mu, "M L ", hr^{-1}))) + # see demo(plotmath)
  theme(axis.text.x=element_text(angle=-70, hjust=0), text=element_text(size=fs), axis.text=element_text(size=fs))
print(hyd_rate_comp_arg_leu)
#ggsave("plots/PEEC and TN Final/for_pub_2016_12_21_peec_and_tn_arg_vs_leu.tiff", hyd_rate_comp_arg_leu, height=2.5, width=7.08, units="in", dpi=300, compression= 'lzw')

#for Gly vs leu
hyd_rate_comp_Gly_leu <- ggplot(v0_wide, aes(x=LeuAP, y=GlyAP)) +
  #geom_pointrange(aes(ymin=GlyAP-GlyAP.sd, ymax=GlyAP+GlyAP.sd)) +
  geom_point() +
  geom_errorbar(aes(ymin=GlyAP-GlyAP.sd, ymax=GlyAP+GlyAP.sd)) +
  geom_errorbarh(aes(xmin=LeuAP-LeuAP.sd, xmax=LeuAP+LeuAP.sd)) +
  geom_smooth(colour="black", method="lm", se=TRUE) +
  scale_x_log10() +
  scale_y_log10() +
  ylab(expression(paste("log" [10], " GlyAP", ", ", mu, "M L ", hr^{-1}))) + # see demo(plotmath)
  xlab(expression(paste("log" [10], " LeuAP", ", ", mu, "M L ", hr^{-1}))) + # see demo(plotmath)
  theme(axis.text.x=element_text(angle=-70, hjust=0), text=element_text(size=fs), axis.text=element_text(size=fs))
#facet_wrap(~when, scales="free_y") +
print(hyd_rate_comp_Gly_leu)
#ggsave("plots/PEEC and TN Final/for_pub_2016_12_21_peec_and_tn_gly_vs_leu.tiff", hyd_rate_comp_Gly_leu, height=2.5, width=7.08, units="in", dpi=300, compression= 'lzw')

#for GlyGlyArg vs. Leu
hyd_rate_comp_GlyGlyArg_leu <- ggplot(v0_wide, aes(x=LeuAP, y=Trypsin)) +
  geom_point() +
  #geom_pointrange(aes(ymin=GlyGlyArgAP-GlyGlyArgAP.sd, ymax=GlyGlyArgAP+GlyGlyArgAP.sd)) +
  geom_errorbar(aes(ymin=Trypsin-Trypsin.sd, ymax=Trypsin+Trypsin.sd)) +
  geom_errorbarh(aes(xmin=LeuAP-LeuAP.sd, xmax=LeuAP+LeuAP.sd)) +
  geom_smooth(colour="black", method="lm", se=TRUE) +
  scale_x_log10() +
  scale_y_log10() +
  ylab(expression(paste("log" [10], " Trypsin", ", ", mu, "M L ", hr^{-1}))) + # see demo(plotmath)
  xlab(expression(paste("log" [10], " LeuAP", ", ", mu, "M L ", hr^{-1}))) + # see demo(plotmath)
  theme(axis.text.x=element_text(angle=-70, hjust=0), text=element_text(size=fs), axis.text=element_text(size=fs))
print(hyd_rate_comp_GlyGlyArg_leu)
#ggsave("plots/PEEC and TN Final/2016_03_11_glyglyarg_vs_leu.tiff", hyd_rate_comp_GlyGlyArg_leu, height=2.5, width=7.08, units="in", dpi=300, compression= 'lzw')

#for Pyr vs. Leu
hyd_rate_comp_Pyr_leu <- ggplot(v0_wide, aes(x=LeuAP, y=PyrAP)) +
  geom_point() +
  #geom_pointrange(aes(ymin=GlyGlyArgAP-GlyGlyArgAP.sd, ymax=GlyGlyArgAP+GlyGlyArgAP.sd)) +
  geom_errorbar(aes(ymin=PyrAP-PyrAP.sd, ymax=PyrAP+PyrAP.sd)) +
  geom_errorbarh(aes(xmin=LeuAP-LeuAP.sd, xmax=LeuAP+LeuAP.sd)) +
  geom_smooth(colour="black", method="lm", se=TRUE) +
  scale_x_log10() +
  scale_y_log10() +
  ylab(expression(paste("log" [10], " PyrAP", ", ", mu, "M L ", hr^{-1}))) + # see demo(plotmath)
  xlab(expression(paste("log" [10], " LeuAP", ", ", mu, "M L ", hr^{-1}))) + # see demo(plotmath)
  theme(axis.text.x=element_text(angle=-70, hjust=0), text=element_text(size=fs), axis.text=element_text(size=fs))
print(hyd_rate_comp_Pyr_leu)

#create a figure that has the above 3 plots in one

#going to try using the package  gridExtra package cowplot
# library("gridExtra")
# library("cowplot")

p_all_hydrolysis <- plot_grid(hyd_rate_comp_arg_leu, hyd_rate_comp_Gly_leu, hyd_rate_comp_GlyGlyArg_leu, 
                              hyd_rate_comp_Pyr_leu, NULL, labels=c("A", "B", "C", "D"), nrow=2, ncol=2)
#save_plot("plots/PEEC and TN Final/2017_03_27_all_hydrolysis_comp_v3.tiff", p_all_hydrolysis, ncol = 2, nrow = 2)
print(p_all_hydrolysis)

#same plot without pyrleu comp
p_all_hydrolysis_2 <- plot_grid(hyd_rate_comp_arg_leu, hyd_rate_comp_Gly_leu, hyd_rate_comp_GlyGlyArg_leu, 
                                NULL, labels=c("A", "B", "C", ""), nrow=2, ncol=2)
#save_plot("plots/PEEC and TN Final/2017_08_22_all_hydrolysis_comp_v4.tiff", p_all_hydrolysis_2, ncol = 2, nrow = 2)
print(p_all_hydrolysis)

# save_plot("plots/PEEC and TN Final/2017_03_27_all_hydrolysis_comp_v2.tiff", p_all_hydrolysis, 
#           base_height=2.5, base_width=7.08,  units="in", dpi=400, compression= 'lzw')
save_plot("plots/scatter_plots.tiff", p_all_hydrolysis, 
          base_height=5, base_width=6,  units="in", dpi=400, compression= 'lzw')

# Do propoer linear model analysis
# #Arg vs. Leu
# mArgLeu <- lm(ArgAP ~ LeuAP, data=v0_wide)
# summary(mArgLeu)
# # plot(mArgLeu)

m_log_ArgLeu <- lm(log10(ArgAP) ~ log10(LeuAP), data=v0_wide)  ## Added log to normalize the plots
summary(m_log_ArgLeu)
#plot(m_log_ArgLeu)

# #Gly vs. Leu
# mGlyLeu <- lm(GlyAP ~ LeuAP, data=v0_wide)
# summary(mGlyLeu)
# plot(mGlyLeu) # plot function: returns homoscedasticity tests, Hmm, this lm looks pretty abysmal

m_log_GlyLeu <- lm(log10(GlyAP) ~ log10(LeuAP), data=v0_wide)
summary(m_log_GlyLeu)
#plot(m_log_GlyLeu)

#GGR vs. Leu
# mGGRLeu <- lm(Trypsin ~ LeuAP, data=v0_wide)
# summary(mGGRLeu)
# plot(mGGRLeu)

m_log_GGRLeu <- lm(log10(Trypsin) ~ log10(LeuAP), data=v0_wide)  ## Added log to normalize the plots
summary(m_log_GGRLeu)
# plot(m_log_GGRLeu)

#pyr vs. leu
mPyrLeu <- lm(log10(PyrAP) ~ log10(LeuAP), data=v0_wide)
summary(mPyrLeu)

#interquartile range calculations for trypsin to leu-AMC

#mGGRLeu <- lm(Trypsin ~ LeuAP, data=v0_wide)
# min(v0_wide$LeuAP)
# max(v0_wide$LeuAP) 
# quantile(v0_wide$LeuAP, c(0.25, 0.5, 0.75))
# 
# min(v0_wide$Trypsin)
# max(v0_wide$Trypsin) 
# quantile(v0_wide$Trypsin, c(0.25, 0.5, 0.75))

v0_wide$tryp.ratio <- v0_wide$Trypsin / v0_wide$LeuAP
ggplot(v0_wide, aes(x=tryp.ratio)) + 
  #geom_histogram(bins=10) +
  geom_density() +
  geom_rug(sides="t") +
  scale_x_log10() + 
  annotation_logticks() + 
  xlab(expression(paste(v[0], ",trypsin ", "/ ", v[0], ",LeuAP")))
ggsave("plots/S2_tryp_ratio_distribution.tiff", height=3, width=3, units="in", dpi=300, compression="lzw")
#ratio_df <- subset(v0_wide, select = c("site", "fluorimeter", "when", "LeuAP", "Trypsin", "PyrAP"))

v0_wide$pyr.ratio <- v0_wide$PyrAP / v0_wide$LeuAP
v0_wide$arg.ratio <- v0_wide$ArgAP / v0_wide$LeuAP
v0_wide$gly.ratio <- v0_wide$GlyAP / v0_wide$LeuAP


ratio_df <- melt(v0_wide[ , c("site", "fluorimeter", "when", "tryp.ratio", "pyr.ratio", "arg.ratio", "gly.ratio")], id.vars=c("site", "fluorimeter", "when"), value.name="ratio")
ratio_df <- ratio_df[ratio_df$ratio > 0 & !is.na(ratio_df$ratio), ]


library(cvequality)
ratio_df_no_pyr <- subset(ratio_df, variable != "pyr.ratio")
ratio_df_no_pyr$log.ratio <- log10(ratio_df_no_pyr$ratio)

#with(ratio_df_no_pyr, asymptotic_test(log.ratio, variable))
#asymptotic_test(ratio_df_no_pyr$log.ratio, ratio_df_no_pyr$variable)


p_ratio_dist <- ggplot(ratio_df, aes(x=ratio, fill=variable)) + 
  geom_density(alpha=0.3) + 
  scale_x_log10()
print(p_ratio_dist)
# ggsave("PEEC_Final/Final_plots/PEEC and TN Final/supplemental/S2_tryp_ratio_distribution.png", height=3, width=3, units="in", dpi=300)
# ggplot(ratio_df, aes(x=variable, y=log10(ratio))) + 
#  geom_boxplot() + 
#  geom_point(position=position_jitter(width=0.3)) 
  

# ratio_df$tryp.ratio <- ratio_df$Trypsin / ratio_df$LeuAP
# ggplot(ratio_df, aes(x=tryp.ratio)) + 
#   #geom_histogram(bins=10) +
#   geom_density() +
#   geom_rug(sides="t") +
#   scale_x_log10() + 
#   annotation_logticks() + 
#   xlab(expression(paste(v[0], ",trypsin ", "/ ", v[0], ",LeuAP")))
# ggsave("PEEC_Final/Final_plots/PEEC and TN Final/supplemental/S2_ratio_distribution.png", height=3, width=3, units="in", dpi=300)

ratios_asymptotic_test <- with(ratio_df, mslr_test(nr=1e4, ratio, variable))

geoMean(ratio_df$tryp.ratio)
max(ratio_df$tryp.ratio)
min(ratio_df$tryp.ratio)
quantile(ratio_df$ratio, c(0.25, 0.5, 0.75))



m2 <- lm(LeuAP ~ Trypsin, data=df.LeuGGR)
summary(m2)
den <- density(df.LeuGGR$LeuAP)
den2 <- density(df.LeuGGR$Trypsin)
plot(den)
plot(den2)

quantile(df.LeuGGR)

p_density <- ggplot(mGGRLeu, aes(x=v0_mean)) +
  geom_density()
print(p_density)
#####################################
#finish this 
#make a density plot of trypsin to leu
v0_trypsin_leu <- data.frame(v0_wide$site, v0_wide$fluorimeter, v0_wide$when, v0_wide$Trypsin, v0_wide$LeuAP, v0_wide$Trypsin.sd, v0_wide$LeuAP.sd)
names(v0_trypsin_leu) <- c("site", "fluorimeter", "when", "Trypsin", "LeuAP", "Trypsin.sd", "LeuAP.sd")
p_density <- ggplot(v0_trypsin_leu, aes(x=v0_mean)) +
  geom_density()
print(p_density)

#Mapping 
#install.packages("ggmap")
#install.packages("RColorBrewer")
#install.packages("ggrepel")
library(ggplot2)
library(ggrepel)
library(ggmap)
library(RColorBrewer)

#Generate Large map
#big_map <- get_googlemap(center=c(-83.921761, 35.959802), zoom=10, maptype =  "terrain")
big_map <- get_googlemap(center= c(-80.365531, 38.365531), zoom=6, maptype =  "terrain")


p_big <- ggmap(big_map) + 
  #geom_polygon(data=map_data("state"), aes(x=long, y=lat, group=group), col="black", fill=NA) + 
  #coord_map(xlim = c(-80, -68), ylim=c(32.5, 42.5)) +
  theme(text=element_text(size=8))
print(p_big)

#Add a detailed pennsylvania map
PAbase <- get_googlemap(center=c(-74.9149, 41.17067), zoom=13, maptype = "terrain")

p_PAbase <- ggmap(PAbase) +
  theme(text=element_text(size=9))
print(p_PAbase)

# Load data
coords <- read.csv("data/2015_06_03_peec_sample_sites.csv")
coords_2 <- read.csv("data/2015_06_01_sampling_sites.csv")

p_whole_both <- ggmap(p_big) +
  geom_point(data=coords, aes(x=long, y=lat), colour="black") +
  geom_point(data=coords_2, aes(x=long, y=lat), colour="red") 
print(p_whole_both)
#ggsave("plots/PEEC and TN Final/p_whole_both.tiff", p_whole_both, height=2.5, width=7.08, units="in", dpi=300, compression='lzw')

#Generate a detailed TN map
p_TN_map <- ggmap(big_map) + 
  geom_point(data=coords, aes(x=long, y=lat), colour="red", zoom=8) + 
  ggrepel::geom_label_repel(data=coords, aes(x=long, y=lat, label=site.init), size=3, colour="black", fill="white") 
print(p_TN_map)
#ggsave("plots/PEEC and TN Final/TN_map_two.tiff", p_TN_map, height=2.5, width=7.08, units="in", dpi=300, compression= 'lzw')

p_PA_map <- ggmap(PAbase) +
  geom_point(data=coords_2, aes(x=long, y=lat), colour="red", zoom=8) +
  ggrepel::geom_label_repel(data=coords_2, aes(x=long, y=lat, label=site.init), size=3, colour="black", fill="white")
print(p_PA_map)
#ggsave("plots/PEEC and TN Final/PA_map_final.tiff", p_PA_map, height=2.5, width=7.08, units="in", dpi=300, compression= 'lzw')

#create a figure that has the above 2 maps in one plot
#load gridExtra and cowplot packages
library("gridExtra")
library("cowplot")

#Make a plot with gridextra that has the TN and PA map with lables, this will be the base for the less detailed map
bottom_row_map <- plot_grid(p_TN_map, p_PA_map, labels = c('TN', 'PA'), ncol = 1, align = 'h', rel_widths = c(1, 1))
print(bottom_row_map)

#Generate a plot that has all three maps in one image
p_maps_cow_top <- plot_grid(p_whole_both, bottom_row_map,
                                      NULL, labels=c('A', ''), ncol=1, rel_heights = c(1,1))
print(p_maps_cow_top)
#save_plot("plots/PEEC and TN Final/2017_08_28_finished_map.tiff", p_maps_cow_top, ncol = 2, nrow = 2)


p_all_maps <- plot_grid(p_TN_map, p_PA_map, NULL,
                        nrow=1, ncol=2)

print(p_all_maps)
#save_plot("plots/PEEC and TN Final/2017_03_27_all_maps.tiff", p_all_maps, 
          base_height=2.5, base_width=7.08,  units="in", dpi=400, compression= 'lzw')

