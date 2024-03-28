### All the libraries used to arrange the data and statistical analysis
library('RColorBrewer')
library('viridis')
library('ggthemes')
library('pals')
library("yarrr")
library("rcartocolor")

library('extrafont')
library('scales')
library('psych')
library('ggsci')
library('gratia')

library('tools')
library('knitr')
library('dplyr')
library('dtplyr')
library('readxl')
library('writexl')
library('reshape')
library('stringr')
library('data.table')

library('grid')
library('lattice')
library('sciplot')
library('plotrix')
library('ggplot2')
library('ggridges')
library('ggh4x')
library('ggtext')
library('ggpubr')

library('ggforce')
library('rstatix')
library('Matrix')
library('tidyr')
library('tidyverse')

library('tibble')
library('purrr')
library('readr')
library('forcats')

library('MASS')
library('lme4')
library('afex')
library('bestNormalize')
library('pomp')
library('Lahman')
library('zoo')
library('car')
library('minpack.lm')

library('nlme')
library('emmeans')
library('nlraa')
library('mgcv')


### Minimum packages necessary for the modeling
library('data.table')
library('ggplot2')
library('nlme')
library('nlraa')
library('emmeans')
library('mgcv')



### Fit a generalized additive model (GAM) to the data
### Visualize the fitness of the model
### Store the model and fitness outcome

# Read in a data file
# The test data set is stored in the figshare.com. Enter "10.6084/m9.figshare.25467865" at figshare.com to download the test data sets (.csv) named below
gam.dt <- as.data.table(read.csv(file=paste0(getwd(), "/Data_Rsums_mc2r_7.5_2min.csv"), header=T, sep=",", strip.white=T, stringsAsFactors=T))
# gam.dt <- as.data.table(read.csv(file=paste0(getwd(), "/Data_Rsums_mc2r_7.5_7.5min.csv"), header=T, sep=",", strip.white=T, stringsAsFactors=T))

# Make sure the explantory variable is a factor, not an ordered factor
gam.dt[, `:=` (anid = factor(anid, ordered=F),
               genoF = factor(geno, ordered=F))]

# Build a model
theK <- 130       # Define the k value. For the dark-light repeat assays, a value between 100 and 130 usually works
fm0 <- gam(rsums ~ genoF + s(Time, k=theK, by=genoF) + s(anid, bs="re"), data=gam.dt)   # Takes less than a minute

# Print out the summary of the model and the model fitness
capture.output(
      cat("GAM modeling summary: \n\n"),
      cat("theK =", theK, "\n"),
      summary(fm0),
      cat("\n\n\n\n\n"),
      
      cat("gam.check() summary: \n\n"),
      gam.check(fm0),
      cat("\n\n"),
      
      file=paste0(getwd(), "01_GAM_summary.txt"), type="output")

# Print out the graphs of the model fitness
gamck.type1 <- "deviance"  ## "pearson" & "response" are other available choices
gamck.resid <- residuals(fm0, type = gamck.type1)
gamck.linpred <- napredict(fm0$na.action, fm0$linear.predictors)
gamck.observed.y <- napredict(fm0$na.action, fm0$y)

canvasSize <- 10
pdf.options(pointsize=6, useKerning=T, compress=T)
pdf(paste0(getwd(), "/02_GAMcheck_Graphs_fm0.pdf"), width = canvasSize, height = (canvasSize*0.8))
qq.gam(fm0, rep = 0, level = 0.9, type = gamck.type1, rl.col = 2, rep.col = "gray80")
hist(gamck.resid, xlab = "Residuals", main = "Histogram of residuals")
plot(gamck.linpred, gamck.resid, main = "Resids vs. linear pred.", xlab = "linear predictor", ylab = "residuals")
plot(fitted(fm0), gamck.observed.y, xlab = "Fitted Values", ylab = "Response", main = "Response vs. Fitted Values")
dev.off()

# Extract the predicted values from the model, fm0
gam.dt$prd <- predict(fm0, exclude = "s(anid)")       # Takes about a minute

# Store the predicted values
write.table(gam.dt, file=paste0(getwd(), "/03_gamPredict_fm0.csv"), append=F, sep=",", row.names=F, col.names=T)

# Plot the model. A rough plotting before ggpplot
pdf.options(pointsize=6, useKerning=T, compress=T)
pdf(paste0(getwd(), "/04_plot_fm0.pdf"), width = canvasSize, height = (canvasSize*0.8))
plot(fm0, rug=F)
dev.off()

# Similar to the model fitness graphs. This provides better graphics. 
pdf(paste0(getwd(), "/05_drawAppraise_fm0.pdf"), width = canvasSize, height = (canvasSize*0.8))
print(draw(fm0, residuals=T, rug=F))
print(appraise(fm0))
dev.off()



### Visualize the predited values with the raw data

# Sort the genotype so the display is in the order of WT, HT, and HM
gam.dt[, geno := fct_relevel(geno, "WT", "HT", "HM")]

# Produce the graph
gg1 <- ggplot(data=gam.dt, aes(x=Time/60, y=rsums, color=geno)) +
      geom_point(size=0.2, alpha=0.02) +              # Visualize the raw data as points
      geom_line(aes(y=prd), linewidth=1.5) +          # Visualize the predicted values of the model as lines
      xlab("Time (mins)") +                           # x axis label
      ylab("Distance moved (mm/min)") +               # y axis label
      theme_bw() +                                    # Simplest theme for the ggplot
      scale_color_aaas()                              # Color schemes for the journal Science

# Print the ggplot
pdf(file=paste0(getwd(), "/06_gamPlot_.pdf"), width=10, height=8, pointsize=6)
print(gg1)
dev.off()


# Pairwise comparison at a specific time point
# At time = 1100 
pw.means <- emmeans(fm0, ~ genoF, at = list(Time = 1100))
print(pw.means)
pw.est <- pairs(emmeans(fm0, ~ genoF, at = list(Time = 1100)))
print(pw.est)
plot(pairs(emmeans(fm0, ~ genoF, at = list(Time = 1100)))) + 
      geom_vline(xintercept = 0)

# At time = 900
pw.means <- emmeans(fm0, ~ genoF, at = list(Time = 900))
print(pw.means)
pw.est <- pairs(emmeans(fm0, ~ genoF, at = list(Time = 900)))
print(pw.est)
plot(pairs(emmeans(fm0, ~ genoF, at = list(Time = 900)))) + 
      geom_vline(xintercept = 0)

# Loop the pairwise comparison for the entire time of an experiment
# This can take 20 - 50 minutes depending on the experimental condition and computing power
# Do not run the scripts below unless you need the significance tested for the whole experimental window
pw.means <- list(); length(pw.means) <- length(unique(gam.dt[, Time]))
pw.est <- list(); length(pw.est) <- length(unique(gam.dt[, Time]))

for (idx in unique(gam.dt[, Time])) {
      pw.means[[i]] <- emmeans(fm0, ~genoF, at=list(Time=idx))                # emmeans object
      pw.est[[i]] <- pairs(emmeans(fm0, ~genoF, at=list(Time=idx)))           # emmeans object
      
      pw.means.dt <- rbind(pw.means.dt, as.data.table(pw.means[[i]])[, Time := idx])      # Convert the output to a data table, so you can save them later
      pw.est.dt <- rbind(pw.est.dt, as.data.table(pw.est[[i]])[, Time := idx])
      i<-i+1
}

write.table(pw.means.dt, file=paste0(getwd(), "/07_emmeans_means.csv"), append=F, sep=",", col.names=T, row.names=F)
write.table(pw.est.dt, file=paste0(getwd(), "/07_emmeans_estimate.csv"), append=F, sep=",", col.names=T, row.names=F)





