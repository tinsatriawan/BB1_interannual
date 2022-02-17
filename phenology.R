# data wrangling 
source("data_wrangling.R")

library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(smooth)
library(zoo)
library(modelbased)
library(phenofit)
library(tidyverse)
library(ggpubr)


## growing season division 

INPUT <- check_input(t = as.Date(weekly.GF$floor_date), 
                     y = (weekly.GF$GPP_f_RF))

# INPUT <- check_input(t = as.Date(weekly.nGF$floor_date), 
#                      y = -(weekly.nGF$co2gC))
iters <- 2
lambda <- 100
wFUN <- wTSM
r_max <- 0.2 
r_min <- 0.0
r_minPeakHeight <- 0.05
calendarYear <- FALSE
nptperyear <- 52

brks2 <- season_mov(INPUT, 
                    rFUN = "smooth_wWHIT",
                    wFUN = wFUN,
                    iters = iters, wmin = 0.1,
                    IsOptim_lambda = FALSE,
                    lambda = lambda, nf = 3, frame = frame,
                    maxExtendMonth = 12,
                    r_max = r_max, r_min = r_min,
                    r_minPeakHeight = r_minPeakHeight,
                    calendarYear = calendarYear,
                    # ...,
                    # IsPlot.vc = FALSE,
                    # plotdat = INPUT, print = TRUE,
                    # titlestr = "")
                    IsPlot = FALSE, IsPlot.OnlyBad = FALSE,
                    minpeakdistance = 1/3 * nptperyear, 
                    MaxPeaksPerYear = 3,
                    MaxTroughsPerYear = 4,
                    ypeak_min = 0.08)

plot_season(INPUT, brks2)


## curve fitting
fit <- curvefits(INPUT, brks2, 
                 options = list(
                 methods = c("AG", "Zhang", "Beck", "Elmore"), #,"klos",, 'Gu'
                 wFUN = wBisquare,
                 nextend = 2, 
                 maxExtendMonth = 5, minExtendMonth = 1, 
                 minPercValid = 0.2,
                 print = FALSE,
                 iters = 2,
                 # use.rough = TRUE,
                 verbose = FALSE))


d_fit <- get_fitting(fit)
d_gof <- get_GOF(fit)

p <- plot_curvefits(d_fit, brks2, title = NULL, cex = 1.5, ylab = "co2", angle = 0)
grid::grid.newpage()
grid::grid.draw(p)

l_param <- get_param(fit)
# pheno_plot<- plot_phenofit(d_fit, brks2, NULL)

pheno_plot<- plot_phenofit(obj = list(INPUT = INPUT, fit = fit, seasons = brks2),
                           type = "pheno", 
                           methods = c("Zhang", "AG", "Beck"),
                           title.ylab = "c_flux", title.xlab = "Time",
                           theme = coord_cartesian(xlim = c(as.Date(min(weekly.nGF$floor_date)), 
                                                            as.Date(max(weekly.nGF$floor_date)))))




##  Extract phenology
TRS <- c(0.1, 0.2, 0.6)
l_pheno <- get_pheno(fit, wFUN = wBisquare, TRS = TRS, IsPlot = FALSE) # %>% map(~melt_list(., "meth"))
print(l_pheno$doy$Zhang)
pheno <- l_pheno$doy %>% melt_list("meth")


zhang <- l_pheno$doy$Zhang

pheno_index <- data.frame(CUP = zhang$DER.SOS )

l_pheno <- get_pheno(fit[1:7], wFUN = wBisquare, TRS = TRS, IsPlot = T) # %>% map(~melt_list(., "meth"))

