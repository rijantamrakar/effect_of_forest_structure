# Author: Rijan Tamrakar
# R scripts developed for PhD data analysis
# Supervisors: Alexander Knohl, Mark Rayment, Fernando Moyano

# Contains different functions to calculate average values

#------------------------------------------------------------------------------#
rm(list = ls())                                                                 # removing old files if any to clean up the memory
#------------------------------------------------------------------------------#
#### Loading required packages ####
#------------------------------------------------------------------------------#

if('greenbrown' %in% installed.packages()[,"Package"] == FALSE) {
        install.packages("C:/Users/Rijan/Documents/phd_work/r-scripts/greenbrown_2.4.3.tar.gz", repos = NULL, type="source")
}
library(greenbrown) 

# list of required packages
required.packages <- c("dplyr",                                                 # package is good for data manipulation
                       "lubridate",                                             # package for working with date objects
                       "zoo",                                                   # package for working with date objects
                       "gdata",                                                 # package for data manipulation
                       'data.table',                                            # package for data importing
                       'lmodel2',                                               # package for applying lmodel2 and permutation of models
                       'QuantPsyc',                                             # package for regression models and testing and beta values                                       
                       'reshape2',                                              # package for data change from long into wide and vice versa
                       'leaps',                                                 # package for best model selection using subset method
                       'bestglm',                                               # package for best model selection using subset method
                       'psych',                                                 # package for multivariate analysis and scale construction using factor analysis, PCA, cluster analysis and reliability analysis. 
                       'pastecs',                                               # package for Regulation, decomposition and analysis of space-time series
                       'car',                                                   # package for helping regression model with data
                       'Hmisc',                                                 # package for getting lagged value in timeseries data
                       'plotrix',                                               # for drawing table in graph
                       'ggplot2',
                       "fpp",                                                   # time series decomposition
                       "xts",
                       'forecast',
                       'strucchange',
                       'Kendall',
                       'bfast',
                       'phenopix',
                       'ncdf4',
                       'chron',
                       'RColorBrewer',
                       'lattice',
                       'matrixStats' 
)  
# checking if any packages are not installed
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
# Installing uninstalled packages
if(length(new.packages)) {
        install.packages(new.packages)
} 
# loading all required packages
lapply(required.packages, require, character.only = TRUE)

#------------------------------------------------------------------------------#
#### Timeseries flux plot annual ####
#------------------------------------------------------------------------------#

timeseries_avg <- function(FileName,                                            # Filepath to the data
                            Problematic.Year,                                   # List of years when data are not good or unavailable
                            Site,                                               # flux site name
                            timestep)                                           # can be daily, weekly, monthly, annual
        {                      
        
        data           <- read.csv(FileName)
        data           <- dplyr::select(data, 'TIMESTAMP', 'NEE_VUT_REF', 'RECO_DT_VUT_REF', 'GPP_DT_VUT_REF')
        names(data)    <- c('timestamp', 'nee', 'reco', 'gpp')
        data[data == -9999] <- NA
        data$nep       <- -data$nee
        data           <- dplyr::filter(data, !timestamp %in% Problematic.Year)
        
        
        
        ReqdVar    <- c('nep', 'reco', 'gpp')
        average    <- colMeans2(as.matrix(data[,ReqdVar]))
        stdDev     <- colSds(as.matrix(data[,ReqdVar])) 
        CoefVar    <- round(100 * stdDev/average, 1)
        
        results    <- data.frame(ReqdVar, average, stdDev, CoefVar)
        results$no.of.year <- nrow(results)
        results$site <- Site
        results$timestep <- timestep
        
        
        
        # plot the data if asked 
        if(PlotData == 'YES') {
                maxVal <- max(data$nep, data$reco, data$gpp, na.rm = T)
                minVal <- min(data$nep, data$reco, data$gpp, na.rm = T)
                
                with(data, plot(timestamp, 
                                nep, 
                                xlab = '', 
                                ylab = '',
                                typ = 'l', 
                                xlim = c(min(data$timestamp) - 2, max(data$timestamp) + 1), 
                                ylim = c(minVal, maxVal)
                )
                )
                with(data, points(timestamp, nep, pch = 15))
                
                with(data, points(timestamp, gpp, pch = 16))
                with(data, lines(timestamp, gpp, pch = 16))
                
                with(data, points(timestamp, reco, pch = 17))
                with(data, lines(timestamp, reco, pch = 17))
                
                text(min(data$timestamp) + 1, minVal + 50, Site)
                text(min(data$timestamp) -1.5, average[1], 'NEP', cex = 0.7)
                text(min(data$timestamp) -1.5, average[2], 'Reco', cex = 0.7)
                text(min(data$timestamp) -1.5, average[3], 'GPP', cex = 0.7)
                
                mtext(expression( "CO"[2]~ "["~ g ~ C ~ m^{-2}~ yr^{-1}~ ']'), side = 2, line = 2.2, cex = 0.8)
        }
        
        return(results)
}

sitelistwithstructuredata <- substr(list.files('main_analysis/strc.data/final_data', full.names = F), 1, 6)
dir.data <- 'main_analysis/flux.data/ann'
dir.out   <- 'main_analysis/flux.data/annual_flux_all_sites'
data.problem <- read.csv('main_analysis/flux.data/Data_with_problem.csv')

for(j in 1:length(sitelistwithstructuredata)) {
        Site <- sitelistwithstructuredata[j]
        if(length(FileName) == 1) {
                data.problem.for.site <- data.problem %>%
                        filter(site == Site) %>%
                        filter(dont_use_annual == 1)
                
                Problematic.Year <- data.problem.for.site$Not_good_year
                FileName <- list.files(dir.data, pattern = Site, full.names = T)
                
                timeseries_avg(FileName, site)     
        }
}





