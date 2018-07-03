# Author: Rijan Tamrakar
# R scripts developed for PhD data analysis
# Supervisors: Alexander Knohl, Mark Rayment, Fernando Moyano

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
#### Timeseries flux plot annual function ####
#------------------------------------------------------------------------------#

timeseries_avg <- function(FileName,                                            # Filepath to the data
                           Problematic.Year,                                    # List of years when data are not good or unavailable
                           Site,                                                # flux site name
                           timestep,                                            # can be daily, weekly, monthly, annual
                           PlotData)
{                      
        
        data           <- read.csv(FileName)                                    # reading data
        data           <- dplyr::select(data, 'TIMESTAMP', 
                                        'NEE_VUT_REF', 
                                        'RECO_DT_VUT_REF', 
                                        'GPP_DT_VUT_REF')                       # selecting the carbon dioxide flux columns
        names(data)    <- c('timestamp', 'nee', 'reco', 'gpp')                  # renaming columns
        data[data == -9999] <- NA                                               # converting -9999 into NA
        data$nep       <- -data$nee                                             # getting nep
        
        data2           <- dplyr::filter(data, !timestamp %in% Problematic.Year) # removing problematic year
        
        ReqdVar    <- c('nep', 'reco', 'gpp')                                   # names of variable for result output
        average    <- round(colMeans2(as.matrix(data2[,ReqdVar]), na.rm = T), 1)                      # calculating average using MatixStats package
        stdDev     <- round(colSds(as.matrix(data2[,ReqdVar]), na.rm = T), 1)                         # calculating standard deviation using MatixStats package
        CoefVar    <- round(100 * stdDev/abs(average), 1)                            # calculating coefficient of variation using MatixStats package
        
        results    <- data.frame(ReqdVar, average, stdDev, CoefVar)             # creating a dataframe of results
        results$no.of.year <- nrow(data2)                                        # no of years with data available
        results$site <- Site                                                    # Site for which statistics are calculated
        results$timestep <- timestep                                            # timesteps (daily, weekly, annual, ) 
        
        # plot the data if asked 
        if(PlotData == 'YES') {
                maxY <- max(data$nep, data$reco, data$gpp, na.rm = T)           # maximum value for y - axis
                minY <- min(data$nep, data$reco, data$gpp, na.rm = T)           # minimum value for y - axis
                
                with(data, plot(timestamp, 
                                nep, 
                                xlab = '', 
                                ylab = '',
                                typ = 'l', 
                                xlim = c(min(data$timestamp) - 2, max(data$timestamp) + 1), 
                                ylim = c(minY, maxY)
                )
                )
                with(data, points(timestamp, nep, pch = 15))
                
                with(data, points(timestamp, gpp, pch = 16))
                with(data, lines(timestamp, gpp, pch = 16))
                
                with(data, points(timestamp, reco, pch = 17))
                with(data, lines(timestamp, reco, pch = 17))
                
                text(min(data$timestamp) + 1, minY + 50, Site)
                text(min(data$timestamp) -1.5, average[1], 'NEP', cex = 0.7)
                text(min(data$timestamp) -1.5, average[2], 'Reco', cex = 0.7)
                text(min(data$timestamp) -1.5, average[3], 'GPP', cex = 0.7)
                mtext(expression( "CO"[2]~ "["~ g ~ C ~ m^{-2}~ yr^{-1}~ ']'), 
                      side = 2, 
                      line = 2.2, 
                      cex = 0.8)
        } # end of the plotting
        return(results)
}

#------------------------------------------------------------------------------#
#### Timeseries flux plot calculation ####
#------------------------------------------------------------------------------#
sitelistwithstructuredata <- substr(list.files('main_analysis/strc.data/final_data', full.names = F), 1, 6)
dir.data     <- 'main_analysis/flux.data/ann'
dir.out      <- 'main_analysis/flux.data/annual_flux_all_sites'
data.problem <- read.csv('main_analysis/flux.data/Data_with_problem.csv')
AllResults   <- NULL
dir.out      <- 'main_analysis/flux.data/annual_flux_all_sites'
for(j in 1:length(sitelistwithstructuredata)) {
        Site <- sitelistwithstructuredata[j]
        FileName <- list.files(dir.data, pattern = Site, full.names = T)
        if(length(FileName) == 1) {
                data.problem.for.site <- data.problem %>%
                        filter(site == Site) %>%
                        filter(dont_use_annual == 1)
                
                Problematic.Year <- data.problem.for.site$Not_good_year
                PlotData <- 'NO'
                timestep <- 'ann'
                
                Res <- timeseries_avg(FileName, Problematic.Year, Site, timestep, PlotData)  
                AllResults <- rbind(AllResults, Res)
        }
}

dir.out <- 'main_analysis/flux.data/annual_flux_all_sites'
write.csv(AllResults, file.path(dir.out, 'annualco2fluxesstatistics.csv'), row.names = F)

#------------------------------------------------------------------------------#
#### Corelation with structural index at daily scale  ####
#------------------------------------------------------------------------------#
stru.file <- 'main_analysis/strc.data/structural.index.one.value.csv'
anom.file <- 'main_analysis/anomaly_detection/daily/daily_anomalies.csv'

stru.data <- read.csv(stru.file)
anom.data <- read.csv(anom.file)
Var       <- c('gpp') #, 'gpp', 'reco')
YVar      <- c("HigherNegativeTail995", 
               "areaUnderHigherNegativeTail995",
               "daysUnderHigherNegativeTail995",
               "PerCentdaysUnderHigherNegativeTail995",
               "HigherNegativeTail100",
               "areaUnderHigherNegativeTail100",
               "daysUnderHigherNegativeTail100",
               "PerCentdaysUnderHigherNegativeTail100")

XVar       <- c("CVDBH",
                "meanDBH",
                "noTreeperHa",
                "scaleWeibull",
                "sdDBH",
                "ShannonSize",
                "ShannonSps",
                "shapeWeibull",
                "SimpsonSize",
                "SimpsonSps",
                "TBAm2ha",
                "TotalSpeciesRichnessBA95"
)

resultsCorr <- NULL
for(i in 1:length(Var)) {
        anom.data.var <- dplyr::filter(anom.data, var == var[i])
        data <- merge(anom.data.var, stru.data, by = 'site')
        
        
        
        for(j in 1:length(YVar)) {
                par(oma=c(0,4,0.5,0.5), mar = c(4,0,0,0))
                par(mfrow=c(3,4))
                for(k in 1:length(XVar)) {
                        dataYXvar <- data[, c(YVar[j], XVar[k])]
                        names(dataYXvar) <- c('YVar', 'XVar')
                        corrResults <- cor.test(dataYXvar$YVar, dataYXvar$XVar)
                        pVal <- corrResults$p.value
                        corrValue <- cor(dataYXvar$YVar, dataYXvar$XVar)
                        
                        resultsCorr <- rbind(resultsCorr, c('Var' = Var[i],
                                                            'YVar' = YVar[j],
                                                            'XVar' = XVar[k],
                                                            'corrVal' = corrValue,
                                                            'pVal' = pVal))
                        
                        with(dataYXvar, plot(XVar, YVar, yaxt = 'n', xlab = ''))
                        mtext(XVar[k], outer = F, side = 1,  line = 2.2, cex = 0.8)
                        
                }
                mtext(YVar[j], outer = T, side = 2,  line = 2.2, cex = 0.8)
        }
        
        
}


#------------------------------------------------------------------------------#
#### Corelation with structural index at annual scale  ####
#------------------------------------------------------------------------------#
stru.file <- 'main_analysis/strc.data/structural.index.one.value.csv'
anom.file <- 'main_analysis/flux.data/annual_flux_all_sites/annualco2fluxesstatistics.csv'

stru.data <- read.csv(stru.file)
anom.data <- read.csv(anom.file)
Var       <- c('nep', 'gpp', 'reco')
YVar      <- c("average", 
               "CoefVar")

XVar       <- c("CVDBH",
                "meanDBH",
                "sdDBH",
                "ShannonSize",
                "ShannonSps",
                "shapeWeibull",
                "SimpsonSize",
                "SimpsonSps",
                "TotalSpeciesRichnessBA95"
)

resultsCorr <- NULL
for(i in 1:length(Var)) {
        anom.data.var <- dplyr::filter(anom.data, ReqdVar == Var[i])
        data <- merge(anom.data.var, stru.data, by = 'site')
        
        for(j in 1:length(YVar)) {
                par(oma=c(0,4,0.5,0.5), mar = c(4,0,0,0))
                par(mfrow=c(3,3))
                for(k in 1:length(XVar)) {
                        dataYXvar <- data[, c(YVar[j], XVar[k])]
                        names(dataYXvar) <- c('YVar', 'XVar')
                        minY <- min(dataYXvar$YVar)
                        maxY <- max(dataYXvar$YVar)
                        
                        corrResults <- cor.test(dataYXvar$YVar, dataYXvar$XVar)
                        pVal <- corrResults$p.value
                        corrValue <- cor(dataYXvar$YVar, dataYXvar$XVar)
                        
                        resultsCorr <- rbind(resultsCorr, c('Var' = Var[i],
                                                            'YVar' = YVar[j],
                                                            'XVar' = XVar[k],
                                                            'corrVal' = corrValue,
                                                            'pVal' = pVal))
                        
                        with(dataYXvar, plot(XVar, 
                                             YVar, 
                                             yaxt = 'n', 
                                             xlab = '', 
                                             ylab = '')
                        )
                        if(k %in% c(1,4,7)) axis(2, at = round(seq(minY, maxY, by = (maxY - minY)/4)), labels = T)
                        
                        mtext(XVar[k], outer = F, side = 1,  line = 2.2, cex = 0.8)
                        
                }
                mtext(paste(YVar[j], Var[i]), outer = T, side = 2,  line = 2.2, cex = 0.8)
        }
}

