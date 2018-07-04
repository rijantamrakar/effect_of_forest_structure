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


################################################################################
################################################################################
#------------------------------------------------------------------------------#
#### Preparing list of sites for analysis  ####
#------------------------------------------------------------------------------#
final_data <- read.csv('Communication_with_PIs/finallist.csv')
list_of_sites <- as.character(final_data$Site.CODE)

available_data <- data.frame(sn = character(length(list_of_sites)))
available_data$site <- list_of_sites
available_data <- available_data['site']

#------------------------------------------------------------------------------#
#### Prepare tier1 available data  ####
#------------------------------------------------------------------------------#
# what data is available in fluxnet 2015 tier 1 data
dir.data      <- 'Communication_with_PIs/fluxnet2015_tier1'
list_fn2015   <- list.files(dir.data, pattern = 'zip')
list_of_sites <- substr(list_fn2015, 5,10)

tier1_sy            <- as.numeric(substr(list_fn2015, 32, 35))
tier1_ey            <- as.numeric(substr(list_fn2015, 37, 40))

available_data_tier <- data.frame(sn = character(length(list_of_sites)))
available_data_tier$site <- list_of_sites
available_data_tier <- available_data_tier['site']
available_data_tier$tier1_sy   <- tier1_sy
available_data_tier$tier1_ey   <- tier1_ey

# merge with available data
available_data <- merge(available_data, 
                        available_data_tier, 
                        by = 'site', 
                        all.x = T)

available_data$tier1_ny <- available_data$tier1_ey - available_data$tier1_sy

#------------------------------------------------------------------------------#
#### Prepare tier2 available data  ####
#------------------------------------------------------------------------------#
# what data is available in fluxnet 2015 tier 1 data
dir.data      <- 'Communication_with_PIs/fluxnet2015_tier2'
list_fn2015   <- list.files(dir.data, pattern = 'zip')
list_of_sites <- substr(list_fn2015, 5,10)

tier1_sy            <- as.numeric(substr(list_fn2015, 32, 35))
tier1_ey            <- as.numeric(substr(list_fn2015, 37, 40))

available_data_tier <- data.frame(sn = character(length(list_of_sites)))
available_data_tier$site <- list_of_sites
available_data_tier <- available_data_tier['site']
available_data_tier$tier2_sy   <- tier1_sy
available_data_tier$tier2_ey   <- tier1_ey

# merge with available data
available_data <- merge(available_data, 
                        available_data_tier, 
                        by = 'site', 
                        all.x = T)

available_data$tier2_ny <- available_data$tier2_ey - available_data$tier2_sy

#------------------------------------------------------------------------------#
#### Merge with final data  ####
#------------------------------------------------------------------------------#
final_data <- merge(final_data, available_data, by.x = 'Site.CODE', by.y = 'site')

#------------------------------------------------------------------------------#
#### Add data from fernandez  ####
#------------------------------------------------------------------------------#
fernandex_data <- read.csv('Fernandez_paper_2017/trend_data_Fernandez_2017.csv')
fernandex_data <- dplyr::select(fernandex_data, -Forest, -FT, -Y)
header_fernand <- names(fernandex_data)
header_fernand <- paste0('F_', header_fernand)
names(fernandex_data) <- header_fernand
final_data <- merge(final_data, fernandex_data, by.x = 'Site.CODE', by.y = 'F_Code', all.x = T)

#------------------------------------------------------------------------------#
#### Add data from talie  ####
#------------------------------------------------------------------------------#
musavi_data <- read.csv('Talie_paper_2017/Cv_gpp_Talie_add.csv')
musavi_data <- dplyr::select(musavi_data, -X, -PFT, -ClimateGroup, -LAT, -LONG)

header_musavi <- names(musavi_data)
header_musavi <- paste0('M_', header_musavi)
names(musavi_data) <- header_musavi
final_data <- merge(final_data, musavi_data, by.x = 'Site.CODE', by.y = 'M_site.code', all.x = T)

write.csv(final_data, file.path(paste0('Communication_with_PIs/finallist.', format(Sys.time(), "%Y%m%d"), '.csv')), row.names = F)
################################################################################
################################################################################

#------------------------------------------------------------------------------#
# preparing required daily Fluxnet2015 data ####
#------------------------------------------------------------------------------#
reqd.var <- c("TIMESTAMP", 
              "TA_F", 
              "SW_IN_F", 
              "VPD_F", 
              'NEE_VUT_REF',
              'RECO_DT_VUT_REF', 
              'GPP_DT_VUT_REF')

dir.data   <- 'flux_data/Unzipped_Fluxnet2015_database_tier_2'
filetype   <- 'FULLSET_DD'
timestep   <- 'day'
dir.out    <- "main_analysis/flux.data"
FinalSites <- read.csv('main_analysis/finallist_20180626.csv') 
ListSites  <- as.character(FinalSites$site)

SitesFluxnet2015 <- list.files(dir.data)
nodata <- NULL
for(i in 1:length(ListSites)) {
        if(length(grep(ListSites[i], SitesFluxnet2015)) == 1) {
                
                cat(paste('Reading data for site', ListSites[i]), '\n')
                site_Filename<- list.files(file.path(dir.data, ListSites[i]), pattern = filetype, full.names = T)
                
                dataheader   <- names(read.table(site_Filename, nrow = 1, sep = ',', header = T))
                col.var      <- which(dataheader %in% reqd.var)
                
                site_Data    <- data.table::fread(site_Filename, select = col.var)
                write.csv(site_Data, file.path(dir.out, timestep, paste0('SEL_FLX_', ListSites[i], '_FLUXNET2015_', filetype, '.csv')), row.names = F)
                
        } else {nodata <- append(nodata, ListSites[i]) }
}

#------------------------------------------------------------------------------#
# preparing required annual Fluxnet data ####
#------------------------------------------------------------------------------#
reqd.var <- c("TIMESTAMP", 
              "TA_F", 
              "SW_IN_F", 
              "VPD_F", 
              'NEE_VUT_REF',
              'RECO_DT_VUT_REF', 
              'GPP_DT_VUT_REF')

dir.data   <- 'flux_data/Unzipped_Fluxnet2015_database_tier_2'
filetype   <- 'FULLSET_YY'
timestep   <- 'ann'
dir.out    <- "main_analysis/flux.data"
FinalSites <- read.csv('main_analysis/finallist_20180626.csv') 
ListSites  <- as.character(FinalSites$site)

SitesFluxnet2015 <- list.files(dir.data)
nodata <- NULL
for(i in 1:length(ListSites)) {
        if(length(grep(ListSites[i], SitesFluxnet2015)) == 1) {
                
                cat(paste('Reading data for site', ListSites[i]), '\n')
                site_Filename<- list.files(file.path(dir.data, ListSites[i]), pattern = filetype, full.names = T)
                
                dataheader   <- names(read.table(site_Filename, nrow = 1, sep = ',', header = T))
                col.var      <- which(dataheader %in% reqd.var)
                
                site_Data    <- data.table::fread(site_Filename, select = col.var)
                write.csv(site_Data, file.path(dir.out, timestep, paste0('SEL_FLX_', ListSites[i], '_FLUXNET2015_', filetype, '.csv')), row.names = F)
                
        } else {nodata <- append(nodata, ListSites[i]) }
}

#------------------------------------------------------------------------------#
# Plotting half-hourly fluxnet2015 data ####
#------------------------------------------------------------------------------#
reqd.var <- c("TIMESTAMP_START", 
              'NEE_VUT_REF',
              "NEE_CUT_REF_QC" ,
              'RECO_DT_VUT_REF', 
              'GPP_DT_VUT_REF')

dir.data <- 'flux_data/Unzipped_Fluxnet2015_database_tier_2'
filetype <- 'FULLSET_HH'
timestep <- 'HH'
dir.out  <- "main_analysis/flux.data"
FinalSites <- read.csv('main_analysis/finallist_20180626.csv') 
ListSites  <- as.character(FinalSites$site)

dir.out <- 'flux_data/fluxdata_compare'
pdf(file.path(dir.out, 'halfhourly.plots.pdf'), compress = T)
par(mfrow = c(4,3))
par(mar = c(4,4,0.5,0.5))
for(i in 1:length(ListSites)) {
        cat(paste('plotting for', ListSites[i]), '\n')
        site <- ListSites[i]
        dir.site   <- list.files(dir.data, pattern = site, full.names = T)
        if(length(dir.site) == 1) {
                filename   <- list.files(dir.site, pattern = filetype, full.names = T)
                dataheader <- names(read.table(filename, nrow = 1, sep = ',', header = T))
                col.var    <- which(dataheader %in% reqd.var)
                data       <- data.table::fread(filename, select = col.var)
                names(data)<- c('timestamp', 'NEE', 'NEE_QC', 'RECO', 'GPP')
                data$time  <- ymd_hm(data$timestamp)
                data$year  <- year(data$time)
                
                Years <- unique(data$year)
                for(j in 1:length(Years)) {
                        cat(paste('plotting for', ListSites[i], 'Year', Years[j]), '\n')
                        with(data, plot(time[year == Years[j]], NEE[year == Years[j]], ylab = 'NEP', xlab = Years[j], pch = '.'))
                        legend('bottomleft', paste(site, 'NEP'), bty = 'n')
                        with(data, plot(time[year == Years[j]], GPP[year == Years[j]], ylab = 'GPP', xlab = Years[j], pch = '.'))
                        legend('topleft', paste('GPP'), bty = 'n')
                        with(data, plot(time[year == Years[j]], GPP[year == Years[j]], ylab = 'RECO', xlab = Years[j], pch = '.'))
                        legend('topleft', paste('RECO'), bty = 'n')
                }     
        }
        
}
dev.off()

#------------------------------------------------------------------------------#
# Plotting daily fluxnet2015 data ####
#------------------------------------------------------------------------------#
reqd.var <- c("TIMESTAMP", 
              'NEE_VUT_REF',
              "NEE_CUT_REF_QC" ,
              'RECO_DT_VUT_REF', 
              'GPP_DT_VUT_REF')

dir.data   <- 'flux_data/Unzipped_Fluxnet2015_database_tier_2'
filetype   <- 'FULLSET_DD'
timestep   <- 'DD'
dir.out    <- "main_analysis/flux.data"
FinalSites <- read.csv('main_analysis/finallist_20180626.csv') 
ListSites  <- as.character(FinalSites$site)

dir.out <- 'flux_data/fluxdata_compare'
pdf(file.path(dir.out, 'daily.plots.pdf'), compress = T)
par(mfrow = c(4,3),
    mar = c(4,4,0.5,0.5))

for(i in 1:length(ListSites)) {
        cat(paste('plotting for', ListSites[i]), '\n')
        site <- ListSites[i]
        dir.site   <- list.files(dir.data, pattern = site, full.names = T)
        if(length(dir.site) == 1) {
                filename   <- list.files(dir.site, pattern = filetype, full.names = T)
                dataheader <- names(read.table(filename, nrow = 1, sep = ',', header = T))
                col.var    <- which(dataheader %in% reqd.var)
                data       <- data.table::fread(filename, select = col.var)
                names(data)<- c('timestamp', 'NEE', 'NEE_QC', 'RECO', 'GPP')
                data$time  <- ymd(data$timestamp)
                data$year  <- year(data$time)
                
                Years <- unique(data$year)
                for(j in 1:length(Years)) {
                        cat(paste('plotting for', ListSites[i], 'Year', Years[j]), '\n')
                        with(data, plot(time[year == Years[j]], NEE[year == Years[j]], ylab = 'NEP', xlab = Years[j], pch = '.'))
                        legend('bottomleft', paste(site, 'NEP'), bty = 'n')
                        with(data, plot(time[year == Years[j]], GPP[year == Years[j]], ylab = 'GPP', xlab = Years[j], pch = '.'))
                        legend('topleft', paste('GPP'), bty = 'n')
                        with(data, plot(time[year == Years[j]], GPP[year == Years[j]], ylab = 'RECO', xlab = Years[j], pch = '.'))
                        legend('topleft', paste('RECO'), bty = 'n')
                }     
        }
}
dev.off()

#------------------------------------------------------------------------------#
#### Function for plotting annual time series #### can be modified to for smaller timescales easily
#------------------------------------------------------------------------------#
timeseries_plot <- function(FileName, site) {
        data           <- read.csv(FileName)
        data           <- dplyr::select(data, 'TIMESTAMP', 'NEE_VUT_REF', 'RECO_DT_VUT_REF', 'GPP_DT_VUT_REF')
        names(data)    <- c('timestamp', 'nee', 'reco', 'gpp')
        data[data == -9999] <- NA
        data$nep       <- -data$nee
        
        maxVal <- max(data$nep, data$reco, data$gpp, na.rm = T)
        minVal <- min(data$nep, data$reco, data$gpp, na.rm = T)
        
        with(data, plot(timestamp, 
                        nep, 
                        xlab = '', 
                        ylab = '',
                        typ = 'l', 
                        xlim = c(min(data$timestamp) -1, max(data$timestamp) + 1), 
                        ylim = c(minVal, maxVal)
        )
        )
        with(data, points(timestamp, nep, pch = 15))
        
        with(data, points(timestamp, gpp, pch = 16))
        with(data, lines(timestamp, gpp, pch = 16))
        
        with(data, points(timestamp, reco, pch = 17))
        with(data, lines(timestamp, reco, pch = 17))
        
        text(min(data$timestamp) + 1, minVal + 50, site)
        
        mtext(expression( "CO"[2]~ "["~ g ~ C ~ m^{-2}~ yr^{-1}~ ']'), side = 2, line = 2.2, cex = 0.8)
}

#------------------------------------------------------------------------------#
#### PLotting annual flux timeseries ####
#------------------------------------------------------------------------------#
sitelistwithstructuredata <- substr(list.files('main_analysis/strc.data/final_data', full.names = F), 1, 6)
dir.data <- 'main_analysis/flux.data/ann'
dir.out   <- 'main_analysis/flux.data/annual_flux_all_sites'

pdf(file.path(dir.out, 'annual.plots.pdf'))
par(mfrow = c(3,1))
par(mar = c(4,4,0.5,0.5))
for(j in 1:length(sitelistwithstructuredata)) {
        site <- sitelistwithstructuredata[j]
        FileName <- list.files(dir.data, pattern = site, full.names = T)
        if(length(FileName) == 1) {
                timeseries_plot(FileName, site)     
        }
}
dev.off()

#------------------------------------------------------------------------------#
#### PI flux_timeseries preparation function ####
#------------------------------------------------------------------------------#
# could be different for different site
flux_timeseries_annual <- function(FileName, 
                                   dataline
                                   ){
        head.data   <- read.delim(FileName, nrows=1, skip = dataline, sep=',', header=T)
        head.name   <- names(head.data)
        
        col.var     <- which(names(head.data) %in% reqd.var)
        data        <- fread(FileName, select = col.var, skip = dataline, header = T)
        data[data == -999] <- NA
        
        data <- as.data.frame(data)
        
        dailydata <- data %>%
                group_by(year, day) %>%
                summarise(nep  = mean(nee_co2_filled, na.rm = T), 
                          gpp  = mean(gpp, na.rm = T),
                          reco = mean(reco, na.rm = T),
                          non_na_count = sum(is.na(nep))) %>%
                mutate(nep = -nep*1.0368, 
                       gpp = gpp*1.0368,
                       reco = reco*1.0368)
        
        anndata <- dailydata %>%
                group_by(year) %>%
                summarise(nep = sum(nep, na.rm = T),
                          gpp = sum(gpp, na.rm = T),
                          reco = sum(reco, na.rm = T), 
                          na_count = sum(non_na_count))
        
        return(list(dailydata, anndata))
        
} # end of the function flux_timeseries

#------------------------------------------------------------------------------#
#### flux_timeseries for US-WCr ####
#------------------------------------------------------------------------------#
## annual timeseries of WCr ###
dir.data  <- 'Communication_with_PIs/US-Wcr/flux'
ListFiles <- list.files(dir.data, full.names = T)
reqd.var  <- c('year', 'day', 'nee_co2_filled', 'reco', 'gpp')
site      <- 'US-WCr' 

DataDaily <- NULL
DataAnn   <- NULL

for(i in 1:length(ListFiles)) {
        cat('data processing for', ListFiles[i], '\n')
        FileName <- ListFiles[i]
        dataline <- 126
        datareturn <- flux_timeseries_annual (FileName, dataline)
        
        DataDaily <- rbind(DataDaily, datareturn[[1]])
        DataAnn   <- rbind(DataAnn, datareturn[[2]])
}

dir.out <-'flux_data/PI_flux_data'
write.csv(DataDaily, file.path(dir.out, paste0(site,'_PI_Flux_data_daily.csv')), row.names = F)
write.csv(DataAnn, file.path(dir.out, paste0(site,'_PI_Flux_data_ann.csv')), row.names = F)

#------------------------------------------------------------------------------#
#### flux_timeseries for US-Syv ####
#------------------------------------------------------------------------------#
## annual timeseries of US-Syv ###
dir.data  <- 'Communication_with_PIs/US-Syv/flux'
ListFiles <- list.files(dir.data, full.names = T)
reqd.var  <- c('year', 'day', 'nee_co2_filled', 'reco', 'gpp')
site      <- 'US-Syv' 

DataDaily <- NULL
DataAnn   <- NULL

for(i in 1:length(ListFiles)) {
        cat('data processing for', ListFiles[i], '\n')
        FileName <- ListFiles[i]
        dataline <- 114
        datareturn <- flux_timeseries_annual (FileName, dataline)
        
        DataDaily <- rbind(DataDaily, datareturn[[1]])
        DataAnn   <- rbind(DataAnn, datareturn[[2]])
}

dir.out <-'flux_data/PI_flux_data'
write.csv(DataDaily, file.path(dir.out, paste0(site,'_PI_Flux_data_daily.csv')), row.names = F)
write.csv(DataAnn, file.path(dir.out, paste0(site,'_PI_Flux_data_ann.csv')), row.names = F)

#------------------------------------------------------------------------------#
#### flux_timeseries for US-Pfa ####
#------------------------------------------------------------------------------#
dir.data  <- 'Communication_with_PIs/US-Pfa/flux'
ListFiles <- list.files(dir.data, full.names = T)
reqd.var  <- c('year', 'day', 'nee_co2_filled', 'reco', 'gpp')
site      <- 'US-PFa' 

DataDaily <- NULL
DataAnn   <- NULL

for(i in 1:length(ListFiles)) {
        cat('data processing for', ListFiles[i], '\n')
        FileName <- ListFiles[i]
        dataline <- 124
        datareturn <- flux_timeseries_annual (FileName, dataline)
        
        DataDaily <- rbind(DataDaily, datareturn[[1]])
        DataAnn   <- rbind(DataAnn, datareturn[[2]])
}

dir.out <-'flux_data/PI_flux_data'
write.csv(DataDaily, file.path(dir.out, paste0(site,'_PI_Flux_data_daily.csv')), row.names = F)
write.csv(DataAnn, file.path(dir.out, paste0(site,'_PI_Flux_data_ann.csv')), row.names = F)

#------------------------------------------------------------------------------#
#### flux_timeseries for US-Ha1 ####
#------------------------------------------------------------------------------#
FileName <- 'Communication_with_PIs/US-Ha/US-Ha1-10m10m/HF_9215_filled.txt'
site     <- 'US-Ha1'
head.data   <- read.delim(FileName, nrows=1, sep='\t', header=T)
head.name   <- names(head.data)

reqd.var <- c('Year.Year', 'DoY.Day', 'nee.e.6mol.m2.s', 'Resp.e.e.6mol.m2.s', 'gee.e.6mol.m2.s')

col.var     <- which(names(head.data) %in% reqd.var)
data        <- fread(FileName, select = col.var, header = T)

data <- as.data.frame(data)
names(data) <-  c('year', 'day', 'nee_co2_filled', 'reco', 'gpp')

dailydata <- data %>%
        group_by(year, day) %>%
        summarise(nep  = mean(nee_co2_filled, na.rm = T), 
                  gpp  = mean(gpp, na.rm = T),
                  reco = mean(reco, na.rm = T),
                  non_na_count = sum(is.na(nep))) %>%
        mutate(nep = -nep*1.0368, 
               gpp = gpp*1.0368,
               reco = reco*1.0368)

anndata <- dailydata %>%
        group_by(year) %>%
        summarise(nep = sum(nep, na.rm = T),
                  gpp = sum(gpp, na.rm = T),
                  reco = sum(reco, na.rm = T), 
                  na_count = sum(non_na_count))

dir.out <-'flux_data/PI_flux_data'
write.csv(dailydata, file.path(dir.out, paste0(site,'_PI_Flux_data_daily.csv')), row.names = F)
write.csv(anndata, file.path(dir.out, paste0(site,'_PI_Flux_data_ann.csv')), row.names = F)

#------------------------------------------------------------------------------#
#### flux_timeseries for DE-Hai  ####
#------------------------------------------------------------------------------#
FileName <- 'Communication_with_PIs/DE-Hai/daily_flux.csv'
site     <- 'DE-Hai'
data     <- read.csv(FileName)
data$timestamp <- dmy(data$timestamp)
data$year <- year(data$timestamp)
data$day  <- yday(data$timestamp)
data$non_na_count <- 0
DataDaily <- dplyr::select(data, year, day, nep, gpp, reco, non_na_count)
dir.out <-'flux_data/PI_flux_data'
write.csv(DataDaily, file.path(dir.out, paste0(site,'_PI_Flux_data_daily.csv')), row.names = F)

#------------------------------------------------------------------------------#
#### flux_timeseries for DE-Lnf  ####
#------------------------------------------------------------------------------#
FileName <- 'Communication_with_PIs/DE-Lnf/daily_flux.csv'
site     <- 'DE-Lnf'
data     <- read.csv(FileName)
data$timestamp <- dmy(data$timestamp)
data$year <- year(data$timestamp)
data$day  <- yday(data$timestamp)
data$non_na_count <- 0
DataDaily <- dplyr::select(data, year, day, nep, gpp, reco, non_na_count)
dir.out <-'flux_data/PI_flux_data'
write.csv(DataDaily, file.path(dir.out, paste0(site,'_PI_Flux_data_daily.csv')), row.names = F)

#------------------------------------------------------------------------------#
#### compare annual PI_flux series with Fluxnet ####
#------------------------------------------------------------------------------#
{
dir.FluxNet2015 <- 'main_analysis/flux.data/ann'
dir.PIFlux      <- 'flux_data/PI_flux_data'

ListSites <- substr(list.files(dir.PIFlux, pattern = 'ann'), 1, 6)

var <- c('nep', 'gpp', 'reco')
dir.out <- 'flux_data/fluxdata_compare'
pdf(file.path(dir.out, 'annual_compare.plots.pdf'))
par(mfrow = c(3,3))
par(mar = c(4,4,0.5,0.5))
for(i in 1:length(ListSites)) {
        File.FluxNet2015 <- list.files(dir.FluxNet2015, pattern = ListSites[i], full.names = T)
        data.FluxNet2015 <- read.csv(File.FluxNet2015)
        data.FluxNet2015 <- dplyr::select(data.FluxNet2015, 
                                          year = TIMESTAMP,
                                          nep = NEE_VUT_REF,
                                          reco = RECO_DT_VUT_REF,
                                          gpp = GPP_DT_VUT_REF)
        data.FluxNet2015 <- dplyr::mutate(data.FluxNet2015, nep = -nep)
        
        File.PIFlux     <- file.path(dir.PIFlux, paste0(ListSites[i], '_PI_Flux_data_ann.csv'))
        data.PIFlux     <- read.csv(File.PIFlux)
        
        data <- merge(data.PIFlux, data.FluxNet2015, by = 'year')
        
        data <- dplyr::filter(data, na_count == 0)
        
        for(j in 1:length(var)) {
                dataVar <- data[ ,c(paste0(var[j], '.x'), paste0(var[j], '.y'))]
                names(dataVar) <- c('PIFlux', 'FluxNet2015')
                
                lmData <- lmodel2(FluxNet2015~PIFlux, dataVar, nperm = 99)
                slope  <- round(lmData$regression.results[2,3], 2)
                intercept <- round(lmData$regression.results[2,2], 2)
                
                with(dataVar, plot(PIFlux, FluxNet2015, 
                                   xlab = 'PI', 
                                   ylab = 'FluxNet2015', 
                                   ylim = c(min(FluxNet2015, PIFlux), 
                                            max(FluxNet2015, PIFlux)), 
                                   xlim = c(min(FluxNet2015, PIFlux), 
                                            max(FluxNet2015, PIFlux))
                ))
                
                legend('topleft',paste(ListSites[i], var[j]) , bty = 'n')
                abline(intercept, slope, lty = 2)
                abline(0, 1, lty = 2, col = 2)
                
                legend('bottomright', paste('Y = ', intercept, '+', slope, 'X'), bty = 'n')
        }
}
dev.off()
}

#------------------------------------------------------------------------------#
#### compare daily PI_flux series with Fluxnet ####
#------------------------------------------------------------------------------#
dir.FluxNet2015 <- 'main_analysis/flux.data/day'
dir.PIFlux      <- 'flux_data/PI_flux_data'

ListSites <- substr(list.files(dir.PIFlux, pattern = 'daily'), 1, 6)

var <- c('nep', 'gpp', 'reco')
dir.out <- 'flux_data/fluxdata_compare'
pdf(file.path(dir.out, 'daily_compare.plots.pdf'))
par(mfrow = c(3,3))
par(mar = c(4,4,0.5,0.5))
for(i in 1:length(ListSites)) {
        File.FluxNet2015 <- list.files(dir.FluxNet2015, pattern = ListSites[i], full.names = T)
        data.FluxNet2015 <- read.csv(File.FluxNet2015)
        data.FluxNet2015 <- dplyr::select(data.FluxNet2015, 
                                          TIMESTAMP,
                                          nep = NEE_VUT_REF,
                                          reco = RECO_DT_VUT_REF,
                                          gpp = GPP_DT_VUT_REF)
        data.FluxNet2015[data.FluxNet2015 == -9999] <- NA
        data.FluxNet2015 <- na.omit(data.FluxNet2015)
        data.FluxNet2015 <- dplyr::mutate(data.FluxNet2015, nep = -nep)
        data.FluxNet2015$timestamp <- ymd(data.FluxNet2015$TIMESTAMP)
        data.FluxNet2015$year <- year(data.FluxNet2015$timestamp)
        data.FluxNet2015$day  <- yday(data.FluxNet2015$timestamp)
        
        
        File.PIFlux     <- file.path(dir.PIFlux, paste0(ListSites[i], '_PI_Flux_data_daily.csv'))
        data.PIFlux     <- read.csv(File.PIFlux)
        
        data <- merge(data.PIFlux, data.FluxNet2015, by = c('year', 'day'))
        
        data <- dplyr::filter(data, non_na_count == 0)
        
        for(j in 1:length(var)) {
                dataVar <- data[ ,c(paste0(var[j], '.x'), paste0(var[j], '.y'))]
                names(dataVar) <- c('PIFlux', 'FluxNet2015')
                
                lmData <- lmodel2(FluxNet2015~PIFlux, dataVar, nperm = 99)
                slope  <- round(lmData$regression.results[2,3], 2)
                intercept <- round(lmData$regression.results[2,2], 2)
                
                with(dataVar, plot(PIFlux, FluxNet2015, 
                                   xlab = 'PI', 
                                   ylab = 'FluxNet2015', 
                                   ylim = c(min(FluxNet2015, PIFlux), 
                                            max(FluxNet2015, PIFlux)), 
                                   xlim = c(min(FluxNet2015, PIFlux), 
                                            max(FluxNet2015, PIFlux)),
                                   pch = '.'
                ))
                
                legend('topleft',paste(ListSites[i], var[j]) , bty = 'n')
                abline(intercept, slope, lty = 2)
                abline(0, 1, lty = 2, col = 2)
                
                legend('bottomright', 
                       paste('Y = ', intercept, '+', slope, 'X'), 
                       bty = 'n')
        }
}
dev.off()

#------------------------------------------------------------------------------#
#### Function for calculating annual averages for CO2 fluxes ####
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
        
        data2      <- dplyr::filter(data, !timestamp %in% Problematic.Year) # removing problematic year
        
        ReqdVar    <- c('nep', 'reco', 'gpp')                                   # names of variable for result output
        average    <- round(colMeans2(as.matrix(data2[,ReqdVar]), na.rm = T), 1)# calculating average using MatixStats package
        stdDev     <- round(colSds(as.matrix(data2[,ReqdVar]), na.rm = T), 1)   # calculating standard deviation using MatixStats package
        CoefVar    <- round(100 * stdDev/abs(average), 1)                       # calculating coefficient of variation using MatixStats package
        results    <- data.frame(ReqdVar, average, stdDev, CoefVar)             # creating a dataframe of results
        results$no.of.year <- nrow(data2)                                       # no of years with data available
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
#### calculating average values for carbon dioxide fluxes ####
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
#### Calculating corelation between structural index and co2 fluxes at daily scale  ####
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
#### Calculating corelation between structural index and co2 fluxes at annual scale  ####
#------------------------------------------------------------------------------#
{
stru.file <- 'main_analysis/strc.data/structural.index.one.value.csv'
anom.file <- 'main_analysis/flux.data/annual_flux_all_sites/annualco2fluxesstatistics.csv'
dir.out   <- 'main_analysis/CorrelationBetweenAnoAndStruc'

stru.data <- read.csv(stru.file)
anom.data <- read.csv(anom.file)
Var       <- c('nep', 'gpp', 'reco')
YVar      <- c("average", 
               "CoefVar")

XVar       <- c("meanDBH",
                "CVDBH",
                "ShannonSize",
                "ShannonSps",
                "shapeWeibull",
                "TotalSpeciesRichnessBA95")

resultsCorr <- NULL
for(i in 1:length(Var)) {
        anom.data.var <- dplyr::filter(anom.data, ReqdVar == Var[i])
        data <- merge(anom.data.var, stru.data, by = 'site')
        
        for(j in 1:length(YVar)) {
                PlotName <- paste(Var[i], 
                                  YVar[j], 
                                  'Correlation.structural.Index.annual.jpeg', 
                                  sep = '.')
                jpeg(file=file.path(dir.out, PlotName),
                     width  = 90,
                     height = 160, 
                     units  = 'mm', 
                     res    = 300
                     )
                par(oma = c(0,4,0.5,0.5), 
                    mar = c(4,0,0,0),
                    mfrow=c(3,2)
                    )
                
                for(k in 1:length(XVar)) {
                        dataYXvar <- data[, c(YVar[j], XVar[k])]
                        names(dataYXvar) <- c('YVar', 'XVar')
                        minY <- min(dataYXvar$YVar)
                        maxY <- max(dataYXvar$YVar)
                        
                        corrResults <- cor.test(dataYXvar$YVar, dataYXvar$XVar)
                        pVal        <- corrResults$p.value
                        corrValue   <- cor(dataYXvar$YVar, dataYXvar$XVar)
                        
                        resultsCorr <- rbind(resultsCorr, 
                                             c('Var' = Var[i],
                                               'YVar' = YVar[j],
                                               'XVar' = XVar[k],
                                               'corrVal' = corrValue,
                                               'Rsq' = round(corrValue^2, 2), 
                                               'pVal' = pVal))
                        
                        with(dataYXvar, plot(XVar, 
                                             YVar, 
                                             yaxt = 'n', 
                                             xlab = '', 
                                             ylab = '')
                        )
                        if(k %in% c(1,3,5)) 
                                axis(2, 
                                     at = round(seq(minY, maxY, by = (maxY - minY)/4)), 
                                     labels = T)
                        
                        mtext(XVar[k], outer = F, side = 1,  line = 2.2, cex = 0.8)
                        legend('topright', paste('corr =', round(corrValue, 2)), bty = 'n')
                        
                } # end for all dependent variable
                mtext(paste(YVar[j], Var[i]), outer = T, side = 2,  line = 2.2, cex = 0.8)
                dev.off()
        } # end for all independent variables Yvar
} # end for all carbon dioxide variables
FileName <- paste(Var[i], 
                  YVar[j], 
                  'Correlation.structural.Index.annual.csv', 
                  sep = '.')
write.csv(resultsCorr, file.path(dir.out, FileName), row.names = F)
}

#------------------------------------------------------------------------------#
#### Function: daily timeseries decomposition ####
#------------------------------------------------------------------------------#
ts_decompositon_daily <- function(FileName, 
                                  Site, 
                                  Problematic.Year 
                                  ) {
        data           <- read.csv(FileName)
        data           <- dplyr::select(data, 
                                        'TIMESTAMP', 
                                        'NEE_VUT_REF', 
                                        'RECO_DT_VUT_REF', 
                                        'GPP_DT_VUT_REF'
                                        )
        
        names(data)    <- c('timestamp', 'nee', 'reco', 'gpp')
        data$timestamp <- ymd(data$timestamp, tz = 'GMT')
        data$year      <- year(data$timestamp)
        
        data$nee       <- with(data, ifelse(year %in% Problematic.Year, NA, data$nee))
        data$reco      <- with(data, ifelse(year %in% Problematic.Year, NA, data$reco))
        data$gpp       <- with(data, ifelse(year %in% Problematic.Year, NA, data$gpp))
        
        data$doy       <- yday(data$timestamp)
        data           <- dplyr::arrange(data, year, doy)
        data[data == -9999] <- NA
        data$nep       <- -data$nee
        
        # doing data analysis for growing season
        data <- dplyr::filter(data, doy > 120, doy < 301)
        startyear <- data[1, 'year']
        endyear   <- data[nrow(data), 'year']
        Freq = 180
        var <- c('nep', 'gpp', 'reco') 
        DataReturn <- data.frame('year' = rep(c(startyear:endyear), each = Freq), 
                                 'doy' = rep(c(121:300), length(unique(c(startyear:endyear))))
        )
        for(i in 1:length(var)) {
                dataVar <- data[,c('year', 'doy', var[i])] 
                names(dataVar) <- c('year', 'doy', 'var')
                
                #--------------------------------------------------------------#
                # calculate trend based using annual averages
                datats <- ts(dataVar$var, start=c(startyear,1,1), frequency=Freq)
                trd  <- try(Trend(datats, method="AAT"))
                pval <- trd$pval 
                
                if(pval <= 0.1) {
                        DataTrend <- data.frame('year' = rep(c(startyear:endyear), 1))
                        DataTrend$trend <- trd$trend 
                        dataVar <- merge(dataVar, DataTrend, by = 'year')
                } else {
                        dataVar$trend <- mean(dataVar$var, na.rm = T)
                }
                
                # detrending the data 
                dataVar$detrendVar <- dataVar$var-dataVar$trend
                # calculating seasonal cycle after removing trend
                seasonalcycle <- dataVar %>%
                        group_by(doy) %>%
                        summarise(seasonalVar = mean(detrendVar, na.rm = T))
                # merging seasonal data with main dataset
                dataVar<- merge(dataVar, 
                                seasonalcycle, 
                                by='doy', 
                                all = T)
                
                # random component calculated after trend and seasonal cycle is removed
                dataVar$random     <- dataVar$var - dataVar$trend - dataVar$seasonalVar
                # z score of random component
                dataVar$zrandom    <- (dataVar$random - mean(dataVar$random, na.rm = T))/sd(dataVar$random, na.rm = T)
                # absolute z score of random component
                dataVar$abszrandom <- abs(dataVar$zrandom)
                
                #--------------------------------------------------------------#
                # seasonal cycle two doesn't consider the trend in data
                seasonalcycle2 <- dataVar %>%
                        group_by(doy) %>%
                        summarise(seasonalVar2 = mean(var, na.rm = T))
                
                dataVar <- merge(dataVar, seasonalcycle2, by='doy', all = T)
                dataVar$random2 <- dataVar$var - dataVar$seasonalVar2
                dataVar$zrandom2 <- (dataVar$random2 - mean(dataVar$random2, na.rm = T))/sd(dataVar$random2, na.rm = T)
                dataVar$abszrandom2 <- abs(dataVar$zrandom2)
                #--------------------------------------------------------------#
                names(dataVar) <- c('doy', 
                                    'year', 
                                    var[i], 
                                    paste0(var[i], 'trend'), 
                                    paste0(var[i], 'detrended'),
                                    paste0(var[i], 'seasonal'),
                                    paste0(var[i], 'random'),
                                    paste0(var[i], 'zrandom'),
                                    paste0(var[i], 'abszrandom'),
                                    paste0(var[i], 'seasonal2'),
                                    paste0(var[i], 'random2'),
                                    paste0(var[i], 'zrandom2'),
                                    paste0(var[i], 'abszrandom2')
                )
                
                DataReturn <- merge(DataReturn, dataVar, by = c('year', 'doy'))
        } # end for analysis for ith co2 component
        DataReturn[is.na(DataReturn)] <- -9999
        DataReturn$site <- Site
        return(DataReturn)
} # end of the funtion

#------------------------------------------------------------------------------#
#### Daily ts decomposition ####
#------------------------------------------------------------------------------#
sitelistwithstructuredata <- substr(list.files('main_analysis/strc.data/final_data', full.names = F), 1, 6)

dir.data  <- 'main_analysis/flux.data/day'
dir.out   <- 'main_analysis/anomaly_detection/daily/ts_decomposed_data'

data.problem <- read.csv('main_analysis/flux.data/Data_with_problem.csv')

for(j in 1:length(sitelistwithstructuredata)) {
        Site <- sitelistwithstructuredata[j]
        FileName <- list.files(dir.data, pattern = Site, full.names = T)
        if(length(FileName) == 1) {
                data.problem.for.site <- data.problem %>%
                        filter(site == Site) %>%
                        filter(dont_use_annual == 1)
                Problematic.Year <- data.problem.for.site$Not_good_year
                site.data <- ts_decompositon_daily (FileName, Site, Problematic.Year) 
                write.csv(site.data, file.path(dir.out, paste0(Site, '_flux_ts_decomposed.csv')), row.names = F)
        }
}

#------------------------------------------------------------------------------#
#### Anomaly calculation function ####
#------------------------------------------------------------------------------#
anomaly_calculation <- function(FileName, site, randomtype, Plot, Freq, dir.plot) {
        data <- read.csv(FileName)
        data[data == -9999] <- NA
        var  <- c('nep', 'gpp', 'reco')
        ReturnResults <- NULL
        for(i in 1:length(var)) {
                cat('calculating extreme data of', var[i], 'for site', site, '\n')
                
                if(randomtype == 'TrendRemoved') {
                        reqdVar <- c('year', 'doy', paste0(var[i], c('', 'trend', 'seasonal', 'random', 'zrandom', 'abszrandom')))
                } else {
                        reqdVar <- c('year', 'doy', paste0(var[i], c('', 'trend', 'seasonal2', 'random2', 'zrandom2', 'abszrandom2')))
                }
                
                dataVar     <- data[ ,reqdVar]
                names(dataVar) <- c('year', 'doy', 'var', 'trend', 'seasonal', 'random', 'zrandom', 'abszrandom')
                
                if(randomtype != 'TrendRemoved') {
                      dataVar$trend <- mean(dataVar$var, na.rm = T)
                } 
                
                dataVar2 <- na.omit(dataVar)
                
                if(var[i] == 'reco') {
                        maxPositiveZscore  <- abs(min(dataVar2$zrandom, na.rm = T))
                        quantile_995z      <- as.numeric(quantile(abs(dataVar2$zrandom[dataVar2$zrandom < 0]), 0.995, na.rm = T))
                        HigherNegativeTail995 <- -quantile_995z + max(dataVar2$zrandom, na.rm = T)
                        areaUnderHigherNegativeTail995 <- 100* sum(abs(dataVar2$zrandom[dataVar2$zrandom > quantile_995z]))/sum(dataVar2$abszrandom)
                        daysUnderHigherNegativeTail995 <- length(dataVar2$zrandom[dataVar2$zrandom > quantile_995z])
                        PerCentdaysUnderHigherNegativeTail995 <- 100* length(dataVar2$zrandom[dataVar2$zrandom > quantile_995z])/nrow(dataVar2)
                        
                        HigherNegativeTail100 <- -maxPositiveZscore + max(dataVar2$zrandom)
                        areaUnderHigherNegativeTail100 <- 100* sum(abs(dataVar2$zrandom[dataVar2$zrandom > maxPositiveZscore]))/sum(dataVar2$abszrandom)
                        daysUnderHigherNegativeTail100 <- length(dataVar2$zrandom[dataVar2$zrandom > maxPositiveZscore])
                        PerCentdaysUnderHigherNegativeTail100 <- 100* length(dataVar2$zrandom[dataVar2$zrandom > maxPositiveZscore])/nrow(dataVar2)
                } else {
                        maxPositiveZscore <- max(dataVar2$zrandom, na.rm = T)
                        quantile_995z      <- as.numeric(quantile(dataVar2$zrandom[dataVar2$zrandom > 0], 0.995, na.rm = T))
                        HigherNegativeTail995 <- -quantile_995z - min(dataVar2$zrandom)
                        areaUnderHigherNegativeTail995 <- 100* sum(abs(dataVar2$zrandom[dataVar2$zrandom < -quantile_995z]))/sum(dataVar2$abszrandom)
                        daysUnderHigherNegativeTail995 <- length(dataVar2$zrandom[dataVar2$zrandom < -quantile_995z])
                        PerCentdaysUnderHigherNegativeTail995 <- 100 * length(dataVar2$zrandom[dataVar2$zrandom < -quantile_995z])/nrow(dataVar2)
                        
                        HigherNegativeTail100 <- -maxPositiveZscore - min(dataVar2$zrandom)
                        areaUnderHigherNegativeTail100 <- 100* sum(abs(dataVar2$zrandom[dataVar2$zrandom < -maxPositiveZscore]))/sum(dataVar2$abszrandom)
                        daysUnderHigherNegativeTail100 <- length(dataVar2$zrandom[dataVar2$zrandom < -maxPositiveZscore])
                        PerCentdaysUnderHigherNegativeTail100 <- 100 * length(dataVar2$zrandom[dataVar2$zrandom < -maxPositiveZscore]) / nrow(dataVar2)
                }
                
                ReturnResults <- rbind(ReturnResults, c('site' = site, 
                                                        'var' = var[i],
                                                        'HigherNegativeTail995' = HigherNegativeTail995,
                                                        'areaUnderHigherNegativeTail995' = areaUnderHigherNegativeTail995,
                                                        'daysUnderHigherNegativeTail995' = daysUnderHigherNegativeTail995,
                                                        'PerCentdaysUnderHigherNegativeTail995' = PerCentdaysUnderHigherNegativeTail995,
                                                        'HigherNegativeTail100' = HigherNegativeTail100,
                                                        'areaUnderHigherNegativeTail100' = areaUnderHigherNegativeTail100,
                                                        'daysUnderHigherNegativeTail100' = daysUnderHigherNegativeTail100,
                                                        'PerCentdaysUnderHigherNegativeTail100' = PerCentdaysUnderHigherNegativeTail100
                                                        )
                                       ) # end of rbind
                
                # plot decomposed data
                if(Plot == 'YES') {
                        cat(paste('creating plot of', var[i], 'for site', site), '\n')
                        plot.ts(ts(dataVar$var, start=c(startyear,1,1), frequency=Freq), lty = 2, xaxt = 'n')
                        if (var[i] == 'nep') mtext(expression(NEP ~ "["~g ~ C ~ m^{-2} ~ wk^{-1} ~ ']'), side = 2, outer = F, cex =0.7, line = 2.2, col = "black")
                        if (var[i] == 'gpp') mtext(expression(GPP ~ "["~g ~ C ~ m^{-2} ~ wk^{-1} ~ ']'), side = 2, outer = F, cex =0.7, line = 2.2, col = "black")
                        if (var[i] == 'reco') mtext(expression(Reco ~ "["~g ~ C ~ m^{-2} ~ wk^{-1} ~ ']'), side = 2, outer = F, cex =0.7, line = 2.2, col = "black")
                        lines(ts(dataVar$trend, start=c(startyear,1,1), frequency=Freq), col = 2)
                        legend('topleft', paste(site, var[i], randomtype), bty = 'n')
                        
                        plot.ts(ts(dataVar$seasonal, start=c(startyear,1,1), frequency=Freq), lty = 1, col=1, xaxt = 'n')
                        if (var[i] == 'nep') mtext(expression(NEP[s] ~ "["~g ~ C ~ m^{-2} ~ wk^{-1} ~ ']'), side = 2, outer = F, cex =0.7, line = 2.2, col = "black")
                        if (var[i] == 'gpp') mtext(expression(GPP[s] ~ "["~g ~ C ~ m^{-2} ~ wk^{-1} ~ ']'), side = 2, outer = F, cex =0.7, line = 2.2, col = "black")
                        if (var[i] == 'reco') mtext(expression(Reco[s] ~ "["~g ~ C ~ m^{-2} ~ wk^{-1} ~ ']'), side = 2, outer = F, cex =0.7, line = 2.2, col = "black")
                        
                        plot.ts(ts(dataVar$random, start=c(startyear,1,1), frequency=Freq), col=1, typ = 'p', pch = '.')
                        if (var[i] == 'nep') mtext(expression(NEP[r] ~ "["~g ~ C ~ m^{-2} ~ wk^{-1} ~ ']'), side = 2, outer = F, cex =0.7, line = 2.2, col = "black")
                        if (var[i] == 'gpp') mtext(expression(GPP[r] ~ "["~g ~ C ~ m^{-2} ~ wk^{-1} ~ ']'), side = 2, outer = F, cex =0.7, line = 2.2, col = "black")
                        if (var[i] == 'reco') mtext(expression(Reco[r] ~ "["~g ~ C ~ m^{-2} ~ wk^{-1} ~ ']'), side = 2, outer = F, cex =0.7, line = 2.2, col = "black")
                        abline(h=0, col = 2)
                        
                        if(var[i] == 'reco') {
                                maxPositiveZscore <- -min(dataVar$random, na.rm = T)
                                quantile_995      <- as.numeric(quantile(abs(dataVar$random[dataVar$random < 0]), 0.995, na.rm = T))
                                dataRandom        <- ts(dataVar$random, start=c(startyear,1,1), frequency=Freq)
                                dataRandom995     <- ifelse(dataRandom < quantile_995, NA, dataRandom)
                                dataRandommax     <- ifelse(dataRandom < maxPositiveZscore, NA, dataRandom)   
                        } else {
                                maxPositiveZscore <- max(dataVar$random, na.rm = T)
                                quantile_995      <- as.numeric(quantile(dataVar$random[dataVar$random > 0], 0.995, na.rm = T))
                                dataRandom        <- ts(dataVar$random, start=c(startyear,1,1), frequency=Freq)
                                dataRandom995     <- ifelse(dataRandom > -quantile_995, NA, dataRandom)
                                dataRandommax     <- ifelse(dataRandom > -maxPositiveZscore, NA, dataRandom)
                        }
                        
                        points(dataRandom995, col=3, typ = 'p')
                        points(dataRandommax, col=2, typ = 'p')
                } # end of plot when requested 
        } # end of anomaly calculation for ith co2 var
        return(ReturnResults)
}

#------------------------------------------------------------------------------#
#### Daily Anomaly calculation  ####
#------------------------------------------------------------------------------#
dir.data <- 'main_analysis/anomaly_detection/daily/ts_decomposed_data'
sitelist <- substr(list.files(dir.data, full.names = F), 1, 6)
dir.plot <- 'main_analysis/anomaly_detection/daily/figures'
dir.out  <- 'main_analysis/anomaly_detection/daily'
Plot <- 'YES'
Freq <- 180
resultsAll <- NULL

{
randomtype <- 'TrendRemoved' # 'TrendKept'
pdf(file.path(dir.out, 'dailyTrendRemovedAnomaly.plots.pdf'))
par(mfrow=c(3,1), oma=c(4,4,0.5,0.5), mar = c(0,0,0,0))

for(j in 1:length(sitelist)) {
        site <- sitelist[j]
        FileName <- list.files(dir.data, pattern = site, full.names = T)
        if(length(FileName) == 1) {
                
                site.data <- anomaly_calculation (FileName, 
                                                  site, 
                                                  randomtype,
                                                  Plot, 
                                                  Freq, 
                                                  dir.plot)
                resultsAll<- rbind(resultsAll, site.data)
        }
}
dev.off()
write.csv(resultsAll, file.path(dir.out, 'daily_anomalies_Trend_Removed.csv'), row.names = F)
}

{
        randomtype <- 'TrendKept' # 'TrendRemoved'
        pdf(file.path(dir.out, 'dailyTrendKeptAnomaly.plots.pdf'))
        par(mfrow=c(3,1), oma=c(4,4,0.5,0.5), mar = c(0,0,0,0))
        
        for(j in 1:length(sitelist)) {
                site <- sitelist[j]
                FileName <- list.files(dir.data, pattern = site, full.names = T)
                if(length(FileName) == 1) {
                        
                        site.data <- anomaly_calculation (FileName, 
                                                          site, 
                                                          randomtype,
                                                          Plot, 
                                                          Freq, 
                                                          dir.plot)
                        resultsAll<- rbind(resultsAll, site.data)
                }
        }
        dev.off()
        write.csv(resultsAll, file.path(dir.out, 'daily_anomalies_Trend_Kept.csv'), row.names = F)
}
