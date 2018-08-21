FileName <- 'main_analysis/finallist_Short_20180712ver2.csv'  
data     <- read.csv(FileName)
sitelistwithstructuredata <- sort(as.character(unique(data$site)))

#------------------------------------------------------------------------------#
#### calculating average values for carbon dioxide fluxes ####
#------------------------------------------------------------------------------#
dir.data     <- 'main_analysis/flux.data/ann'
dir.out      <- 'main_analysis/flux.data/annual_flux_all_sites'
data.problem <- read.csv('main_analysis/flux.data/Data_with_problem.csv')
AllResults   <- NULL
dir.out      <- 'main_analysis/flux.data/annual_flux_all_sites'

pdf(file.path(dir.out, 'annual.plots.pdf'))
par(mfrow = c(3,1))
par(mar = c(4,4,0.5,0.5))
for(j in 1:length(sitelistwithstructuredata)) {
        Site <- sitelistwithstructuredata[j]
        FileName <- list.files(dir.data, pattern = Site, full.names = T)
        if(length(FileName) == 1) {
                data.problem.for.site <- data.problem %>%
                        filter(site == Site) %>%
                        filter(dont_use_annual == 1)
                
                Problematic.Year <- data.problem.for.site$Not_good_year
                PlotData <- 'YES'
                timestep <- 'ann'
                
                Res <- timeseries_avg(FileName, Problematic.Year, Site, timestep, PlotData)  
                AllResults <- rbind(AllResults, Res)
        }
}
dev.off()
dir.out <- 'main_analysis/flux.data/annual_flux_all_sites'
write.csv(AllResults, file.path(dir.out, 'annualco2fluxesstatistics.csv'), row.names = F)


#------------------------------------------------------------------------------#
#### combine with other data ####
#------------------------------------------------------------------------------#
annual_data <- read.csv( 'main_analysis/flux.data/annual_flux_all_sites/annualco2fluxesstatistics.csv')
annual_data <- dplyr::select(annual_data, -timestep, -no.of.year)
annual_data <- reshape(annual_data, idvar = "site", timevar = "ReqdVar", direction = "wide")

data <- merge(data, annual_data, by = 'site', all = T)
data$age2 <- log(data$Age)

#------------------------------------------------------------------------------#
#### Calculating corelation between structural index and gppsat  ####
#------------------------------------------------------------------------------#
{
        dir.out   <- 'main_analysis/CorrelationBetweenAnoAndStruc'
        pdf(paste0(dir.out, '/annualEverything.pdf'), width = 6, height = 9, compress = T)
        par(oma = c(0,4,3,3), 
            mar = c(4,0,0,0),
            mfrow=c(3,3)
        )
        
        YVar      <- c('average.nep',
                       'stdDev.nep',
                       'CoefVar.nep',
                       'average.gpp',
                       'stdDev.gpp',
                       'CoefVar.gpp',
                       'average.reco',
                       'stdDev.reco',
                       'CoefVar.reco',
                       "AvgAmax", 
                       "SdAmax",
                       'CvAmax',
                       'AvgGPPsat',
                       'SdGPPsat',
                       'CvGPPsat')
        
        XVar       <- c('MAT',
                        'MAP',
                        'age2',
                        "TBAm2ha",
                        "meanDBH",
                        "CVDBH",
                        "ShannonSize",
                        "shapeWeibull",
                        "TotalSpeciesRichnessBA")
        
        resultsCorr <- NULL
        forest_type <- as.character(unique(data$IGBP2))
        
        for(j in 1:length(YVar)) {
                
                for(k in 1:length(XVar)) {
                        dataYXvar <- data[, c('IGBP2', 'Lat', YVar[j], XVar[k])]
                        names(dataYXvar) <- c('IGBP2', 'Lat','YVar', 'XVar')
                        dataYXvar <- na.omit(dataYXvar)
                        minY <- min(dataYXvar$YVar, na.rm = T)
                        maxY <- max(dataYXvar$YVar, na.rm = T)
                        
                        corrResults <- cor.test(dataYXvar$YVar, dataYXvar$XVar)
                        pVal        <- round(corrResults$p.value, 3)
                        corrValue   <- round(cor(dataYXvar$YVar, dataYXvar$XVar), 3)
                        Rsq         <- round(corrValue^2, 2)
                        
                        resultsCorr <- rbind(resultsCorr, 
                                             c('YVar' = YVar[j],
                                               'XVar' = XVar[k],
                                               'corrVal' = corrValue,
                                               'Rsq' = Rsq, 
                                               'pVal' = pVal,
                                               'FT' = 'all'))
                        
                        with(dataYXvar, plot(XVar, 
                                             YVar, 
                                             yaxt = 'n', 
                                             xlab = '', 
                                             ylab = '', 
                                             pch  = 16)
                        )
                        
                        axis(2, 
                             at = round(seq(minY, maxY, by = (maxY - minY)/4), 2), 
                             labels = F)
                        
                        if(k %in% c(1,4,7)) 
                                axis(2, 
                                     at = round(seq(minY, maxY, by = (maxY - minY)/4), 2), 
                                     labels = T)
                        
                        mtext(XVar[k], outer = F, side = 1,  line = 2.2, cex = 0.8)
                        legend('topright', paste('corr =', round(corrValue, 2), '(p =', pVal, ')'), bty = 'n')
                        
                        for(i in 1:length(forest_type)) {
                                with(dataYXvar, points(XVar[IGBP2 == forest_type[i]], 
                                                       YVar[IGBP2 == forest_type[i]], 
                                                       col = i,
                                                       pch = 16)
                                     ) 
                                
                                
                                dataYXvar_FT <- dplyr::filter(dataYXvar, IGBP2 == forest_type[i])
                                
                                corrResults <- cor.test(dataYXvar_FT$YVar, dataYXvar_FT$XVar)
                                pVal        <- round(corrResults$p.value, 3)
                                corrValue   <- round(cor(dataYXvar_FT$YVar, dataYXvar_FT$XVar), 3)
                                Rsq         <- round(corrValue^2, 2)
                                
                                resultsCorr <- rbind(resultsCorr, 
                                                     c('YVar' = YVar[j],
                                                       'XVar' = XVar[k],
                                                       'corrVal' = corrValue,
                                                       'Rsq' = Rsq, 
                                                       'pVal' = pVal,
                                                       'FT' = forest_type[i]))
                        }
                        if(k == 1) legend('topleft', forest_type, pch = 16, col = 1:length(forest_type), bty = 'n', ncol = 1)
                
                } # end for all dependent variable
                mtext(paste(YVar[j]), outer = T, side = 2,  line = 2.2, cex = 0.8)
                
        } # end for all independent variables Yvar
        
        dev.off()
        FileName <- 'Correlation.structural.Index.Allsat.csv'
        write.csv(resultsCorr, file.path(dir.out, FileName), row.names = F)
}


#------------------------------------------------------------------------------#
# How does structural data look like ####
#------------------------------------------------------------------------------#
data2 <- dplyr::select(data,
                                       MAT,
                                       MAP,
                                       Age,
                                       meanDBH, 
                                       sdDBH, 
                                       CVDBH, 
                                       noTreeperHa, 
                                       TBAm2ha, 
                                       TotalSpeciesRichnessBA,
                                       ShannonSps,
                                       ShannonSize,
                                       scaleWeibull,
                                       shapeWeibull
)


dir.out <- 'main_analysis/strc.data'
pdf(file.path(dir.out, 'independentvariablesAnnual.pdf'), compress = T, width = 10, height = 10)
par(mfrow = c(1,1), mar = c(4,4,0.5,0.5))
chart.Correlation(data2, histogram=TRUE, pch=19, method = 'peason')
dev.off()



par(mfrow = c(5, 1), oma = c(4,4,0.5, 0.5), mar = c(0,0,0,0))
# plot Age mean with respect to their forest types
forest_type <- as.character(unique(data$IGBP2))
with(data, plot(Age, meanDBH))
for(i in 1:length(forest_type)) {
        with(data, points(Age[IGBP == forest_type[i]], meanDBH[IGBP == forest_type[i]], col = i))  
}

with(data, plot(Lat, meanDBH))
with(data, points(Lat[ Age < 20 ], meanDBH[Age < 20 ], col = 2))
with(data, points(Lat[ Age >= 20 & Age < 40 ], meanDBH[ Age >= 20 & Age < 40 ], col = 3))


with(data, plot(meanDBH, AvgAmax))
with(data, plot(Age, AvgAmax))

with(data, plot(meanDBH, AvgGPPsat))
with(data, plot(Age, AvgGPPsat))



with(data, plot(meanDBH, F_NEPtrend))
with(data, plot(meanDBH, F_GPPtrend))
with(data, plot(meanDBH, F_Rtrend))

with(data, plot(Age, F_NEPtrend))
with(data, plot(Age, F_GPPtrend))
with(data, plot(Age, F_Rtrend))

with(data, plot(TBAm2ha, F_NEPtrend))
for(i in 1:length(forest_type)) {
        with(data, points(TBAm2ha[IGBP == forest_type[i]], F_NEPtrend[IGBP == forest_type[i]], col = i))  
}
with(data, plot(TBAm2ha, F_GPPtrend))
for(i in 1:length(forest_type)) {
        with(data, points(TBAm2ha[IGBP == forest_type[i]], F_GPPtrend[IGBP == forest_type[i]], col = i))  
}
with(data, plot(TBAm2ha, F_Rtrend))
for(i in 1:length(forest_type)) {
        with(data, points(TBAm2ha[IGBP == forest_type[i]], F_Rtrend[IGBP == forest_type[i]], col = i))  
}




with(data, plot(Age, AvgGPPsat))


View(data)


########### dummy data ###################
data <- read.csv('paper_work/presentation/dummy_data.csv')
jpeg(file='paper_work/presentation/dummy_figure.jpeg', width=120,height=100, units='mm', res=300)
par(mfrow = c(1,1), mar = c(2,4,0.5,0.5))
with(data, plot(year, NEE, ylim = c(100, 600), ylab = '', pch = 16))
with(data, lines(year, NEE))
with(data, points(year[year == 2003], NEE[year == 2003], col = 'Red', pch = 16))
with(data, points(year[year == 2004], NEE[year == 2004], col = 'Blue', pch = 16))
average = mean(data$NEE[data$year != 2003])
abline(h = average, lty = 2)
mtext(expression( NEE~ "["~ g ~ C ~ m^{-2}~ yr^{-1}~ ']'), 
      side = 2, 
      line = 2.2, 
      cex = 1)
abline(h = 200, lty = 2, col = 'Red')
abline(h = data$NEE[data$year == 2004], lty = 2, col = 'Blue')
text(2003, average - 20, 'Yn (average NEE excluding drought year)')
text(2002.5, data$NEE[data$year == 2004] - 20, 'Ye+1 (NEE after drought year)', col = 'Blue')
text(2002.5, 180, 'Ye (NEE during drought year)', col = 'Red')
dev.off()
