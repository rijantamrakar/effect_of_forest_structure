#----------------------------------------------------------------------#
# Reading eddypro datafile
data.file  <- 'flux_data/Unzipped_Fluxnet2015_database_tier_2/DE-Hai/FLX_DE-Hai_FLUXNET2015_FULLSET_HH_2000-2012_1-3.csv'
head.data  <- read.delim(data.file, sep = ',', header = F, nrow = 2)  
head.data  <- unname(unlist(head.data[1,]))  

reqd.var   <- c('TIMESTAMP_END', 
                'TA_F', 
                'SW_IN_F',
                'VPD_F', 
                'P_F', 
                'RH', 
                'PPFD_IN', 
                'PPFD_DIF', 
                'TS_F_MDS_1', 
                'TS_F_MDS_2',
                'NIGHT',
                'NEE_VUT_REF',
                'NEE_VUT_REF_QC')

col.var    <- which(head.data %in% reqd.var)

data  <- data.table::fread(data.file, select = col.var, skip = 1, header = F)
names(data) <- reqd.var
data[data == -9999] <- NA
data$NEP <- -data$NEE_VUT_REF
# parse the date time
data$timestamp <- ymd_hm(data$TIMESTAMP_END)
data$doy       <- yday(data$timestamp)
data$year      <- year(data$timestamp)
data <- filter(data, doy > 100, doy < 301)

# start working year by year and week by week
Years <- unique(data$year)

Year = 2002
DataYear <- dplyr::filter(data, year== Year)

#------------------------------------------------------------------------------#
Periods = seq(120, 280, 1)
# night time respiration
REcoLloydTaylor <- function(Rb, Eo, Tref, Tair, To) {
        (Rb * exp(Eo * (1/(Tref - To) - 1/(Tair - To))))
}

AmaxVPD <- function(pmax, k, VPD) {
        pmax * exp(-k * (VPD - 10))
}

CostFunction <- function(Yobs, Ymod) {
        DiffbetweenModandMeas <- Yobs - Ymod
        squaredSumofDiff <- sum(DiffbetweenModandMeas^2)
        VarianceMeas <- (sd(Yobs))
        cost <- round(squaredSumofDiff/VarianceMeas, 4)
        return(cost)
}


# constants
To   <- 227.13 
Tref <- 273.15 + 10

dir.out <- 'main_analysis/anomaly_detection/half-hourly/NightTimeRespiration'
pdf(file.path(dir.out, 'DE-Hai_2002.pdf'), compress = T)
par(mfrow = c(2,2), mar = c(4,4,0.5,0.5))
Eo.all <- NULL
for(i in 1:length(Periods)) {
        Data.Night <- dplyr::filter(DataYear, SW_IN_F <= 20, doy >= Periods[i]-7, doy <= Periods[i] + 7 ) 
        # get original data only
        Data.Night <- dplyr::filter(Data.Night, NEE_VUT_REF_QC == 0)
        
        print(paste('Day', Periods[i]))
        if(nrow(Data.Night) > 10) {
                # Reco
                REco <- Data.Night$NEE_VUT_REF
                Tair <- 273.15 + Data.Night$TA_F
                meanNEE <- mean(Data.Night$NEE_VUT_REF)
                InitialEoList  <- c(100, 50, 200)
                InitialRbList  <- c(meanNEE, meanNEE/2, 2*meanNEE)
                Eo.Res.i.Parm <- NULL
                for (j in 1:3){
                        for(k in 1:3) {
                                print(paste('Day', Periods[i], 'PARM', j, k))
                                MReco   <- try(nls(REco ~ REcoLloydTaylor (Rb, Eo, Tref, Tair, To), 
                                                   start = c(Rb = InitialRbList[j], 
                                                             Eo = InitialEoList[k]))
                                )
                                if(class(MReco)!="try-error") {
                                        summaryRecoModel = summary(MReco)
                                        Eodiff = summaryRecoModel$coefficients['Eo', 'Estimate']
                                        if(Eodiff >= 50 & Eodiff <= 400) {
                                                Yobs  <- REco 
                                                Ymod  <- predict(MReco)
                                                lmodObsMod <- summary(lm(Yobs~Ymod))
                                                Modelcost <- CostFunction(Yobs, Ymod)  
                                                
                                                Eo.Res.i.Parm <- rbind(Eo.Res.i.Parm, c('Eo' = Eodiff, 
                                                                                        'Cost' = Modelcost,
                                                                                        'R2' = round(lmodObsMod$r.squared, 2),
                                                                                        'n' = nrow(Data.Night) 
                                                                                        )
                                                )
                                                
                                                if(j == 1 & k == 1) {
                                                        plot(REco~Tair, 
                                                             ylab = 'Reco', 
                                                             xlab = 'Tair [deg K]')
                                                        points(predict(MReco)~Tair, col = 2)
                                                        leg <- paste('doy =',
                                                                     Periods[i],
                                                                     ', Rsq =',
                                                                     round(lmodObsMod$r.squared, 2),
                                                                     ', n =',
                                                                     nrow(Data.Night)
                                                        )
                                                        legend('topleft', legend= leg, bty = 'n')
                                                }
                                        } # end of if the Eo 
                                } # end of try error
                        } # end of k
                } # end of j
                
                if(length(Eo.Res.i.Parm) > 0) {
                        Eo.Res.i.Parm <- as.data.frame(Eo.Res.i.Parm)
                        Eo.Res.i.Parm <- Eo.Res.i.Parm[!duplicated(Eo.Res.i.Parm$Cost), ]
                        minCost       <- min(Eo.Res.i.Parm$Cost)
                        Eo.final      <- Eo.Res.i.Parm$Eo[Eo.Res.i.Parm$Cost == minCost]
                        R2            <- Eo.Res.i.Parm$R2[Eo.Res.i.Parm$Cost == minCost]
                        n             <- Eo.Res.i.Parm$n[Eo.Res.i.Parm$Cost == minCost]
                        error         <- 0
                        
                        if(Eo.final < 50) {
                                Eo.final <-  Eo.all[i-1, 'Eo']   
                                R2       <- 0
                                n        <- 0
                                error    <- 1
                        }
                        
                        if(Eo.final > 400) {
                                Eo.final <-  Eo.all[i-1, 'Eo']   
                                R2       <- 0
                                n        <- 0
                                error    <- 2
                        }
                        
                } else {
                        Eo.final <-  Eo.all[i-1, 'Eo']   
                        R2       <- 0
                        n        <- 0
                        error    <- 4
                }
                
                Eo.all <- rbind(Eo.all, c('doy'   = Periods[i], 
                                          'Eo'    = Eo.final, 
                                          'R2'    = R2, 
                                          'n'     = n,
                                          'error' = error
                                          )
                                )
                
        } else {                                                                # if no of data is greater than 10
                Eo.final <-  Eo.all[i-1, 'Eo'] 
                R2       <- 0
                n        <- 0
                error    <- 5
                
                Eo.all <- rbind(Eo.all, c('doy'   = Periods[i], 
                                          'Eo'    = Eo.final, 
                                          'R2'    = R2, 
                                          'n'     = n,
                                          'error' = error
                )
                )
        }  
}
dev.off()

Eo.all <- as.data.frame(Eo.all)

dir.out <- 'main_analysis/anomaly_detection/half-hourly/NightTimeRespiration'
pdf(file.path(dir.out, 'DE-Hai_Eo.pdf'), compress = T)
par(mfrow = c(1,1))
par(mar = c(0, 4, 0, 4), oma = c(2, 0, 0.5, 0.5)) 
with(Eo.all, plot(doy, Eo))      
# Plot R2
par(new = T)
with(Eo.all, plot(doy, R2, axes=F, xlab=NA, ylab=NA, cex=0.85, typ = 'l', ylim = c(0,0.30), xlim = c(120,280)))
axis(side = 4)
mtext(side = 4, line = 2.7, expression(R^2 ~ '[ ]'), cex = 0.85)
dev.off()


Eo.all.Filtered <- filter(Eo.all, R2 >0.1)
Eo <- Eo.mean <- mean(Eo.all.Filtered$Eo)


dir.out <- 'main_analysis/anomaly_detection/half-hourly/NightTimeRespiration'
pdf(file.path(dir.out, 'DE-Hai_2002_LRC.pdf'), compress = T)
par(mfrow = c(2,2), mar = c(4,4,0.5,0.5))
res.all <- NULL
for(i in 1:length(Periods)) {
        Data.Day <- dplyr::filter(DataYear, NIGHT == 0, doy >= Periods[i]-3, doy <= Periods[i] + 3 ) 
        Data.Day <- dplyr::filter(Data.Day, NEE_VUT_REF_QC == 0)
        if(nrow(Data.Day) > 20) {
                VPD     <- Data.Day$VPD_F
                SW_IN_F <- Data.Day$SW_IN_F
                NEP     <- -Data.Day$NEE_VUT_REF
                Tair    <- Data.Day$TA_F + 273.15
                
                InitialEoList  <- c(100, 50, 200)
                InitialRbList  <- c(meanNEE, meanNEE/2, 2*meanNEE)
                
                cc<-try(nls(NEP~((alpha * AmaxVPD (pmax, k, VPD) * SW_IN_F)/(alpha * SW_IN_F + AmaxVPD (pmax, k, VPD)) -  REcoLloydTaylor (Rb, Eo, Tref, Tair, To)),
                            start = list("pmax" = 20,"alpha"= 0.04, 'k' = 0.03, 'Rb' = 5),
                            lower = list("pmax" = 0, "alpha"= 0,    'k' = 0, 'Rb' = 0),
                            upper = list("pmax" = 200, "alpha"= 0.22,    'k' = 100, 'Rb' = 20),
                            algorithm = 'port',
                            trace=F
                ))
                
                if(class(cc)!="try-error") {
                        lmodObsMod = lmodel2(predict(cc)~NEP)
                        rsqModel <- round(lmodObsMod$rsquare, 2)
                        
                        summaryGPPmodel = summary(cc)
                        pmax1 = summaryGPPmodel$coefficients['pmax', 'Estimate']
                        alpha1 = summaryGPPmodel$coefficients['alpha', 'Estimate']
                        k1 = summaryGPPmodel$coefficients['k', 'Estimate']
                        Ro1 = summaryGPPmodel$coefficients['Rb', 'Estimate']
                        
                        res.all  <-rbind(res.all, c("Day"   = Periods[i], 
                                                    "n"     = length(resid(cc)),
                                                    "error" = 0,
                                                    "r2 "   = rsqModel,
                                                    'pmax'  = pmax1,
                                                    'alpha' = alpha1,
                                                    'k'     = k1,
                                                    'Rb'    = Ro1,
                                                    'Eo'    = Eo
                        )
                        )
                        
                        plot(NEP~SW_IN_F, xlim = c(0, 2000), ylim = c(-10,50))
                        points(predict(cc)~SW_IN_F, col = '2')
                        leg <- paste('doy =',
                                     Periods[i],
                                     ', Rsq =',
                                     rsqModel,
                                     ', n =',
                                     nrow(Data.Day)
                        )
                        legend('topleft', legend= leg, bty = 'n')
                } else {                                                 # end of if if the class of cc is not "Error"
                        res.all  <-rbind(res.all, c("Day"   = Periods[i], 
                                                    "n"     = 0,
                                                    "error" = 1,
                                                    "r2 "   = rsqModel,
                                                    'pmax'  = pmax1,
                                                    'alpha' = alpha1,
                                                    'k'     = k1,
                                                    'Rb'    = Ro1,
                                                    'Eo'    = Eo
                        )
                        )
                        
                        
                }
        } else {
                res.all  <-rbind(res.all, c("Day"   = Periods[i], 
                                            "n"     = 0,
                                            "error" = 2,
                                            "r2 "   = rsqModel,
                                            'pmax'  = pmax1,
                                            'alpha' = alpha1,
                                            'k'     = k1,
                                            'Rb'    = Ro1,
                                            'Eo'    = Eo
                )
                )  
        }     
}

dev.off()        
        





# get daytime data and
data <- dplyr::filter(data, NIGHT == 0)
# get original data only
data <- dplyr::filter(data, NEE_VUT_REF_QC == 0)







for(i in 1:length(k)) {
        
        datashort <- filter(data, doy >= k[i] & doy <= k[i] + 6)
        
        cc<-try(nlrq(NEP~((alpha * amax * SW_IN_F)/(alpha * SW_IN_F + amax) + Rd),
                    start = list("amax"=15,"alpha"=0.04,"Rd"=-1.5),
                    tau = 0.90,
                    trace=F, 
                    method="L-BFGS-B",
                    data = datashort))
        
        cc<-try(nlrq(NEP~((alpha * (pmax * exp(-k * (VPD_F - 10))) * SW_IN_F)/(alpha * SW_IN_F + (pmax * exp(-k * (VPD_F - 10)))) - Rd),
                     start = list("pmax" = 21,"alpha"= 0.04,"Rd"=3.9, 'k' =0.03),
                     tau = 0.50,
                     trace=F, 
                     method="L-BFGS-B",
                     data = datashort))
        
        with(datashort, plot(NEP~SW_IN_F, ylim = c(-20, 40)))
        points(datashort[,"SW_IN_F"],predict(cc),col=8)
        
        datashort$ResidSW = resid(cc)
        
        with(datashort, plot(ResidSW~TA_F))
        summary(with(datashort, lm(ResidSW~TA_F, data = datashort)))
        
        with(datashort, plot(ResidSW~VPD_F))
        with(datashort, plot(ResidSW~SW_IN_F))
        
}

plot(NEP~SW_IN_F, ylim = c(-20, 20), xlim = c(0, 2000))
points(predict(cc)~SW_IN_F, col = 'Red')

plot(lm, ylim = c(-40, 40), xlim = c(-40, 40), 
     ylab = 'predicted NEE',
     xlab = 'observed NEE')
abline(a = 0, b = 1, lwd = 2)