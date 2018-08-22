# SPEI ##

data(wichita)
View(wichita)
names(wichita)
summary(wichita)
wichita$PET <- thornthwaite(wichita$TMED, 37.6475)
wichita$BAL <- wichita$PRCP-wichita$PET

# Convert to a ts (time series) object for convenience
wichita <- ts(wichita[,-c(1,2)], end=c(2011,10), frequency=12)
plot(wichita)

spei1 <- spei(wichita[,'BAL'], 1)
summary(spei1)
names(spei1)
spei1$call
spei1$fitted
spei1$coefficients
# Plot spei object
par(mfrow=c(2,1))
plot(spei1, main='Wichita, SPEI-1')



HainichData <- 'C:/Users/Rijan/Documents/phd_work/data_analysis/analysed_data/correlation_diff_timescales/month/De.Ha.month.csv'
HainichData <- read.csv(HainichData)
HainichData <- dplyr::select(HainichData, year, month, Tair, rain_mm)
HainichData <- mutate(HainichData, pet = thornthwaite(Tair, 51.0794084),
                      bal = rain_mm - pet)
HainichData <- ts(HainichData[,-c(1,2)], start=c(2000,1), frequency=12)
plot(HainichData)

spei1 <- spei(HainichData[,'bal'], 1)
summary(spei1)
names(spei1)
spei1$call
spei1$fitted
spei1$coefficients
# Plot spei object
par(mfrow=c(2,1))
plot(spei1, main='Wichita, SPEI-1')

data.spei <- ts_df(spei1$fitted)
