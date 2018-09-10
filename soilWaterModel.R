# Daily soil water model using simple bucket method # 

# Data required
# timeseries of mean air temperature, max air temperature, min air temperature, global radiation, precipitation


# Parameters required
# 
#dfafdasf


# preparing data for soil water modelling 
FileName  <- 'WAI_calc/FLX_DE-Hai_FLUXNET2015_FULLSET_DD_2000-2012_1-3.csv'
reqd.var  <- c('TIMESTAMP', 'TA_F', 'SW_IN_F', 'P_F')
head.data <- read.delim(FileName, nrows=1, sep=',', header=T)
head.name <- names(head.data)

col.var   <- which(names(head.data) %in% reqd.var)
data      <- fread(FileName, select = col.var, header = T)
names(data) <- c('date', 'Ta', 'Rg', 'Ppt')
data$date   <- ymd(data$date)

data.fapar <- fortify.zoo(mults)
data.fapar$Index <- ymd(data.fapar$Index)

data.final <- merge(data, 
                    data.fapar[,c('Index', 'fapar_f')], 
                    by.x = 'date',
                    by.y = 'Index')

model_SM <- function(data){
        
        #### Extracting site details ####
        filepath_site <- file.path (directory.in, 'site_detail.txt')            # file path for reading the file with site characteristic
        site_detail   <- read.table(filepath_site, sep='\t', header=T)          # read the file with site detail 
        site_detail   <- filter (site_detail, forest==forest_site)              # details of selected site
        lat           <- site_detail$latitude                                   # latitude of selected site
        lon           <- site_detail$longitude                                  # longitude of selected site
        ele           <- site_detail$elevation                                  # eleveation of selected site
        max_lai       <- site_detail$lai                                        # LAI of selected site
        s_d           <- site_detail$Depth                                      # soil depth of selected site
        FC            <- site_detail$FC                                         # extractable water of selected site
        WP            <- site_detail$WP                                         # wilting point of selected site
        EW            <- FC - WP
        
        #### extracting growing season details####
        filepath_gs <- file.path (directory.in, 'gsse.txt')
        gs          <- read.table(filepath_gs, sep='\t', header=T)
        
        #### constant parameters ####
        lambda        <- 2.45                                                   # Latent heat of vaporization (MJ/kg)
        pressure      <- 101.3 * ((293-(0.0065*ele))/293)^5.26                  # Atmospheric pressure (kPa)
        gamma         <- 0.000665 * pressure                                    # psychrometric constant (kPa/oC)
        
        #### soil moisture####
        filelist      <- list.files(directory.in, pattern=forest_site)          # list of files of meteo data of particular site
        year          <- as.numeric(substr (filelist, 3, 6))                    # a vector containing list of site and year
        
        
        for (i in year) {
                filepath <- file.path (directory.in, paste(forest_site, i, '_meteo.txt', sep=''))
                dataset  <- read.table(filepath, sep='\t', header=T)
                dataset  <- dataset %>%
                        mutate(delta = 4098*0.6108*exp(17.27*tair/(tair+237.3))* (tair+237.3)^-2) %>%      # slope of saturation vapor pressure curve (kPa/K)
                        mutate(pet   = ((delta*rn)+(gamma*6.43*(1+(0.536*u))*vpd))/(lambda*(delta+gamma))) # penman equation using Shuttleworth (1993) formulation (mm/day)
                #### LAI values for the year####
                bud     <- filter (gs, site==forest_site, year==i )               
                gss     <- bud$gss
                gse     <- bud$gse                
                dataset$lai <- with(dataset, ifelse(DoY < gss | DoY > gse, 0,
                                                    ifelse(DoY >= gss & DoY <= gss+30, (max_lai*DoY/30)- (max_lai*gss/30),
                                                           ifelse(DoY > gss+30 & DoY<= gse-30, max_lai, 
                                                                  (-max_lai*DoY/30)+(max_lai*gse/30)))))
                dataset$gc  <- dataset$lai/
                        
                        ####for soil water modelling####                
                DoY     <- dataset$DoY
                nt      <- length(DoY)
                p       <- dataset$precip
                lai     <- dataset$lai
                pet     <- dataset$pet
                SM      <- FC
                out     <- matrix (ncol=10,nrow=nt)    
                
                
                for (j in 1:nt){
                        REW <- (SM-WP)/EW
                        
                        # Calculation of transpiration
                        if (REW >= 0.4) {
                                if (lai[j]>=6) {
                                        tr <- 0.75 * pet[j]      
                                } else {
                                        tr <- 0.125 * lai[j] *pet [j]    
                                }
                        } else {
                                if (lai[j]>=6) {
                                        tr <- (1.5*REW + 0.15)*pet[j]                                       
                                } else {
                                        tr <- (((0.125*lai[j]-0.15)*REW/0.4) + 0.15)*pet[j]
                                }
                        }
                        
                        # calculation of interception
                        if(p[j] < 1) {
                                in_rain <- p[j]
                        } else {
                                il      <- 1- (exp(-0.5*lai[j]))               ### intercepted light based on beer lambert law, using extinction coefficient of 0.5
                                # in_rain <-  max(min(p[j] - exp (0.186 + (0.0027*il)+(0.229*p[j])+(-0.0043*p[j]*p[j])), 3), 0)                     
                                in_rain <-  max(, 1)    
                        }
                        
                        t_SM    <- SM + p[j]-tr-in_rain
                        
                        
                        if (t_SM>FC){
                                dr <- t_SM - FC
                        } else {
                                dr <- 0
                        } 
                        
                        out[j, ] <- c(i, DoY[j], p[j], lai[j], pet[j], tr, in_rain, dr, SM, REW)
                        
                        SM <- (SM + p[j] - tr - in_rain - dr)
                        
                        colnames(out) <- c("year","DoY","precip", 'lai', 'pet', 'tr', 'in', 'dr', 'sm', 'rew')
                        
                }
                write.table(out, file.path(directory.out, paste(forest_site,'_', i, '_sm_PMmodel.txt', sep='')), sep='\t', row.names=F)                           
                
        }
        
}
