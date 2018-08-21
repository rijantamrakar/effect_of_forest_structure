### preparing a list of sites and completing data ###

#### Structural data ####
DataStru   <- read.csv ('main_analysis/strc.data/structural.index.one.value.csv')
#### GPP sat data ####
DataGPPsat <- read.csv('main_analysis/GPPsat/Annual_extractions_GPP1000_Amax.csv')
DataGPPsatAvg  <- DataGPPsat %>%
        group_by(Site_code) %>%
        summarise(AvgAmax   = mean(Amax, na.rm = T),
                  SdAmax    =  sd(Amax, na.rm = T),
                  AvgGPPsat = mean(GPP1000, na.rm = T),
                  SdGPPsat  = sd(GPP1000, na.rm = T),
                  NoYears   = n()
        ) %>%
        mutate(CvAmax  = SdAmax/AvgAmax,
               CvGPPsat = SdGPPsat/AvgGPPsat
        )

# Combined GPPsat and structural data #
DataCombined <- merge(DataStru, DataGPPsatAvg, by.x = 'site', by.y = 'Site_code', all = T)

# other site details
DataSiteDetails <- read.csv('main_analysis/finallist_20180626.csv')
DataSiteDetails <- dplyr::select(DataSiteDetails, 
                                 -Name, 
                                 -Country, 
                                 -PI, 
                                 -email,
                                 -X20180621)

# combine
DataCombined <- merge(DataCombined, DataSiteDetails, by = 'site', all.y = T)
write.csv(DataCombined, 'main_analysis/finallist_Short_20180712ver2.csv', row.names = F)
