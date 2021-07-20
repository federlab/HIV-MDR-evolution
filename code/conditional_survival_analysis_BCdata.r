
treatmentTypes <- c("3TC", "NRTI", "NNRTI", "PI", "total")

cond_surv <- foreach(treatmentType = treatmentTypes, .combine = "rbind")%do%{

    print(paste0("Analyzing ", treatmentType))

    #read in the data - 
    #NOTE: this data must be acquired via the BC Centre for Excellence in HIV/AIDS
    km <- read.table(paste0("BC_KM_", treatmentType, ".csv"), sep = ",", header = TRUE)

    #I think we should just be selecting the columns STRATUM, follow_X_res, and Censor here
    km <- km[c(1, 3, 4)]

    #standardize column names
    colnames(km) <- c('STRATUM', "follow_res", "Censor")

    #For this survival analysis, we need our censor events to be encoded by 0s
    km <- cbind(km, Event = 1 - km$Censor)

    #We only want to keep individuals from strata 1 and 2
    km.toanalyze <- km[km$STRATUM == 1 | km$STRATUM == 2, ]

    #For each year of follow up
    conditional_survivals <- foreach(k = 0:9, .combine = "rbind")%do%{

                          #This should fit the conditional survival probability for year k + 1
                          #conditional only on individuals with no events up to year k
                          #the time1 = follow_res, event1 = 0 means that there are no additional 
                          #events we're tracking beyond the initial endpoint
        cond_surv_year_k <- summary(survCOND(survCS(
            time1 = follow_res, event1 = Event, 
            Stime = follow_res, event = Event) ~ 1, 
            x = k, y = k+1,
            data = km.toanalyze, method = "KMW"))
                          
        #Let's also keep track of the sample sizes at the start of each year
        n <- nrow(km.toanalyze[km.toanalyze$follow_res >= k, ])

        return(c(as.numeric(cond_surv_year_k), n)) 

    }

    colnames(conditional_survivals) <- 
        c("Year", "point.estimate", "95_conf_low", "95_conf_high", "n")

    km.toanalyze.with.dummy <- km.toanalyze
    km.toanalyze.with.dummy$dummy = 0
    surv_object <- Surv(time = km.toanalyze.with.dummy$follow_res, 
                        event = km.toanalyze.with.dummy$Event)
    fit1 <- survfit(surv_object ~ dummy, data = km.toanalyze.with.dummy)
    km.est <- (((summary(fit1, times = 0:10)))[[6]])[-1]
    
    tbl_df(conditional_survivals) %>% mutate(km = km.est, 
                                             treatmentType = treatmentType)# %>% 
#    mutate(year_kp1_over_year_k = ifelse(is.na(lag(km)), km, km/lag(km))) %>% 
#    select(Year, point.estimate, year_kp1_over_year_k)
#   This commented code can also be used to compute a rough estimate of the conditional
#   survival probability from the km estimates.     

}

write.table(cond_surv, "../processed/conditional_survival_analysis.csv", 
            sep = ",", col.names = TRUE, 
            row.names = FALSE, quote = FALSE)


