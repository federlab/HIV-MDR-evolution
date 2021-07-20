
runTime <- function(combs, maxDoseT, missp, Rscale, doseType = "random", Ne = 1.5*10^5){


    mutinf <- return_mutinf()
    druginf <- return_druginf()
    
    hours <- list()
    hours[[combs[1]]] <- 0
    hours[[combs[2]]] <- 0
    hours[[combs[3]]] <- 0

    #There are three dose missing types:
    #These patterns are NOT explored in this MS, but I will include them anyway

    #1) Totally random, with dose taking probability p
    if(doseType == "random"){
        missedDoses <- runif(maxDoseT, 0, 1) < missp

	#Each year, they have a missp probability of a week long treatment interruption (i.e., 336 time periods)
	missedYears <- runif(10, 0, 1) < missp
	#Each year is 48 * 7 * 4 * 12 = 16128 time periods
	for(i in 1:10){
	      if(i == TRUE){
	      	   #choose a time within the year
		   missedDoseStart <- 16128 * (i - 1) + sample(1:16128, 1)
		   missedDoses[missedDoseStart:(missedDoseStart + 336 * 2)] <- 0
	      }
	}
	
        for(i in 2:maxDoseT){
            if(missedDoses[i] == TRUE){
                for(D in combs){           
                    hours[[D]] <- append(hours[[D]], tail(hours[[D]], n = 1) + 0.5)
                }
            }else{
                for(D in combs){
                    di <- druginf[[which(names(druginf) == D)]]
                    dosing <- di[which(names(di) == "dosing")]
                    hour <- tail(hours[[D]], n = 1)
                    if(i %% (24/dosing) == 0){
                        hour = 0
                    }else{
                        hour = hour + 0.5
                    }
                    hours[[D]] <- append(hours[[D]], hour)
                }
            }
        }

    }

    #Weekends off, weekdays on (Five on, Two off)
    if(doseType == "FOTO"){

        #Weekday index?
        #240 half hour chunks in 5 days
        #96 half hour chunks in a weekend
        missedDoses <- rep(c(rep(FALSE, 240), rep(TRUE, 96)),
                           ceiling(maxDoseT/(240+96)))[1:maxDoseT]

        for(i in 2:maxDoseT){
            if(missedDoses[i] == TRUE){
                for(D in combs){           
                    hours[[D]] <- append(hours[[D]], tail(hours[[D]], n = 1) + 0.5)
                }
            }else{
                for(D in combs){
                    di <- druginf[[which(names(druginf) == D)]]
                    dosing <- di[which(names(di) == "dosing")]
                    hour <- tail(hours[[D]], n = 1)
                    if(i %% (24/dosing) == 0){
                        hour = 0
                    }else{
                        hour = hour + 0.5
                    }
                    hours[[D]] <- append(hours[[D]], hour)
                }
            }
        }

    }

    #The dosing approach that maximizes the time at which resistant strains have higher R than
    # sensitive strains
    if(doseType == "resistantAdvantage"){

          dosingPeriods <- foreach(D = combs, .combine = "rbind")%do%{

            di <- druginf[[which(names(druginf) == D)]]
            mi <- mutinf[[which(names(druginf) == D)]]

            IC_50 <- di[which(names(di) == "IC_50")]
            hl <- di[which(names(di) == "hl")]
            m <- di[which(names(di) == "slope")]
            C_max <- di[which(names(di) == "c_max")]
            dosing <- di[which(names(di) == "dosing")]

            s <- as.numeric(mi[which(names(mi) == "s")])
            mu <- as.numeric(mi[which(names(mi) == "mu")])
            rho <- as.numeric(mi[which(names(mi) == "rho")])
            sigma <- as.numeric((mi[which(names(mi) == "sigma")]))

            mut <- 1
            checkPeriod <- 1000

            #clunky
            relHours <- bind_rows(bind_cols(t =  1:checkPeriod, nhours = seq(0.5, checkPeriod/2, by  = 0.5),
                                mut = rep(1, checkPeriod)),
                      bind_cols(t =  1:checkPeriod, nhours = seq(0.5, checkPeriod/2, by  = 0.5),
                                mut = rep(0, checkPeriod)))   %>% 
                      mutate(drug = D) %>%
                      mutate( R = (1 - mut * s)/(1 + (decay_instant(nhours, hl, C_max )/(IC_50 * (1 - mut*(1 - rho))))^(m * (1 + mut * (sigma)))))


            worstT <- (relHours %>% spread(mut, R) %>%
                 filter(`1` <= `0`) %>% slice(n = 1))$t

            missedDoses <- rep(c(rep(TRUE, worstT - 1), rep(FALSE, 1)),
                           ceiling(maxDoseT/worstT))[1:maxDoseT]


            for(i in 2:maxDoseT){

                if(missedDoses[i] == TRUE){
                    hours[[D]] <- append(hours[[D]], tail(hours[[D]], n = 1) + 0.5)

                }else{
                    hour = 0
                    hours[[D]] <- append(hours[[D]], hour)
                    
                }
            }
        }
    }


    Rinf <- foreach(D = combs, .combine = "rbind")%do%{

         di <- druginf[[which(names(druginf) == D)]]
         mi <- mutinf[[which(names(druginf) == D)]]

         IC_50 <- di[which(names(di) == "IC_50")]
         hl <- di[which(names(di) == "hl")]
         m <- di[which(names(di) == "slope")]
         C_max <- di[which(names(di) == "c_max")]
         dosing <- di[which(names(di) == "dosing")]

         s <- as.numeric(mi[which(names(mi) == "s")])
         mu <- as.numeric(mi[which(names(mi) == "mu")])
         rho <- as.numeric(mi[which(names(mi) == "rho")])
         sigma <- as.numeric((mi[which(names(mi) == "sigma")]))
       

        mut <- 1
                                        #clunky
         bind_rows(bind_cols(t =  1:maxDoseT, nhours = hours[[D]], mut = rep(1, maxDoseT)),
              bind_cols(t =  1:maxDoseT, nhours = hours[[D]], mut = rep(0, maxDoseT)))   %>% 
              mutate(drug = D) %>%
              mutate( R = (1 - mut * s)/(1 + (decay_instant(nhours, hl, C_max )/(IC_50 * (1 - mut*(1 - rho))))^(m * (1 + mut * (sigma)))))

    }

    toTrack <- tbl_df(expand.grid(A = c(1, 0), B = c(1, 0), C = c(1, 0)))
    names(toTrack) = combs
    druginds <- toTrack %>% mutate(ind = 1:n()) %>% gather(drug, mut, -ind)

   #Now, we go through and compute R at each of these timepoints for each of the resistance profiles
    allRs <- left_join(druginds, Rinf, by = c("drug", "mut")) %>%
        group_by(ind, t) %>% summarize(R = Rscale * prod(R), .groups = 'drop')

    arsplit <- allRs %>% group_split(t)

    N <- Ne
    A_i <- 3000 #3000
    lambda <- Ne * 1 * 10/(10 - 1)
    dy <- 1
    dt <- .02
    mu <- .00001

    mutinf <- foreach(D = combs, .combine = 'rbind')%do%{

        mi <- mutinf[[which(names(mutinf) == D)]]
        s <- mi[which(names(mi) == "s")]
        mu <- mi[which(names(mi) == "mu")]
        rho <- mi[which(names(mi) == "rho")]
        sigma <- mi[which(names(mi) == "sigma")]

        return(c(D, mu, s, sigma, rho))
    }

    colnames(mutinf) <- c("drug", "mu", "s", "sigma", "rho")
    mutinf <- tbl_df(mutinf)

    init.n <- left_join(druginds, mutinf, by = "drug") %>%
        mutate(mu_over_s = as.numeric(mu)/as.numeric(s) ) %>%
            group_by(ind) %>% mutate(mu_over_s = ifelse(mut == 1, mu_over_s, 1)) %>%
                summarize(lam = prod(mu_over_s), .groups = "drop")  %>%
                    mutate(i.N = round(lam*N))

    tmp <- init.n %>% filter(ind < 8)
    init.n <- init.n %>% mutate(i.N = ifelse(ind == 8, N - sum(tmp$i.N), i.N)) %>%
        mutate(lam = ifelse(ind == 8, 1 - sum(tmp$lam), lam))

                                        #new cells
    ys <- matrix(0, nrow = 8, ncol = maxDoseT + 1)
    ys[,1] <- init.n$i.N

    toMerge <- toTrack %>% mutate(ind = 1:n()) 
    names(toMerge) <- c("m1", "m2", "m3", "ind")

    ref <- bind_rows(
        bind_cols(ind = 0, name = "m1"), 
        bind_cols(ind = 0, name = "m2"),
        bind_cols(ind = 0, name = "m3"))


    start_time <- Sys.time()
    i <- 1

    res <- init.n$lam
    while(i < maxDoseT & ! stopCondition(i, ys[,i]) ){ 

        emergingCells <- rpois(8, A_i * dt * res)

        rs <- arsplit[[i]]$R

        yvals <- ys[,i]
        
        denom <- lambda + sum(( dy * rs * yvals))
        rates <- (yvals*dy*lambda*rs*dt)/denom

        newys <- rpois(8, rates)
        
        mutFrom.s1 <- rbinom(8, newys, as.numeric(mutinf$mu[1]))
        mutFrom.s2 <- rbinom(8, newys, as.numeric(mutinf$mu[2]))
        mutFrom.s3 <- rbinom(8, newys, as.numeric(mutinf$mu[3]))

        mutFrom <- rep(0, 8)
        mutTo <- rep(0, 8)

        for(ind in 1:8){

            if(newys[ind] > 0){
                m1 <- sample(1:newys[ind], mutFrom.s1[ind], replace = FALSE)
                m2 <- sample(1:newys[ind], mutFrom.s2[ind], replace = FALSE)
                m3 <- sample(1:newys[ind], mutFrom.s3[ind], replace = FALSE)

                base <- toTrack[ind,]
                
                mutsTo <- bind_rows(bind_rows(
                    bind_cols(ind = m1, name = rep("m1", length(m1))), 
                    bind_cols(ind = m2, name = rep("m2", length(m2))),
                    bind_cols(ind = m3, name = rep("m3", length(m3)))) %>% 
                        mutate(i = 1), ref %>% mutate(i =0)) %>%
                            spread(name, i) %>% group_by(m1, m2, m3) %>% 
                                summarize(n = n(), .groups = "drop")  %>% 
                                    filter(!(m1 == 0 & m2 == 0 & m3 == 0)) %>%
                                        mutate(m1 = ifelse(is.na(m1), 0, m1), 
                                               m2 = ifelse(is.na(m2), 0, m2),
                                               m3 = ifelse(is.na(m3), 0, m3))
               
                #So, we actually want something different here - we want to mutate from
                #based on the index we are considering
                #Here, ind = 5, or 110. Our mutation is the same (second position), 
                #so we should go to 100, instead of 000 -> 010
                
                myInd <- ind
                flips <- toMerge %>% filter(ind == myInd)

                mutsTo <- mutsTo %>% mutate(m1 =  abs(flips$m1 - m1),
                                            m2 =  abs(flips$m2 - m2),
                                            m3 =  abs(flips$m3 - m3))
                
                mutFrom[ind] <- sum(mutsTo$n)
                mutTo <- mutTo + (left_join(toMerge, mutsTo, by = c("m1", "m2", "m3")) %>% 
                                      mutate(n = ifelse(is.na(n), 0, n)))$n
                   
            }

        }
        
        ys[,i + 1] <- rbinom(8, emergingCells + newys + yvals - mutFrom + mutTo, exp(-dy * dt))

        i <- i + 1

    }

    end_time <- Sys.time()
    print(end_time - start_time)
    
    conc_by_drug <- Rinf

    R_by_geno <- left_join(druginds, Rinf, by = c("drug", "mut")) %>%
        group_by(ind, t) %>% summarize(R = Rscale * prod(R), .groups = 'drop')

    return(list(f = ys, drugs = conc_by_drug, R = R_by_geno))
}

stopCondition <- function(t, ysi){
    
    #every 6 months
    if(t %% (2*24*7*4*3) == 0  ){

        if(sum(ysi[1:8]) >= 5000){
            return(TRUE)
        }
        return(FALSE)
    }
    return(FALSE)

}


##############################################
################ Code to run #################
##############################################


#Relates to parallelization - feel free to comment
ncores <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
registerDoParallel(ncores)

toTrack <- tbl_df(expand.grid(A = c(1, 0), B = c(1, 0), C = c(1, 0)))
druginds <- toTrack %>% mutate(ind = 1:n()) %>% gather(drug, mut, -ind)
indinf <- toTrack %>% mutate(comb = paste0(A, B, C)) %>%
    mutate(ind = 1:8)  %>% select(comb, ind)

#2 periods/hour, 24 hours/day, 7 days/week, 4 weeks/month, 12 months/year * 10 years
maxDoseT <- 2 * 24 * 7 * 4 * 12 * 10
ind <- 1

dirp <- "../dat/time_model/"
print(paste0("making directory: ", dirp))

system(paste0("mkdir ", dirp))

#This refers to the random nature at which doses are missed
treatmentType <- "random"
combs <- c("3TC", "D4T", "NFV")

#These nested loops draw from our cluster set up
#As a test to see if things are working properly, I'd recommend lowering
#these numbers
foreach(ind = 1:20)%dopar%{
    foreach(ind = 1:75 )%do%{

        R00 <- 10   

        Ne <- floor(10^rnorm(1, 5.2, .5))
        while(Ne < 10^4 | Ne > 10^6){
            Ne <- floor(10^rnorm(1, 5.2, .5))
        }
        rseed <- sample(1:100000, 1)

        #Adherence
        pv <- runif(1, 0, 1)

        stringToPrint <- paste(paste0(combs, collapse = "-"),
                               treatmentType,
                               pv,
                               R00,
                               ind, Ne, rseed, sep  = "_")

        dat <- runTime(combs = combs, maxDoseT = maxDoseT, missp = pv,
                      Rscale = R00, doseType = treatmentType, Ne = Ne)

        relInds <- c(1, which(1:ncol(dat$f) %% (2 * 24 * 7 * 4 ) == 0))

        write.table(dat$f[,relInds], paste0(dirp, stringToPrint, "_f.csv"),
                    quote = FALSE, sep = ",",
                    row.names = FALSE, col.names = FALSE)

    }
}


		      







