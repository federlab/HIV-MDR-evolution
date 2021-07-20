
runSpace <- function(t, combs = c("3TC",  "AZT",  "NFV"), Rscale = 10, 
	    		write = FALSE, dirn = NULL, Ne = Ne, rseed ){

#For each potential virus in N
#multinomially sample from the 8 potential compartments
#    rseed <- sample(1:10000, 1)

    set.seed(rseed)
   
    ord <- expand.grid(A = c(1, 0), B = c(1, 0), C = c(1, 0))
    or <- ord %>% mutate(geno = paste0(A, B, C)) %>%
        mutate(geno = factor(geno, levels =
                                 c("000",
                                   "100", "010", "001",
                                   "110", "101", "011",
                                   "111" ))) %>% 
                                       arrange(geno)

    or <- (or %>% mutate(f = 0) %>% 
            mutate(f = ifelse(geno == "000", 0.01, 
                       ifelse(geno == "110", 0.03, 
                       ifelse(geno == "111", 1, f)))))

    probs <- or$f

    allvals <-  Ne * probs

    sancSize <- allvals[1]

    toTrack <- tbl_df(expand.grid(A = c(1, 0), B = c(1, 0), C = c(1, 0), comp = c("0", "A", "B", "C", "AB", "AC", "BC", "ABC")))
    initNs <- bind_cols( comp = c("0", "A", "B", "C", "AB", "AC", "BC", "ABC"), N = allvals)

    #Compute epsilon values:
    druginf <- return_druginf()

    #These are fitness costs
    mutinf <- return_mutinf()
    
    druginf <- list(druginf$`3TC`,
                    druginf$`D4T`,
                    druginf$`NFV`,
                    druginf$`LPV`)

#    names(druginf) <- c('3TC', 'D4T', 'NFV')
#    names(druginf) <- c('3TC', 'D4T', 'LPV')
    names(druginf) <- c('3TC', 'D4T', 'NFV', 'LPV')
    

    druginf[['3TC']] <- c(druginf[['3TC']], csf_ratio = ((0.05 + 1.14)/2)/((4.3 + 8.7)/2) )
    druginf[['D4T']] <- c(druginf[['D4T']], csf_ratio = ((0.2 + 0.36)/2)/((3.3 + 6.4)/2) )
    druginf[['NFV']] <- c(druginf[['NFV']], csf_ratio = ((0 + 0.012)/2)/((5.6 + 8.45)/2) )
    druginf[['LPV']] <- c(druginf[['LPV']], csf_ratio = 0.003) #((16.75))/((67945)) )

    drugEps <- foreach(drug = names(druginf), .combine = "rbind")%do%{

        R_00 <- 1

        di <- druginf[[drug]]

    #This is ugly, but whatever for now
        IC_50 <- di[which(names(di) == "IC_50")]
        hl <- di[which(names(di) == "hl")]
        m <- di[which(names(di) == "slope")]
        C_max <- di[which(names(di) == "c_max")]
        dosing <- di[which(names(di) == "dosing")]
        csf_ratio <- di[which(names(di) == "csf_ratio")]

        mi <- mutinf[[drug]]
        s <- mi[which(names(mi) == "s")]
        mu <- mi[which(names(mi) == "mu")]
        rho <- mi[which(names(mi) == "rho")]
        sigma <- mi[which(names(mi) == "sigma")]
        
        decay_int <- function(t, hl, C_max){
            -hl * C_max * 2^(-t/hl)/log(2) -  -hl * C_max * 2^(-0/hl)/log(2)
        }

        #Integrate over one dosage period
        av_concentration <- decay_int(24/dosing, hl, C_max)/(24/dosing)

        D <- lseq(10^-3, 10^5, length.out = 18)

        computeRs <- tbl_df( expand.grid( c(0, 1), c(av_concentration, av_concentration * csf_ratio, 0)))
        names(computeRs) <- c("mut", "conc")

        computeRs <- computeRs %>% mutate(comp = c("ABC", "ABC", "AB", "AB", "0", "0" ))

        computeRs %>% mutate( R = (1 - mut * s)/(1 + (conc/(IC_50 * (1 - mut*(1 - rho))))^(m * (1 + mut * (sigma)))), drug = drug)

    }

    comb_match <- tibble(combs)
    names(comb_match) <- "drug"

#Under this entry, let's just have it so that the first thing is A, the second thing is B
#etc...
    ord.drugs <- left_join(comb_match, drugEps, by = "drug")

    eps <- as.numeric(ord.drugs$R)

    A_eps <- (eps[1])
    B_eps <- (eps[2])
    C_eps <- (eps[3])


#Where did I get these computations?
    A_s<- mutinf[[combs[1]]]['s'] 
    B_s<- mutinf[[combs[2]]]['s'] 
    C_s<- mutinf[[combs[3]]]['s'] 
#I need to compute these the same way I do in the time sims

    relMuts <- c(mutinf[[combs[1]]]['mu'],  
                 mutinf[[combs[2]]]['mu'], 
                 mutinf[[combs[3]]]['mu'])

   #migration and death rates
    death <- 1

    cullNs <- initNs %>% rename(cullN = N)

    toTrack <- toTrack %>% mutate(comp = paste(comp))

    drugInfToMerge <- drugEps %>% filter(is.element(drug, combs)) %>% 
                              mutate(drug = ifelse(drug == "3TC", "A", 
                              ifelse(drug == "D4T", "B", "C"))) %>%
                                  mutate(R.A = ifelse(drug == "A", R, 1),
                              R.B = ifelse(drug == "B", R, 1), 
                              R.C = ifelse(drug == "C", R, 1)) %>% 
                                  spread( drug, mut)


   tmp_p <-  left_join(left_join(left_join(left_join(toTrack, initNs), 
              drugInfToMerge %>% select(comp, A, R.A), by = c('comp', 'A')),  
              drugInfToMerge %>% select(comp, B, R.B), by = c('comp', 'B')), 
              drugInfToMerge %>% select(comp, C, R.C), by = c('comp', 'C')) %>%
         mutate(N = ifelse((A == 0 & B == 0 & C == 0 & comp == "0"), N, 0)) %>%
         mutate(R = Rscale * R.A * R.B * R.C) %>% 
                  mutate(ind = 1:n()) %>%
                  group_by(comp) %>% 
                  mutate(numtype = cur_group_id()) %>% ungroup() %>%
                  mutate(geno = paste0(A, B, C)) %>% group_by(geno) %>%
                  mutate(numgeno = cur_group_id()) %>% ungroup() %>% 
                      mutate(R = ifelse(is.na(R), 0, R))

    plasma_inds <- which(tmp_p$comp == "ABC")

#Ok, so, tmp_p tracks each of the 8 genotypes in each of the 8 compartments. Good

    tautot <- 0

    R <- tmp_p$R
    N <- tmp_p$N

    types <- as.matrix(tmp_p %>% dplyr::select(A, B, C,  numtype, numgeno, ind))
    denom <- nrow(tmp_p)

    ordType <- unique(tmp_p$numtype)

   indList <- foreach(i = 1:8)%do%{
       which(types[,'numtype'] == i)
   }
   
    cullNs <- cullNs %>% group_by(comp) %>% 
    	   mutate(numtype = cur_group_id()) %>% ungroup() %>%
#    	   arrange(factor(numtype, levels = ordType))
        arrange(numtype)

   cn.mat <- c(as.matrix(cullNs[,'cullN']))
  
   initNs <- left_join(tmp_p %>% select(comp, numtype) %>% unique(),
                       initNs, by = "comp") %>%
             arrange(factor(numtype, levels = ordType))


    migPs <- c(initNs$N)/sum(initNs$N)
    migPs[migPs > 0] <- 1/3

    daymarker <- 0

#Our dataframe is 
    datSto <- matrix(data = rep(NA, denom * 28*12*10 * 3 + 64 * 3), ncol = 3)

    k <- 1

    newday <- 1

    while(k < t & !stopCondition_space(daymarker, N, plasma_inds, types, newday)){

        newday <- 0

        #These computations need to be done each time, because N is changing
        Rtot <- R*N
        d <- death*N
        mig <- mrate*N

        alphas <- cumsum(c(Rtot, d, mig))
        if(max(alphas) == 0){
            print(N)
        }
        tautot <- tautot + rexp(1, max(alphas))
        alphas <- alphas/max(alphas)

        indToDec <- min(which(runif(1, 0, 1) < alphas)) - 1
#eventType gives us the k
        eventType <- floor(indToDec/denom - .0000000001)
#remainder gives us the i,j
        remainder <- indToDec %% denom + 1 #(eventType > 0)
        
    #Division
        if(eventType == 0){

            toMut <- runif(3) < relMuts
            if(sum(toMut) > 0){

                toFlip <- types[remainder,]
                indToFlip <- which(toMut > 0) #sample(1:3, 1)
                toFlip[indToFlip] <- 1 - toFlip[indToFlip]
                
                addInd <- which((types[,1] == toFlip[1] &  
                                     (types[,2] == toFlip[2] &
                                      types[,3] == toFlip[3])) &
                                      types[,4] == toFlip[4])

                N[addInd] <- N[addInd] + 1
           
            }else{
                
                N[remainder] <- N[remainder] + 1

            }

        #Death
        }else if(eventType == 1){

            N[remainder] <- N[remainder] - 1
       
    #Migration
        }else if(eventType == 2){

            drawNewComp <- (remainder %% 8 + 8 * sample(0:7 + (remainder == 8), 1, prob = migPs))

            while(drawNewComp == remainder){
                 drawNewComp <- (remainder %% 8 + 8 * sample(0:7 + (remainder == 8), 1, prob = migPs))
            }


            N[remainder] <- N[remainder] -1 
            N[drawNewComp] <- N[drawNewComp] + 1

        }


        #Could we make this a faster step by testing the most common culling?
        if(N[8] > sancSize){
            toCull <- 1

        }else{

            psize <- unlist(lapply(indList, function(x){sum(N[x]) }))
            toCull <- which(cn.mat - psize  < 0)

        }
              
        if(length(toCull) > 0){
            relInds <-  indList[[toCull]]
            N[relInds] <- N[relInds] -  rmultinom(1, 1, N[relInds])

        }

        if(floor(tautot)  >= daymarker){

            daymarker <- daymarker + 1

            newday <- 1

            datSto[((denom * daymarker) + 1):(denom*(daymarker + 1)),] <- cbind(N, daymarker, tautot)

	    if(write == TRUE & daymarker %% 28 == 0 ){

	        print(paste0("Recording Day ",daymarker))
                write.table(cbind(N, k, daymarker), dirn,
                        quote = FALSE, sep = ",",
                        row.names = FALSE, col.names = FALSE, append = TRUE)
                }
        }

        k <- k + 1

    }

    colnames(datSto) <- c("N", "k", "tautot")
    datSto <- tbl_df(datSto)
    runEvs <- datSto %>% filter(!is.na(N))

}



stopCondition_space <- function(daymarker, N, plasma_inds, types, newday){

    if(newday == 1){

    #10 years
        if(daymarker >= (7*4*12*10)){
            return(TRUE)
	}

    #every 3 months
        if(daymarker %% (7*4*3) == 0 & daymarker > 0 ){

            #Note - here I have population size > 200, but I check for the more restrictive
            #150 condition that we report in the main text elsewhere

            if(sum(N[plasma_inds]) >= 200){

                return(TRUE)
            }
            return(FALSE)
        }
    }
    
    return(FALSE)

}



##############################################
################ Code to run #################
##############################################


ncores <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
registerDoParallel(ncores)

dirp <-  "../dat/space_model/"
system(paste0("mkdir ", dirp))

ind <- 1
mrate <-0.1
combs <- c("3TC",  "D4T",  "NFV")
R00 <- 10

stringToPrint <- paste(paste0(combs, collapse = "-"),
                       R00,
                       ind, sep  = "_")


foreach(i = 1:20)%dopar%{

    foreach(iter = 1:75)%do%{

      Ne <- floor(10^rnorm(1, 5.2, .5))
      while(Ne < 10^4 | Ne > 10^6){
          Ne <- floor(10^rnorm(1, 5.2, .5))
      }
      rseed <- sample(1:100000, 1)
      toWrite <-  paste0(dirp, "fast_",stringToPrint, "_m",mrate, "_Ne", Ne, "_seed",rseed,".csv")

      dat<-runSpace(500000000000, combs, Rscale = 10, write = FALSE,
              dirn = toWrite, Ne = Ne, rseed = rseed)

      write.table(dat, toWrite,
                        quote = FALSE, sep = ",",
                        row.names = FALSE, col.names = FALSE, append = TRUE)

  }
}





#This tests the case of a narrower population size distribution
dirp <-  "../dat/space_model_narrow/"
system(paste0("mkdir ", dirp))

stringToPrint <- paste(paste0(combs, collapse = "-"),
                       R00,
                       ind, sep  = "_")

foreach(i = 1:20)%dopar%{
    foreach(iter = 1:75)%do%{

      Ne <- floor(10^rnorm(1, 5.2, .05))
      rseed <- sample(1:100000, 1)

      toWrite <-  paste0(dirp, "fast_",stringToPrint, "_m",mrate, "_Ne", Ne, "_seed",rseed,".csv")

      dat<-runSpace(500000000000, combs, Rscale = 10, write = FALSE,
              dirn = toWrite, Ne = Ne, rseed = rseed)

      write.table(dat, toWrite,
                        quote = FALSE, sep = ",",
                        row.names = FALSE, col.names = FALSE, append = TRUE)

  }
}
