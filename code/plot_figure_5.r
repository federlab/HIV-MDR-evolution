#Plot data from spatial and temporal model runs
red <- rgb(222, 39, 47, 255, maxColorValue = 255)
yellow <-  rgb(254, 229, 157, 255, maxColorValue = 255)
lblue <- rgb(168, 206, 225, 255, maxColorValue = 255)
dblue <- rgb(43, 122, 176, 255, maxColorValue = 255)                
cols <- c(red, yellow, lblue, dblue, dblue, dblue, "black")
names(cols) <- c("NNRTI", "PI", "NRTI", "M184VI", "3TC", "3TC/FTC", "total")

numMutCols <- c(rgb(239, 239, 239, 255, maxColorValue = 255), 
                rgb(189, 189, 189, 255, maxColorValue = 255),
                rgb(129, 129, 129, 255, maxColorValue = 255),
                rgb(0, 0, 0, 255, maxColorValue = 255)
                )
names(numMutCols) <- c("0", "1", "2", "3")

indinf <- tbl_df(expand.grid(1:0, 1:0, 1:0)) %>% 
        mutate(comb = paste0(Var1, Var2, Var3)) %>% 
        select(comb) %>% mutate(ind = 1:8)


toTrack <- tbl_df(expand.grid(A = c(1, 0), B = c(1, 0), C = c(1, 0), comp = c("0", "A", "B", "C", "AB", "AC", "BC", "ABC")))
initNs <- bind_cols( comp = c("0", "A", "B", "C", "AB", "AC", "BC", "ABC"), N = 0)
tmp_p <- left_join(toTrack, initNs) %>% mutate(geno = paste0(A, B, C))

fancy_scientific <- function(l) {
#Modified from: https://groups.google.com/forum/#!topic/ggplot2/a_xhMoQyxZ4
     # turn in to character string in scientific notation
     l <- format(l, scientific = TRUE)
     # quote the part before the exponent to keep all the digits
     l <- gsub("^(.*)e", "'\\1'e", l)
     # turn the 'e+' into plotmath format
     l <- gsub("e", "%*%10^", l)
     # return this as an expression
     l <- gsub("0e\\+00","0",l)
     l <- gsub("\\+", "", l)
     parse(text=l)
 }


readInTime <- function(infiles){

    id.time  <- foreach(fi = infiles, .combine = "rbind")%do%{

        indat <- tbl_df(read.table(paste0(dirp, fi), sep = ","))

        toP <- indat %>% mutate(ind= 1:8) %>% gather(time, N, -ind) %>% 
            mutate(time = as.numeric(gsub("V", "", time))) %>% filter(time > 1)
        toP <- left_join(indinf,toP, by = "ind") %>% arrange(time)

        fi_split <- strsplit(fi, split = "_")[[1]]
        adh <- as.numeric(fi_split[3])
        R00 <- as.numeric(fi_split[4])
        Ne <-  as.numeric(fi_split[6])

#Let's find the first time that the population size goes above the detectable threshold
        firstTimeDetect <- toP %>% group_by(time) %>% summarize(totalN = sum(N)) %>% 
            filter(totalN > detect, time %% 3 == 0) %>% slice(n = 1)

        if(nrow(firstTimeDetect) > 0){

#We now want to find which mutations are detectable at 15% frequency
            callFreqs <- toP %>% filter(time == firstTimeDetect$time) %>% 
                separate(comb, into = c("_", "l1", "l2", "l3"), sep = "") %>% 
                    select(-'_', -ind, -time)  %>% mutate(f = N/sum(N))

            trialDat <- bind_cols(callFreqs %>% filter(l1 == 1) %>% summarize(f1 = sum(f)), 
                                  callFreqs %>% filter(l2 == 1) %>% summarize(f2 = sum(f)),
                                  callFreqs %>% filter(l3 == 1) %>% summarize(f3 = sum(f))) %>% 
                        mutate(numDetectableMuts = (f1 > 0.15) +(f2 > 0.15) + (f3 > 0.15), 
                                             time = firstTimeDetect$time, VF = TRUE) 
            
        }else{

            trialDat <- t(as.matrix(c(NA, NA, NA, NA, max(toP$time), FALSE)))
            colnames(trialDat) <- c("f1", "f2", "f3", "numDetectableMuts", "time", "VF")
            trialDat <- tbl_df(trialDat)
        }

        trialDat %>% mutate(adh = adh, R00 = R00, Ne = Ne)  %>% 
            mutate(time = time/12)

    }

    return(id.time)
}



readInSpace <- function(infiles){

    id.space <- foreach(i = 1:length(infiles), .combine = "rbind")%do%{

        fi <- infiles[i]
        indat <- tbl_df(read.table(paste0(dirp, fi), sep = ",", skip = 0))
        names(indat) <- c('N', 'td', 't')
        genoinf <- tmp_p %>% select(comp, geno)
        id <- indat %>% group_by(t) %>% 
            mutate(ind = 1:n(), comp = floor((ind - 1)/8), ind = ind%%8) %>% 
                select(-td) %>% ungroup()
        id <- id %>% group_by(t) %>% 
            mutate(comp = genoinf$comp, geno= genoinf$geno) %>% ungroup() %>% 
                mutate(t = floor(t) + 1)

        firstTimeDetect <- id %>% filter(comp == "ABC") %>% 
            group_by(t) %>% summarize(totalN = sum(N)) %>% 
            filter(totalN > detect, t %% (28 * 3) == 0) %>% slice(n = 1)

        fi_split <- strsplit(fi, split = "_")[[1]]
        m <- as.numeric(gsub("m", "", fi_split[5]))
        R00 <- as.numeric(fi_split[3])
        Ne <-  as.numeric(gsub("Ne", "", fi_split[6]))

        if(nrow(firstTimeDetect) > 0){

        #We now want to find which mutations are detectable at 15% frequency
            callFreqs <- id  %>% filter(comp == "ABC") %>% 
                filter(t == firstTimeDetect$t) %>% 
                separate(geno, into = c("_", "l1", "l2", "l3"), sep = "") %>% 
                select(-'_', -ind, -t)  %>% mutate(f = N/sum(N))
        
            trialDat <- bind_cols(callFreqs %>% filter(l1 == 1) %>% summarize(f1 = sum(f)), 
                  callFreqs %>% filter(l2 == 1) %>% summarize(f2 = sum(f)),
                  callFreqs %>% filter(l3 == 1) %>% summarize(f3 = sum(f))) %>% 
                  mutate(numDetectableMuts = (f1 > 0.15) +(f2 > 0.15) + (f3 > 0.15), 
                         t = firstTimeDetect$t, VF = TRUE) 

        }else{

            trialDat <- t(as.matrix(c(NA, NA, NA, NA, max(id$t), FALSE)))
            colnames(trialDat) <- c("f1", "f2", "f3", "numDetectableMuts", "t", "VF")
            trialDat <- tbl_df(trialDat)
        }

    #Let's put everything in terms of 48 week years
        trialDat %>% mutate(m = m, R00 = R00, Ne = Ne) %>% mutate(t = t / (4 * 7 * 12) )

    }

    return(id.space)

}



plot_sims <- function(id, justConds = FALSE){

    condSurvs <- foreach(i = 1:3, .combine = "rbind")%do%{

        loc <- id[, c(i, 5, 9)]
        names(loc) <- c('f', 't', 'Ne')
        loc <- loc %>% mutate(event = ifelse(f > 0.15, 1, 0)) %>% #select(t, event) %>% 
            mutate(event = ifelse(is.na(event), 0, event))

        conditional_survivals <- foreach(k = 0:9, .combine = "rbind")%do%{

                          #This should fit the conditional survival probability for year k + 1
                          #conditional only on individuals with no events up to year k
                          #the time1 = follow_res, event1 = 0 means that there are no additional 
                          #events we're tracking beyond the initial endpoint
            
                          cond_surv_year_k <- summary(survCOND(survCS(
                                                  time1 = t, event1 = event, 
                                                  Stime = t, event = event) ~ 1, 
                                           x = k, y = k+1,
                                           data = loc, method = "KMW"))

                          return(c(as.numeric(cond_surv_year_k))) 

                      }

        colnames(conditional_survivals) <- 
            c("Year", "point.estimate", "X95_conf_low", "X95_conf_high")
    
        conditional_survivals <- tbl_df(conditional_survivals) %>% mutate(treatmentType = i)

    }

    p1 <- condSurvs %>% mutate(point.estimate = 1 - point.estimate, 
                        X95_conf_low = 1 - X95_conf_low, 
                        X95_conf_high = 1 - X95_conf_high) %>% 
          filter(treatmentType != "total") %>%
          mutate(treatmentType = ifelse(treatmentType == 1, "3TC", 
                                 ifelse(treatmentType == 2, "NRTI",
                                 ifelse(treatmentType == 3, "PI", NA)))) %>%     
          ggplot() + 
              geom_point(aes(x =  Year, y = point.estimate, col = treatmentType)) + 
              geom_line(aes(x = Year , y = point.estimate, col = treatmentType, 
                            group = treatmentType)) + 
              geom_segment(aes(x = Year , xend = Year, y = X95_conf_low, 
                               yend = X95_conf_high, col = treatmentType, 
                               group = treatmentType))  +
              scale_color_manual(values = cols) + labs(x = "Year of treatment",
                                                        y = "Probability of resistance") + 
              theme_classic() + theme(legend.position="none") +
              theme(legend.title = element_blank()) + 
              scale_x_continuous(breaks = seq(0, 10, by = 2), labels = seq(0, 10, by = 2)) 

    if(justConds){
        return(condSurvs)
    }


    p2 <- id %>% filter(VF == 1) %>% rowwise() %>% 
    ggplot() +
        geom_bar(aes(x = 1,  fill = factor(numDetectableMuts)),
                 col = "white",  lwd = 0.05) + 
        theme_void() + 
        scale_fill_manual(values = numMutCols) +
        labs(x = "", y = "Proportion", fill = "Number\nof DRMs") +
        coord_flip() +
        theme( strip.background = element_blank(), strip.text.y = element_blank()) +
        theme(legend.position="none")

    p3 <- id %>% filter(VF == 1 & numDetectableMuts == 1) %>% 
        mutate(drugtype = ifelse(f1 > 0.15, "3TC", 
                                      ifelse(f2 > 0.15, "NRTI", "PI"))) %>% 
        ggplot() + 
        geom_bar(aes(x = 1,fill = drugtype), col = "white",
             width = 1, lwd = .2) + 
        theme_minimal() + coord_flip() +
        theme(strip.background = element_blank(), strip.text.y = element_blank() ) + 
        labs(y = "Proportion", fill = "Mutation type") +
        theme(axis.title.y = element_blank()) + 
        scale_fill_manual(values = cols)  +theme_void() + 
        theme(legend.position="none") + theme(legend.title = element_blank()) 

    return( list(plot_grid(p1, p2, p3, rel_heights = c(60, 10, 10), nrow = 3), condSurvs))

}


numsamp <- 1500

#Read in temporal heterogeneity data

detect <- 5000
dirv <- "time_model/"
dirp <- (paste0("../dat/", dirv))
fis <- list.files(dirp)
fis <- fis[file.info(paste0(dirp, fis))$size > 0]
set.seed(934119)
fis <- sample(fis , numsamp)
NFV.time <- readInTime(fis)


#Read in spatial heterogeneity data
detect  <- 150
dirv <- "space_model/"
dirp <- (paste0("../dat/", dirv))
dirfiles <- list.files(dirp)
dirfiles <- dirfiles[file.info(paste0(dirp, dirfiles))$size > 0]
set.seed(934119)
dirfiles <- sample(list.files(dirp) , numsamp)
NFV.space <- readInSpace(dirfiles)


dirv <- "space_model_narrow/"
dirp <- (paste0("../dat/", dirv))
dirfiles <- list.files(dirp)
dirfiles <- dirfiles[file.info(paste0(dirp, dirfiles))$size > 0]
set.seed(934119)
dirfiles <- sample(list.files(dirp) , numsamp)
NFV.space.narrowNe <- readInSpace(dirfiles)



#Under the temporal model, 
#given that a patient doesn't fail in the first year, 
#what is the probaability that they fail?
NFV.time %>% filter(time >= 1) %>% 
    summarize(suppressedFor10Years = mean(is.na(numDetectableMuts)))


#IQRs for narrow and normal dist'n
NFV.space %>% summarize(quantile(Ne, 0.25), 
                       quantile(Ne, 0.75)
                       )

NFV.space.narrowNe %>% summarize(quantile(Ne, 0.25), 
                       quantile(Ne, 0.75)
                       )



nfv.t <- plot_sims(NFV.time)
nfv.s <- plot_sims(NFV.space)

space.narrow <- plot_sims(NFV.space.narrowNe)

plotToPrint <- plot_grid(nfv.t[[1]], nfv.s[[1]], nrow = 1)

jpeg('../../HIV_MDR/figures/F5_data.jpg', width = 6.5, height = 3, units = "in", res = 300)
print(plotToPrint)
dev.off()

pdf('../../HIV_MDR/figures/F5_supp2.pdf', height = 4, width = 5)
space.narrow[[2]] %>% mutate(point.estimate = 1 - point.estimate, 
                        X95_conf_low = 1 - X95_conf_low, 
                        X95_conf_high = 1 - X95_conf_high) %>% 
                            filter(treatmentType != "total") %>%
                                mutate(treatmentType = ifelse(treatmentType == 1, "3TC", 
                                                       ifelse(treatmentType == 2, "NRTI",
                                                              ifelse(treatmentType == 3, "PI", NA)))) %>%     ggplot() + 
    geom_point(aes(x =  Year, y = point.estimate, col = treatmentType)) + 
    geom_line(aes(x = Year , y = point.estimate, col = treatmentType, group = treatmentType)) + 
    geom_segment(aes(x = Year , xend = Year, y = X95_conf_low, yend = X95_conf_high, col = treatmentType, group = treatmentType))  + scale_color_manual(values = cols) + labs(x = "Year of treatment", y = "Yearly probability of resistance evolution") + theme_classic() + theme(legend.position="none") + theme(legend.title = element_blank()) + 
        scale_x_continuous(breaks = seq(0, 10, by = 2), labels = seq(0, 10, by = 2)) 
dev.off()



#We would also like to do some sort of comparison of patient end points on Ne
pdf('../../HIV_MDR/figures/F5_supp1.pdf', height = 4, width =5)
NFV.space %>% ggplot() + geom_jitter(aes(x = t, y = Ne), width = 0.05, alpha = 0.25) + 
        scale_x_continuous(breaks = seq(0, 10, by = 2), labels = seq(0, 10, by = 2)) + 
        scale_y_log10(labels=fancy_scientific) + 
            theme_classic() + labs(x = "Years of treatment until clinical endpoint reached",
                                   y = "Plasma population size")
dev.off()

