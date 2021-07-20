
toPlot <- tbl_df(read.table("../processed/summarized_first_mutations.txt", stringsAsFactors = FALSE, header = TRUE))

ListMutRates <- read.table(file = "../processed/ListMutRatesPerDrug.csv", stringsAsFactors = FALSE, header = TRUE, sep = ",")


#Treatments to test
treatsToTest <- toPlot %>% 
    filter(regCats == "M184VI-NRTI-PI"  | regCats == "M184VI-NNRTI-NRTI") %>% 
    filter(Reg != "3TC+ABC+AZT+EFV") %>% arrange(Reg)


PIRegs <- treatsToTest %>% filter(regCats == "M184VI-NRTI-PI") %>% 
    summarize(Reg = unique(Reg))

NNRTIRegs <- treatsToTest %>% filter(regCats == "M184VI-NNRTI-NRTI") %>% 
    summarize(Reg = unique(Reg))

PI.pvals <- foreach(PI = PIRegs$Reg, .combine = "rbind")%do%{

    #For PIs, we're testing the hypothesis - M184V happens first
    checkReg <- treatsToTest  %>% filter(Reg == PI) %>%
        group_by(M184V = drugtype == "M184VI") %>% 
        summarize(count = sum(count), .groups = "drop") %>% arrange(M184V)

    splitPIName <- strsplit(PI, split = "\\+")[[1]]
    #sum the mutation rates of the non-3TC/FTC drugs

    namesToMerge = tbl_df(splitPIName)
    names(namesToMerge) <- c("ListDrugs")

    relMuts <- left_join(namesToMerge, ListMutRates, by = "ListDrugs") %>% 
        mutate(M184V = grepl("TC", ListDrugs)) %>% 
        group_by(M184V) %>% 
        summarize(relMuts = sum(ListMutTargetRate), .groups = "drop") %>% 
        mutate(relMuts = relMuts/sum(relMuts))

    bt = binom.test((checkReg %>% filter(M184V == TRUE))$count, 
               sum(checkReg$count),
               p = (relMuts %>% filter(M184V == TRUE))$relMuts, alternative = "greater")

    c(PI, "3TC/FTC first", bt$statistic, bt$parameter,  round(bt$null.value, 3), signif(bt$p.value, 3))

}

colnames(PI.pvals)[c(1, 2, 6)] <- c("reg", "category", "pval")
PI.pvals <- tbl_df(PI.pvals) %>% 
    mutate(pval = as.numeric(pval)) #%>%
#        filter(!grepl("ABC", reg))
#Here, I exclude ABC regimens because it's not clear what the hypothesis here
#should be

PI.toprint <- PI.pvals %>% 
    arrange(pval) %>% mutate(rank = 1:n()) %>% 
    mutate(BHval = rank/n() * 0.05) %>% 
    mutate(sig = pval < BHval)


NNRTI.pvals <- foreach(NNRTI = NNRTIRegs$Reg, .combine = "rbind")%do%{

    #For NNRTIs, we're testing the hypothesis - NNRTI happens first
    checkReg <- treatsToTest  %>% filter(Reg == NNRTI) %>%
        group_by(NNRTI = drugtype == "NNRTI") %>% 
        summarize(count = sum(count), .groups = "drop") %>% arrange(NNRTI)

    splitNNRTIName <- strsplit(NNRTI, split = "\\+")[[1]]
    #sum the mutation rates of the NNRTI drugs

    namesToMerge = tbl_df(splitNNRTIName)
    names(namesToMerge) <- c("ListDrugs")

    relMuts <- left_join(namesToMerge, ListMutRates, by = "ListDrugs") %>% 
        mutate(NNRTI = is.element(ListDrugs, c("EFV", "NVP"))) %>% 
        group_by(NNRTI) %>% 
        summarize(relMuts = sum(ListMutTargetRate), .groups = "drop") %>% 
        mutate(relMuts = relMuts/sum(relMuts))

    bt = binom.test((checkReg %>% filter(NNRTI == TRUE))$count, 
               sum(checkReg$count),
               p = (relMuts %>% filter(NNRTI == TRUE))$relMuts, alternative = "greater")

    c(NNRTI, "NNRTI first", bt$statistic, bt$parameter,  round(bt$null.value, 3), signif(bt$p.value, 3))
}



colnames(NNRTI.pvals)[c(1, 2, 6)] <- c("reg", "category", "pval")
NNRTI.pvals <- tbl_df(NNRTI.pvals) %>% 
    mutate(pval = as.numeric(pval)) 
#Here, I do NOT exclude ABC regimens because both ABC and 3TC are being counted in the same
#category


NNRTI.toprint <- NNRTI.pvals %>% 
    arrange(pval) %>% mutate(rank = 1:n()) %>% 
    mutate(BHval = rank/n() * 0.05) %>% 
    mutate(sig = pval < BHval)  

write.csv(NNRTI.toprint, "../tables/NNRTI_pvals.csv")
write.csv(PI.toprint, "../tables/PI_pvals.csv")







