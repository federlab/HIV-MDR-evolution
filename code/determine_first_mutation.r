byMut <- tbl_df(read.table("../processed/processed_mutations.txt", stringsAsFactors = FALSE, header = TRUE))
byDrug <- tbl_df(read.table("../processed/processed_drugs.txt", stringsAsFactors = FALSE, header = TRUE))

takeFirst <- function(x){return(x[1,])}


tmp <- byMut %>% group_by(PtID, Reg, Weeks) %>% 
    filter(!is.na(id_full)) %>% filter(id > 0) %>%
        group_by(PtID, MedlineID, Weeks, mut, AA_full) %>%
        do(takeFirst(.)) %>% group_by(PtID, Reg, Weeks) %>% 
    mutate(fs = sum(id_full == 1 , na.rm = TRUE), 
           ps = sum(id_full == 2 , na.rm = TRUE), 
           s = fs + ps) %>% ungroup()


#This is an example of sweep of I and V at RT184 
#This is a soft sweep at 184, and counted as full
tmp %>% filter(PtID == 6011)%>% select(-MedlineID)

#Now, what about actual partial sweeps
tmp %>% filter(PtID == 81805) %>% select(-MedlineID)
#Here, RT100 and RT184 are counted as partial sweeps, 
#because the ambiguous forms match non-resistant alleles



#Here we would like to pull out all patients with exactly one DRM
sings <- tmp %>% group_by(PtID, Reg) %>%
    #which can be accomplished by having exactly one full sweep
    filter((fs == 1 &  id_full == 1) | 
#    filter((s == 1 &  id_full == 1) | 
    #OR, by having only a partial sweep
           (id_full == 2 & fs == 0 & ps == 1 )) %>%
    #some patients have multiple timepoints - we want to take the first one
    #the fulfills this condition
    arrange(PtID, Weeks) %>% do(takeFirst(.)) %>% ungroup()


mutkey <- byDrug %>% select(mut, type) %>% distinct()

toPlot <- left_join(sings, mutkey, by = "mut") %>% 
    group_by(Reg) %>% 
        summarize(NNRTI = sum(type == "NNRTI"), 
                  NRTI = sum(type == "NRTI"),
                  PI = sum(type == "PI"), 
                  M184VI = sum(type == "M184VI"),
                  n = n(), .groups = "drop") %>% 
    ungroup() %>% gather(drugtype, count, -Reg, -n) %>% filter(n >= 10) 

#What are non-M184V/I mutations that confer resistance to 3TC or FTC?
byDrug %>% filter(reg == "3TC")
byDrug %>% filter(reg == "FTC")

#We are going to test that mutations that confer resistance to 3TC/FTC
#occur more often than mutations that do not
#For the purpose of plotting, we are only going to plot M184V
#but let's figure out how many instances 
other3tcs <- left_join(sings, mutkey, by = "mut")  %>% 
    filter(mut == "RT65" | mut == "RT151") %>% group_by(Reg, mut) %>%
        summarize(n())

#There really are not very many cases here that will affect our analysis
other3tcs %>% rowwise() %>% 
    mutate(numDrugs = length(strsplit(Reg, "\\+")[[1]])) %>% 
    filter(numDrugs == 3) %>% ungroup()

# 3TC+AZT+TDF is a triple NRTI therapy - we don't analyze it further
# 3TC+EFV+TDF/3TC+D4T+EFV are NNRTI-based therapies,
#   so it's not important to distingiush 
#   NRTI-only resistance mutations from NRTI/3TC resistance mutations


#For the purpose of plotting, we assign ABC and DDI to be in the NRTI class
# instead of the M184VI class (like FTC and 3TC). In practice, DDI does not appear
# in any triple drug therapies we examine. ABC appears in two cases which we mark 
# separately. 
drugKey <- byDrug %>% select(reg, type) %>% distinct() %>% 
    group_by(reg) %>% arrange(type) %>% summarize(type = type[1], .groups = "drop") %>%
    ungroup() %>%
    mutate(type = ifelse(type == "M184VI" & grepl('ABC|DDI', reg), "NRTI", type))

#Check that the correct drugs are applied to the correct category
byDrug %>% select(reg, type) %>% distinct() %>% 
    group_by(reg) %>% arrange(type) %>% summarize(type = type[1]) %>%
        ungroup()


#These are regularized category types
regCatsToAdd <- unlist(lapply(strsplit(toPlot$Reg, "\\+"), function(x){
    tmp = tbl_df(x)
    names(tmp) <- "reg"
return(paste(unique(sort(left_join(tmp, drugKey, by = "reg")$type)), collapse = "-"))
}))

#Order the factors for plotting purposes
regCatsToAdd <- factor(regCatsToAdd,
          levels = c(
              "NRTI", 
              "M184VI", 
              "M184VI-NRTI",
              "M184VI-NNRTI-NRTI",
              "M184VI-NNRTI",
              "NNRTI-NRTI", 
              "NNRTI", 
              "M184VI-NNRTI-NRTI-PI",
              "NNRTI-NRTI-PI",
              "NNRTI-PI",
              "NRTI-PI", 
              "M184VI-PI", 
              "M184VI-NRTI-PI" 
          ))

toPlot <- toPlot %>% mutate(regCats = regCatsToAdd)




write.table(toPlot, "../processed/summarized_first_mutations.txt", col.names = TRUE, row.names = FALSE)



simpleexpand <- function(x){
    toRet <- foreach(i = 1:(x$count), .combine = "rbind")%do%{
        x
    }
    return(toRet)
}

tmp_plot<- toPlot %>% group_by(Reg, drugtype) %>% filter(count > 0) %>% 
    arrange(count, drugtype) %>%
    do(simpleexpand(.)) %>% group_by(Reg) %>% mutate(c = (1)/n)

orderedFacts <- c(
                  (toPlot %>% group_by(regCats) %>% mutate(count = count/n) %>% 
    spread(drugtype, count) %>% arrange(regCats, PI,NNRTI, M184VI, NRTI  ) )$Reg, 
"Expectation")

plot_noexp <- tmp_plot %>% 
    filter(is.element(regCats, c( "M184VI-NRTI-PI",  "M184VI-NNRTI-NRTI"))) %>%
    ungroup() %>% 
        mutate(Reg = factor(Reg, levels = orderedFacts)) %>%
filter(Reg != "3TC+ABC+AZT+EFV" ) %>% 
#I'd also like to organize the drugs based on 3TC/FTC first
#Then other NRTI
#Then NNRTI or PI last
     mutate(Ap = ifelse(grepl("3TC", Reg), "3TC", "FTC")) %>%
     mutate(Bp = ifelse(grepl("AZT", Reg), "AZT", 
                 ifelse(grepl("D4T", Reg), "D4T", 
                 ifelse(grepl("ABC", Reg), "ABC", 
                 ifelse(grepl("TDF", Reg), "TDF", NA))))) %>% 
     mutate(Cp = ifelse(grepl("EFV", Reg), "EFV", 
                 ifelse(grepl("LPV", Reg), "LPV", 
                 ifelse(grepl("NFV", Reg), "NFV", 
                 ifelse(grepl("NVP", Reg), "NVP", 
                 ifelse(grepl("IDV", Reg), "IDV", 
                 ifelse(grepl("SQV", Reg), "SQV", NA)))))))  %>% 
mutate(Reg = paste(Ap, Bp, Cp, sep = "-"))


ordReg <- (plot_noexp %>% ungroup() %>%  select(Reg, n) %>%
    unique() %>% 
    separate(Reg, into = c("A", "B", "C")) %>%
        arrange(desc(C), desc(A), desc(B)) %>% 
            mutate(Reg = paste0(paste(A, B, C, sep = "-"), " (n=", n, ")")))$Reg

Cplot <- plot_noexp %>% mutate(Reg = factor(paste0(Reg, " (n=", n, ")"), levels = ordReg))





numMuts <- byMut %>% group_by(PtID, Reg, Weeks) %>% 
    filter(!is.na(id_full)) %>% 
        group_by(PtID, MedlineID, Weeks, mut, AA_full) %>%
        do(takeFirst(.)) %>% group_by(PtID, Reg, Weeks) %>% 
    mutate(fs = sum(id_full == 1 , na.rm = TRUE), 
           ps = sum(id_full == 2 , na.rm = TRUE), 
           s = fs + ps) %>%  group_by(PtID, Reg)  %>% 
               mutate(hasMut = id_full >= 1) %>% summarize(numMuts = sum(hasMut))

regCatsToAdd <- unlist(lapply(strsplit(numMuts$Reg, "\\+"), function(x){
    tmp = tbl_df(x)
    names(tmp) <- "reg"
return(paste(unique(sort(left_join(tmp, drugKey, by = "reg")$type)), collapse = "-"))
}))


numMuts <- numMuts %>% ungroup() %>% mutate(regCats = factor(regCatsToAdd,
          levels = c(
              "NRTI", 
              "M184VI", 
              "M184VI-NRTI",
              "M184VI-NNRTI-NRTI",
              "M184VI-NNRTI",
              "NNRTI-NRTI", 
              "NNRTI", 
              "M184VI-NNRTI-NRTI-PI",
              "NNRTI-NRTI-PI",
              "NNRTI-PI",
              "NRTI-PI", 
              "M184VI-PI", 
              "M184VI-NRTI-PI" 
          )))

numMutsPlot <- numMuts %>% group_by(PtID, Reg, regCats) %>% 
    filter(is.element(regCats, c( "M184VI-NRTI-PI",  "M184VI-NNRTI-NRTI"))) %>%
        mutate(numMuts = ifelse(numMuts >= 3, "3+", paste(numMuts))) %>% 
    group_by(Reg, regCats) %>% mutate(y = (1/n()), n = n()) %>% 
        filter(n > 20) %>% filter(Reg != "3TC+ABC+AZT+EFV" ) %>% 
#I'd also like to organize the drugs based on 3TC/FTC first
#Then other NRTI
#Then NNRTI or PI last
#     separate(Reg, into = c("A", "B", "C")) %>%
     mutate(Ap = ifelse(grepl("3TC", Reg), "3TC", "FTC")) %>%
     mutate(Bp = ifelse(grepl("AZT", Reg), "AZT", 
                 ifelse(grepl("D4T", Reg), "D4T", 
                 ifelse(grepl("ABC", Reg), "ABC", 
                 ifelse(grepl("TDF", Reg), "TDF", NA))))) %>% 
     mutate(Cp = ifelse(grepl("EFV", Reg), "EFV", 
                 ifelse(grepl("LPV", Reg), "LPV", 
                 ifelse(grepl("NFV", Reg), "NFV", 
                 ifelse(grepl("NVP", Reg), "NVP", 
                 ifelse(grepl("IDV", Reg), "IDV", 
                 ifelse(grepl("SQV", Reg), "SQV", NA)))))))  %>% 
mutate(Reg = paste(Ap, Bp, Cp, sep = "-"))

ordReg <- (numMutsPlot %>% ungroup() %>%  select(Reg, n) %>%
    unique() %>% 
    separate(Reg, into = c("A", "B", "C")) %>%
        arrange(desc(C), desc(A), desc(B)) %>% 
            mutate(Reg = paste0(paste(A, B, C, sep = "-"), " (n=", n, ")")))$Reg

Bplot <- numMutsPlot %>% 
    mutate(Reg = paste0(Reg, " (n=", n, ")")) %>% 
        mutate(Reg = factor(Reg, levels = ordReg))

write.table(Bplot, "../processed/data_for_figure_2b.txt", col.names = TRUE, row.names = FALSE)
write.table(Cplot, "../processed/data_for_figure_2c.txt", col.names = TRUE, row.names = FALSE)
