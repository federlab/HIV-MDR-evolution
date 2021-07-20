
#This inputs the formatted data from the processing step. 
#We will use this information to call DRMs
proc <- tbl_df(read.table("../processed/stanford_hivdb_long.txt", header = TRUE, stringsAsFactors = FALSE))

mut.raw <- tbl_df(read.csv('../dat/who-muts.csv', sep = ",", header = FALSE, 
                           stringsAsFactors = FALSE))

names(mut.raw) <- c("pos", "AA", "reg", "gene", "type")


#Ok, so we would ultimately like to format for each drug
#the list of all possible mutations that confer resistance to that drug

NRTIs <- "AZT|3TC|ABC|DDI|D4T|FTC|TDF"

drugsToCheck <- c("AZT", unique((mut.raw %>% 
                                     filter(reg != 'NRTI'))$reg))

allDrugs <- unique(unlist(strsplit(proc$Reg, "\\+")))
excDrugs <- allDrugs[!is.element(allDrugs, drugsToCheck)]

# This includes only patients who have drug combinations consisting only of
# drugs we have information to check
proc.p <- proc %>% filter(!grepl(paste(excDrugs, collapse = "|"), Reg))


byDrug <- foreach(drug = drugsToCheck, .combine = "rbind")%do%{

    isNRTI <-  grepl(NRTIs, drug)
    mut.raw %>% filter(ifelse(isNRTI & reg == "NRTI" ,1, 0) |
                           ifelse(grepl(drug, reg), 1, 0)) %>%
                select(pos, AA, gene, type) %>% mutate(reg = drug) 
}

byDrug <- byDrug %>% mutate(mut = paste0(gene, pos))

allmutations <- tbl_df(unique((byDrug %>% mutate(mut = paste0(gene, pos)) )$mut))
names(allmutations) <- "mut"


#Ok, I want these added as columns
#Then, the goal will be to list one of four options:
#1: the mutation confers resistance and is present
#2: the mutation confers resistance and is partially present
#3: the mutation confers resistance and is absent
#4: the mutation does not confer resistance or is not sampled

proc.long <- proc.p %>% 
    select(PtID, MedlineID, Reg, Weeks, matches("RT[0-9]+"), matches("PR[0-9]+")) %>%
    gather("mut", "id", -PtID, -Reg, -Weeks, -MedlineID)  %>%
    group_by(PtID, MedlineID, Reg, Weeks, mut) %>% 
    summarize(id = paste0(unique(id), collapse = ""))

proc.long <- proc.long %>% spread(mut, id)


assessMuts <- function(seqVal, resVal){

    if(sum(is.na(seqVal) ) > 0 | sum(seqVal == "NA") > 0){ return(NA) }

    resVal <- strsplit(resVal, split = "")[[1]]
    seqVal<- unique(strsplit(seqVal, split = "")[[1]])
    
    #if all of the seqVals are present in resVals, count it as
    # resistant
    if(sum(!is.element(seqVal, resVal)) == 0){
        return(1)
    }

    #If none of the seqVals are present in resVals, count it as
    #partial
    if(sum(is.element(seqVal, resVal)) == 0){
        return(0)
    }

    return(2)

}






byMut <- foreach(rownum = 1:nrow(proc.long), .combine = "rbind")%do%{

    if(rownum %% 100 == 0){ print(rownum) }

    rowCheck <- proc.long[rownum,]

    relmuts <- byDrug %>% 
        filter(grepl(gsub("\\+", "|", rowCheck$Reg), reg)) %>% 
            select(-reg) %>% distinct()
   
    #check for duplicates
    dups <- relmuts$mut[duplicated(relmuts$mut)]
    if(length(dups) > 0){
        toAdd <- foreach(i = 1:length(dups), .combine = "rbind")%do%{
            combmut <- paste0(unique(unlist(strsplit((relmuts %>% filter(mut == dups[i]))$AA, split = ""))), collapse = "")
            (relmuts %>% filter(mut == dups[i]))[1,] %>% mutate(AA = combmut)
        }
        relmuts <- bind_rows(
            relmuts %>% filter(!grepl(paste(dups, collapse = "|"), mut)),
            toAdd)
    }


    relIdent<- tbl_df(foreach(mutrow = 1:nrow(relmuts), .combine = 'rbind')%do%{

        rowinf <- relmuts[mutrow,]
        seqVal <- paste(rowCheck %>% ungroup %>% select(rowinf$mut) )
        resVal <- rowinf$AA

        overall_ident <- assessMuts(seqVal, resVal)
        if(nchar(resVal) == 1){
            sub_ident <- (t(c(rowinf$mut, resVal, overall_ident, resVal, overall_ident)))
            colnames(sub_ident) <- paste0("V", 1:5)
            sub_ident <-  tbl_df(sub_ident)
        }else{
            sub_ident <- foreach(res = strsplit(resVal, "")[[1]], 
                             .combine = "rbind")%do%{
               c(rowinf$mut, res, assessMuts(seqVal, res), resVal, overall_ident)
           }
            colnames(sub_ident) <- paste0("V", 1:5)
            sub_ident <-  tbl_df(sub_ident)
        }
        return(sub_ident)
        
    })

    names(relIdent) <- c("mut", "AA", "id", "AA_full", "id_full")

    left_join(allmutations, relIdent, by = "mut") %>% 
        mutate(PtID = rowCheck$PtID, MedlineID = rowCheck$MedlineID,
               Reg = rowCheck$Reg, Weeks = rowCheck$Weeks)

}



write.table(byMut, "../processed/processed_mutations.txt", col.names = TRUE, row.names = FALSE)
write.table(byDrug, "../processed/processed_drugs.txt", col.names = TRUE, row.names = FALSE)




