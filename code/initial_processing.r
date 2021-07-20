
dat <- tbl_df(read.table("../dat/dataset.01.24.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t"))

maxRTlastAA <- max(dat$RTLastAA)
maxPRlastAA <- as.numeric(max((dat %>% filter((PRLastAA) != "NULL"))$PRLastAA ))


#Let's start by fixing the treatment names, 
#because that will allow us to throw away all the 
#patients with useless treatments. 

processedNames <- unlist(lapply(strsplit(dat$RegimenName, split = "\\+"),
                         function(x){paste0(sort(trimws(x)),collapse = "+")}))
dat <- dat %>% mutate(Reg = processedNames)
#Filter any patients where we don't know the specific treatment
dat <- dat %>% filter(!grepl("RTI|PI|known|aAPA", Reg))
dat <- dat %>% mutate(RTSeq = toupper(RTSeq),
                      PRSeq = toupper(PRSeq))

ambs <- function(nuc){
        if(grepl("A|T|C|G", nuc)){return(nuc)}
        if(nuc == "R"){ return(c("A", "G")) }
        if(nuc == "Y"){ return(c("C", "T")) }
        if(nuc == "S"){ return(c("G", "C")) }
        if(nuc == "W"){ return(c("A", "T")) }
        if(nuc == "K"){ return(c("G", "T")) }
        if(nuc == "M"){ return(c("A", "C")) }
        if(nuc == "B"){ return(c("C", "G", "T")) }
        if(nuc == "D"){ return(c("A", "G", "T")) }
        if(nuc == "H"){ return(c("A", "C", "T")) }
        if(nuc == "V"){ return(c("A", "C", "G")) }
        if(nuc == "N"){ return(c("A", "C", "G", "T")) }
        return(NA)
    }       

nucPossibilities <- function(nucvect){

    return(expand.grid(ambs(nucvect[1]), ambs(nucvect[2]), 
                       ambs(nucvect[3]), stringsAsFactors = FALSE))


    
}

process <- function(x){

    #This function is just going to translate the initial data into workable alleles

    x <- tbl_df(x)

    #Alignment of RT
    RTseq <- strsplit(x$RTSeq, split = "")[[1]]

    internalCounter <- 0

    #Check each position
    RTalign <- foreach(i = 1:maxRTlastAA, .combine = "rbind")%do%{

        if(i < x$RTFirstAA | i > x$RTLastAA){
            return(c(i, NA))
        }
        
        internalCounter <- internalCounter + 1
        toTranslate <- RTseq[(internalCounter*3 - 2):(internalCounter * 3)]

        #Here, if toTranslate is not just composed of ATCG, 
        #I want to handle this as an ambiguous read
        #Simply expand the vector 
        transOpts <- nucPossibilities(toTranslate)

        allOpts <- foreach(k = 1:nrow(transOpts), .combine = paste0)%do%{
            if(sum(is.na(transOpts[k,])) > 0){
                return(NA)                
            }else{
                return(translate(as.matrix(transOpts[k,])))
            }
            
        }

        return(c(i, allOpts))
    }


    colnames(RTalign) <- c("pos",  "AA")

    RTalignp <- tbl_df(RTalign) %>% mutate(pos = as.numeric(pos)) %>% 
        mutate(PtID = x$PtID) %>% mutate(pos = paste0("RT", pos)) %>% 
        spread(pos, AA)


        #Alignment of RT
    PRseq <- strsplit(x$PRSeq, split = "")[[1]]

    if(is.na(x$PRFirstAA) | x$PRFirstAA == "NULL" ){

        PRalign <- tbl_df(1:maxPRlastAA) %>% mutate(pos = NA)

    }else{

        prfirst <- as.numeric(x$PRFirstAA)
        prlast <- as.numeric(x$PRLastAA)
    
        internalCounter <- 0
        PRalign <- foreach(i = 1:maxPRlastAA, .combine = "rbind")%do%{

            if(i < prfirst | i > prlast){
                return(c(i, NA))
            }
        
            internalCounter <- internalCounter + 1
            toTranslate <- PRseq[(internalCounter*3 - 2):(internalCounter * 3)]

        #Here, if toTranslate is not just composed of ATCG, 
        #I want to handle this as an ambiguous read

        #Simply expand the vector 
            transOpts <- nucPossibilities(toTranslate)
            allOpts <- foreach(k = 1:nrow(transOpts), .combine = paste0)%do%{
                if(sum(is.na(transOpts[k,])) > 0){
                    return(NA)                
                }else{
                    return(translate(as.matrix(transOpts[k,])))
                }
            }

            return(c(i, allOpts))
        }
    }
    colnames(PRalign) <- c("pos",  "AA")
    PRalignp <- tbl_df(PRalign) %>% mutate(pos = as.numeric(pos)) %>%
        mutate(PtID = x$PtID) %>% mutate(pos = paste0("PR", pos)) %>% 
            spread(pos, AA)

    return(left_join(left_join(x, RTalignp, by = "PtID"),
                      PRalignp, by = "PtID") %>%
            select(-RTFirstAA, -RTLastAA, -RTSeq, -PRFirstAA, -PRLastAA, -PRSeq))

}


proc <- dat %>% rowwise() %>% do(process(.))

write.table(proc, "../processed/stanford_hivdb_long.txt", col.names = TRUE, row.names = FALSE)


