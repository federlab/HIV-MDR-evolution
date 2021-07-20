
#get all possible codons and amino acids
allcodons<-c(); allAA<-c()
for (i in c("a","c","g","t")){
  for (j in c("a","c","g","t")){
    for (k in c("a","c","g","t")){
      codon=paste(i,j,k,collapse ="",sep="")
      allcodons<-c(allcodons,codon)
      allAA<-c(allAA,translate(c(i,j,k)))
    }
  }
}

#Read consensus data (from Los Alamos Database)
consensusfasta<-read.dna("../dat/HIV1_CON_2004_POL_DNA.fasta", format = "fasta",as.character=TRUE)	
#where is the start of POL? 
polstart=regexpr("cctca",paste(consensusfasta[which(row.names(consensusfasta)=="CONSENSUS_B"),],collapse=""))[1]
#consensusB_Pol contains B sequence for first 1500 nucleotides of Pol 
consensusB_Pol<-consensusfasta[which(row.names(consensusfasta)=="CONSENSUS_B"), polstart:(polstart+1500)]

#read file with mutation rates (source Abrams paper
#Abram ME, Ferris AL, Shao W, Alvord WG, Hughes SH. Nature, position, and frequency of mutations made in a single cycle of HIV-1 replication. J Virol. 2010; 84(19):9864–9878. https://doi.org/10.1128/ JVI.00915-10 PMID: 20660205)
read.csv("../dat/HIVMutRates.csv")->mutrates

#Read the file with the WHO list of mutations
read.csv(file = "../dat/who-muts.csv", header = FALSE)->WHOmutations
names(WHOmutations)<-c("pos", "mut", "drug", "gene", "class")

#Now I want to add the WT codon (from consensus B) for each position. 
for (i in 1:nrow(WHOmutations)){
  if (WHOmutations$gene[i]=="RT"){AApos = WHOmutations$pos[i]+99}
  if (WHOmutations$gene[i]=="PR"){AApos = WHOmutations$pos[i]}
  WHOmutations$wt[i]<-translate(consensusB_Pol)[AApos]
  WHOmutations$wtcodon[i]<-paste(consensusB_Pol[(3*AApos-2):(3*AApos)],collapse = "")
}

#Now, for each of the AAs in column "pos", I want a separate row. 
#create an empty new data frame
WHO_long<-WHOmutations[1,] #copy existing dataframe
WHO_long<-WHO_long[-1,] #Make the dataframe empty
WHO_long$mut<-as.character(WHO_long$mut) #the original list of levels isn't enough, so we change to character
#go through all rows of WHOmutations and put the long form into WHO_long
for (i in 1:nrow(WHOmutations)){
  print(WHOmutations$pos[i])
  #read the number of letters in WHOmutations$pos[i], store in num_res_AA
  num_res_AA = nchar(as.character(WHOmutations$mut[i]))
  #for each of the num_res_AA AAs, create a new row in a new dataframe, copying everything from WHOmutations[i,], but only letter j from WHOmutations$pos[i] 
  for (j in 1:num_res_AA){
    WHO_long[nrow(WHO_long)+1,]<-WHOmutations[i,] #create new row in WHO_long
    relevantAA<-substr(as.character(WHOmutations$mut[i]), j,j) #determine the relevant mutant AA
    print(relevantAA)
    WHO_long$mut[nrow(WHO_long)]<-relevantAA #put the relevant AA in the right row and column
  }
}

#Now I need to look at every line and calculate the mutation rate to get to the resistance mutation! 
WHO_long$resmutrate<-0
#Go through all rows of WHO_long
for (i in 1:nrow(WHO_long)){
  print(" ")
  print(paste("row=",i))
  #for each row get the resistant mut and get the indeces in allAA for that AA. 
  #There could be more than one index because of the redundancy of the genetic code
  whichAA<-which(allAA%in%as.character(WHO_long$mut[i]))
  #and get all possible codons that make the resistant mut
  possibleREScodons<-allcodons[whichAA]
  print(possibleREScodons)
  #so now, for each of these codons, I have to determine if we can reach it in one step and what the mut. rate is for that step. 
  wtcodon=WHO_long$wtcodon[i]; resmutrate=0
  print(wtcodon) 
  #let's determine all possible one step mutants from the WT
  for (j in 1:3){ #for each position in the wt codon, 
    print(paste("nuc pos =", j))
    if (substr(wtcodon, j, j) =="a"){muts<-c("g","t","c")}
    if (substr(wtcodon, j, j) =="c"){muts<-c("t","a","g")}
    if (substr(wtcodon, j, j) =="g"){muts<-c("a","c","t")}
    if (substr(wtcodon, j, j) =="t"){muts<-c("c","a","g")}
    for (l in muts){ 
      #create each one-step mutant for all three possible muts 
      mutcodon = wtcodon; substr(mutcodon, j, j) <- l #determine a one step mutant
      if (mutcodon %in% possibleREScodons){ 
        #if the one-step mutant is actually a relevant resistant codon, then we need to determine the mutation rate for that one-step mutant
        #next few lines to determine mutation type
        WTnuc=toupper(substr(wtcodon, j, j));if (WTnuc=="T") WTnuc="U" 
        MUTnuc=toupper(l); if (MUTnuc=="T") MUTnuc="U"
        mutationtype=paste(WTnuc,MUTnuc,collapse = "",sep="")
        print(mutationtype)
        #get the relevant mutation rate from the mutrates df (from source X)
        relevant_one_step_mutrate = mutrates$Probability[which(mutrates$Nucleotide.substitution==mutationtype)]
        resmutrate = resmutrate + relevant_one_step_mutrate
      }}
  }
  WHO_long$resmutrate[i]<-resmutrate
}

ListDrugs<-unique(WHO_long$drug)
ListMutTargetRate<-c()
for (d in ListDrugs){
#  print(d)
  ListMutTargetRate<-c(ListMutTargetRate,sum(WHO_long$resmutrate[WHO_long$drug==d]))
}

ListMutRates<-data.frame(ListDrugs, ListMutTargetRate)

#write.csv(x = ListMutRates, file = "ListMutRatesPerDrug.csv")
#pdf(width = 12, height = 6, file="targetMutationRate.pdf")
#barplot(ListMutTargetRate, names.arg = ListDrugs, main = "Total mutation rate towards resistance")
#dev.off()

#source("RResistanceMutations.r")
AllMuts <- WHO_long

#Which are the resistant codons?
for (i in 1:nrow(AllMuts)){
  resAA<-as.vector(strsplit(AllMuts$mut[i],split=""))
  rescodons<-allcodons[which(allAA%in%resAA[[1]])] #I need to change this to include all possible mutations

  #from the WT, all 1-step mutants (i include the wt because  easier to program)
  onestepmuts<-c();wtcodon=AllMuts$wtcodon[i]
  resmutrates=c()
  for (j in 1:3){
    if (substr(wtcodon, j, j) =="a"){muts<-c("g","t","c")}
    if (substr(wtcodon, j, j) =="c"){muts<-c("t","a","g")}
    if (substr(wtcodon, j, j) =="g"){muts<-c("a","c","t")}
    if (substr(wtcodon, j, j) =="t"){muts<-c("c","a","g")}
    for (l in muts){
      mutcodon = wtcodon
      substr(mutcodon, j, j) <- l
      onestepmuts<-c(onestepmuts, mutcodon)
      x=toupper(substr(wtcodon, j, j));if (x=="T") x="U"
      y=toupper(l); if (y=="T") y="U"
      mutationtype=paste(x,y,collapse = "",sep="")
      resmutrates=c(resmutrates,mutrates$Probability[which(mutrates$Nucleotide.substitution==mutationtype)])
    }}

  #get the highest mut rate for a one step mut that leads to a res codon!
  whichonesteps<-which(onestepmuts%in%rescodons)
  onestepmuts<-onestepmuts[whichonesteps]
  resmutrates<-resmutrates[whichonesteps]

  if (length(whichonesteps)>0){
    m=sum(resmutrates) #resmutrates[which.max(resmutrates)]
    o=onestepmuts[which.max(resmutrates)] #this takes the most likely mutation
    AllMuts$mutrate[i]<-m
    AllMuts$rescodon[i]<-o
  }
}

#AF June 16
CummMutRateNNRTI = 1- prod(1-AllMuts$mutrate[which(AllMuts$class=="NNRTI")])
CummMutRateNRTI = 1- prod(1-AllMuts$mutrate[which(AllMuts$class=="NRTI")])
CummMutRate184 = 1- prod(1-AllMuts$mutrate[which(AllMuts$class=="M184VI")])
CummMutRatePI = 1- prod(1-AllMuts$mutrate[which(AllMuts$class=="PI")])

#AF June 16
tbl_df(AllMuts) %>% group_by(drug) %>% summarize(mutrate = sum(resmutrate))
write.csv(x = ListMutRates, file = "../processed/ListMutRatesPerDrug.csv")

ListMutRates

#For NNRTI treatments
ProbNNRTIFirst<-CummMutRateNNRTI/(CummMutRateNNRTI+CummMutRateNRTI+CummMutRate184)
ProbNRTIFirst<-CummMutRateNRTI/(CummMutRateNNRTI+CummMutRateNRTI+CummMutRate184)
Prob184First<-CummMutRate184/(CummMutRateNNRTI+CummMutRateNRTI+CummMutRate184)

Pred<-data.frame(NNRTI_based = c(0,ProbNNRTIFirst,ProbNRTIFirst,Prob184First), PI_based = rep(0,4),
                 row.names=c("PI","NNRTI","NRTI","184"))

#For PI treatments
ProbPIFirst<-CummMutRatePI/(CummMutRatePI+CummMutRateNRTI+CummMutRate184)
ProbNRTIFirst<-CummMutRateNRTI/(CummMutRatePI+CummMutRateNRTI+CummMutRate184)
Prob184First<-CummMutRate184/(CummMutRatePI+CummMutRateNRTI+CummMutRate184)
Pred[c(3,4,1),2]<-c(ProbNRTIFirst,Prob184First,ProbPIFirst)

pdf("whichMutationFirst.pdf",width=10,height=8)
barplot(as.matrix(Pred),col=c("yellow","red","lightblue","blue"),
        horiz=TRUE, cex.axis=2, main = "Which resistance comes first? Predictions",cex.main=2,cex.lab=1.6, cex.names=1.6)
#legend(x=0.2,y=0.2,legend=row.names(Pred))

text(0.14,0.65,"NNRTI first",cex=1.5)
text(0.58,0.65,"NRTI first",cex=1.5)
text(0.92,0.65,"184 first",,cex=1.5)

text(0.2,1.9,"PI first",cex=1.5)
text(0.65,1.9,"NRTI first",cex=1.5)
text(0.935,1.9,"184 first",,cex=1.5)
dev.off()

