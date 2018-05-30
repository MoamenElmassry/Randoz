###################################################Tamiru###Salmonella

#Location of working files
cd /panfs/pan1.be-md.ncbi.nlm.nih.gov/product_manager_research_projects/

#Blast database location
/net/frosty/vol/blast/db/disk.blast/newest_blast/blast/nucl_db

#To pull genome accession from a protein id for one protein or multiple in a file
esearch -query YP_001949744.1 -db protein | elink -target nuccore | efetch -format acc
for i in `cat prot_ids`;do esearch -query $i -db protein | elink -target nuccore | efetch -format acc >> acc.foo; done

#To pull protein fasta using accession
esearch -query YP_001949744.1 -db protein | efetch -format fasta

#Retrieve proteins for all phages using their accessions
for i in `cat phage_accns`;
do wget `esearch -db assembly -query $i | efetch -format docsum | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" '{print $0"/"$NF"_protein.faa.gz"}'`; done

#extraction
gunzip *gz

#concatenate proteins from phages
cat *.faa > typhi_ptyphii_accns.faa

#make d blst database from the proteins
makeblastdb -in typhi_ptyphii_accns.faa -dbtype prot -parse_seqids

#Blast the proteins in the typhi or the top 2 typhi agaisnt datasbe of all proteins in all phages
deltablast -db typhi_ptyphii_accns.faa -query GCF_000890795.1_ViralProj64609_protein.faa -outfmt 6 -evalue .00000001 -out NC_015271_proteins_into_typhi_acc_list -max_target_seqs 5000 &                  

#Retrieve DNA sequences for proteins accession IDs
First, extract acc numbers from NCBI with suitable filters
for i in `cat betaglucuronidase_accLIST.seq`;do efetch -db protein -format fasta_cds_na -id $i >> BG_seq.fasta ; done

#Retrieve DNA sequences for proteins accession IDs directly
esearch -db protein -query "rgpa" | efilter -source refseq -organism bacteria | efetch -db protein -format fasta_cds_na

#Retrieve protein sequences for proteins accession IDs directly
esearch -db protein -query "rgpa" | efilter -source refseq -organism bacteria | efetch -db protein -format fasta

#extract DNA seuqnce for uids from a file
 esearch -db protein -query prot | efilter -source refseq -organism human | efetch -format fasta_cds_na > human_metab_enz.fasta

#get protein acc from uid
esearch -db protein -query prot | efetch -db protein -format acc > prot.acc 



############################################################MEGABLAST

#This is the database of the refseq proteins realted to bacteria as of 5/7/2018
makeblastdb -in refseqbacterialprot.fasta -dbtype prot -parse_seqids

#To find bacterial proteins similar to that of humans that are related to drug metabolism
deltablast -db refseqbacterialprot.fasta -query human_protein.fasta -outfmt 6 -evalue .0000001 -out humanVSbacteria -max_target_seqs 5000 -num_threads 8 



#make a blast database
makeblastdb -in my_reference.fa -out my_reference -parse_seqids -dbtype nucl

#Megablast against SRA samples
magicblast -sra SRR6865476,SRR1031317,SRR1031519,SRR1031758,SRR1031317,SRR1765556,SRR2155334,SRR1818231,SRR5548988,SRR1031278,SRR1031478,SRR6845410,SRR6845425,SRR640383,SRR1031585,SRR1031153,SRR1031817,SRR6442942,SRR1031849,SRR638774,SRR640417,SRR639509,SRR5826792,SRR5826659,SRR5826679 -db my_reference -paired -num_threads 8 -outfmt tabular -no_unaligned -splice F > metagenome_blast_pred_bacterial

esearch -db sra -query "Simvastatin metagenome" |  efetch -format docsum | xtract -pattern Runs -ACC @acc  -element "&ACC"


transcriptome
magicblast -sra  -db betagluc -paired -num_threads 10 -outfmt tabular -no_unaligned -reftype transcriptome > betaglu_females





##########################
library(plyr)
#Read the table from magic-blast in
tab=read.table("betaglu_females",header=F)
#fixing column names
col_names=c("query.acc.","reference.acc.","% identity","not used","not.used",
            "not.used","query.start","query.end","reference.start","reference.end","not.used",
            "not.used","score","query.strand","reference.strand","query.length","BTOP",
            "num.placements","not.used","compartment","left.overhang","right.overhang",
            "mate.reference","mate.ref..start","composite.score")
colnames(tab)=col_names
#remove duplicate hits for eahc sequence
tab$seq=sub("..$","",tab$query.acc.)
tab2=tab[!duplicated(tab[,"seq"]),]
#clean SRA names to aggregate them
tab2$SRA=sub("\\..*","",tab2$seq)
#tab2$length=abs(tab2$query.start-tab2$query.end)
#tab2=tab2[tab2$length>=50,]
#tab2=tab2[tab2$X..identity>=90,]
#count freq of each enzyme for SRA
df <- count(tab2, c('SRA','reference.acc.'))
#save the clean table
# write.csv(df,"tablegene.csv")
#fuction for retrieving SRA MBS size
library(rentrez)
getSize <- function(ids) {size_mega <- c()
  for(i in ids) {term =  paste0(i, "[ACCN]")
    run = entrez_search(db = "sra", term = term)
    exp_descrip = entrez_summary(db = "sra", id = run[[1]])
    x = exp_descrip$run
    ## extract range defined by "total bases" and "load_done"
    size = substr(x, start = regexpr("total_bases=", x)[[1]][1] + attr(regexpr("total_bases=", x), "match.length"),stop = regexpr(" load_done", x)[[1]][1])
    ### some clean up
    size = gsub('\"', "", size, fixed = TRUE)
    ## converto to numeric and Mb
    size = as.numeric(size)/1e6
    size_mega <- c(size_mega, size)}
size_mega}
#retrieve the size of SRAs
df$MBS=getSize(df$SRA)
#Noprmalize using the size of each SRA
df$count=(df$freq*mean(df$MBS)/df$MBS)
#Transform to log10
df$trans_count=log2(df$count)
#Divide in tertiles or percentiles
df$group <- as.numeric(cut(df$trans_count, 4))

#Plotting boxplot of the pecentiles
g <- ggplot(subset(df, group %in% c(1, 2, 3)), aes(group, trans_count))
g + geom_violin(aes(fill=factor(group)))+ scale_fill_brewer(palette="Set1")+guides(fill=FALSE)+
   labs(x="Female Gut Metagenome",y="log2 of Beta-Glucuronidase Hits")+theme_classic(base_size = 30)+
   geom_point(color="black", position = "jitter")+scale_x_continuous(breaks=c(1,2,3),labels=c("Low","Medium","High"))

#Plotting boxplot of the pecentiles
g <- ggplot(df, aes(group, trans_count))
g + geom_violin(aes(fill=factor(group)))+ scale_fill_brewer(palette="Set1")+guides(fill=FALSE)+
   labs(x="Female Gut Metagenome",y="log2 of Beta-Glucuronidase Hits")+theme_classic(base_size = 30)+
   geom_point()+scale_x_continuous(breaks=c(1,2,3,4),labels=c("Low","Medium","High","Very High"))

#Plotting boxplot of the pecentiles
df$group <- as.numeric(cut(df$trans_count, 3))
g <- ggplot(df, aes(group, trans_count))
g + geom_boxplot(aes(fill=factor(group)))+ scale_fill_brewer(palette="Set1")+guides(fill=FALSE)+
   labs(x="Female Gut Metagenome",y="log10 of Beta-Glucuronidase Hits")+theme_classic(base_size = 30)+
   geom_point()+scale_x_continuous(breaks=c(1, 2,3),labels=c("Low", "Medium","High"))

#Plotting boxplot of the pecentiles
df$group <- as.numeric(cut(df$trans_count, 2))
g <- ggplot(df, aes(group, trans_count))
g + geom_boxplot(aes(fill=factor(group)))+ scale_fill_brewer(palette="Set1")+guides(fill=FALSE)+
   labs(x="Female Gut Metagenome",y="log10 of Beta-Glucuronidase Hits")+theme_classic(base_size = 30)+
   geom_point()+scale_x_continuous(breaks=c(1,2),labels=c("Low", "High"))

#######################################################
tab=read.csv("MBS.csv",header=T)
NeedMBS5=merge(NeedMBS, tab, all.x = TRUE)
#############################################CODE for getting MBS for SRA# Jose experimenting
entrez_dbs()
entrez_db_summary(db = "sra")
entrez_db_searchable(db = "sra")
run = entrez_search(db = "sra", term = "ERR2009629[ACCN]")

str(run)
run[[1]]

run = entrez_summary(db = "sra", id = run[[1]])
run[1:5]

faa  = entrez_fetch(db = "sra", id = run[[1]], rettype = "xml", parsed =TRUE)
faa

rootNode <- XML:: xmlRoot(faa)
rootNode[1]

data = XML::xmlSApply(rootNode,function(x) XML::xmlSApply(x, XML::xmlValue))
cd.catalog <- data.frame(t(data),row.names=NULL)

XML::xpathSApply(rootNode, "serica", XML::xmlValue)

mlSApply(rootNode,function(x) XML::xmlSApply(x, XML::xmlValue))

## 'getSize' takes a vector containing NCBI SRA accessions and returns and vector 
## with the run size (Mb) 
# Dependencies
library(rentrez)


getSize <- function(ids) {
  
  size_mega <- c()
  
  for(i in ids) {
    
    term =  paste0(i, "[ACCN]")
    run = entrez_search(db = "sra", term = term)
    exp_descrip = entrez_summary(db = "sra", id = run[[1]])
    x = exp_descrip$run
    
    ## extract range defined by "total bases" and "load_done"
    size = substr(x, start = regexpr("total_bases=", x)[[1]][1] + attr(regexpr("total_bases=", x), "match.length"),
                  stop = regexpr(" load_done", x)[[1]][1])
    ### some clean up
    size = gsub('\"', "", size, fixed = TRUE)
    
    
    ## converto to numeric and Mb
    size = as.numeric(size)/1e6
    
    
    size_mega <- c(size_mega, size)
    
  }
  
  size_mega
  
}

## Usage 
ids <- c("DRR071071", "ERR1912953")
getSize(ids)

#Another way to do it#####################################################
library(rentrez)
library(stringr)

getSize <- function(ids) {
  
  sizes <- c()
  
  for(i in ids) {
    
    term =  paste0(i, "[ACCN]")
    run = entrez_search(db = "sra", term = term)
    exp_descrip = entrez_summary(db = "sra", id = run[[1]])
    x = exp_descrip$run
    
    size = str_sub(x, start =str_locate(x, "total_bases=")[2]+2, 
                   end =str_locate(x, "load_done")[1]-3)
    
    size_mega = as.numeric(size)/1e6
    
    sizes <- c(sizes, size_mega)
    
  }
  
  sizes
  
}




## Dependencies
library(rentrez)

## Define function
get_CDDacc <- function(prot_vector) {
   
   # initializes a vector
   acc = vector(mode = "character", length = length(seq(prot_vector)))
   
   # loop over the protein ids list 
   for(i in seq(prot_vector)) {
      
      # NCBI Eutils through rentrez pckg.
      protein <-  entrez_summary(db="protein", id= prot_vector[i])
      p_links = entrez_link(dbfrom='protein', id=protein$uid, db='cdd')
      cdd = p_links$links$protein_cdd_concise_2
      
      # if CDD
      if(! is.null(cdd)) {
         region = entrez_summary(db="cdd", id= cdd)
         acc[i] = paste(region$accession,"_",region$abstract)
      } else {acc[i] = 0}
      
   }
   
   acc  
   
}

## Usage
xps = c("YP_001949746.1", "YP_001949758", "YP_001949763", "YP_001949766", "YP_001949773", "YP_001949775", "YP_001949776", "YP_001949796", "YP_001949797", "YP_001949798", "YP_004306652", "YP_004306653", "YP_004306654", "YP_004306656", "YP_004306662", "YP_004306669", "YP_004306672", "YP_004306675", "YP_004306679", "YP_004306680", "YP_004306696")
y=as.data.frame(xps)
y$CDD=get_CDDacc(y$xps)
library(plyr)
dff <- count(y, 'CDD')
dff <- dff[order(-dff$freq),]
