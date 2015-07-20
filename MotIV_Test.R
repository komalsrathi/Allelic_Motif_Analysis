#######################################  multiple TF motifs ####################################### 
require('ggplot2')
require('tidyr')
require('dplyr')
library(MotifDb)
library(seqLogo)
library(MotIV) 
library(BSgenome.Hsapiens.UCSC.hg19)
library(plyr)
library(doMC)
registerDoMC(4)

# setwd
setwd('projects/allelic_motif_analysis')

# input data & take +/- 10 bp from snp position
snps = read.table('cisetlSNPs.csv', skip = 1, sep=',', header=T)
colnames(snps) = c('ID','Chr','Pos','Minor','Major')
snps$start = snps$Pos-10
snps$end = snps$Pos+10
snps = unique(snps)

# load genome 
genome <- BSgenome.Hsapiens.UCSC.hg19

# MotifDB only Homo Sapiens
mdb.human <- MotifDb[grep("Hsapiens",values(MotifDb)$organism),] 

# import transfac data, get vertebrate motifs & convert to proportions 
motifs.transfac <- readPWMfile(file = 'meme_matrix_files/MotIV_transfac.txt')
motifs.transfac <- sapply(X = motifs.transfac, FUN = function(x) x/rowSums(x))
motifs.transfac <- do.call(list, rapply(motifs.transfac, function(x) ifelse(is.nan(x),0,x), how="replace"))
motifs.transfac <- motifs.transfac[grep('^V_',names(motifs.transfac))] # vertebrates only

# cardiac specific transcription factors only
# motifs.transfac <- motifs.transfac[grep('GATA4|MEF2A|MEF2C|MITF|TBX5_|TBX1_|NKX25',names(motifs.transfac))]

# jaspar data
jaspar.t <- sapply(X = jaspar, FUN = function(x) x/rowSums(x))
jaspar.t <- do.call(list, rapply(jaspar.t, function(x) ifelse(is.nan(x),0,x), how="replace"))

myfunc <- function(x,seq)
{
  hits <- matchPWM(pwm = x,subject = seq,with.score=T,min.score = .5)
  scores <- mcols(hits)$score
  start <- start(hits)
  end <- end(hits)
  res <- data.frame(score=scores,start=start,end=end)
  res <- res[which.max(res$score),]
  return(res)
}

# wrapper for myfunc
runAllelicMotif <- function(x, motifdb)
{
  print(x)
  # get sequence & its complement
  ref <- getSeq(genome,x$Chr,start=x$start,end=x$end,strand='+')
  ref.comp <- reverseComplement(ref)
  ref.res <- ldply(.data = motifdb, .fun = function(i) myfunc(i,ref))
  ref.comp.res <- ldply(.data = motifdb, .fun = function(i) myfunc(i,ref.comp))
  
  # alter sequence 
  alt <- replaceLetterAt(x = ref, at = x$end-x$Pos, letter = as.character(x$Major))
  alt.comp <- reverseComplement(alt)
  alt.res <- ldply(.data = motifdb, .fun = function(i) myfunc(i,alt))
  alt.comp.res <- ldply(.data = motifdb, .fun = function(i) myfunc(i,alt.comp))
  
  # merge
  res1 = merge(ref.res, alt.res, by=c('.id','start','end'), all=TRUE, suffixes = c('ref','alt'))
  res2 = merge(ref.comp.res, alt.comp.res, by=c('.id','start','end'), all=TRUE, suffixes = c('ref.comp','alt.comp'))
  res1$score.diff = res1$scoreref-res1$scorealt
  res2$score.comp.diff = res2$scoreref.comp-res2$scorealt.comp
  res = merge(res1, res2, by=c('.id','start','end'), all=TRUE)
  res = res[complete.cases(res),] 
  return(res)
}

# use mdb.human database
mdb.human.res = ddply(.data = snps, .variables = 'ID', .fun = function(x) runAllelicMotif(x, mdb.human))

# use transfac database
transfac.res = ddply(.data = snps, .variables = 'ID', .fun = function(x) runAllelicMotif(x, motifs.transfac))

# finding max altered TF
transfac.res$score.diff = abs(transfac.res$score.diff)
tmp = transfac.res[,1:7]
tmp = tmp[-which(tmp$score.diff==0),]
tt = ddply(tmp, .(ID), summarise, x = max(score.diff))
tmp = merge(tmp, tt, by.x=c('ID','score.diff'),by.y=c('ID','x'))
tmp = tmp[order(tmp$ID,tmp$score.diff),]
count_(x = tmp,vars = '.id',sort = TRUE)
