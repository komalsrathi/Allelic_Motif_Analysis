library(devtools)
library(atSNP)

load('atSNP.RData')

# ---------------------------------------- load motif library -------------------------------------- #
# transpose the motif database elements
mdb.human.t = sapply(X = mdb.human, FUN = function(x) t(x))
motifs.transfac.t = sapply(X = motifs.transfac, FUN = function(x) t(x))
jaspar.t <- sapply(X = jaspar.t, FUN = function(x) t(x))

# ---------------------------------------- load motif library -------------------------------------- #


# ---------------------------------------- load snp data from snp info file ------------------------------------------- #
snp_tbl = read.table('cisetlSNPs.csv',skip = 1,sep=',',header=T)
colnames(snp_tbl) = c('snpid','chr','snp','a1','a2')
snp_tbl = snp_tbl[,c(1,4,5,2,3)]

write.table(snp_tbl, file = "./atSNP_test/test_snp_file.txt", row.names = FALSE, quote = FALSE)

snpInfo <- LoadSNPData("./atSNP_test/test_snp_file.txt", genome.lib = "BSgenome.Hsapiens.UCSC.hg19", half.window.size = 30, default.par = T, mutation = F)
ncol(snpInfo$sequence) == nrow(snp_tbl)
snpInfo$rsid.rm
# ---------------------------------------- load snp data from snp info file ------------------------------------------- #

# ---------------------------------------- load snp data from rsid file ------------------------------------------- #
library(SNPlocs.Hsapiens.dbSNP.20120608)
snp.chr <- read.table('cis_eSNPS_LV-DNASe')
snp.chr <- as.character(snp.chr$V1)
snpInfo = LoadSNPData(snpids=snp.chr,genome.lib="BSgenome.Hsapiens.UCSC.hg19",
                      snp.lib='SNPlocs.Hsapiens.dbSNP.20120608',
                      half.window.size=30,default.par=T,mutation=F)
# ---------------------------------------- load snp data from rsid file ------------------------------------------- #

# ------------------------------- compute affinity scores ------------------------------------------- #
# motif_library <- mdb.human.t
motif_library <- motifs.transfac.t

atsnp.scores <- ComputeMotifScore(motif_library, snpInfo, ncores=20)

# ------------------------------- compute affinity scores ------------------------------------------- #

# ----------------------------------- compute p values ---------------------------------------------- #
atsnp.result <- ComputePValues(motif.lib = motif_library, snp.info = snpInfo, motif.scores = atsnp.scores$motif.scores, ncores = 20)
atsnp.result <- atsnp.result[order(pval_rank),list(snpid,motif,pval_ref,pval_snp,pval_rank)] # order by pval_rank (sig. of affinity changes)

atsnp.result[pval_rank<=0.1,] # pval_rank threshold
atsnp.result[,pval_rank_bh := p.adjust(pval_rank, method='BH')] # fdr adjustment
atsnp.result.sub = atsnp.result[pval_rank<=0.05,] #threshold on fdr adjusted pvalue
atsnp.result.sub.count = count_(x = atsnp.result.sub,vars = "motif",sort = T) # count motif occurances

atsnp.result[,pval_rank_bh := p.adjust(pval_rank, method='BH'),by = motif] # by motif if lot of snps
atsnp.result[,pval_rank_bh := p.adjust(pval_rank, method='BH'),by = snpid] # by snpid if many motifs
# ----------------------------------- compute p values ---------------------------------------------- #

# ----------------------------------- additional analysis ---------------------------------------------- #
match_result <- MatchSubsequence(snp.tbl = atsnp.scores$snp.tbl,
                                 motif.scores = atsnp.result,
                                 motif.lib = motif_library,
                                 snpids = c("rs10910078", "rs4486391"),
                                 motifs = names(motif_library)[1:2],
                                 ncores = 2)
match_result[, list(snpid, motif, IUPAC, ref_match_seq, snp_match_seq)]


## ----include=TRUE,eval=TRUE, echo=TRUE,fig.align="center",dpi=600,fig.width=6,fig.height=6----
plotMotifMatch(snp.tbl = atsnp.scores$snp.tbl,
               motif.scores = atsnp.scores$motif.scores,
               snpid = atsnp.scores$snp.tbl$snpid[1],
               motif.lib = motif_library,
               motif = atsnp.scores$motif.scores$motif[1])
## ----eval=TRUE,echo=FALSE,results="markup",cache=FALSE----------------------------------

# ----------------------------------- additional analysis ---------------------------------------------- #

# get the motifs that have minimum pvalue for each rsid
library(plyr)
atsnp.result.sub = as.data.frame(atsnp.result.sub)
tt = ddply(.data = atsnp.result.sub,.variables =  .(snpid), summarise, x = min(pval_rank))
tt = merge(tt, atsnp.result.sub, by.x=c('snpid','x'), by.y=c('snpid','pval_rank'))
tt = tt[order(tt$snpid,tt$x),]
tt.count = count_(x = tt,vars = 'motif',sort = TRUE)

# compute ranks for atsnp.result
ranks = ddply(atsnp.result,.(motif),transform,Order = rank(pval_rank))
ranks.mat <- dcast(data = ranks, formula = motif~snpid, value.var = 'Order')
rownames(ranks.mat) <- ranks.mat$motif
ranks.mat <- ranks.mat[,-1]
ranks.mat$mean = rowMeans(ranks.mat)
ranks.mat$min <- apply(ranks.mat,1,FUN = min)
ranks.mat$max <- apply(ranks.mat,1,FUN = max)
ranks.mat$sum <- apply(ranks.mat,1,FUN = sum)
# motif.ranks <- ranks.mat[,5270:5273]
