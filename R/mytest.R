mytest <- function(hgu95av2.db, hgu95av2ENTREZID, GeneAnswers) {
  library(hgu95av2.db)
  sample_probes <-sample(ls(hgu95av2ENTREZID),40)
  sample_genes <-sapply(sample_probes,function(x)
    hgu95av2ENTREZID[[x]])
  sim<-mgeneSim(sample_genes,ont="MF",organism="human",measure="Wang")
  d <- semData('org.Hs.eg.db', ont="MF")
  require(GeneAnswers)
  
  data('humanExpr')
  data('humanGeneInput')
  y <- geneAnswersBuilder(humanGeneInput, 'org.Hs.eg.db', 
                          categoryType='KEGG', testType='hyperG', 
                          pvalueT=0.1, geneExpressionProfile=humanExpr, 
                          verbose=FALSE)
  yy <- y@enrichmentInfo
  
  require(clusterProfiler)
  x <- enrichKEGG(humanGeneInput$GeneID, pvalueCutoff=0.2, 
                  qvalueCutoff=0.2, minGSSize=1)
  xx <- summary(x)
  
  id <- sub("hsa", "", xx$ID)
  idx <- id %in% rownames(yy)
  
  p.clusterProfiler <- xx$pvalue[idx]
  p.GeneAnswers <- yy[id[idx],]$"p value"
  
  cor(p.clusterProfiler, p.GeneAnswers)
  
  p.clusterProfiler-p.GeneAnswers
  
  require(DOSE)
  require(clusterProfiler)
  data(geneList)
  deg <- names(geneList)[abs(geneList)>2]
  
  gda <- read.delim("/media/H_driver/Aimin_project/DO_data/all_gene_disease_associations.tsv")
  dim(gda)
  head(gda)
  
  disease2gene=gda[, c("diseaseId", "geneId")]
  disease2name=gda[, c("diseaseId", "diseaseName")]
  
  x = enricher(deg, TERM2GENE=disease2gene, TERM2NAME=disease2name)
  head(summary(x))
  barplot(x)
  
  y = GSEA(geneList, TERM2GENE=disease2gene, TERM2NAME=disease2name)
  head(y)
  gseaplot(y, "umls:C0003872")
  
  library(AnnotationHub)
  hub <- AnnotationHub()
  
  query(hub, "Cricetulus")
  query(hub, "human")
  
  Cgriseus <- hub[["AH48061"]]
  sample_gene <- sample(keys(Cgriseus), 100)
  str(sample_gene)
  
  library(clusterProfiler)
  
  sample_test <- enrichGO(sample_gene, OrgDb=Cgriseus, pvalueCutoff=1, qvalueCutoff=1)
  head(summary(sample_test))
  
  
  library(org.Hs.eg.db)
  data(geneList)
  gene <- names(geneList)[abs(geneList) > 2]
  gene.df <- bitr(gene, fromType = "ENTREZID", 
                  toType = c("ENSEMBL", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
  
  head(gene.df)
  
  sink("test_ok_HT.txt")
  
  ego.HT <- enrichGO(gene          = gene,
                  universe      = rownames(Example.Go.adjusted.by.exon[[2]]),
                  keytype = "ENSEMBL",
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 1,
                  qvalueCutoff  = 0.8,GOFromGOSeq=Example.Go.adjusted.by.exon,whichway = "HT")
  
  
  
  dim(ego.HT@result)
  
  sink()
  enrichMap(ego.HT, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
  
  sink("test_ok_WHT.txt")
  
  ego.WHT <- enrichGO(gene          = gene,
                     universe      = rownames(Example.Go.adjusted.by.exon[[2]]),
                     keytype = "ENSEMBL",
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 1,
                     qvalueCutoff  = 0.8,GOFromGOSeq=Example.Go.adjusted.by.exon,whichway = "WHT")
  
  
  
  dim(ego.WHT@result)
 
  
  sink()
  enrichMap(ego.WHT, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
  
  
  
  
  
  
  enrichMap(ego, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
  cnetplot(ego)
  plotGOgraph(ego)
  
  save.image(file="test_out.RData")
  
  savehistory(file="test_out.Rhistory")
  
  qExtID2TermID.df <- data.frame(extID=rep(names(ego$qExtID2TermID),times=lapply(ego$qExtID2TermID, length)),termID=ego$qTermID)
  
  qExtID2TermID.df <- unique(qExtID2TermID.df)
  
  qTermID2ExtID <- with(qExtID2TermID.df,split(as.character(extID), as.character(termID)))
  
  head(summary(ego))
  
  ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                   OrgDb         = org.Hs.eg.db,
                   keytype       = 'ENSEMBL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
  head(summary(ego2))
  
  ego3 <- enrichGO(gene         = gene.df$SYMBOL,
                   OrgDb         = org.Hs.eg.db,
                   keytype       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
  head(summary(ego3))
  
  ego <- setReadable(ego, OrgDb = org.Hs.eg.db)
  ego2 <- setReadable(ego2, OrgDb = org.Hs.eg.db)
  head(summary(ego), n=3)
  
  
  head(summary(ego2), n=3)
  
  
  
  gsecc <- gseGO(geneList=geneList, ont="CC", OrgDb=org.Hs.eg.db, verbose=F)
  head(summary(gsecc))
  
  gseaplot(gsecc, geneSetID="GO:0000779")
  
  gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
  c5 <- read.gmt(gmtfile)
  egmt <- enricher(gene, TERM2GENE=c5)
  head(summary(egmt))
  
  egmt <- setReadable(egmt, OrgDb=org.Hs.eg.db, keytype="ENTREZID")
  head(summary(egmt))
  
  gsegmt <- GSEA(geneList, TERM2GENE=c5, verbose=F)
  head(summary(gsegmt))
  
  enrichMap(gsegmt, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
}
