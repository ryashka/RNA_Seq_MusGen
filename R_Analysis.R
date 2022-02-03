#Make column names legible
convertInputColNames = function(fc){
  require(stringr)
  colnames(fc)[6:ncol(fc)] = str_split_fixed(colnames(fc),"\\.",Inf)[,ncol(str_split_fixed(colnames(fc),"\\.",Inf))-1][6:ncol(fc)]
  cols = as.data.frame(str_split_fixed(colnames(fc),"_",Inf)[,1:2][6:ncol(fc),])
  cols$namesOfCols = paste(cols$V1, cols$V2, sep="_")
  colnames(fc)[6:ncol(fc)] = cols$namesOfCols
  return(fc)
}
#Make row names legible
convertInputRowNames = function(fc){
  source("~/Documents/Programming/programs/ScriptsAndFunctions/UsefulFunctions.R")
  return(add_MGI(fc))
}
makeSampleTable = function(fc, condition1cols){
  justCounts = fc[,6:29]
  samples = colnames(justCounts)
  condition1 = 
}
#To generate TPMs of B cells for original Mus81 and Gen1
analyze_GSE72018 = function(){
  bm = read.delim("~/Downloads/GSE72018_Bcell_featureCounts.tsv")
  bm = addLengths(bm,1)
  row.names(bm) = bm[,1]
  bm = add_MGI(bm)
  bm = bm[!is.na(bm$Length),]
  cts = bm[,2:(ncol(bm)-2)]
  tpm_matrix = countsToTPM(cts,bm)
  head(tpm_matrix)
  bm[,2:(ncol(bm)-2)] = tpm_matrix
  #write.csv(bm,"~/Documents/Programming/Keith_seq/analyzed_data/GSE72018_tpm.csv")
}
#This script does the complete analysis

scriptForAnalysis = function(){
  #Import relevant functions
  source("~/Documents/Programming/programs/ScriptsAndFunctions/UsefulFunctions.R")
  fc = read.delim("~/Documents/Programming/Keith_seq/processed_data/output_files/featureCounts_files/combinedFeatureCounts.txt", row.names = 1, skip = 1)
  fc = convertInputColNames(fc)
  WTcols = 24:29
  GKOcols = 12:17
  MKOcols = 18:23
  DKOcols = 6:11
  colsM = cbind(WTcols,GKOcols, MKOcols, DKOcols)
  colnames(colsM) = c("WT", "GKO", "MKO", "DKO")
  genesOfInterest = c("Aicda", "Mus81", "Gen1")
  writeOutputTo = "~/Documents/Programming/Keith_seq/analyzed_data/210902/"
  #Generate all comparisons of pairs of genotypes (e.g. DKO vs WT)
  for (i in 1:4){
    for (j in 1:4){
      comparisonName = paste0(colnames(colsM)[i],"vs",colnames(colsM)[j])
      deseq_ij = deseq_and_volcano(fc, colsM[,i], colsM[,j], writeOutputTo, comparisonName, genesOfInterest)
      #go_and_kegg(passed_deseq = deseq_ij, output_directory = writeOutputTo, fileNames = comparisonName)
      #double_go(passed_deseq = deseq_ij, output_directory = writeOutputTo, fileNames = comparisonName)
    }
  }
  #GSEA and GO analyses
  fc_tpm_2 = fc
  fc_tpm_2$ENSEMBL = str_split_fixed(row.names(fc_tpm_2),'\\.',Inf)[,1]
  human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset="mmusculus_gene_ensembl")
  mmus_ids <- getBM(attributes=c('ensembl_gene_id','external_gene_name'),filters = 'ensembl_gene_id',
                    values = fc_tpm_2$ENSEMBL, mart = mouse)
  library(tibble)
  library(dplyr)
  library(tidyr)
  mmus_ids <- as_tibble(mmus_ids)
  colnames(mmus_ids)[1] = "ENSEMBL"
  fc_tpm_2 = inner_join(fc_tpm_2, mmus_ids)
  fc_tpm_2 = as.data.frame(fc_tpm_2)
  fc_tpm_2[,8:ncol(fc_tpm_2)-2] = countsToTPM(fc_tpm_2[,8:ncol(fc_tpm_2)-2],fc_tpm_2)
  enriched.pathway = as.data.frame(cbind(c("HALLMARK_P53_PATHWAY",
                   "HALLMARK_G2M_CHECKPOINT",
                   "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                   "HALLMARK_APOPTOSIS",
                   "HALLMARK_DNA_REPAIR",
                   "HALLMARK_MYC_TARGETS_V1", 
                   "HALLMARK_MYC_TARGETS_V2",
                   "HALLMARK_E2F_TARGETS"),
                   c(TRUE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE)))
  colnames(enriched.pathway) = c("Pathway","Enriched")
  head(fc_tpm_2)
  library(RColorBrewer)
  library(ComplexHeatmap)
  row.names(fc_tpm_2) = make.unique(fc_tpm_2$external_gene_name)
  #Generate unique GSEA plots
  for (i in 1:8){
    pathway_diff_genes = pathways.hallmark %>% 
      enframe("pathway", "SYMBOL") %>% 
      unnest() %>% 
      inner_join(res, by="SYMBOL") %>%
      filter(pathway == enriched.pathway[i,1])
    pathway_diff_genes = pathway_diff_genes[!duplicated(pathway_diff_genes$ENSEMBL),]
    pathway_diff_genes = pathway_diff_genes[!is.na(pathway_diff_genes$padj) & pathway_diff_genes$padj < 0.05,]
    if(enriched.pathway[i,2]){
      if(i == 1){
        pathway_diff_genes = pathway_diff_genes[pathway_diff_genes$log2FoldChange > 1,]
      } else {
        pathway_diff_genes = pathway_diff_genes[pathway_diff_genes$log2FoldChange > 0.5,]
      }
      coul <- colorRampPalette(brewer.pal(8, "YlOrRd"))(25)
    } else {
      pathway_diff_genes = pathway_diff_genes[pathway_diff_genes$log2FoldChange < -0.2,]
      coul <- colorRampPalette(rev(brewer.pal(8, "YlGnBu")))(25)
    }
    mat = t(scale(t(as.matrix(fc_tpm_2[c(fc_tpm_2$ENSEMBL %in% pathway_diff_genes$ENSEMBL),cols]))))
    pdf(paste0(writeOutputTo,enriched.pathway[i,1],"_heatmap.pdf"),
               width = 5, height = 5)
    #heatmap(as.matrix(fc_tpm_2[c(fc_tpm_2$ENSEMBL %in% pathway_diff_genes$ENSEMBL),cols]),
    #        Rowv = NA, Colv = NA, scale = "row", col = coul)
    print(Heatmap(mat, col = coul, cluster_rows = FALSE, cluster_columns = FALSE, show_column_dend = FALSE,
            column_title = enriched.pathway[i,1],
            heatmap_legend_param = list(
              title = "Z_Score",
              col = coul, 
              legend_height = unit(6, "cm")),))
    dev.off()
  }
  ensembl.Human = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id",
                         values = deseq_output$ENSEMBL, mart = mouse,
                         attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=T)
  colnames(ensembl.Human) = c("ENSEMBL","ensembl_gene_id")
  deseq_for_fgsea = inner_join(ensembl.Human,as.data.frame(deseq_output))[,c("ensembl_gene_id","stat")]
  deseq_ranks = deseq_for_fgsea$stat
  names(deseq_ranks) = deseq_for_fgsea$ensembl_gene_id
  #Adaptation of fgsea for plots
  for (i in 1:8){
    enriched_hallmark = as.data.frame(pathways.hallmark[[enriched.pathway[i,1]]])
    colnames(enriched_hallmark)[1] = "hgnc_symbol"
    enriched_human = getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),filters = 'hgnc_symbol',
                      values = enriched_hallmark$hgnc_symbol, mart = human)
    pathway = enriched_human$ensembl_gene_id
    stats = deseq_ranks
    gseaParam = 1
    ticksSize = 0.2
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                            returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    diff <- (max(tops) - min(bottoms))/8
    x = y = NULL
    if(enriched.pathway[i,2]){
      g = ggplot(toPlot, aes(x = x, y = y)) +
        geom_point(color = "green", size = 0.1) +
        geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") +
        geom_hline(yintercept = 0, colour = "black") + geom_line(color = "green") + theme_bw() +
        geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2),
                     size = ticksSize) +
        theme(panel.border = element_blank(), panel.grid.minor = element_blank()) +
        theme(plot.title = element_text(hjust = 0.5)) +
        labs(x = "Rank in Ordered Gene List", y = "Enrichment Score") +
        labs(title=enriched.pathway[i,1])
    } else {
      g = ggplot(toPlot, aes(x = x, y = y)) +
        geom_point(color = "green", size = 0.1) +
        geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") +
        geom_hline(yintercept = 0, colour = "black") + geom_line(color = "green") + theme_bw() +
        geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2),
                     size = ticksSize) +
        theme(panel.border = element_blank(), panel.grid.minor = element_blank()) +
        theme(plot.title = element_text(hjust = 0.5)) +
        labs(x = "Rank in Ordered Gene List", y = "Enrichment Score") +
        labs(title=enriched.pathway[i,1])
    }
    
    pdf(paste0(writeOutputTo,enriched.pathway[i,1],"_fgsea.pdf"),
        width = 6, height = 3)
    print(g)
    dev.off()
  }
  }
#Original sequencing data from Weifeng
#Generate TPM files for LI-stim cells
wf_tpm = function(){
  require("biomaRt")
  wf_data = read.delim("~/Downloads/GSE128321_RawCts_byTx_CHD4KO.txt")
  library(plyr)
  wf_gene = ddply(wf_data,"gene.name",numcolwise(sum))
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  colnames(wf_gene)[1] = "mgi_symbol"
  annotLookup <- getBM(
    mart = mouse,
    attributes = c(
      'mgi_symbol',
      'ensembl_gene_id'),
    uniqueRows = TRUE)
  wf_gene_ens = join(wf_gene,annotLookup)
  wf_gene_ens = wf_gene_ens[!is.na(wf_gene_ens$ensembl_gene_id),]
  
  wf_gene_ens = addLengths(wf_gene_ens,14)
  wf_gene_ens = wf_gene_ens[!is.na(wf_gene_ens$Length),]
  tpm_matrix = countsToTPM(wf_gene_ens[,3:(ncol(wf_gene_ens)-2)],wf_gene_ens)
  head(tpm_matrix)
  wf_gene_ens[,3:(ncol(wf_gene_ens)-2)] = tpm_matrix
  write.csv(wf_gene_ens,"~/Documents/Programming/Keith_seq/analyzed_data/GSE128321.csv")
}
