convertRowsToMGI = function(cts){
  require(stringr)
  G_list = read.csv("~/Documents/Programming/programs/ScriptsAndFunctions/Ensembl_To_Mgi_Symbol.csv",row.names =  1)
  genesInOrder <- as.data.frame(row.names(cts))
  genesInOrder$order = 1:nrow(genesInOrder)
  genesInOrder$EnsemblID = str_split_fixed(genesInOrder$`row.names(cts)`,'\\.',2)[,1]
  genesInOrder = merge(genesInOrder,G_list[!duplicated(G_list$ensembl_gene_id),],
                       by.x = "EnsemblID",by.y="ensembl_gene_id",all.x=TRUE)
  genesInOrder = genesInOrder[order(genesInOrder$order),]
  0 %in% match(genesInOrder$`row.names(cts)`,row.names(cts),nomatch=0)
  genesInOrder_2 = genesInOrder
  genesInOrder_2$EnsemblID = as.character(genesInOrder_2$EnsemblID)
  genesInOrder_2$mgi_symbol_character = as.character(genesInOrder_2$external_gene_name)
  genesInOrder_2$mgi_symbol_character = ifelse(genesInOrder_2$mgi_symbol_character == "",
                                               genesInOrder_2$ensembl_id,genesInOrder_2$mgi_symbol_character)
  genesInOrder_2$mgi_symbol_character = ifelse(is.na(genesInOrder_2$mgi_symbol_character),
                                               genesInOrder_2$EnsemblID,genesInOrder_2$mgi_symbol_character)
  genesInOrder_3 = genesInOrder_2[,c("order","mgi_symbol_character","EnsemblID")]
  rownames(cts) = make.unique(genesInOrder_2$mgi_symbol_character)
  cts$ENSEMBL = genesInOrder_3$EnsemblID
  return(cts)
}
add_Entrez = function(df, drop_empty_rows){
  df$MGI = row.names(df)
  entrez_names = bitr(df$ENSEMBL, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Mm.eg.db")
  if(!drop_empty_rows){
    df = merge(as.data.frame(df),entrez_names[!duplicated(df$ENSEMBL),],
          by.x = "ENSEMBL",by.y="ENSEMBL", all.x=TRUE)
  } else {
    df = merge(as.data.frame(df),entrez_names[!duplicated(df$ENSEMBL),],
               by.x = "ENSEMBL",by.y="ENSEMBL")
  }
  return(df)
}
deseq_cts = function(fc, groupACols, groupBCols, saveResOrderedTo, nameOutputCSV){
  #FC matrix should have row names in ENSEMBL format
  countMatrix = fc[,c(groupACols,groupBCols)]
  groupA = rep("GroupA",length(groupACols))
  groupB = rep("GroupB",length(groupBCols))
  sampleTypes = c(groupA,groupB)
  
  library(reshape2)
  sampleTable = melt(data.frame(colnames(fc[,c(groupACols,groupBCols)]),sampleTypes))
  colnames(sampleTable) <- c("sampleName","sampleType")
  sampleTable$sampleName = factor(sampleTable$sampleName)
  sampleTable$sampleType = factor(sampleTable$sampleType)
  library("DESeq2")
  rownames(sampleTable) = sampleTable$sampleName
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                                colData = sampleTable,
                                design = ~ sampleType)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  dds = DESeq(dds)
  res = results(dds)
  resOrdered <- convertRowsToMGI(res[order(res$padj), ])
  write.csv(resOrdered, paste0(saveResOrderedTo,nameOutputCSV,".csv"))
  return(resOrdered)
}
deseq_volcano <- function(deseq_resOrdered, saveVolcanoTo, nameOutputVolcano, genesOfInterest){
  df = as.data.frame(deseq_resOrdered)
  require(dplyr)
  require(ggplot2)
  require(ggrepel)
  top_genes = rownames(slice_min(df[df$log2FoldChange > 0,], n=10, order_by=padj ))
  bot_genes = rownames(slice_min(df[df$log2FoldChange < 0,], n=10, order_by=padj ))
  
  df$shaped_genes = ifelse(rownames(df) %in% top_genes,"Top_10","Other_gene")
  df$shaped_genes = ifelse(rownames(df) %in% bot_genes,"Bot_10",df$shaped_genes)
  
  df$log2FCOffAxisPositive = ifelse(df$log2FoldChange > 2,TRUE,FALSE)
  df$log2FCOffAxisNegative = ifelse(df$log2FoldChange < -2,TRUE,FALSE)
  df$log2FCOffAxis = ifelse(df$log2FCOffAxisPositive | df$log2FCOffAxisPositive,TRUE,FALSE)
  df$log2FCWithOffAxisValues = ifelse(df$log2FCOffAxisPositive,
                                      2,
                                      df$log2FoldChange)
  df$log2FCWithOffAxisValues = ifelse(df$log2FCOffAxisNegative,
                                      -2,
                                      df$log2FCWithOffAxisValues)
  df$padjOffAxis = ifelse(df$padj > 0.1,TRUE,FALSE)
  df$padjWithOffAxisValues = ifelse(df$padjOffAxis,
                                      0.1,
                                      df$padj)
  df$OffAxis = ifelse(df$log2FCOffAxis | df$padjOffAxis,TRUE,FALSE)
  df$shaped_genes = ifelse(df$OffAxis, "Off_Axis", df$shaped_genes)
  
  p53_genes = read.csv("~/Documents/Programming/reference_files/p53_Gene_List.csv")
  g2m_genes = read.csv("~/Documents/Programming/reference_files/G2M_Gene_List.csv")
  ifna_genes = read.csv("~/Documents/Programming/reference_files/Ifna_Gene_List.csv")
  apoptosis_genes = read.csv("~/Documents/Programming/reference_files/Apoptosis_Gene_List.csv")
  DNArepair_genes = read.csv("~/Documents/Programming/reference_files/DNA_repair_Gene_List.csv")
  
  # df$colored_genes = ifelse(df$ENSEMBL %in% ifna_genes$ENSEMBL, "Ifna", "pink")
  # df$colored_genes = ifelse(df$ENSEMBL %in% g2m_genes$ENSEMBL, "G2M", df$colored_genes)
  # df$colored_genes = ifelse(df$ENSEMBL %in% p53_genes$ENSEMBL, "p53", df$colored_genes)
  # df$colored_genes = ifelse(df$padj > 0.05, "pink", df$colored_genes)
  # df$colored_genes = ifelse(df$log2FoldChange < 1, "pink", df$colored_genes)
  # df$colored_genes = ifelse(rownames(df) %in% genesOfInterest, "Keith_Gene", df$colored_genes)

  df$colored_genes = ifelse(rownames(df) %in% top_genes, "top_gene", "pink")
  df$colored_genes = ifelse(rownames(df) %in% bot_genes, "bot_gene", df$colored_genes)
  df$colored_genes = ifelse(rownames(df) %in% genesOfInterest, "Keith_Gene", df$colored_genes)
  
  #p53_top_genes = rownames(slice_min(df[df$padj < 0.05 & df$ENSEMBL %in% p53_genes$ENSEMBL,], n=10, order_by=-log2FoldChange ))
  #g2m_top_genes = rownames(slice_min(df[df$padj < 0.05 & df$ENSEMBL %in% g2m_genes$ENSEMBL,], n=5, order_by=log2FoldChange ))
  #ifna_top_genes = rownames(slice_min(df[df$padj < 0.05 & df$ENSEMBL %in% ifna_genes$ENSEMBL,], n=8, order_by=-log2FoldChange ))
   

   df$labeled_genes = ifelse(rownames(df) %in% genesOfInterest,"Interesting_Gene", "Other_Gene")
   df$labeled_genes = ifelse(rownames(df) %in% top_genes,"Top10_Gene", df$labeled_genes)
   df$labeled_genes = ifelse(rownames(df) %in% bot_genes,"Bot10_Gene", df$labeled_genes)
  # df$labeled_genes = ifelse(rownames(df) %in% mycv1_top_genes,"Mycv1_Gene", df$labeled_genes)
  # df$labeled_genes = ifelse(rownames(df) %in% mycv2_top_genes,"Mycv2_Gene", df$labeled_genes)
  # df$labeled_genes = ifelse(rownames(df) %in% ifna_top_genes,"Ifna_Gene", df$labeled_genes)
  # df$labeled_genes = ifelse(rownames(df) %in% g2m_top_genes,"G2M_Gene", df$labeled_genes)
  # df$labeled_genes = ifelse(rownames(df) %in% p53_top_genes,"p53_Gene", df$labeled_genes)
  
  gg_volcano = ggplot(df, aes(x = log2FCWithOffAxisValues, y = -log10(padjWithOffAxisValues), 
                 color = colored_genes, shape = shaped_genes)) + 
    geom_point() + theme_bw() + 
    scale_y_continuous(trans = "log10") +
    xlim(-2,2) +
    #coord_trans(y="log10") +
    #geom_vline(xintercept = -2, linetype="dotted", color = "grey28", size = 1.5) +
    #geom_vline(xintercept = 2, linetype="dotted", color = "grey28", size = 1.5) +
    geom_hline(yintercept = -log10(0.05), linetype="dotted", color = "grey28", size = 1.5) +
    scale_color_manual(values=c("forestgreen","red","grey50","purple")) +
    scale_shape_manual(values=c(19, 17, 19, 19)) +
    geom_text_repel(aes(label = ifelse(as.character(`labeled_genes`) != "Other_Gene",
                                       as.character(rownames(df)),'')),
                    nudge_x = 0.2, nudge_y = 0.1,
                    force = 20, force_pull = 3,
                    segment.size = 0.1, segment.alpha = 1, 
                    point.padding = NA, box.padding = 1,
                    size = 4, max.overlaps = Inf)
  ggsave(file = paste0(saveVolcanoTo,nameOutputVolcano,".pdf"), 
         plot = gg_volcano,
         device = "pdf")
}
deseq_and_volcano = function(fc, groupACols, groupBCols, saveTo, nameBoth, genesOfInterest){
  deseq_resOrdered = deseq_cts(fc, groupACols, groupBCols, saveTo, nameBoth)
  deseq_volcano(deseq_resOrdered, saveTo, nameBoth, genesOfInterest)
  return(deseq_resOrdered)
}

go_and_kegg = function(input_directory = "~/Documents/Programming/", passed_deseq = NULL, output_directory, fileNames){
  if(is.null(passed_deseq)){
    passed_deseq = read.csv(paste0(input_directory,"/deseq_output.csv"), row.names = 1)
  }
  double_go(passed_deseq = passed_deseq, output_directory = output_directory, fileNames = fileNames)
  double_kegg(passed_deseq = passed_deseq, output_directory = output_directory, fileNames = fileNames)
}
double_go = function(input_directory = "~/Documents/Programming/", passed_deseq = NULL, output_directory, fileNames){
  if(is.null(passed_deseq)){
    passed_deseq = read.csv(paste0(input_directory,"/deseq_output.csv"), row.names = 1)
  }
  if(passed_deseq$padj[1] < 0.05){
    pd = passed_deseq[!is.na(passed_deseq$padj),]
    significant_deseq = pd[pd$padj < 0.05,]
    comparisonName = paste0(fileNames,"_p<0.05")
    go_cp(passed_deseq = significant_deseq, output_directory = output_directory, fileNames = comparisonName)
  }
  go_cp(passed_deseq = passed_deseq, output_directory = output_directory, fileNames = fileNames)
}
go_cp = function(input_directory = "~/Documents/Programming/", passed_deseq = NULL, output_directory, fileNames){
  require(clusterProfiler)
  require(enrichplot)
  require(ggplot2)
  require(icesTAF)
  organism = "org.Mm.eg.db"
  if(is.null(passed_deseq)){
    passed_deseq = read.csv(paste0(input_directory,"/deseq_output.csv"), row.names = 1)
  }
  library(organism, character.only = TRUE)
  original_gene_list <- passed_deseq$log2FoldChange
  names(original_gene_list) <- passed_deseq$ENSEMBL
  gene_list = na.omit(original_gene_list)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  #Maybe go back here and switch gseaMultilevel by removing nperm
  gse = gseGO(geneList=gene_list, 
              ont ="ALL", 
              keyType = "ENSEMBL", 
              minGSSize = 3,
              maxGSSize = 800,
              nPermSimple = 10000,
              pvalueCutoff = 0.05, 
              verbose = TRUE, 
              OrgDb = org.Mm.eg.db
              )
  
  require(DOSE)
  require(ggnewscale)
  mkdir(paste0(output_directory,"/GO_plots/"))
  mkdir(paste0(output_directory,"/GSEA_plots/"))
  go_directory = paste0(output_directory,"/GO_plots/")
  gsea_directory = paste0(output_directory,"/GSEA_plots/")
  go_dotplot = dotplot(gse, showCategory=10, split=".sign", x = "Count") + facet_grid(.~.sign)
  ggsave(file = paste0(go_directory,fileNames,"_go_dotplot.pdf"),
         plot = go_dotplot,
         device = "pdf")
  #emap = emapplot(gse, showCategory = 10)
  go_cnet = cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
  ggsave(file = paste0(go_directory,fileNames,"_go_cnetplot.pdf"),
         plot = go_cnet,
         device = "pdf")
  go_ridge = ridgeplot(gse) + labs(x = "enrichment distribution")
  ggsave(file = paste0(go_directory,fileNames,"_go_ridgeplot.pdf"),
         plot = go_ridge,
         device = "pdf")
  for(i in 1:10){
    single_gsea = gseaplot(gse, by = "all", title = gse$Description[i], geneSetID = i)
    ggsave(file = paste0(gsea_directory,fileNames,"_gsea_set_",i),
           plot = single_gsea,
           device = "pdf")
  }
}
double_kegg = function(input_directory = "~/Documents/Programming/", passed_deseq = NULL, output_directory, fileNames){
  if(is.null(passed_deseq)){
    passed_deseq = read.csv(paste0(input_directory,"/deseq_output.csv"), row.names = 1)
  }
  if(passed_deseq$padj[1] < 0.05){
    pd = passed_deseq[!is.na(passed_deseq$padj),]
    significant_deseq = pd[pd$padj < 0.05,]
    comparisonName = paste0(fileNames,"_p<0.05")
    kegg_cp(passed_deseq = significant_deseq, output_directory = output_directory, fileNames = comparisonName)
  }
  kegg_cp(passed_deseq = passed_deseq, output_directory = output_directory, fileNames = fileNames)
}

kegg_cp = function(input_directory = "~/Documents/Programming/", passed_deseq = NULL, output_directory, fileNames){
  require(clusterProfiler)
  require(enrichplot)
  require(ggplot2)
  require(icesTAF)
  organism = "org.Mm.eg.db"
  if(is.null(passed_deseq)){
    passed_deseq = read.csv(paste0(input_directory,"/deseq_output.csv"), row.names = 1)
  }
  
  #library(organism, character.only = TRUE)
  original_gene_list <- passed_deseq$log2FoldChange
  ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Mm.eg.db") # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
  dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
  
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  df2 = passed_deseq[passed_deseq$X %in% dedup_ids$ENSEMBL,]
  
  # Create a new column in df2 with the corresponding ENTREZ IDs
  df2$Y = dedup_ids$ENTREZID
  
  # Create a vector of the gene unuiverse
  kegg_gene_list <- df2$log2FoldChange
  
  # Name vector with ENTREZ ids
  names(kegg_gene_list) <- df2$Y
  
  # omit any NA values 
  kegg_gene_list<-na.omit(kegg_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
  
  # From significant results, we want to filter on log2fold change
  kegg_genes <- kegg_sig_genes_df$log2FoldChange
  
  # Name the vector with the CONVERTED ID!
  names(kegg_genes) <- kegg_sig_genes_df$Y
  
  # omit NA values
  kegg_genes <- na.omit(kegg_genes)
  
  # filter on log2fold change (PARAMETER)
  kegg_genes <- names(kegg_genes)[abs(kegg_genes) > 2]
  
  #Maybe go back here and switch gseaMultilevel by removing nperm
  gse = gseKEGG(geneList=kegg_genes, 
              ont ="ALL", 
              keyType = "ncbi-geneid", 
              minGSSize = 3,
              maxGSSize = 800,
              nPermSimple = 10000,
              pvalueCutoff = 0.05, 
              verbose = TRUE, 
              OrgDb = org.Mm.eg.db
  )
  
  require(DOSE)
  require(ggnewscale)
  mkdir(paste0(output_directory,"/KEGG_plots/"))
  mkdir(paste0(output_directory,"KEGG/GSEA_plots/"))
  go_directory = paste0(output_directory,"/KEGG_plots/")
  gsea_directory = paste0(output_directory,"KEGG_plots/GSEA_plots/")
  go_dotplot = dotplot(gse, showCategory=10, split=".sign", x = "Count") + facet_grid(.~.sign)
  ggsave(file = paste0(go_directory,fileNames,"_kegg_dotplot.pdf"),
         plot = go_dotplot,
         device = "pdf")
  #emap = emapplot(gse, showCategory = 10)
  go_cnet = cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
  ggsave(file = paste0(go_directory,fileNames,"_kegg_cnetplot.pdf"),
         plot = go_cnet,
         device = "pdf")
  go_ridge = ridgeplot(gse) + labs(x = "enrichment distribution")
  ggsave(file = paste0(go_directory,fileNames,"_kegg_ridgeplot.pdf"),
         plot = go_ridge,
         device = "pdf")
  for(i in 1:10){
    single_gsea = gseaplot(gse, by = "all", title = gse$Description[i], geneSetID = i)
    ggsave(file = paste0(gsea_directory,fileNames,"_gsea_set_",i),
           plot = single_gsea,
           device = "pdf")
  }
}
makeSampleTable = function(fc){
  justCounts = fc[,6:ncol(fc)]
  sampleTable = as.data.frame(colnames(justCounts))
  colnames(sampleTable) = "Sample"
  sampleTable$sampleType = as.factor(str_split_fixed(sampleTable$Sample,"\\_",Inf)[,1])
  return(sampleTable)
}
makePCAplot = function(fc){
  library("DESeq2")
  library("PCAtools")
  sampleTable = makeSampleTable(fc)
  fc = convertRowsToMGI(fc)
  rownames(sampleTable) = sampleTable$Sample
  sampleTable$Gen1 = c(rep("GenKO",12),rep("GenW",12))
  sampleTable$Mus81 = c(rep("MusKO",6),rep("MusW",6),rep("MusKO",6),rep("MusW",6))
  testTable = sampleTable[1:12,]
  ddsfm <- DESeqDataSetFromMatrix(countData = fc[,7:ncol(fc)-1],
                                colData = sampleTable,
                                design = ~ sampleType)
  keep <- rowSums(counts(ddsfm)) >= 10
  ddsfm <- ddsfm[keep,]
  dds = DESeq(ddsfm)
  vst <- assay(vst(dds))
  p <- pca(vst, metadata = testTable, removeVar = 0.1)
  screeplot(p, axisLabSize = 18, titleLabSize = 22)
  biplot(p)
  biplot(p, showLoadings = TRUE, labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
  rld <- rlog(dds, blind=TRUE)
  se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                             colData=colData(dds))
  plotPCA( DESeqTransform( se ) , intgroup=c("sampleType"), ntop = 400)
  library(pcaExplorer)
  pcaExplorer(dds = ddsfm)
  ddsfm_pca <- prcomp(as.data.frame(t(fc[,7:ncol(fc)-1])))
  df_out <- as.data.frame(ddsfm_pca$x)
  df_out$group <- str_split_fixed(as.character(colnames(ddsfm)), "_", Inf)[,1]
  head(df_out)
  df_out_r <- as.data.frame(ddsfm_pca$rotation)
  df_out_r$feature <- row.names(df_out_r)
  
  
  
  p<-ggplot(df_out_r,aes(x=PC1,y=PC3,color=feature ))
  p<-p+geom_point()
  p
  
  plot(ddsfm_pca$x[,1], ddsfm_pca$x[,2])
  head(iris)
  df <- iris
  df <- as.data.frame(iris)
  row.names(df) <- paste(df$Species, row.names(df), sep="_") 
  df$Species <- NULL
  df_pca <- prcomp(df)
  plot(df_pca$x[,1], df_pca$x[,2])
  
  head(df)
}
countsToTPM <- function(cts,countMatrix) {
  rpk = (cts/countMatrix$Length)
  tpm_cts = sweep(rpk,2,colSums(rpk),`/`)*1e6
  return(tpm_cts)
}
goana_biased = function(deseq_output){
  library(edgeR)
  library(kableExtra)
  with_Entrez = add_Entrez(deseq_output, FALSE)
  sig_deseq = with_Entrez[!is.na(with_Entrez$pvalue) & with_Entrez$pvalue < 0.05,]
  dedup_sig = sig_deseq[!duplicated(sig_deseq[c("ENTREZID")]),]
  dedup_sig = dedup_sig[!is.na(dedup_sig$ENTREZID),]
  row.names(dedup_sig) = dedup_sig$ENTREZID
  GO = goana(dedup_sig$ENTREZID,coef = "padj",species = "Mm", trend = TRUE)
  topGO(GO, ontology = "BP", number = 10) %>%
    kable(caption = "Top enriched or depleted GO-terms.") %>%
    kable_styling(latex_options = "hold_position")
}
convert_mouse_human = function(gene_set, human_to_mouse, ensembl_to_ensembl){
  require("biomaRt")
  if(ensembl_to_ensembl){
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    if(human_to_mouse){
      ensembl.Mouse = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = gene_set$ensembl.ID , mart = human, attributesL = c("ensembl_gene_id"), martL = mouse, uniqueRows=T)
      return(ensembl.Mouse)
      } else {
        ensembl.Human = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = gene_set$ENSEMBL , mart = mouse, attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=T)
        return(ensembl.Human)
      }
  }
}
fgsea_analysis = function(deseq_output, saveDirectory, enriched.pathway = NULL
)
  {
  
  library(tidyverse)
  library(msigdbr)
  library(fgsea)
  library(ggtext)
  library(glue)
  ensembl_human = convert_mouse_human(deseq_output,FALSE,ensembl_to_ensembl = TRUE)
  
  library(org.Hs.eg.db)
  ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                      key=ensembl_human$Gene.stable.ID.1, 
                                      columns="SYMBOL",
                                      keytype="ENSEMBL")
  ens2symbol <- as_tibble(ens2symbol)
  ens2symbol = inner_join(ensembl_human, ens2symbol, by=c("Gene.stable.ID.1"="ENSEMBL"))
  res = as_tibble(deseq_output)
  res <- inner_join(res, ens2symbol, by=c("ENSEMBL"="Gene.stable.ID"))
  
  res2 <- res %>% 
    dplyr::select(SYMBOL, stat) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(SYMBOL) %>% 
    summarize(stat = mean(stat))
  ranks <- deframe(res2)
  pathways.hallmark <- gmtPathways("~/Downloads/h.all.v7.4.symbols.gmt")
  pathways.hallmark %>% 
    head() %>% 
    lapply(head)
  fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
    arrange(padj) %>% 
    DT::datatable()
  fgseaResWithoutLE = fgseaResTidy[,1:7]
  write.csv(fgseaResWithoutLE, paste0(saveDirectory,"/fgseaResWithoutLE.csv"))
  fgseaResTidySig = fgseaResTidy[fgseaResTidy$padj<0.05,]
  fgseaResTidySig$SpecialPathway = as.factor(fgseaResTidySig$pathway)
  if(is.null(enriched.pathway)){
    enriched.pathway = as.data.frame(fgseaResTidySig[,c("pathway","NES")])
  }
  bold.labels <- ifelse(fgseaResTidySig$SpecialPathway %in% enriched.pathway[i,1], yes = "bold", no = "plain")
  
  
  hallmarkPlot = ggplot(fgseaResTidySig, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=NES<0)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Hallmark pathways NES from GSEA") +
    theme(axis.text.y = element_text(face = rev(bold.labels)))
  
  ggsave(file.path(saveDirectory,"Hallmark.pdf"), hallmarkPlot)
  for (i in 1:nrow(enriched.pathway)){
    pathway_diff_genes = pathways.hallmark %>% 
      enframe("pathway", "SYMBOL") %>% 
      unnest(cols = c(SYMBOL)) %>% 
      inner_join(res, by="SYMBOL") %>%
      filter(pathway == enriched.pathway[i,1])
    pathway_diff_genes = pathway_diff_genes[!duplicated(pathway_diff_genes$ENSEMBL),]
    pathway_diff_genes = pathway_diff_genes[!is.na(pathway_diff_genes$padj) & pathway_diff_genes$padj < 0.05,]
    G_list = read.csv("~/Documents/Programming/Rahul_Sequencing/Project_10863_B/R_files/Ensembl_To_Mgi_Symbol.csv",row.names =  1)
    colnames(G_list) = c("ENSEMBL","MGI")
    pathway_diff_genes = left_join(pathway_diff_genes,G_list[!duplicated(G_list$ENSEMBL),])
    write.csv(pathway_diff_genes,paste0(saveDirectory,"/",enriched.pathway[i,1],".csv"))
  } 
  
}
