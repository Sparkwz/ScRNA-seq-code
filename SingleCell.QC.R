if (!dir.exists("./Outdata")){
  dir.create("./Outdata")
}

if (!dir.exists("./Outplot")){
  dir.create("./Outplot")
}

#### 1.QC caculate
sc.QC.caculate = function(merged_seurat,org = "human",cellcycle = T){
  ###Runs Seurat for QC calculate
  ##Inputs:
  #merged_seurat = Seruat object
  #org = "human" or "mouse"
  
  ##Outputs:
  #merged_seurat with
  #log10GenesPerUMI
  #Mitochondrial ratio, percent.mt
  #Ribosomal ratio, percent.ribo
  #Erythrocyte ratio, percent.hb
  #Cell cycle, S.Score and G2M.Score
  
if(org == "human"){
  ### 1.1 Number of genes detected per UMI
  # Add number of genes per UMI for each cell to metadata
  merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
  
  ### 1.2 QC index calculation
  ## 1.2.1 Mitochondrial ratio
  mito_genes=rownames(merged_seurat)[grep("^MT-", rownames(merged_seurat))]
  print("mito_genes:")
  print(mito_genes)
  merged_seurat <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-",col.name = "percent.mt")
  
  ## 1.2.2 Ribosomal ratio
  ribo_genes=rownames(merged_seurat)[grep("^RP[SL]", rownames(merged_seurat),ignore.case = T)]
  print("ribo_genes:")
  print(ribo_genes)
  
  merged_seurat <- PercentageFeatureSet(merged_seurat,"^RP[SL]",col.name = "percent.ribo")
  print(summary(merged_seurat@meta.data$percent.ribo))
  
  ## 1.2.3 Erythrocyte ratio
  hb_genes <- rownames(merged_seurat)[grep("^HB[^(P)]", rownames(merged_seurat),ignore.case = T)]
  print("hb_genes")
  print(hb_genes)
  merged_seurat=PercentageFeatureSet(merged_seurat, "^HB[^(P)]", col.name = "percent.hb")
  print(summary(merged_seurat@meta.data$percent.hb))
  
  if(cellcycle){
  ## 1.2.4 Cell cycle
  # human cell cycle genes
  s.genes=Seurat::cc.genes.updated.2019$s.genes
  print(s.genes)
  g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
  
  merged_seurat=CellCycleScoring(object = merged_seurat, 
                                 s.features = s.genes, 
                                 g2m.features = g2m.genes, 
                                 set.ident = TRUE)
  }
}
if(org == "mouse"){
    ### 1.1 Number of genes detected per UMI
    # Add number of genes per UMI for each cell to metadata
    merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
    
    ### 1.2 QC index calculation
    ## 1.2.1 Mitochondrial Ratio
    mito_genes=rownames(merged_seurat)[grep("^mt-", rownames(merged_seurat))]
    print("mito_genes:")
    print(mito_genes)
    merged_seurat <- PercentageFeatureSet(object = merged_seurat, pattern = "^mt-",col.name = "percent.mt")
    
    ## 1.2.2 Ribosomal ratio
    ribo_genes=rownames(merged_seurat)[grep("^Rp[sl]", rownames(merged_seurat),ignore.case = T)]
    print("ribo_genes:")
    print(ribo_genes)
    
    merged_seurat <- PercentageFeatureSet(merged_seurat,"^Rp[sl]",col.name = "percent.ribo")
    print(summary(merged_seurat@meta.data$percent.ribo))
    
    ## 1.2.3 Erythrocyte ratio
    hb_genes <- rownames(merged_seurat)[grep("^Hb[^(p)]", rownames(merged_seurat),ignore.case = T)]
    print("hb_genes")
    print(hb_genes)
    merged_seurat=PercentageFeatureSet(merged_seurat, "^Hb[^(p)]", col.name = "percent.hb")
    print(summary(merged_seurat@meta.data$percent.hb))
    
    if(cellcycle){
    ## 1.2.4 Cell cycle
    load("~/ref_annotation_Geneset/7.CellCycleRdata/scRNAseq_cc.genes.updated.2019_human_mouse.Rdata")
    s.genes=mus_s.genes
    g2m.genes=mus_g2m.genes
    
    merged_seurat=CellCycleScoring(object = merged_seurat, 
                                   s.features = s.genes, 
                                   g2m.features = g2m.genes, 
                                   set.ident = TRUE)
    }
  }
return(merged_seurat)
}

#### 2.QC plot
sc.QC.plot = function(metadata = merged_seurat@meta.data,group=NA,legend = "right"){
  ### refer to https://github.com/hbctraining/scRNA-seq_online/blob/master/lessons/04_SC_quality_control.md
  # Add cell IDs to metadata
  mytheme <- theme(legend.position = legend )
  
  metadata$cells <- rownames(metadata)
  if(!is.na(group[1])){
    metadata$group = group
  }
  # Rename columns
  metadata <- metadata %>%
    dplyr::rename(seq_folder = orig.ident,
                  nUMI = nCount_RNA,
                  nGene = nFeature_RNA)
  
  # Visualize the number of cell counts per sample
  p1 = metadata %>% 
    ggplot(aes(x=group, fill=group)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells")+mytheme
  
  # UMI counts (transcripts) per cell and Genes detected per cell
  # Visualize the number UMIs/transcripts per cell
  p2 = metadata %>% 
    ggplot(aes(color=group, x=nUMI, fill= group)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)+mytheme
  
  # Visualize the distribution of genes detected per cell via histogram
  p3 = metadata %>% 
    ggplot(aes(color=group, x=nGene, fill= group)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 300)+mytheme
  
  # Visualize the distribution of genes detected per cell via boxplot
  p4 = metadata %>% 
    ggplot(aes(x=group, y=log10(nGene), fill=group)) + 
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells vs NGenes")+theme(legend.position = "none")
  
  # Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
  p5 = metadata %>% 
    ggplot(aes(x=nUMI, y=nGene, color=percent.mt)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(~group)
  
  # Visualize the distribution of mitochondrial gene expression detected per cell
  p6 = metadata %>% 
    ggplot(aes(color=group, x=percent.mt, fill=group)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 20)
  
  # 计算核糖体基因比例
  p7 = metadata %>% 
    ggplot(aes(color=group, x=percent.ribo, fill=group)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 3)+mytheme
  
  #  计算红血细胞基因比例
  p8 = metadata %>% 
    ggplot(aes(color=group, x=percent.hb, fill=group)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 0.1)+mytheme
  
  # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
  p9 = metadata %>%
    ggplot(aes(x=log10GenesPerUMI, color = group, fill=group)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)+mytheme
  
  p.QC = (p1 | p2 | p3) / (p4 | p5 | p6) / (p7 | p8 | p9)  
  return(p.QC)
}

#### 3. Check for cell cycle 
scCellCycle <- function(seurat.qc = merged_sce,
                        org = "human",
                        num.cores= 4,
                        num.mem = 10,
                        save.dir = "./Outplot/Step3.check_data"){
  ###Runs Seurat for QC calculate
  ##Inputs:
  #seurat.qc = Seruat object
  #org = "human" or "mouse"
  ##Outputs:
  #cell cycle check pdf
  
  library(Seurat)
  library(dplyr)
  library(future)
  library(ggplot2)
  library(patchwork)
  plan("multisession", workers = num.cores)
  options(future.globals.maxSize = num.mem * 1024^3)
  if(org %in% c("human","mouse")){
    if(org == "human"){
      s.genes=Seurat::cc.genes.updated.2019$s.genes
      g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
    }else{
      load("~/ref_annotation_Geneset/7.CellCycleRdata/scRNAseq_cc.genes.updated.2019_human_mouse.Rdata")
      s.genes=mus_s.genes
      g2m.genes=mus_g2m.genes
    }
  }else{stop('org = "human" or "mouse"')}

  seurat.qc <- NormalizeData(seurat.qc,
                            normalization.method = "LogNormalize", 
                            scale.factor = 10000) %>%
    FindVariableFeatures(selection.method = "vst",
                         nfeatures = 2000,
                         verbose = FALSE) %>%
    ScaleData() %>%
    RunPCA()

  p1 <- VlnPlot(seurat.qc, features = c("PCNA","TOP2A","MCM6","MKI67"), group.by = "Phase", ncol = 4)
  p2 <- VlnPlot(seurat.qc, features = c("S.Score", "G2M.Score"), group.by = "orig.ident", 
             ncol = 2, pt.size = 0.1)
  p3 <- DimPlot(seurat.qc, reduction = "pca",group.by= "Phase")
  p4 <- DimPlot(seurat.qc,reduction = "pca", group.by= "Phase",split.by = "Phase")
  
  # Running a PCA on cell cycle genes
  seurat.qc <- RunPCA(seurat.qc, features = c(s.genes, g2m.genes))
  p5 = DimPlot(seurat.qc,
               reduction = "pca",
               group.by= "Phase")
  
  p6 = DimPlot(seurat.qc,
               reduction = "pca",
               group.by= "Phase",
               split.by = "Phase")
  
  p.cc = wrap_plots(ncol = 1,p1,p2,(p3 | p4),(p5 | p6))
  ggsave(p.cc, filename = paste0(save.dir,"/Step3.cellcycle.plot.before.scale.PDF"),height = 10,width = 8)

### scale cell cycle score
  seurat.qc$CC.Difference <- seurat.qc$S.Score - seurat.qc$G2M.Score
  seurat.qc <- ScaleData(seurat.qc, vars.to.regress = "CC.Difference")%>%
    RunPCA( features = c(s.genes, g2m.genes))
  p7 = DimPlot(seurat.qc,
               reduction = "pca",
               group.by= "Phase")
  
  p8 = DimPlot(seurat.qc,
               reduction = "pca",
               group.by= "Phase",
               split.by = "Phase")
  p.cc.scale = wrap_plots(ncol = 1,p7,p8)
  ggsave(p.cc.scale, filename = paste0(save.dir,"/Step3.cellcycle.plot.after.scale.PDF"),height = 5,width = 5)
  return(seurat.qc)
}
