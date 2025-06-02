library(ggplot2)
library(ggpubr)
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

# 自定义函数，用于计算均值和 SEM
mean_sem <- function(x) {
  return(data.frame(y = mean(x), ymin = mean(x) - sd(x) / sqrt(length(x)), ymax = mean(x) + sd(x) / sqrt(length(x))))
}

batch_cor <- function(gene = gene, exprSet = exprSet, rownames = gene.list,method="sp"){
  library(future.apply)
  plan("multisession", workers = 2)
  plan()
  #设置可用的内存
  options(future.globals.maxSize = 10 * 1024^3)
  y <- as.numeric(exprSet[gene,])
  rownames <- rownames[!rownames%in%gene]
  do.call(rbind,future_lapply(rownames, function(x){
    dd  <- cor.test(as.numeric(exprSet[x,]),y,method=method)
    data.frame(gene=gene,mRNAs=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}

if(T){
  text.size = 8
  text.angle = 45
  text.hjust = 1
  legend.position = "right"
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"), 
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                   #axis.line = element_line(color = "black"),
                   #axis.ticks = element_line(color = "black"),
                   #panel.grid.minor.y = element_blank(),
                   #panel.grid.minor.x = element_blank(),
                   panel.grid=element_blank(), # 去网格线
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size)
                   # panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")
                   # strip.background = element_rect(color="black",size= 1, linetype="solid") # fill="#FC4E07", 
  )
}

#### 0.QC caculate
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

#### 1.百分比条图、箱式图
plot.clusters.group = function (data = seurat_data,
                                Group,
                                legend.position = "top",
                                celltype.id,widths = c(3,1),
                                log =TRUE,
                                order=T,
                                text.size = 12,
                                legend.title = "Group",
                                color = 1,ylab = "Cell proportion (%)",
                                xlab = "",cell.counts = T){ 
  ## take an integrated Seurat object, plot distributions over orig.ident
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  library(paletteer)
  mytheme = theme(plot.title = element_text(size = text.size,color="black",hjust = 0.5),
                  axis.title = element_text(size = c(text.size-2),color ="black"), 
                  axis.text = element_text(size=c(text.size-2),color = "black"),
                  #axis.line = element_line(color = "black"),
                  #axis.ticks = element_line(color = "black"),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  # panel.grid=element_blank(), # 去网格线
                  # legend.position = "none",
                  legend.text = element_text(size= c(text.size-4)),
                  legend.title= element_text(size= c(text.size-4)),
                  # axis.text.x = element_text(angle = 45, hjust=1, vjust=1)
  )
  
  count_table <- table(data@meta.data[,Group], data@meta.data[,celltype.id])
  count_mtx <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx <- melt(count_mtx)
  melt_mtx$cluster <- melt_mtx$cluster
  
  cluster_size <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
  
  if(!is.factor(data@meta.data[,Group])){
    data@meta.data[,Group] = as.factor(data@meta.data[,Group])
    cluster_size$cluster <- factor(cluster_size$cluster,levels = levels(data@meta.data[,Group]))
    melt_mtx$cluster <- factor(melt_mtx$cluster,levels = levels(data@meta.data[,Group]))
  }
  if("0" %in% cluster_size$cluster){
    sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
  }else{
    sorted_labels <- paste(cluster_size$cluster[order(cluster_size$value)])
  }
  if(order){
    cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
    melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
  }else{
    cluster_size$cluster <- factor(cluster_size$cluster,levels = levels(data@meta.data[,Group]))
    melt_mtx$cluster <- factor(melt_mtx$cluster,levels = levels(data@meta.data[,Group]))
  }
  colnames(melt_mtx)[2] <- "dataset"
  
  if(log){
    p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per cluster") + ylab("") + mytheme
  }else{
    p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + xlab("Cells per cluster") + ylab("") + mytheme
  }
  
  if(color==1){
    if(length(unique(melt_mtx$dataset)) < 21){
      p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
        geom_bar(position="fill", stat="identity",) + theme_bw()  + 
        scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.9))+
        ylab(ylab) + xlab(xlab) + 
        theme(legend.position=legend.position) + guides(fill = guide_legend(title = legend.title)) +
        scale_y_continuous(labels =function(x) x * 100, expand = c(0.01, 0.01))+ mytheme
    }else{
      warning("The color limit is <21")
      p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
        geom_bar(position="fill", stat="identity") + theme_bw()  + 
        ylab(ylab) + xlab(xlab) +
        theme(legend.position=legend.position) + guides(fill = guide_legend(title = legend.title)) +
        scale_y_continuous(labels =function(x) x * 100, expand = c(0.01, 0.01))+ mytheme
    }
  }
  ###########################
  if(color==2){
    if(length(unique(melt_mtx$dataset)) < 9){
      p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
        geom_bar(position="fill", stat="identity",) + theme_bw() + 
        scale_fill_brewer(palette = "Set2")+
        ylab(ylab) + xlab(xlab) + 
        theme(legend.position=legend.position) + guides(fill = guide_legend(title = legend.title)) +
        scale_y_continuous(labels =function(x) x * 100, expand = c(0.01, 0.01))+ mytheme
    }else{
      warning("The color limit is <9")
      p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
        geom_bar(position="fill", stat="identity") + theme_bw()  + 
        ylab(ylab) + xlab(xlab) +
        theme(legend.position=legend.position) + guides(fill = guide_legend(title = legend.title)) +
        scale_y_continuous(labels =function(x) x * 100, expand = c(0.01, 0.01))+ mytheme
    }
  }
  if(cell.counts){
    p2 = wrap_plots(ncol = 2,p2 + coord_flip(),p1,widths = widths)
  }else{p2}
  return(p2)
} 

# 箱式图
#测试数据
scCellBoxplot <- function(meta.data,
                          sampleID, group, celltype,filter.cutoff = 20,p.method = "kruskal.test",
                          p.label="p.format",p.size = 2.8,ncol = 6,
                          text.size = 10,...){
    library(dplyr)
    library(ggpubr)
    if(T){
        mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                         axis.ticks = element_line(color = "black"),
                         axis.title = element_text(size = text.size,color ="black"), 
                         axis.text = element_text(size=text.size,color = "black"),
                         panel.grid=element_blank(), # 去网格线
                         legend.text = element_text(size= text.size),
                         legend.title= element_text(size= text.size)
        )
    }
    
    meta.tb = meta.data
    meta.tb$Celltype = meta.tb[,celltype]
    meta.tb$Group = meta.tb[,group]
    meta.tb$SampleID = meta.tb[,sampleID]
    keep.sample <- meta.tb %>%
        group_by(SampleID) %>%
        dplyr::summarize(count = n()) %>%
        filter(count > filter.cutoff) %>%
        pull(SampleID)
    
    meta.tb <- meta.tb %>%
        filter(SampleID %in% keep.sample)
    
    dat.plot <- meta.tb %>%
        mutate(Celltype = as.character(Celltype),
               Group = as.character(Group),
               SampleID = as.character(SampleID))
    
    dat.plot <- suppressMessages(dat.plot %>%
        group_by(SampleID, Celltype, Group) %>%
        dplyr::summarize(N = n()))
    
    sum.data <- suppressMessages(dat.plot %>%
        group_by(SampleID, Group) %>%
        dplyr::summarize(NTotal = sum(N)))
    
    dat.plot = suppressMessages(left_join(dat.plot, sum.data))
    dat.plot$freq = dat.plot$N / dat.plot$NTotal *100 
    head(dat.plot)
    
    if(is.factor(meta.data$Group)){
        dat.plot$Group = factor(dat.plot$Group, levels = levels(meta.data$Group))
    }
    
    p.list <- lapply(unique(sort(dat.plot$Celltype)),function(mcls){
        library(ggpubr)
        library(rstatix)
        p <- ggplot(filter(dat.plot,Celltype== mcls), aes(x = Group, y = freq)) + 
            geom_boxplot(aes(color = Group), outlier.shape = NA) +
            geom_point(aes(color = Group), size = 0.8, alpha = 1, 
                       position = position_jitter(width = 0.2))+
            # stat_summary(fun.data = 'mean_se', geom = "errorbar", color = "red", 
            #              width = 0.4, size = 0.8, position = position_dodge( .9)) +
            labs(title=mcls,x="",y="Frequency (%)")+
            #scale_color_manual(values=c("P"="#E41A1C","N"="#377EB8","T"="#4DAF4A")) +
            #scale_fill_manual(values=c("P"="#E41A1C","N"="#377EB8","T"="#4DAF4A")) +
            # stat_compare_means(label="p.format",size.label = 2.8,
            # comparisons=list(c("CON","H0"),c("H0","H6"),c("CON","H6"))) +
            stat_compare_means(label=p.label,size = p.size, method = p.method) +
            coord_cartesian(clip="off") + 
            theme_classic() + mytheme + 
            theme(legend.position = "none")
        return(p)
    })
    p.prop = wrap_plots(p.list,ncol = ncol)
    return(p.prop)
}
                           
#### 2.UMAP或TSNE可视化
TSNE.UMAP.Plot <- function(object,
                           groupBy,
                           Style=1,
                           plot.title = NA,
                           legend.point.size = 4,
                           reduction = "umap",
                           legend.position = "bottom",
                           label = FALSE,
                           label.size = 4,point.size = 0.5){
  library(tibble)
  library(ggplot2)
  library(paletteer)
  library(ggrepel)
  library(ggsci)
  library(stringr)
  # groupBy = 'celltype'
  # object = sce
  
  if(!is.null(plot.title)){
    if(is.na(plot.title)){
      plot.title = groupBy
    } 
  }
  # (1) 获取非线性降维坐标，可以选择tsne或者umap，前提是存在哦：
  plot_data = object@reductions[[reduction]]@cell.embeddings %>% 
    as.data.frame() %>% 
    dplyr::mutate(Barcode = rownames(.)) %>% 
    dplyr::inner_join(object@meta.data %>% tibble::rownames_to_column('Barcode'), by = 'Barcode') %>% 
    tibble::column_to_rownames('Barcode') %>% 
    dplyr::select(1,2, {{groupBy}}) %>% 
    dplyr::rename(Group = as.name(groupBy))
  colnames(plot_data) = gsub("_"," ",colnames(plot_data))
  
  if(grepl(x = reduction,pattern = '*tsne*' )) {
    colnames(plot_data)[1] = 'tSNE 1'
    colnames(plot_data)[2] = 'tSNE 2'
  } else {
    colnames(plot_data)[1] = 'UMAP 1'
    colnames(plot_data)[2] = 'UMAP 2'
  }
  
  # (2) 生成聚类中心坐标
  centroids = aggregate(as.matrix(plot_data[,c(1,2)]) ~ Group,
                        data = plot_data,
                        FUN = mean)
  
  # (3) 两种风格的绘图
  if(grepl(x = reduction,pattern = '*tsne*' )) {
    x = 'tSNE 1'; y = 'tSNE 2'
    segment.df=data.frame(x=c(-55,-55),
                          xend=c(-30,-55),
                          y=c(-44,-44),
                          yend=c(-44,-19))
  } else {
    x = 'UMAP 1'; y = 'UMAP 2'
    segment.df=data.frame(x=c(-15,-15),
                          xend=c(-5,-15),
                          y=c(-15,-15),
                          yend=c(-15,-5))
  }
  
  if(Style == 1){
    p1 =  DimPlot(object, reduction = reduction,group.by = groupBy,
                  label = label,label.box = T,label.size = label.size,repel = T) + theme_bw() +
      labs( x= colnames(plot_data)[1],y=colnames(plot_data)[2],title = plot.title) + 
      theme(panel.grid=element_blank(), # 去网格线
            plot.title = element_text(size = 12,color="black",hjust = 0.5),
            axis.text.x = element_text(size = 10, color = 'black'),
            axis.text.y = element_text(size = 10, color = 'black'),
            axis.title.x = element_text(size = 10, color = 'black'),
            axis.title.y = element_text(size = 10, color = 'black'),
            axis.ticks = element_line(color = 'black', lineend = 'round'),
            legend.position = legend.position,
            legend.text = element_text(size = 10, color = 'black'),
            legend.title = element_text(size = 10, color = 'black'),
            panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
      guides(color=guide_legend(override.aes = list(size=legend.point.size))) 
    
    if(length(unique(plot_data$Group))<21){
      p1 = p1 + 
        scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 1)) +
        scale_color_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 1))
    }
    return(p1)   
  }
  
  if(Style > 1){
    if(Style == 2){
      p2 = ggplot(data = plot_data, mapping = aes(x = !!as.name(x), y = !!as.name(y))) +
        geom_point(data = centroids,
                   mapping = aes(x = !!as.name(x), y = !!as.name(y),
                                 fill = Group),
                   color = 'white', shape = 21, stroke = 0.5, size = legend.point.size) +
        geom_point(mapping = aes(fill = Group),
                   color = 'white', shape = 21, stroke = 0.5, size = 3,show.legend = F) +
        theme_bw() +
        labs( x= colnames(plot_data)[1],y=colnames(plot_data)[2],title = plot.title) + 
        theme(
          axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 12,color="black",hjust = 0.5),
          axis.text.x = element_text(size = 10, color = 'black'),
          axis.text.y = element_text(size = 10, color = 'black'),
          axis.title.x = element_text(size = 10, color = 'black'),
          axis.title.y = element_text(size = 10, color = 'black'),
          axis.ticks = element_line(color = 'black', lineend = 'round'),
          legend.position = legend.position,
          legend.text = element_text(size = 10, color = 'black'),
          legend.title = element_text(size = 13, color = 'black'),
          panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
    }
    
    if(Style == 3){
      p2 = ggplot(data = plot_data,mapping = aes(x = !!as.name(x), y = !!as.name(y),color=Group))+
        geom_point(mapping = aes(color = Group),size = point.size,show.legend = F) +
        geom_point(data = centroids,
                   mapping = aes(x = !!as.name(x), y = !!as.name(y),
                                 color = Group),alpha = 1,size = 1,show.legend = T)+
        theme_classic() +
        labs( x= colnames(plot_data)[1],y=colnames(plot_data)[2],title = plot.title) +
        tidydr::theme_dr(arrow = grid::arrow(length = unit(0.3, "cm"), type = "open")) +
        # geom_segment(data = segment.df,mapping = aes(x=x,xend=xend,y=y,yend=yend),
        #              arrow = arrow(length=unit(0.3, "cm")),color = "black",size=0.6,show.legend = F) +
        #zfm:这个主题更简洁一点
        coord_cartesian(clip = "off")+
        theme(axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = legend.position,
              plot.title = element_text(size = 12,color="black",hjust = 0.5),
              # axis.title.x.bottom = element_text(hjust = 0.12,size = 12,margin=margin(-5,0,0,0)),
              # axis.title.y.left = element_text(hjust = 0.12,size = 12,margin=margin(0,-8,0,0)),
              axis.title.x.bottom = element_text(hjust = 0.12,size = 12),
              axis.title.y.left = element_text(hjust = 0.12,size = 12),
              axis.title.x = element_text(size = 12, color = 'black'),
              axis.title.y = element_text(size = 12, color = 'black'),
              # legend.text = element_text(size = 10, color = 'black'),
              legend.title = element_text(size = 10, color = 'black'),
              # plot.margin = margin(50,50,50,50),
              legend.background = element_blank())+
        guides(color=guide_legend(override.aes = list(size=legend.point.size))) 
    }
    if(Style == 4){
      p2 = ggplot(data = plot_data,mapping = aes(x = !!as.name(x), y = !!as.name(y),color=Group))+
        geom_point(data = centroids,
                   mapping = aes(x = !!as.name(x), y = !!as.name(y),
                                 fill = Group),
                   color = 'white', shape = 21, stroke = 0.5, size = legend.point.size) +
        geom_point(mapping = aes(fill = Group),
                   color = 'white', shape = 21, stroke = 0.5, size = 3,show.legend = F) +
        theme_classic() +
        labs( x= colnames(plot_data)[1],y=colnames(plot_data)[2],title = plot.title) +
        tidydr::theme_dr(arrow = grid::arrow(length = unit(0.3, "cm"), type = "open")) +
        # geom_segment(data = segment.df,mapping = aes(x=x,xend=xend,y=y,yend=yend),
        #              arrow = arrow(length=unit(0.3, "cm")),color = "black",size=0.6,show.legend = F) +
        #zfm:这个主题更简洁一点
        coord_cartesian(clip = "off")+
        theme(axis.line = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size = 12,color="black",hjust = 0.5),
              legend.position = legend.position,
              # axis.title.x.bottom = element_text(hjust = 0.12,size = 12,margin=margin(-5,0,0,0)),
              # axis.title.y.left = element_text(hjust = 0.12,size = 12,margin=margin(0,-8,0,0)),
              axis.title.x.bottom = element_text(hjust = 0.12,size = 12),
              axis.title.y.left = element_text(hjust = 0.12,size = 12),
              axis.title.x = element_text(size = 12, color = 'black'),
              axis.title.y = element_text(size = 12, color = 'black'),
              # legend.text = element_text(size = 10, color = 'black'),
              legend.title = element_text(size = 10, color = 'black'),
              # plot.margin = margin(50,50,50,50),
              legend.background = element_blank())
    }
    if(label == T){
      if(length(unique(plot_data$Group))<21){
        p3 = p2 + ggrepel::geom_label_repel(data = centroids,
                                            mapping = aes(fill = Group,
                                                          label= Group),
                                            size = label.size,
                                            fontface = 'plain', #bold
                                            color="black", 
                                            # family = 'Arial',
                                            show.legend = FALSE) + 
          scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 1)) +
          scale_color_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 1))
      }else{
        warning("Too many label box (> 20). Show label is off!")
        p3 = p2
        #   p3 = p2 + ggrepel::geom_label_repel(data = centroids,
        #                                       mapping = aes(fill = Group,
        #                                                     label= Group),
        #                                       size = label.size,
        #                                       fontface = 'plain', #bold
        #                                       color="black", 
        #                                       family = 'Arial',
        #                                       show.legend = FALSE)
        #   # scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65)) +
        #   # scale_color_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.65))   
      }
      return(p3)
    }
    if(label == F){
      if(length(unique(plot_data$Group))<21){
        p3 = p2 +
          scale_color_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 1))
        
      }else{
        p3 = p2  
      }
      return(p3)
    }
  }
}  

#### 3.多分组下的细胞比例的比较
plot.time.point = function(celltype = metadata$Celltype,group = metadata$group,
                           sampleID = "",
                           barplot = TRUE,
                           lineplot = FALSE, legend.title="Cell type",color = 1,
                           ylab= "Cellular fraction (%)",
                           xlab ="",plot.title = "",label.y = 100,label.size = 3,
                           my.test.method = "kruskal.test",error.bar = "sem",
                           text.size = 10,
                           text.angle = 45,
                           text.hjust = 1,
                           legend.position = "none",ncol= 8,
                           dot.size = 1,label.dot = TRUE,dot.alpha = 1,...){
  
  ######## R包加载  
  require(ggplot2)
  require(reshape)
  require(ggthemes)
  require(dplyr)
  require(ggpubr)
  if (!requireNamespace("REdaS", quietly = TRUE)){
    install.packages("REdaS") 
  }else(require(REdaS))
  
  ######## 判断因子水平
  if(!is.factor(celltype)){
    celltype = as.factor(celltype)
  }
  
  if(!is.factor(group)){
    group = as.factor(group)
  }
  
  ####### 构建输入数据   
  if(c(sampleID[1] == "")){
    meta.temp = data.frame(Celltype = celltype,Group = group)
  }else{
    meta.temp = data.frame(Celltype = celltype,Group = group,
                           SampleID = sampleID)
  }
  ####### 绘图风格及主题设置
  if(lineplot){
    barplot = FALSE
  }
  if(!barplot){
    lineplot = TRUE
  }
  if(T){
    mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                     axis.ticks = element_line(color = "black"),
                     axis.title = element_text(size = text.size,color ="black"), 
                     axis.text = element_text(size=text.size,color = "black"),
                     axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                     #axis.line = element_line(color = "black"),
                     #axis.ticks = element_line(color = "black"),
                     #panel.grid.minor.y = element_blank(),
                     #panel.grid.minor.x = element_blank(),
                     panel.grid=element_blank(), # 去网格线
                     legend.position = legend.position,
                     legend.text = element_text(size= text.size),
                     legend.title= element_text(size= text.size),
                     panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
                     strip.background = element_rect(color="black",size= 1, linetype="solid") # fill="#FC4E07", 
    )
  }
      
  ############### From my code: Famingzhao  
  # table.temp = metadata[,c("Celltype", "Group","SampleID")]
  if(T){
    warning("Input data including the three variables, method of Xiaoming will be performed!") 
    library(data.table)
    dat.plot = data.table(meta.temp)
    dat.plot <- dat.plot[,.(N=.N),by=c("SampleID","Celltype","Group")]
    sum.data <- dat.plot[,{.(NTotal = sum(.SD$N))},by=c("SampleID","Group")]
    dat.plot = left_join(dat.plot,sum.data)
    dat.plot$freq = dat.plot$N / dat.plot$NTotal *100
    dat.plot$Celltype = factor(dat.plot$Celltype,levels = levels(celltype))
    Data_summary <- summarySE(dat.plot, measurevar="freq", groupvars=c("Group","Celltype"))
    head(Data_summary)
    Data_summary[is.na(Data_summary)]= 0
    ################ Vis
    if(lineplot){
      p1 <- ggplot(dat.plot, aes(x = Group, y = freq,group = Celltype)) + 
        geom_line(data = Data_summary,aes(color=Celltype))+ theme_bw() + 
        facet_wrap(~Celltype, ncol = ncol) + 
        theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1), legend.position="bottom") + 
        labs(y= ylab,x= xlab,title = plot.title) +
        geom_point(data = Data_summary,aes(x= Group, y= freq),pch=19,
                   position=position_dodge(0.5),size= 1,alpha = 1)+
        mytheme + guides(color = guide_legend(title = legend.title)) +
        stat_compare_means(aes(group = Group),
                           label = "p.format",
                           method = my.test.method,
                           # paired = paired,
                           size = label.size,
                           label.y=label.y,
                           label.x = 1.5,
                           hide.ns = T)
      
      if(error.bar == "sem"){
        p1 = p1 +geom_errorbar(data = Data_summary,aes(ymin = freq-se, ymax= freq+se), 
                               width= 0.2, 
                               position= position_dodge(0.5), 
                               color="black",
                               alpha = 0.8)
      }
      if(error.bar == "sd"){
        p1 = p1 + geom_errorbar(data = Data_summary,aes(ymin = freq-sd, ymax= freq+sd), 
                                width= 0.2, 
                                position= position_dodge(0.5), 
                                color="black",
                                alpha = 0.8)
      }
      if(error.bar == "ci"){
        p1 = p1 + geom_errorbar(data = Data_summary,aes(ymin = freq-ci, ymax= freq+ci), 
                                width= 0.2, 
                                position= position_dodge(0.5), 
                                color="black",
                                alpha = 0.8) 
      }
      
      if(label.dot == T) {
        p1 = p1+geom_point(aes(group=SampleID,fill = Group),pch=21,alpha = dot.alpha,
                           position = position_dodge(0.3),size = dot.size)
      }
      if(color == 2){return(
        p1 + scale_color_brewer(palette = "Set2")
      )}else{return(
        p1 + scale_color_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.9))
      )}
    }
    #################
    if(barplot){ 
      p2 <- ggplot(dat.plot, aes(x = Group, y = freq)) + 
        geom_bar(stat="summary",fun=mean,width = NULL,aes(fill=Celltype))+
        # scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
        labs(y= ylab,x= xlab,title = plot.title) + 
        facet_wrap(~Celltype, ncol = ncol) + 
        theme_bw()+ mytheme +
        geom_point(data = Data_summary,aes(x= Group, y= freq),pch=19,
                   position=position_dodge(0.5),size= 1,alpha = 0)+
        theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1), legend.position= "none") +
        stat_compare_means(aes(group = Group),
                           label = "p.format",
                           # paired = paired,
                           method = my.test.method,
                           label.y=label.y,
                           size = label.size,
                           label.x = 1.5,
                           hide.ns = T)
      
      if(error.bar == "sem"){
        p2 = p2 +geom_errorbar(data = Data_summary,aes(ymin = freq-se, ymax= freq+se), 
                               width= 0.2, 
                               position= position_dodge(0.5), 
                               color="black",
                               alpha = 0.8)
      }
      if(error.bar == "sd"){
        p2 = p2 + geom_errorbar(data = Data_summary,aes(ymin = freq-sd, ymax= freq+sd), 
                                width= 0.2, 
                                position= position_dodge(0.5), 
                                color="black",
                                alpha = 0.8)
      }
      if(error.bar == "ci"){
        p2 = p2 + geom_errorbar(data = Data_summary,aes(ymin = freq-ci, ymax= freq+ci), 
                                width= 0.2, 
                                position= position_dodge(0.5), 
                                color="black",
                                alpha = 0.8) 
      }
      
      if(label.dot == T) {
        p2 = p2+geom_point(aes(group=SampleID,fill = Group),pch=21,alpha = dot.alpha,
                           position = position_dodge(0.3),size = dot.size)
      }
      
      if(color == 2){return(
        p2 + scale_fill_brewer(palette = "Set2")
      )}else{return(
        p2 +scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.9))
      )}
    }
  }
  
}

#### 5.半小提琴图
## 包装统计和绘图用的几个函数，代码来自https://gist.github.com/Karel-Kroeze/746685f5613e01ba820a31e57f87ec87
## 绘制误差线用的函数
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE,sd =TRUE,sem = FALSE ) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd

  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#绘制半小提琴图用的函数
GeomSplitViolin <- ggproto(
  "GeomSplitViolin", 
  GeomViolin, 
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- base::transform(data, 
                            xminv = x - violinwidth * (x - xmin), 
                            xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      base:: transform(data, x = if(grp%%2==1) xminv else xmaxv), 
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin", 
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

geom_split_violin <- function (mapping = NULL, 
                               data = NULL, 
                               stat = "ydensity", 
                               position = "identity", ..., 
                               draw_quantiles = NULL, 
                               trim = TRUE, 
                               scale = "area", 
                               na.rm = FALSE, 
                               show.legend = NA, 
                               inherit.aes = TRUE) {
  layer(data = data, 
        mapping = mapping, 
        stat = stat, 
        geom = GeomSplitViolin, 
        position = position, 
        show.legend = show.legend, 
        inherit.aes = inherit.aes, 
        params = list(trim = trim, 
                      scale = scale, 
                      draw_quantiles = draw_quantiles, 
                      na.rm = na.rm, ...)
  )
}

#### 6.选择PCA维数
pcs.select = function(seurat.data = seurat.qc){
  # Determine percent of variation associated with each PC
  pct <- seurat.data[["pca"]]@stdev / sum(seurat.data[["pca"]]@stdev) * 100
  
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # last point where change of % of variation is more than 0.1%.
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  print(paste0("pcs = ",pcs))
  # Create a dataframe with values
  plot_df <- data.frame(pct = pct, 
                        cumu = cumu, 
                        rank = 1:length(pct))
  
  # Elbow plot to visualize 
  p1 = ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
    geom_text() + 
    geom_vline(xintercept = 90, color = "grey") + 
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    theme_bw() 
  print(p1)
  return(pcs)
}

#### 7.gmt文件处理
read.gmt2list = function(filename){
  if(! file.exists(filename)) stop('File ',filename,' not available\n')
  dat=readLines(filename)
  n=length(dat)
  res=list(gs=vector(mode = "list", length = n),set_info = vector(mode = "expression", length = n) )
  names = list(geneset.names=vector(mode = "character", length = n))
  for(i in 1:n){
    s=strsplit(dat[i],'\t')[[1]]
    s = na.omit(dplyr::na_if(x = s,y = ""))
    res$gs[[i]]=s[-c(1:2)]
    names$geneset.names[i]=s[1]
  }
  names(res$gs)=names$geneset.names
  set_info = data.frame(row.names = names$geneset.names,setID = names$geneset.names,
                        name = names$geneset.names, namespace = rep("C2",length(names$geneset.names)),
                        distance = rep(""),length(names$geneset.names))
  res$set_info = set_info
  res
}
# org.Hs.egMsigdbC2REACTOME <- read.gmt2list('D:/Rdata/Code/GSEA_pathway/pca_pathways_ENTREZID.gmt')
# save(org.Hs.egMsigdbC2REACTOME,file = "D:/Rdata/Code/GSEA_pathway/nature_medicine_GSEA/ontology_Rdata/org.Hs.egMsigdbC2REACTOME.Rdata")
write.list2gmt <- function(gs,file){
  sink(file)
  lapply(names(gs), function(i){
    cat( paste(c(i,'tmp',gs[[i]]),collapse='\t') )
    cat('\n')
  })
  sink()
}
# data(gcSample)
# write.list2gmt(gcSample,"./test.gmt")
write.df2gmt <- function(termsID,genelist,file){
  geneset = data.frame(terms = termsID,genelist = genelist)
  name <- unique(geneset$terms)
  description <- name
  names(description) <- name
  genes <- lapply(name, function(name){
    as.vector(geneset[geneset$terms == name,"genelist"])
  })
  names(genes) <- name
  
  setClass("gmtInfo",slots=list(name="vector",description="vector",genes ="list"))
  gmtInfo <- new("gmtInfo",name=name,description=description,genes = genes)
  
  write.gmt1 <- function(filename,gmtInfo){
    if(class(gmtInfo) == "gmtInfo"){
      output <- file(filename, open="wt")
      lapply(1:length(gmtInfo@name),function(x){
        writeLines(paste(c(gmtInfo@description[x], gmtInfo@description[[x]],gmtInfo@genes[[x]]),collapse='\t'),
                   con=output)
      })
      close(output)
    }
  }
  write.gmt1(filename= file,gmtInfo = gmtInfo)
  return(gmtInfo@genes)
}

#### 8. gene set scores
## 8.1 解决as.matrix细胞数量太多的问题
if(F){
  Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){
  int k = z.size() ;
  IntegerMatrix  mat(nrows, ncols);
  for (int i = 0; i < k; i++){
      mat(rp[i],cp[i]) = z[i];
  }
  return mat;
}
' )#[1]
  
  
  as.my.matrix <- function(mat){
    
    row_pos <- mat@i
    col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
    
    tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                    nrows =  mat@Dim[1], ncols = mat@Dim[2])
    
    row.names(tmp) <- mat@Dimnames[[1]]
    colnames(tmp) <- mat@Dimnames[[2]]
    return(tmp)
  }
  # system.time({mymatrix <- as.matrix(mydgmatrix)})#看一下实际运行速度
  ##    user  system elapsed 
  ##   0.238   0.332   0.572
  #不用多说了吧，只用了as.data.frame的1/10，这只是一个测试数据，用真正大型的数据计算时效率会差的更多，并且这种循环的方式可以节省很多的内存。
}

## 8.2 sc.metabolism.Seurat 代码改编
sc.Pathway.Seurat = function (obj, method = "AUCell", imputation = F, ncores = 2, Spatial=F,
                                geneList = "KEGG",assay.names = "pathway",DefaultAssay=F) 
{
  signatures_KEGG_metab <- system.file("data", "KEGG_metabolism_nc.gmt", 
                                       package = "scMetabolism")
  signatures_REACTOME_metab <- system.file("data", "REACTOME_metabolism.gmt", 
                                           package = "scMetabolism")
  
  if (geneList == "KEGG") {
    gmtFile <- signatures_KEGG_metab
    cat("Your choice is: KEGG\n")
  }
  if (geneList == "REACTOME") {
    gmtFile <- signatures_REACTOME_metab
    cat("Your choice is: REACTOME\n")
  }
  ##########    
  if(!geneList %in% c("KEGG","REACTOME")){
    library(GSEABase)
    gmtFile <- geneList
    cat("Custom gene sets by users\n")
    if(! file.exists(gmtFile)) stop('File ',gmtFile,' not available\n')
    geneSets <- getGmt(gmtFile)
  }
  ########## 
  # library(scMetabolism)
  if(Spatial){
    countexp <- obj@assays$Spatial@counts
    countexp <- data.frame(as.matrix(countexp))
  }else{
    countexp <- obj@assays$RNA@counts
    countexp <- data.frame(as.matrix(countexp))
  }
  ##########    
  if (imputation == F) {
    countexp2 <- countexp
  }
  if (imputation == T) {
    library(rsvd)
    cat("Start imputation...\n")
    cat("Citation: George C. Linderman, Jun Zhao, Yuval Kluger. Zero-preserving imputation of scRNA-seq data using low-rank approximation. bioRxiv. doi: https://doi.org/10.1101/397588 \n")
    result.completed <- alra(as.matrix(t(countexp)))
    countexp2 <- result.completed[[3]]
    row.names(countexp2) <- colnames(countexp)
  }
  cat("Start quantify the pathway activity...\n")
  if (method == "VISION") {
    library(VISION)
    n.umi <- colSums(countexp2)
    scaled_counts <- t(t(countexp2)/n.umi) * median(n.umi)
    vis <- Vision(scaled_counts, signatures = gmtFile)
    options(mc.cores = ncores)
    vis <- analyze(vis)
    signature_exp <- data.frame(t(vis@SigScores))
  }
  if (method == "AUCell") {
    library(AUCell)
    library(GSEABase)
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), 
                                           nCores = ncores, plotStats = F)
    geneSets <- getGmt(gmtFile)
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
    signature_exp <- data.frame(getAUC(cells_AUC))
  }
  if (method == "ssGSEA") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile)
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("ssgsea"), 
                    kcdf = c("Poisson"), parallel.sz = ncores)
    signature_exp <- data.frame(gsva_es)
  }
  if (method == "gsva") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile)
    gsva_es <- gsva(as.matrix(countexp2), geneSets, method = c("gsva"), 
                    kcdf = c("Poisson"), parallel.sz = ncores)
    signature_exp <- data.frame(gsva_es)
  }
  colnames(signature_exp) = row.names(obj@meta.data)
  obj[[assay.names]] <- SeuratObject::CreateAssayObject(counts = signature_exp)
  obj <- SeuratObject::SetAssayData(obj, slot = "scale.data",
                                    new.data = as.matrix(signature_exp), assay = assay.names)
  if(DefaultAssay == T){
    DefaultAssay(obj) <- assay.names
  }
  obj
}

DotPlot.Pathway.2 = function (obj, pathway, phenotype, norm = "y") 
{
  input.norm = norm
  input.pathway <- as.character(pathway)
  input.parameter <- phenotype
  metadata <- obj@meta.data
  Pathway.matrix <- obj@assays$Pathway$score
  metadata[, input.parameter] <- as.character(metadata[, input.parameter])
  Pathway.matrix_sub <- t(Pathway.matrix[input.pathway, 
  ])
  gg_table <- c()
  for (i in 1:length(input.pathway)) {
    gg_table <- rbind(gg_table, cbind(metadata[, input.parameter], 
                                      input.pathway[i], Pathway.matrix_sub[, i]))
  }
  gg_table <- data.frame(gg_table)
  gg_table_median <- c()
  input.group.x <- unique(as.character(gg_table[, 1]))
  input.group.y <- unique(as.character(gg_table[, 2]))
  for (x in 1:length(input.group.x)) {
    for (y in 1:length(input.group.y)) {
      gg_table_sub <- subset(gg_table, gg_table[, 1] == 
                               input.group.x[x] & gg_table[, 2] == input.group.y[y])
      gg_table_median <- rbind(gg_table_median, cbind(input.group.x[x], 
                                                      input.group.y[y], median(as.numeric(as.character(gg_table_sub[, 
                                                                                                                    3])))))
    }
  }
  gg_table_median <- data.frame(gg_table_median)
  gg_table_median[, 3] <- as.numeric(as.character(gg_table_median[, 
                                                                  3]))
  gg_table_median_norm <- c()
  input.group.x <- unique(as.character(gg_table[, 1]))
  input.group.y <- unique(as.character(gg_table[, 2]))
  range01 <- function(x) {
    (x - min(x))/(max(x) - min(x))
  }
  if (input.norm == "y") 
    for (y in 1:length(input.group.y)) {
      gg_table_median_sub <- subset(gg_table_median, gg_table_median[, 
                                                                     2] == input.group.y[y])
      norm_value <- range01(as.numeric(as.character(gg_table_median_sub[, 
                                                                        3])))
      gg_table_median_sub[, 3] <- norm_value
      gg_table_median_norm <- rbind(gg_table_median_norm, 
                                    gg_table_median_sub)
    }
  if (input.norm == "x") 
    for (x in 1:length(input.group.x)) {
      gg_table_median_sub <- subset(gg_table_median, gg_table_median[, 
                                                                     1] == input.group.x[x])
      norm_value <- range01(as.numeric(as.character(gg_table_median_sub[, 
                                                                        3])))
      gg_table_median_sub[, 3] <- norm_value
      gg_table_median_norm <- rbind(gg_table_median_norm, 
                                    gg_table_median_sub)
    }
  if (input.norm == "na") 
    gg_table_median_norm <- gg_table_median
  gg_table_median_norm <- data.frame(gg_table_median_norm)
  gg_table_median_norm[, 3] <- as.numeric(as.character(gg_table_median_norm[, 
                                                                            3]))
  library(wesanderson)
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  if(is.factor(pathway)){
    gg_table_median_norm$X2 = factor(gg_table_median_norm$X2 ,levels = levels(pathway))
  }
  
  if(is.factor(obj@meta.data[,phenotype])){
    gg_table_median_norm$X1 = factor(gg_table_median_norm$X1 ,
                                     levels = levels(obj@meta.data[,phenotype]))
  }
  
  ggplot(data = gg_table_median_norm, aes(x = gg_table_median_norm[, 
                                                                   1], y = gg_table_median_norm[, 2], color = gg_table_median_norm[, 
                                                                                                                                   3])) + geom_point(data = gg_table_median_norm, aes(size = gg_table_median_norm[, 
                                                                                                                                                                                                                  3])) + ylab("Metabolic Pathway") + xlab(input.parameter) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, 
                                                  hjust = 1), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
    scale_color_gradientn(colours = pal) + labs(color = "Value", 
                                                size = "Value") + NULL
}

### 9. OR指数 + Ro/e
do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          meta.cluster = cellInfo.tb$meta.cluster,
                          colname.patient = "patient",
                          loc = cellInfo.tb$loc,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0,z.hi.OR=4,z.hi.ROE=2){
  ##input data
  # 
  library("Startrac")
  library(data.table)
  dir.create(dirname(out.prefix),F,T)
  
  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster = as.character(meta.cluster)
  
  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}
  
  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)
  
  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }
  
  startrac.dist <- unclass(Startrac::calTissueDist(cellInfo.tb,byPatient = F,
                                                   colname.cluster="meta.cluster",
                                                   colname.patient = "patient",
                                                   colname.tissue = "loc"))
  startrac.dist <- startrac.dist[,loc.avai.vec]
  
  # cuts <- c(0, 0.8, 1.2,Inf)
  # startrac.dist.bin.values <- factor(c("-", "+/-", "+"),levels=c("-", "+/-", "+"))
  # startrac.dist.bin <- matrix(startrac.dist.bin.values[findInterval(startrac.dist, cuts)],
  #                             ncol=ncol(startrac.dist))
  # colnames(startrac.dist.bin) <- colnames(startrac.dist)
  # rownames(startrac.dist.bin) <- rownames(startrac.dist)
  
  #	sscClust:::plot.matrix.simple(freq.dist.bin,
  #								  col.ht=rev(structure(colorRampPalette(brewer.pal(9,name="Blues"))(10),
  #												   names=0:9 )),
  #								  par.legend=list(labels=rev(sprintf("%s%%~%s%%",10*(0:9),c(10*(1:9),100) )),
  #												  at=0:9),
  #								  out.prefix=sprintf("%s.freq.dist",out.prefix),
  #								  show.number=F,clust.row=T,exp.name=expression(italic(Freq)),
  #								  #palatte=(brewer.pal(n = 7,name = "Blues")),
  #								  #palatte=viridis::viridis(7),
  #								  pdf.width = 4.5, pdf.height = pdf.height)
  
  #	sscClust:::plot.matrix.simple(startrac.dist.bin,
  #								  col.ht=rev(structure(viridis::viridis(3),
  #													   names=levels(startrac.dist.bin.values))),
  #								  out.prefix=sprintf("%s.startrac.dist.bin",out.prefix),
  #								  show.number=F,clust.row=T,exp.name=expression(italic(R)[o/e]),
  #								  pdf.width = pdf.width, pdf.height = pdf.height)
  
  sscVis::plotMatrix.simple(startrac.dist,
                            out.prefix=sprintf("%s.startrac.dist",out.prefix),
                            show.number=F,
                            clust.row=T,
                            #waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(R)[o/e]),
                            z.hi=z.hi.ROE,
                            #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
                            palatte=viridis::viridis(7),
                            pdf.width = pdf.width, pdf.height = pdf.height)
  
  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            #clust.row=T,
                            #par.legend=list(color_bar = "discrete",at=seq(0,4,0.5)),
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=z.hi.OR,
                            #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
                            palatte=viridis::viridis(7),
                            pdf.width = pdf.width, pdf.height = pdf.height)
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx,
                "startrac.dist.ROE" = startrac.dist))
  }else{
    return(OR.dist.mtx)
  }
  
}

test.dist.table <- function(count.dist,min.rowSum=0)
{
  library(plyr)
  count.dist <- count.dist[rowSums(count.dist)>=min.rowSum,,drop=F]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)
  setDT(count.dist.tb,keep.rownames=T)
  count.dist.melt.tb <- melt(count.dist.tb,id.vars="rn")
  colnames(count.dist.melt.tb) <- c("rid","cid","count")
  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col]-this.c
    this.m <- matrix(c(this.c,
                       sum.row[this.row]-this.c,
                       other.col.c,
                       sum(sum.col)-sum.row[this.row]-other.col.c),
                     ncol=2)
    res.test <- fisher.test(this.m)
    data.frame(rid=this.row,
               cid=this.col,
               p.value=res.test$p.value,
               OR=res.test$estimate)
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb,count.dist.melt.ext.tb,
                                  by=c("rid","cid"))
  count.dist.melt.ext.tb = data.table(count.dist.melt.ext.tb)
  count.dist.melt.ext.tb[,adj.p.value:=p.adjust(p.value,"BH")]
  #count.dist.melt.ext.tb[adj.p.value < 0.05,]
  #count.dist.melt.ext.tb[p.value < 0.05 & OR > 0,]
  
  return(count.dist.melt.ext.tb)
  
}

doPlotORHT <- function(OR.all.mtx,out.prefix,note.str="",
                       method.distance="",k=3,
                       do.hclust=T,pdf.width = 5.5, pdf.height = 10)
{
  
  OR.all.mtx.tmp <- OR.all.mtx
  OR.all.mtx.tmp[OR.all.mtx.tmp > 3] <- 3
  
  OR.hclust.row <- NULL
  
  if(do.hclust){
    OR.hclust.row <- run.cutree(OR.all.mtx.tmp,k=k,method.distance=method.distance,method.hclust="ward.D2")
    ##OR.hclust.row <- run.cutreeDynamic(OR.all.mtx,deepSplit=1, minClusterSize=2,method.distance="")
    OR.hclust.row$branch <- dendextend::set(OR.hclust.row$branch,"branches_lwd", 2)
  }else{
    th.OR <- 1.5
    OR.pattern <- OR.all.mtx > th.OR
    colnames(OR.pattern) <- sprintf("bin.%s",colnames(OR.pattern))
    ##OR.ext.tb <- cbind(data.table(cluster.name=rownames(OR.all.mtx)),
    OR.ext.tb <- cbind(data.table(cluster.name=rownames(OR.all.mtx.tmp)),
                       OR.pattern,OR.all.mtx)
    
    OR.tb.list <- list()
    cname.used <- c()
    OR.tb.list[["enrichT"]] <- OR.ext.tb[bin.T==TRUE,][order(-T),]
    cname.used <- c(cname.used,OR.tb.list[["enrichT"]][["cluster.name"]])
    OR.tb.list[["enrichP"]] <- OR.ext.tb[!(cluster.name %in% cname.used) & bin.P==TRUE,][order(-P),]
    cname.used <- c(cname.used,OR.tb.list[["enrichP"]][["cluster.name"]])
    OR.tb.list[["enrichN"]] <- OR.ext.tb[!(cluster.name %in% cname.used) & bin.N==TRUE,][order(-N),]
    cname.used <- c(cname.used,OR.tb.list[["enrichN"]][["cluster.name"]])
    OR.tb.list[["enrichO"]] <- OR.ext.tb[!(cluster.name %in% cname.used),]
    
    OR.tb.order.list <- llply(names(OR.tb.list),function(x){
      #x <- "enrichT"
      x.hclust <- run.cutree(as.matrix(OR.tb.list[[x]][,c("P","N","T")]),
                             k=1,method.distance=method.distance,method.hclust="ward.D2")
      ##OR.tb.list[[x]][rev(x.hclust$hclust$order),]
      OR.tb.list[[x]][(x.hclust$hclust$order),]
      #setorderv()
    })
    names(OR.tb.order.list) <- names(OR.tb.list)
    
    OR.order.mtx <- do.call(rbind,llply(names(OR.tb.order.list),function(x){
      a.mtx <- as.matrix(OR.tb.order.list[[x]][,c("P","N","T")])
      rownames(a.mtx) <- OR.tb.order.list[[x]][["cluster.name"]]
      return(a.mtx)
    }))
    OR.all.mtx <- OR.all.mtx[rownames(OR.order.mtx),]
  }
  
  mapping.col <- g.colSet$meta.cluster
  names(mapping.col) <- mcls2Name[names(mapping.col)]
  #col.row.man <- mapping.col[rownames(OR.all.mtx.tmp)[OR.hclust.row$hclust$order]]
  col.row.man <- mapping.col[rownames(OR.all.mtx)]
  
  sscVis:::plotMatrix.simple(OR.all.mtx,
                             col.ht=circlize::colorRamp2(c(0, 1, 3), viridis::viridis(3)),
                             out.prefix=sprintf("%s.OR.dist.rClust.withDend.a%s",out.prefix,note.str),
                             show.number=F, show.dendrogram=do.hclust,
                             clust.row=if(do.hclust) OR.hclust.row$branch else FALSE,
                             row_dend_width = unit(1.5, "cm"),
                             #par.legend=list(color_bar = "discrete",at=seq(0,4,0.5)),
                             #waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                             exp.name=expression(italic(OR)),
                             z.hi=3,
                             #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
                             #palatte=viridis::viridis(7),
                             par.heatmap=list(cex.row=1.5,
                                              row_names_gp=gpar(col=col.row.man,fontsize=10)),
                             pdf.width = pdf.width, pdf.height = pdf.height)
  
  return(OR.hclust.row)
  
}

gg.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          meta.cluster = cellInfo.tb$meta.cluster,
                          colname.patient = "patient",
                          loc = cellInfo.tb$loc,
                          out.prefix = NULL,
                          pdf.width = 10,
                          pdf.height = 3,
                          cuts =c(0.0000000, 0.0000001, 0.8, 1.5, 2, Inf),
                          bin.label = c("-", "+/-", "+", "++", "+++"),
                          verbose = 1,
                          z.hi.OR = 4,
                          z.hi.ROE = 2,
                          text_color=NULL,
                          text.size = 8,
                          text.angle = 45,
                          text.hjust = 1,
                          legend.position = "right") {
  
  if(T){
    mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                     # axis.ticks = element_line(color = "black"),
                     axis.title = element_text(size = text.size,color ="black"), 
                     axis.text = element_text(size=text.size,color = "black"),
                     axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                     #axis.line = element_line(color = "black"),
                     #axis.ticks = element_line(color = "black"),
                     #panel.grid.minor.y = element_blank(),
                     #panel.grid.minor.x = element_blank(),
                     panel.grid=element_blank(), # 去网格线
                     legend.position = legend.position,
                     legend.text = element_text(size= text.size),
                     legend.title= element_text(size= text.size)
                     # panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")
                     # strip.background = element_rect(color="black",size= 1, linetype="solid") # fill="#FC4E07", 
    )
  }
  
  ### 1.input data
  library("Startrac")
  library(data.table)
  if (!is.null(out.prefix)) {
    dir.create(dirname(out.prefix),F,T)
  }
  
  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster = as.character(meta.cluster)
  
  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}
  
  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  # print(freq.dist.bin)
  
  ### 2.计算OR值
  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }
  
  ### 3.计算Ro/e值
  startrac.dist <- unclass(Startrac::calTissueDist(cellInfo.tb,byPatient = F,
                                                   colname.cluster="meta.cluster",
                                                   colname.patient = "patient",
                                                   colname.tissue = "loc"))
  startrac.dist <- startrac.dist[,loc.avai.vec]
  
  startrac.dist.bin.values <- factor(bin.label,levels=bin.label)
  
  startrac.dist.bin <- matrix(startrac.dist.bin.values[findInterval(startrac.dist, cuts)],
                              ncol=ncol(startrac.dist))
  
  colnames(startrac.dist.bin) <- colnames(startrac.dist)
  rownames(startrac.dist.bin) <- rownames(startrac.dist)
  
  #### 4. 可视化
  ### 4.1 plotMatrix.simple可视化
  # sscVis::plotMatrix.simple(startrac.dist,
  #                           out.prefix=sprintf("%s.startrac.dist",out.prefix),
  #                           show.number=F,
  #                           clust.row=T,
  #                           #waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
  #                           exp.name=expression(italic(R)[o/e]),
  #                           z.hi=2,
  #                           #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
  #                           palatte=viridis::viridis(7),
  #                           pdf.width = 4.5, pdf.height = 4.5)
  
  ### 4.2 pheatmap可视化
  startrac.dist1 = ifelse(startrac.dist>z.hi.ROE,z.hi.ROE,startrac.dist)
  
  invisible({
    plot_1 = pheatmap::pheatmap(
      startrac.dist1,
      color = viridis::viridis(7),
      cluster_rows = T,
      cluster_cols = F,
      fontsize = 10,
      border_color = "white",
      treeheight_col = 10,
      treeheight_row = 10,silent = T
    )
  })
  
  ### 4.3 ggplot2可视化
  plot_data = reshape2::melt(startrac.dist1)
  head(plot_data)
  colnames(plot_data) = c("meta.cluster", "loc", "Roe")
  startrac.dist.bin2 = reshape2::melt(startrac.dist.bin)
  plot_data$bin = startrac.dist.bin2$value %>% factor(levels = bin.label)
  plot_data$meta.cluster = factor(plot_data$meta.cluster,
                                  levels = rev(plot_1$tree_row$labels[plot_1$tree_row$order]))
  head(plot_data)
  
  p3 <- ggplot(data = plot_data, aes(x = loc, y = meta.cluster, fill = Roe)) +
    geom_tile() +
    theme_minimal() +
    scale_x_discrete(expand = c(0.0, 0.0)) +
    scale_y_discrete(expand = c(0.0, 0.0)) +
    scale_fill_gradientn(name = "Ro/e",
                         colours = viridis::viridis(7)) +
    labs(x = NULL, y = NULL, title = "Ro/e") +
    mytheme
  
  suppressWarnings(if(!is.null(text_color)) {
    if(is.null(names(text_color)[1])){
      p3 = p3+theme(axis.text.y = element_text(color = text_color))
    }else{
      p3 = p3+theme(axis.text.y = element_text(color = text_color[levels(plot_data_OR$meta.cluster)]))
    }
  })
  
  colors= c(
    "-" = "white",
    "+/-" = "#fde6ce",
    "+" = "#fcc08b",
    "++" = "#f5904a",
    "+++" = "#e6540d"
  )
  
  labels = c("-" = 0, "+/-" = 0.8, "+" = 1.5, "++" = 2,"+++" = ">2")
  
  ### 4.4 等级可视化
  p4 <- ggplot(data = plot_data, aes(x = loc, y = meta.cluster, fill = bin)) +
    geom_tile() +
    scale_fill_manual(
      name = "Ro/e",
      values =  colors[unique(plot_data$bin)],
      labels = labels[unique(plot_data$bin)]
    ) +
    geom_text(
      aes(label = bin),
      color = "black",
      size = 3,
      show.legend = T
    ) +
    guides(fill = guide_legend(title = "Ro/e",
                               override.aes = list(label = c(
                                 names(colors[unique(plot_data$bin)])
                               )))) + theme_minimal() +
    scale_x_discrete(expand = c(0.0, 0.0)) +
    scale_y_discrete(expand = c(0.0, 0.0)) +
    labs(x = NULL, y = NULL, title = "Ro/e") +
    theme(axis.ticks.y = element_blank(),
          axis.ticks.length.y = unit(0, "cm")) + mytheme
  
  suppressWarnings(if(!is.null(text_color)) {
    if(is.null(names(text_color)[1])){
      p4 = p4+theme(axis.text.y = element_text(color = text_color))
    }else{
      p4 = p4+theme(axis.text.y = element_text(color = text_color[levels(plot_data_OR$meta.cluster)]))
    }
  })
  
  ### 4.5 ggplot2可视化OR指数
  OR.dist.mtx1 = ifelse(OR.dist.mtx > z.hi.OR, z.hi.OR, OR.dist.mtx)
  
  invisible({
    plot_2 = pheatmap::pheatmap(
      OR.dist.mtx1,
      color = viridis::viridis(7),
      cluster_rows = T,
      cluster_cols = F,
      fontsize = 10,
      border_color = "white",
      treeheight_col = 10,
      treeheight_row = 10,silent = T,
    )
  })
  
  plot_data_OR = reshape2::melt(OR.dist.mtx1)
  head(plot_data_OR)
  colnames(plot_data_OR) = c("meta.cluster", "loc", "Roe")
  plot_data_OR$meta.cluster = factor(plot_data_OR$meta.cluster,
                                     levels = rev(plot_2$tree_row$labels[plot_2$tree_row$order]))
  
  p5 <- ggplot(data = plot_data_OR, aes(x = loc, y = meta.cluster, fill = Roe)) +
    geom_tile() +
    theme_minimal() +
    scale_x_discrete(expand = c(0.0, 0.0)) +
    scale_y_discrete(expand = c(0.0, 0.0)) +
    scale_fill_gradientn(name = "OR",
                         colours = viridis::viridis(7)) +
    labs(x = NULL, y = NULL, title = "OR") +
    mytheme
  
  suppressWarnings(if(!is.null(text_color)) {
    if(is.null(names(text_color)[1])){
      p5 = p5+theme(axis.text.y = element_text(color = text_color))
    }else{
      p5 = p5+theme(axis.text.y = element_text(color = text_color[levels(plot_data_OR$meta.cluster)]))
    }
  })
  
  print(wrap_plots(p3,p4,p5))
  
  if (!is.null(out.prefix)) {
    ggsave(filename = paste0(out.prefix,".Roe.ggplot.heatmap.pdf"), 
           height = pdf.height, width = pdf.width)
  }
  
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx,
                "startrac.dist.ROE" = startrac.dist))
  }else{
    return(OR.dist.mtx)
  }
}

#### 10.富集分析
Myenrich <- function(genes, category = c("kegg", "gobp", "self"), gmt.files,
                     geneid = c("SYMBOL", "ENTREZID", "ENSEMBL", "UNIPROT"),
                     pvalueCutoff = 0.05,
                     qvalueCutoff=0.05){
  library(clusterProfiler)
  library(org.Hs.eg.db)
  category <- match.arg(category)
  geneid <- match.arg(geneid)
  # install.packages('R.utils')
  R.utils::setOption( "clusterProfiler.download.method",'auto' )
  if (category == "kegg"){
    if (geneid != "ENTREZID"){
      genes <- clusterProfiler::bitr(genes, fromType = geneid,
                                     toType ="ENTREZID",OrgDb="org.Hs.eg.db") %>% .[, 2] %>% as.character()
    }else{genes <- genes}
    enrich <- enrichKEGG(genes, organism = "hsa",keyType = "kegg",
                         pAdjustMethod = "BH",pvalueCutoff  = pvalueCutoff,
                         qvalueCutoff  = qvalueCutoff, maxGSSize = 5000)
    enrich <- setReadable(enrich, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  }
  if (category == "gobp"){
    if (geneid != "ENTREZID"){
      genes <- clusterProfiler::bitr(genes, fromType = geneid,toType ="ENTREZID",OrgDb="org.Hs.eg.db") %>% .[, 2] %>% as.character()
    }else{genes <- genes}
    enrich <- enrichGO(genes, OrgDb="org.Hs.eg.db",ont= "BP",pAdjustMethod = "BH",
                       pvalueCutoff  = pvalueCutoff,qvalueCutoff  = qvalueCutoff, readable = TRUE)
  }
  if (category == "self"){
    kegmt <- read.gmt(gmt.files) #读gmt文件
    enrich <- enricher(gene=genes,TERM2GENE=kegmt,
                     TERM2NAME = NA)
  }
  return(enrich)
}

#### 11.绘制气泡图
DotPlot_2 <- function(object,
                      features,
                      assay = NULL, 
                      scale = T,
                      dot.range.min = 0,
                      dot.range.max = 3.5,
                      Combine=F,
                      legend.position = "right",
                      label.size = 4,
                      label_widths = 0.1,
                      x.lab = NULL,
                      y.lab = NULL,
                      title = NULL,
                      text.size = 8,
                      text.angle = 90,
                      text.vjust = 0.5,
                      text.hjust = 1,
                      group.by = NULL,
                      color.use = NULL,
                      cols = c("lightgrey", 
                               "blue"),
                      legend.key.size = 0.5,
                      col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                      idents = NULL, split.by = NULL, cluster.idents = FALSE, 
                      scale.by = "radius", scale.min = NA, scale.max = NA,
                      ...
){
  library(dplyr)
  library(ggplot2)
  library(aplot)
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"), 
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust, 
                                              vjust = text.vjust), #,vjust = 0.5
                   panel.grid=element_blank(), # 去网格线
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size)
  )
  
  if(!is.null(group.by)){
    Idents(object) = group.by
  }
  
  if(Combine){
    if(is.null(color.use)){
      color.use <- alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 1)
    }
    df <- data.frame(x = 0, y = levels(object), stringsAsFactors = F )
    df$y <- factor(df$y, levels = df$y )
    p1 <- ggplot(df, aes(x, y, color = factor(y))) +
      geom_point(size = label.size, show.legend = F) +
      scale_color_manual(values = color.use) +
      theme_classic() +
      scale_x_continuous(expand = c(0,0)) + mytheme + 
      theme(
        plot.margin = margin(r=0),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 0,color="white"),
        axis.ticks = element_blank(),
        axis.line = element_blank()
      )
  }
  
  p2 = DotPlot(object = object,cols = cols,
               assay = assay,
               col.min = col.min, col.max = col.max, dot.min = dot.min, dot.scale = dot.scale, 
               idents = idents,split.by = split.by, cluster.idents = cluster.idents, 
               scale.by = scale.by, scale.min = scale.min, scale.max = scale.max,
               features = features,scale = scale,...)+theme_bw()+
    mytheme + labs(x = x.lab,y = y.lab,title = title)  +
    scale_size(range = c(dot.range.min,dot.range.max))+
    theme(legend.key.size = unit(legend.key.size, "cm"))+ 
    guides(size = guide_legend(title = "Per.Exp"),
           color = guide_colorbar(title = paste("Aver.Exp.","Scaled",sep = "\n")))
  if(Combine){
    p3 = p2 %>% insert_left(p1, width=label_widths)
    return(p3)
  }else{
    return(p2) 
  }
}

#### 12. ID转换 mouse to Human
convert_mouse_to_human <- function(gene_list){
  library(dplyr)
  mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  
  return (output)
}
# convert_mouse_to_human("Sca1")

### 13.ID转换
geneID_transform = function(input.data = gset_1, species = "human",
                            fromType = "ENTREZID",toType = c("SYMBOL"), biomaRt = F){
  library(clusterProfiler)
  if(biomaRt){
    (load("~/ref_annotation_Geneset/9.R_code/biomaRt_ensembl_annotation.Rdata"))
    if(fromType == "ENTREZID"){fromType = "entrezgene_id"}
    if(fromType == "ENSEMBL"){fromType = "ensembl_gene_id"}
    if(fromType == "SYMBOL"){fromType = "hgnc_symbol"}
    if(toType == "ENTREZID"){toType = "entrezgene_id"}
    if(toType == "ENSEMBL"){toType = "ensembl_gene_id"}
    if(toType == "SYMBOL"){toType = "hgnc_symbol"}
  }else{
    if(species == "human" | species== "org.Hs.eg.db"){
      library(org.Hs.eg.db) # for human
      species = "org.Hs.eg.db"
    }
    
    if(species == "mouse" | species == "org.Mm.eg.db") {
      library(org.Mm.eg.db) # for mouse
      species = "org.Mm.eg.db"
    }
  }
  
  library(limma)
  # ID 转换
  gene.df <- bitr(row.names(input.data), fromType = fromType, # ENSEMBEL
                  toType = toType,
                  OrgDb = get(species)) # 此处为human的数据集，若为小鼠的，替换为org.Mm.eg.db
  if("entrezgene_id" %in% colnames(gene.df)){
    gene.df$entrezgene_id = as.character(gene.df$entrezgene_id)
  }
  
  if(!is.data.frame(input.data)){input.data = as.data.frame(input.data)}
  data = input.data %>% tibble::rownames_to_column(var = fromType) %>% 
    dplyr::left_join(gene.df) %>% 
    dplyr::select(toType,everything())
  data = na.omit(data)
  
  # 得到表达矩阵
  data = as.data.frame(avereps(data[,-c(1,2)],ID = data[,toType]) )
  return(data)
}

### 14.cellchat
# 14.1 CellChat
computeNetSimilarity_GitHub <- function(object, slot.name = "netP", type = c("functional","structural"), k = NULL, thresh = NULL) {
  type <- match.arg(type)
  prob = methods::slot(object, slot.name)$prob
  if (is.null(k)) {
    if (dim(prob)[3] <= 25) {
      k <- ceiling(sqrt(dim(prob)[3]))
    } else {
      k <- ceiling(sqrt(dim(prob)[3])) + 1
    }
  }
  if (!is.null(thresh)) {
    prob[prob < quantile(c(prob[prob != 0]), thresh)] <- 0
  }
  if (type == "functional") {
    # compute the functional similarity
    D_signalings <- matrix(0, nrow = dim(prob)[3], ncol = dim(prob)[3])
    S2 <- D_signalings; S3 <- D_signalings;
    for (i in 1:(dim(prob)[3]-1)) {
      for (j in (i+1):dim(prob)[3]) {
        Gi <- (prob[ , ,i] > 0)*1
        Gj <- (prob[ , ,j] > 0)*1
        S3[i,j] <- sum(Gi * Gj)/sum(Gi+Gj-Gi*Gj,na.rm=TRUE)
      }
    }
    # define the similarity matrix
    S3[is.na(S3)] <- 0; S3 <- S3 + t(S3); diag(S3) <- 1
    # S_signalings <- S1 *S2
    S_signalings <- S3
  } else if (type == "structural") {
    # compute the structure distance
    D_signalings <- matrix(0, nrow = dim(prob)[3], ncol = dim(prob)[3])
    for (i in 1:(dim(prob)[3]-1)) {
      for (j in (i+1):dim(prob)[3]) {
        Gi <- (prob[ , ,i] > 0)*1
        Gj <- (prob[ , ,j] > 0)*1
        D_signalings[i,j] <- computeNetD_structure(Gi,Gj)
      }
    }
    # define the structure similarity matrix
    D_signalings[is.infinite(D_signalings)] <- 0
    D_signalings[is.na(D_signalings)] <- 0
    D_signalings <- D_signalings + t(D_signalings)
    S_signalings <- 1-D_signalings
  }
  # smooth the similarity matrix using SNN
  SNN <- buildSNN(S_signalings, k = k, prune.SNN = 1/15)
  Similarity <- as.matrix(S_signalings*SNN)
  rownames(Similarity) <- dimnames(prob)[[3]]
  colnames(Similarity) <- dimnames(prob)[[3]]
  comparison <- "single"
  comparison.name <- paste(comparison, collapse = "-")
  if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$matrix)) {
    methods::slot(object, slot.name)$similarity[[type]]$matrix <- NULL
  }
  methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]] <- Similarity
  return(object)
}
computeNetSimilarityPairwise_GitHub <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, k = NULL, thresh = NULL) {
  type <- match.arg(type)
  if (is.null(comparison)) {
    comparison <- 1:length(unique(object@meta$datasets))
  }
  cat("Compute signaling network similarity for datasets", as.character(comparison), '\n')
  comparison.name <- paste(comparison, collapse = "-")
  net <- list()
  signalingAll <- c()
  object.net.nameAll <- c()
  # 1:length(setdiff(names(methods::slot(object, slot.name)), "similarity"))
  for (i in 1:length(comparison)) {
    object.net <- methods::slot(object, slot.name)[[comparison[i]]]
    object.net.name <- names(methods::slot(object, slot.name))[comparison[i]]
    object.net.nameAll <- c(object.net.nameAll, object.net.name)
    net[[i]] = object.net$prob
    signalingAll <- c(signalingAll, paste0(dimnames(net[[i]])[[3]], "--", object.net.name))
    # signalingAll <- c(signalingAll, dimnames(net[[i]])[[3]])
  }
  names(net) <- object.net.nameAll
  net.dim <- sapply(net, dim)[3,]
  nnet <- sum(net.dim)
  position <- cumsum(net.dim); position <- c(0,position)
  if (is.null(k)) {
    if (nnet <= 25) {
      k <- ceiling(sqrt(nnet))
    } else {
      k <- ceiling(sqrt(nnet)) + 1
    }
  }
  if (!is.null(thresh)) {
    for (i in 1:length(net)) {
      neti <- net[[i]]
      neti[neti < quantile(c(neti[neti != 0]), thresh)] <- 0
      net[[i]] <- neti
    }
  }
  if (type == "functional") {
    # compute the functional similarity
    S3 <- matrix(0, nrow = nnet, ncol = nnet)
    for (i in 1:nnet) {
      for (j in 1:nnet) {
        idx.i <- which(position - i >= 0)[1]
        idx.j <- which(position - j >= 0)[1]
        net.i <- net[[idx.i-1]]
        net.j <- net[[idx.j-1]]
        Gi <- (net.i[ , ,i-position[idx.i-1]] > 0)*1
        Gj <- (net.j[ , ,j-position[idx.j-1]] > 0)*1
        S3[i,j] <- sum(Gi * Gj)/sum(Gi+Gj-Gi*Gj,na.rm=TRUE)
      }
    }
    # define the similarity matrix
    S3[is.na(S3)] <- 0;  diag(S3) <- 1
    S_signalings <- S3
  } else if (type == "structural") {
    # compute the structure distance
    D_signalings <- matrix(0, nrow = nnet, ncol = nnet)
    for (i in 1:nnet) {
      for (j in 1:nnet) {
        idx.i <- which(position - i >= 0)[1]
        idx.j <- which(position - j >= 0)[1]
        net.i <- net[[idx.i-1]]
        net.j <- net[[idx.j-1]]
        Gi <- (net.i[ , ,i-position[idx.i-1]] > 0)*1
        Gj <- (net.j[ , ,j-position[idx.j-1]] > 0)*1
        D_signalings[i,j] <- computeNetD_structure(Gi,Gj)
      }
    }
    # define the structure similarity matrix
    D_signalings[is.infinite(D_signalings)] <- 0
    D_signalings[is.na(D_signalings)] <- 0
    S_signalings <- 1-D_signalings
  }
  # smooth the similarity matrix using SNN
  SNN <- buildSNN(S_signalings, k = k, prune.SNN = 1/15)
  Similarity <- as.matrix(S_signalings*SNN)
  rownames(Similarity) <- signalingAll
  colnames(Similarity) <- rownames(Similarity)
  if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$matrix)) {
    methods::slot(object, slot.name)$similarity[[type]]$matrix <- NULL
  }
  # methods::slot(object, slot.name)$similarity[[type]]$matrix <- Similarity
  methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]] <- Similarity
  return(object)
}

### netEmbedding
netEmbedding_GitHub <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, pathway.remove = NULL, k = NULL) {
  if (object@options$mode == "single") {
    comparison <- "single"
    cat("Manifold learning of the signaling networks for a single dataset", '\n')
  } else if (object@options$mode == "merged") {
    if (is.null(comparison)) {
      comparison <- 1:length(unique(object@meta$datasets))
    }
    cat("Manifold learning of the signaling networks for datasets", as.character(comparison), '\n')
  }
  comparison.name <- paste(comparison, collapse = "-")
  Similarity <- methods::slot(object, slot.name)$similarity[[type]]$matrix[[comparison.name]]
  if (is.null(pathway.remove)) {
    pathway.remove <- rownames(Similarity)[which(colSums(Similarity) == 1)]
  }
  if (length(pathway.remove) > 0) {
    pathway.remove.idx <- which(rownames(Similarity) %in% pathway.remove)
    Similarity <- Similarity[-pathway.remove.idx, -pathway.remove.idx]
  }
  if (is.null(k)) {
    k <- ceiling(sqrt(dim(Similarity)[1])) + 1
  }
  options(warn = -1)
  # dimension reduction
  Y <- runUMAP(Similarity, min_dist = 0.3, n_neighbors = k)
  if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$dr)) {
    methods::slot(object, slot.name)$similarity[[type]]$dr <- NULL
  }
  methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]] <- Y
  return(object)
}
netClustering_GitHub <- function(object, slot.name = "netP", type = c("functional","structural"), comparison = NULL, k = NULL, methods = "kmeans", do.plot = TRUE, fig.id = NULL, do.parallel = TRUE, nCores = 4, k.eigen = NULL) {
  type <- match.arg(type)
  if (object@options$mode == "single") {
    comparison <- "single"
    cat("Classification learning of the signaling networks for a single dataset", '\n')
  } else if (object@options$mode == "merged") {
    if (is.null(comparison)) {
      comparison <- 1:length(unique(object@meta$datasets))
    }
    cat("Classification learning of the signaling networks for datasets", as.character(comparison), '\n')
  }
  comparison.name <- paste(comparison, collapse = "-")
  
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  Y[is.na(Y)] <- 0
  data.use <- Y
  if (methods == "kmeans") {
    if (!is.null(k)) {
      clusters = kmeans(data.use,k,nstart=10)$cluster
    } else {
      N <- nrow(data.use)
      kRange <- seq(2,min(N-1, 10),by = 1)
      if (do.parallel) {
        future::plan("multiprocess", workers = nCores)
        options(future.globals.maxSize = 1000 * 1024^2)
      }
      my.sapply <- ifelse(
        test = future::nbrOfWorkers() == 1,
        yes = pbapply::pbsapply,
        no = future.apply::future_sapply
      )
      results = my.sapply(
        X = 1:length(kRange),
        FUN = function(x) {
          idents <- kmeans(data.use,kRange[x],nstart=10)$cluster
          clusIndex <- idents
          #adjMat0 <- as.numeric(outer(clusIndex, clusIndex, FUN = "==")) - outer(1:N, 1:N, "==")
          adjMat0 <- Matrix::Matrix(as.numeric(outer(clusIndex, clusIndex, FUN = "==")), nrow = N, ncol = N)
          return(list(adjMat = adjMat0, ncluster = length(unique(idents))))
        },
        simplify = FALSE
      )
      adjMat <- lapply(results, "[[", 1)
      CM <- Reduce('+', adjMat)/length(kRange)
      res <- computeEigengap(as.matrix(CM))
      numCluster <- res$upper_bound
      clusters = kmeans(data.use,numCluster,nstart=10)$cluster
      if (do.plot) {
        gg <- res$gg.obj
        ggsave(filename= paste0("estimationNumCluster_",fig.id,"_",type,"_dataset_",comparison.name,".pdf"), plot=gg, width = 3.5, height = 3, units = 'in', dpi = 300)
      }
    }
  } else if (methods == "spectral") {
    A <- as.matrix(data.use)
    D <- apply(A, 1, sum)
    L <- diag(D)-A                       # unnormalized version
    L <- diag(D^-0.5)%*%L%*% diag(D^-0.5) # normalized version
    evL <- eigen(L,symmetric=TRUE)  # evL$values is decreasing sorted when symmetric=TRUE
    # pick the first k first k eigenvectors (corresponding k smallest) as data points in spectral space
    plot(rev(evL$values)[1:30])
    Z <- evL$vectors[,(ncol(evL$vectors)-k.eigen1):ncol(evL$vectors)]
    clusters = kmeans(Z,k,nstart=20)$cluster
  }
  if (!is.list(methods::slot(object, slot.name)$similarity[[type]]$group)) {
    methods::slot(object, slot.name)$similarity[[type]]$group <- NULL
  }
  methods::slot(object, slot.name)$similarity[[type]]$group[[comparison.name]] <- clusters
  return(object)
}

# 14.2 词云
computeEnrichmentScore <- function(df, measure = c("ligand", "signaling","LR-pair"), species = c('mouse','human'), color.use = NULL, color.name = "Dark2", n.color = 8,
                                   scale=c(4,.8), min.freq = 0, max.words = 200, random.order = FALSE, rot.per = 0,return.data = FALSE,seed = 1,...) {
  measure <- match.arg(measure)
  species <- match.arg(species)
  LRpairs <- as.character(unique(df$interaction_name))
  ES <- vector(length = length(LRpairs))
  for (i in 1:length(LRpairs)) {
    df.i <- subset(df, interaction_name == LRpairs[i])
    if (length(which(rowSums(is.na(df.i)) > 0)) > 0) {
      df.i <- df.i[-which(rowSums(is.na(df.i)) > 0), ,drop = FALSE]
    }
    ES[i] = mean(abs(df.i$ligand.logFC) * abs(df.i$receptor.logFC) *abs(df.i$ligand.pct.2-df.i$ligand.pct.1)*abs(df.i$receptor.pct.2-df.i$receptor.pct.1))
  }
  if (species == "mouse") {
    CellChatDB <- CellChatDB.mouse
  } else if (species == 'human') {
    CellChatDB <- CellChatDB.human
  }
  df.es <- CellChatDB$interaction[LRpairs, c("ligand",'receptor','pathway_name')]
  df.es$score <- ES
  # summarize the enrichment score
  df.es.ensemble <- df.es %>% group_by(ligand) %>% summarize(total = sum(score))  # avg = mean(score),
  
  set.seed(seed)
  if (is.null(color.use)) {
    color.use <- RColorBrewer::brewer.pal(n.color, color.name)
  }
  
  wordcloud::wordcloud(words = df.es.ensemble$ligand, freq = df.es.ensemble$total, min.freq = min.freq, max.words = max.words,scale=scale,
                       random.order = random.order, rot.per = rot.per, colors = color.use,...)
  if (return.data) {
    return(df.es.ensemble)
  }
}

# 14.3 可视化配受体对
netVisual.V2 = function (object, signaling, signaling.name = NULL, color.use = NULL, 
                         vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, 
                         top = 1, remove.isolate = FALSE, vertex.weight = 1, vertex.weight.max = NULL, 
                         vertex.size.max = NULL, weight.scale = TRUE, edge.weight.max.individual = NULL, 
                         edge.weight.max.aggregate = NULL, edge.width.max = 8, layout = c("circle", 
                                                                                          "hierarchy", "chord"), height = 5, thresh = 0.05, pt.title = 12, 
                         title.space = 6, vertex.label.cex = 0.8, from = NULL, to = NULL, 
                         bidirection = NULL, vertex.size = NULL, out.format = c("svg", 
                                                                                "png"), group = NULL, cell.order = NULL, small.gap = 1, 
                         big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, 
                         legend.pos.x = 20, legend.pos.y = 20, nCol = NULL,out.dir="./NetVisual/", ...) 
{
  if(!exists(out.dir)){dir.create(out.dir)}
  layout <- match.arg(layout)
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents))
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    }
    else {
      vertex.size.max <- 15
    }
  }
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, 
                       key = "pathway_name", matching.exact = T, pair.only = F)
  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  net <- object@net
  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval
  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[, , pairLR.name], 
                                         3, sum) != 0]
  }
  else {
    pairLR.name.use <- pairLR.name[sum(prob[, , pairLR.name]) != 
                                     0]
  }
  if (length(pairLR.name.use) == 0) {
    stop(paste0("There is no significant communication of ", 
                signaling.name))
  }
  else {
    pairLR <- pairLR[pairLR.name.use, ]
  }
  nRow <- length(pairLR.name.use)
  prob <- prob[, , pairLR.name.use]
  pval <- pval[, , pairLR.name.use]
  if (is.null(nCol)) {
    nCol <- min(length(pairLR.name.use), 2)
  }
  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify = "array")
    pval <- replicate(1, pval, simplify = "array")
  }
  if (is.null(edge.weight.max.individual)) {
    edge.weight.max.individual = max(prob)
  }
  prob.sum <- apply(prob, c(1, 2), sum)
  if (is.null(edge.weight.max.aggregate)) {
    edge.weight.max.aggregate = max(prob.sum)
  }
  if (layout == "hierarchy") {
    if (is.element("svg", out.format)) {
      svglite::svglite(file = paste0(out.dir,signaling.name, "_hierarchy_individual.svg"), 
                       width = 8, height = nRow * height)
      par(mfrow = c(nRow, 2), mar = c(5, 4, 4, 2) + 0.1)
      for (i in 1:length(pairLR.name.use)) {
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[, , i]
        netVisual_hierarchy1(prob.i, vertex.receiver = vertex.receiver, 
                             sources.use = sources.use, targets.use = targets.use, 
                             remove.isolate = remove.isolate, top = top, 
                             color.use = color.use, vertex.weight = vertex.weight, 
                             vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                             weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, 
                             edge.width.max = edge.width.max, title.name = signalName_i, 
                             vertex.label.cex = vertex.label.cex, ...)
        netVisual_hierarchy2(prob.i, vertex.receiver = setdiff(1:nrow(prob.i), 
                                                               vertex.receiver), sources.use = sources.use, 
                             targets.use = targets.use, remove.isolate = remove.isolate, 
                             top = top, color.use = color.use, vertex.weight = vertex.weight, 
                             vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                             weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, 
                             edge.width.max = edge.width.max, title.name = signalName_i, 
                             vertex.label.cex = vertex.label.cex, ...)
      }
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(out.dir,signaling.name, "_hierarchy_individual.png"), 
                     width = 8, height = nRow * height, units = "in", 
                     res = 300)
      par(mfrow = c(nRow, 2), mar = c(5, 4, 4, 2) + 0.1)
      for (i in 1:length(pairLR.name.use)) {
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[, , i]
        netVisual_hierarchy1(prob.i, vertex.receiver = vertex.receiver, 
                             sources.use = sources.use, targets.use = targets.use, 
                             remove.isolate = remove.isolate, top = top, 
                             color.use = color.use, vertex.weight = vertex.weight, 
                             vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                             weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, 
                             edge.width.max = edge.width.max, title.name = signalName_i, 
                             vertex.label.cex = vertex.label.cex, ...)
        netVisual_hierarchy2(prob.i, vertex.receiver = setdiff(1:nrow(prob.i), 
                                                               vertex.receiver), sources.use = sources.use, 
                             targets.use = targets.use, remove.isolate = remove.isolate, 
                             top = top, color.use = color.use, vertex.weight = vertex.weight, 
                             vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                             weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, 
                             edge.width.max = edge.width.max, title.name = signalName_i, 
                             vertex.label.cex = vertex.label.cex, ...)
      }
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      grDevices::cairo_pdf(paste0(signaling.name, "_hierarchy_individual.pdf"), 
                           width = 8, height = nRow * height)
      par(mfrow = c(nRow, 2), mar = c(5, 4, 4, 2) + 0.1)
      for (i in 1:length(pairLR.name.use)) {
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[, , i]
        netVisual_hierarchy1(prob.i, vertex.receiver = vertex.receiver, 
                             sources.use = sources.use, targets.use = targets.use, 
                             remove.isolate = remove.isolate, top = top, 
                             color.use = color.use, vertex.weight = vertex.weight, 
                             vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                             weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, 
                             edge.width.max = edge.width.max, title.name = signalName_i, 
                             vertex.label.cex = vertex.label.cex, ...)
        netVisual_hierarchy2(prob.i, vertex.receiver = setdiff(1:nrow(prob.i), 
                                                               vertex.receiver), sources.use = sources.use, 
                             targets.use = targets.use, remove.isolate = remove.isolate, 
                             top = top, color.use = color.use, vertex.weight = vertex.weight, 
                             vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                             weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, 
                             edge.width.max = edge.width.max, title.name = signalName_i, 
                             vertex.label.cex = vertex.label.cex, ...)
      }
      dev.off()
    }
    if (is.element("svg", out.format)) {
      svglite::svglite(file = paste0(out.dir,signaling.name, "_hierarchy_aggregate.svg"), 
                       width = 7, height = 1 * height)
      par(mfrow = c(1, 2), ps = pt.title)
      netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, 
                           sources.use = sources.use, targets.use = targets.use, 
                           remove.isolate = remove.isolate, top = top, 
                           color.use = color.use, vertex.weight = vertex.weight, 
                           vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                           weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, 
                           edge.width.max = edge.width.max, title.name = NULL, 
                           vertex.label.cex = vertex.label.cex, ...)
      netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum), 
                                                               vertex.receiver), sources.use = sources.use, 
                           targets.use = targets.use, remove.isolate = remove.isolate, 
                           top = top, color.use = color.use, vertex.weight = vertex.weight, 
                           vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                           weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, 
                           edge.width.max = edge.width.max, title.name = NULL, 
                           vertex.label.cex = vertex.label.cex, ...)
      graphics::mtext(paste0(signaling.name, " signaling pathway network"), 
                      side = 3, outer = TRUE, cex = 1, line = -title.space)
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(out.dir,signaling.name, "_hierarchy_aggregate.png"), 
                     width = 7, height = 1 * height, units = "in", 
                     res = 300)
      par(mfrow = c(1, 2), ps = pt.title)
      netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, 
                           sources.use = sources.use, targets.use = targets.use, 
                           remove.isolate = remove.isolate, top = top, 
                           color.use = color.use, vertex.weight = vertex.weight, 
                           vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                           weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, 
                           edge.width.max = edge.width.max, title.name = NULL, 
                           vertex.label.cex = vertex.label.cex, ...)
      netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum), 
                                                               vertex.receiver), sources.use = sources.use, 
                           targets.use = targets.use, remove.isolate = remove.isolate, 
                           top = top, color.use = color.use, vertex.weight = vertex.weight, 
                           vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                           weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, 
                           edge.width.max = edge.width.max, title.name = NULL, 
                           vertex.label.cex = vertex.label.cex, ...)
      graphics::mtext(paste0(signaling.name, " signaling pathway network"), 
                      side = 3, outer = TRUE, cex = 1, line = -title.space)
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      grDevices::cairo_pdf(paste0(signaling.name, "_hierarchy_aggregate.pdf"), 
                           width = 7, height = 1 * height)
      par(mfrow = c(1, 2), ps = pt.title)
      netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, 
                           sources.use = sources.use, targets.use = targets.use, 
                           remove.isolate = remove.isolate, top = top, 
                           color.use = color.use, vertex.weight = vertex.weight, 
                           vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                           weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, 
                           edge.width.max = edge.width.max, title.name = NULL, 
                           vertex.label.cex = vertex.label.cex, ...)
      netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum), 
                                                               vertex.receiver), sources.use = sources.use, 
                           targets.use = targets.use, remove.isolate = remove.isolate, 
                           top = top, color.use = color.use, vertex.weight = vertex.weight, 
                           vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                           weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, 
                           edge.width.max = edge.width.max, title.name = NULL, 
                           vertex.label.cex = vertex.label.cex, ...)
      graphics::mtext(paste0(signaling.name, " signaling pathway network"), 
                      side = 3, outer = TRUE, cex = 1, line = -title.space)
      dev.off()
    }
  }
  else if (layout == "circle") {
    if (is.element("svg", out.format)) {
      svglite::svglite(file = paste0(out.dir,signaling.name, "_", 
                                     layout, "_individual.svg"), width = height, 
                       height = nRow * height)
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), 
                    nCol), xpd = TRUE)
      for (i in 1:length(pairLR.name.use)) {
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[, , i]
        netVisual_circle(prob.i, sources.use = sources.use, 
                         targets.use = targets.use, remove.isolate = remove.isolate, 
                         top = top, color.use = color.use, vertex.weight = vertex.weight, 
                         vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                         weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, 
                         edge.width.max = edge.width.max, title.name = signalName_i, 
                         vertex.label.cex = vertex.label.cex, ...)
      }
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(out.dir,signaling.name, "_", layout, 
                            "_individual.png"), width = height, height = nRow * 
                       height, units = "in", res = 300)
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), 
                    nCol), xpd = TRUE)
      for (i in 1:length(pairLR.name.use)) {
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[, , i]
        netVisual_circle(prob.i, sources.use = sources.use, 
                         targets.use = targets.use, remove.isolate = remove.isolate, 
                         top = top, color.use = color.use, vertex.weight = vertex.weight, 
                         vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                         weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, 
                         edge.width.max = edge.width.max, title.name = signalName_i, 
                         vertex.label.cex = vertex.label.cex, ...)
      }
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      grDevices::cairo_pdf(paste0(signaling.name, "_", 
                                  layout, "_individual.pdf"), width = height, 
                           height = nRow * height)
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), 
                    nCol), xpd = TRUE)
      for (i in 1:length(pairLR.name.use)) {
        signalName_i <- pairLR$interaction_name_2[i]
        prob.i <- prob[, , i]
        netVisual_circle(prob.i, sources.use = sources.use, 
                         targets.use = targets.use, remove.isolate = remove.isolate, 
                         top = top, color.use = color.use, vertex.weight = vertex.weight, 
                         vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                         weight.scale = weight.scale, edge.weight.max = edge.weight.max.individual, 
                         edge.width.max = edge.width.max, title.name = signalName_i, 
                         vertex.label.cex = vertex.label.cex, ...)
      }
      dev.off()
    }
    if (is.element("svg", out.format)) {
      svglite(file = paste0(out.dir,signaling.name, "_", layout, 
                            "_aggregate.svg"), width = height, height = 1 * 
                height)
      netVisual_circle(prob.sum, sources.use = sources.use, 
                       targets.use = targets.use, remove.isolate = remove.isolate, 
                       top = top, color.use = color.use, vertex.weight = vertex.weight, 
                       vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                       weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, 
                       edge.width.max = edge.width.max, title.name = paste0(signaling.name, 
                                                                            " signaling pathway network"), vertex.label.cex = vertex.label.cex, 
                       ...)
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(out.dir,signaling.name, "_", layout, 
                            "_aggregate.png"), width = height, height = 1 * 
                       height, units = "in", res = 300)
      netVisual_circle(prob.sum, sources.use = sources.use, 
                       targets.use = targets.use, remove.isolate = remove.isolate, 
                       top = top, color.use = color.use, vertex.weight = vertex.weight, 
                       vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                       weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, 
                       edge.width.max = edge.width.max, title.name = paste0(signaling.name, 
                                                                            " signaling pathway network"), vertex.label.cex = vertex.label.cex, 
                       ...)
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      grDevices::cairo_pdf(paste0(signaling.name, "_", 
                                  layout, "_aggregate.pdf"), width = height, height = 1 * 
                             height)
      netVisual_circle(prob.sum, sources.use = sources.use, 
                       targets.use = targets.use, remove.isolate = remove.isolate, 
                       top = top, color.use = color.use, vertex.weight = vertex.weight, 
                       vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, 
                       weight.scale = weight.scale, edge.weight.max = edge.weight.max.aggregate, 
                       edge.width.max = edge.width.max, title.name = paste0(signaling.name, 
                                                                            " signaling pathway network"), vertex.label.cex = vertex.label.cex, 
                       ...)
      dev.off()
    }
  }
  else if (layout == "chord") {
    if (is.element("svg", out.format)) {
      svglite::svglite(file = paste0(out.dir,signaling.name, "_", 
                                     layout, "_individual.svg"), width = height, 
                       height = nRow * height)
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), 
                    nCol), xpd = TRUE)
      for (i in 1:length(pairLR.name.use)) {
        title.name <- pairLR$interaction_name_2[i]
        net <- prob[, , i]
        netVisual_chord_cell_internal(net, color.use = color.use, 
                                      sources.use = sources.use, targets.use = targets.use, 
                                      remove.isolate = remove.isolate, group = group, 
                                      cell.order = cell.order, lab.cex = vertex.label.cex, 
                                      small.gap = small.gap, big.gap = big.gap, 
                                      scale = scale, reduce = reduce, title.name = title.name, 
                                      show.legend = show.legend, legend.pos.x = legend.pos.x, 
                                      legend.pos.y = legend.pos.y)
      }
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(out.dir,signaling.name, "_", layout, 
                            "_individual.png"), width = height, height = nRow * 
                       height, units = "in", res = 300)
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), 
                    nCol), xpd = TRUE)
      for (i in 1:length(pairLR.name.use)) {
        title.name <- pairLR$interaction_name_2[i]
        net <- prob[, , i]
        netVisual_chord_cell_internal(net, color.use = color.use, 
                                      sources.use = sources.use, targets.use = targets.use, 
                                      remove.isolate = remove.isolate, group = group, 
                                      cell.order = cell.order, lab.cex = vertex.label.cex, 
                                      small.gap = small.gap, big.gap = big.gap, 
                                      scale = scale, reduce = reduce, title.name = title.name, 
                                      show.legend = show.legend, legend.pos.x = legend.pos.x, 
                                      legend.pos.y = legend.pos.y)
      }
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      grDevices::cairo_pdf(paste0(signaling.name, "_", 
                                  layout, "_individual.pdf"), width = height, 
                           height = nRow * height)
      par(mfrow = c(ceiling(length(pairLR.name.use)/nCol), 
                    nCol), xpd = TRUE)
      for (i in 1:length(pairLR.name.use)) {
        title.name <- pairLR$interaction_name_2[i]
        net <- prob[, , i]
        netVisual_chord_cell_internal(net, color.use = color.use, 
                                      sources.use = sources.use, targets.use = targets.use, 
                                      remove.isolate = remove.isolate, group = group, 
                                      cell.order = cell.order, lab.cex = vertex.label.cex, 
                                      small.gap = small.gap, big.gap = big.gap, 
                                      scale = scale, reduce = reduce, title.name = title.name, 
                                      show.legend = show.legend, legend.pos.x = legend.pos.x, 
                                      legend.pos.y = legend.pos.y)
      }
      dev.off()
    }
    if (is.element("svg", out.format)) {
      svglite(file = paste0(out.dir,signaling.name, "_", layout, 
                            "_aggregate.svg"), width = height, height = 1 * 
                height)
      netVisual_chord_cell_internal(prob.sum, color.use = color.use, 
                                    sources.use = sources.use, targets.use = targets.use, 
                                    remove.isolate = remove.isolate, group = group, 
                                    cell.order = cell.order, lab.cex = vertex.label.cex, 
                                    small.gap = small.gap, big.gap = big.gap, scale = scale, 
                                    reduce = reduce, title.name = paste0(signaling.name, 
                                                                         " signaling pathway network"), show.legend = show.legend, 
                                    legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y)
      dev.off()
    }
    if (is.element("png", out.format)) {
      grDevices::png(paste0(out.dir,signaling.name, "_", layout, 
                            "_aggregate.png"), width = height, height = 1 * 
                       height, units = "in", res = 300)
      netVisual_chord_cell_internal(prob.sum, color.use = color.use, 
                                    sources.use = sources.use, targets.use = targets.use, 
                                    remove.isolate = remove.isolate, group = group, 
                                    cell.order = cell.order, lab.cex = vertex.label.cex, 
                                    small.gap = small.gap, big.gap = big.gap, scale = scale, 
                                    reduce = reduce, title.name = paste0(signaling.name, 
                                                                         " signaling pathway network"), show.legend = show.legend, 
                                    legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y)
      dev.off()
    }
    if (is.element("pdf", out.format)) {
      grDevices::cairo_pdf(paste0(signaling.name, "_", 
                                  layout, "_aggregate.pdf"), width = height, height = 1 * 
                             height)
      netVisual_chord_cell_internal(prob.sum, color.use = color.use, 
                                    sources.use = sources.use, targets.use = targets.use, 
                                    remove.isolate = remove.isolate, group = group, 
                                    cell.order = cell.order, lab.cex = vertex.label.cex, 
                                    small.gap = small.gap, big.gap = big.gap, scale = scale, 
                                    reduce = reduce, title.name = paste0(signaling.name, 
                                                                         " signaling pathway network"), show.legend = show.legend, 
                                    legend.pos.x = legend.pos.x, legend.pos.y = legend.pos.y)
      dev.off()
    }
  }
}

#差异性circle plot
netVisual_diffInteraction.V2 <- function (object, comparison = c(1, 2), measure = c("count", 
                                                                                    "weight", 
                                                                                    "count.merged", 
                                                                                    "weight.merged"),
                                          color.use = NULL, 
                                          color.edge = c("#b2182b", "#2166ac"), title.name = NULL, 
                                          sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, 
                                          top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, 
                                          vertex.size.max = 15, vertex.label.cex = 1, vertex.label.color = "black", 
                                          edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6, 
                                          label.edge = FALSE, edge.label.color = "black", edge.label.cex = 0.8, 
                                          edge.curved = 0.2, shape = "circle", layout = in_circle(), 
                                          margin = 0.2, arrow.width = 1, arrow.size = 0.2) 
{
  options(warn = -1)
  measure <- match.arg(measure)
  obj1 <- object@net[[comparison[1]]][[measure]]
  obj2 <- object@net[[comparison[2]]][[measure]]
  net.diff <- obj2 - obj1
  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  }
  else if (measure %in% c("weight", "weight.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                          df.net[["target"]]), sum)
  }
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  net[is.na(net)] = 0
  net[abs(net) < stats::quantile(abs(net), probs = 1 - top)] <- 0
  g <- graph_from_adjacency_matrix(net, mode = "directed", 
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max + 
    5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, 
                       -atan(coords_scale[igraph::V(g), 2]/coords_scale[igraph::V(g), 
                                                                        1]), pi - atan(coords_scale[igraph::V(g), 2]/coords_scale[igraph::V(g), 
                                                                                                                                  1]))
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex <- vertex.label.cex
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$label.cex <- edge.label.cex
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1], 
                               color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, 
                                               alpha.edge)
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * 
      edge.width.max
  }
  else {
    igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
  }
  if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
    igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[, 
                                                                1])] <- loop.angle[edge.start[which(edge.start[, 
                                                                                                               2] == edge.start[, 1]), 1]]
  }
  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x = 1:length(igraph::V(g)), 
                               direction = -1, start = 0)
  label.dist <- vertex.weight/max(vertex.weight) + 2
  plot(g, edge.curved = edge.curved, vertex.shape = shape, 
       layout = coords_scale, margin = margin, vertex.label.dist = label.dist, 
       vertex.label.degree = label.locs, vertex.label.family = "Helvetica", 
       edge.label.family = "Helvetica")
  if (!is.null(title.name)) {
    text(0, 1.5, title.name, cex = 1.1)
  }
  gg <- recordPlot()
  return(gg)
}

#差异性heatmap
netVisual_heatmap.V2 <- function (object, comparison = c(1, 2), measure = c("count", 
                                                                            "weight", 
                                                                            "count.merged", 
                                                                            "weight.merged"), 
                                  signaling = NULL, slot.name = c("netP", "net"), 
                                  color.use = NULL, color.heatmap = c("#2166ac", "#b2182b"), 
                                  title.name = NULL, width = NULL, height = NULL, font.size = 8, 
                                  font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE, 
                                  sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, 
                                  row.show = NULL, col.show = NULL) 
{
  library(ComplexHeatmap)
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (is.list(object@net[[1]])) {
    message("Do heatmap based on a merged object \n")
    obj1 <- object@net[[comparison[1]]][[measure]]
    obj2 <- object@net[[comparison[2]]][[measure]]
    net.diff <- obj2 - obj1
    if (measure%in% c("count","count.merged")) {
      if (is.null(title.name)) {
        title.name = "Differential number of interactions"
      }
    }
    else if (measure%in% c("weight","weight.merged")) {
      if (is.null(title.name)) {
        title.name = "Differential interaction strength"
      }
    }
    legend.name = "Relative values"
  }
  else {
    message("Do heatmap based on a single object \n")
    if (!is.null(signaling)) {
      net.diff <- slot(object, slot.name)$prob[, , signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    }
    else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      }
      else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                          df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
    if (length(idx) > 0) {
      net <- net[-idx, ]
      net <- net[, -idx]
    }
  }
  mat <- net
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat))
  }
  names(color.use) <- colnames(mat)
  if (!is.null(row.show)) {
    mat <- mat[row.show, ]
  }
  if (!is.null(col.show)) {
    mat <- mat[, col.show]
    color.use <- color.use[col.show]
  }
  if (min(mat) < 0) {
    color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)), 
                                   c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                      "\\1", min(mat, na.rm = T))) + 1), 0, round(max(mat, 
                                                                                                                      na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1", 
                                                                                                                                                     max(mat, na.rm = T))) + 1))
  }
  else {
    if (length(color.heatmap) == 3) {
      color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)), 
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 2) {
      color.heatmap.use = colorRamp3(c(min(mat), max(mat)), 
                                     color.heatmap)
    }
    else if (length(color.heatmap) == 1) {
      color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                                 name = color.heatmap))))(100)
    }
    colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*", 
                                                                      "\\1", min(mat, na.rm = T))) + 1), round(max(mat, 
                                                                                                                   na.rm = T), digits = nchar(sub(".*\\.(0*).*", "\\1", 
                                                                                                                                                  max(mat, na.rm = T))) + 1))
  }
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "row", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), 
                                              border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                      show_annotation_name = FALSE)
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), 
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                          show_annotation_name = FALSE)
  if (sum(abs(mat) > 0) == 1) {
    color.heatmap.use = c("white", color.heatmap.use)
  }
  else {
    mat[mat == 0] <- NA
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
                name = legend.name, bottom_annotation = col_annotation, 
                left_annotation = row_annotation, top_annotation = ha2, 
                right_annotation = ha1, cluster_rows = cluster.rows, 
                cluster_columns = cluster.rows, row_names_side = "left", 
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                column_names_gp = gpar(fontsize = font.size), column_title = title.name, 
                column_title_gp = gpar(fontsize = font.size.title), 
                column_names_rot = 90, row_title = "Sources (Sender)", 
                row_title_gp = gpar(fontsize = font.size.title), row_title_rot = 90, 
                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
                                                            fontface = "plain"), title_position = "leftcenter-rot", 
                                            border = NA, legend_height = unit(20, "mm"), labels_gp = gpar(fontsize = 8), 
                                            grid_width = unit(2, "mm")))
  return(ht1)
}

#heatmap
netAnalysis_signalingRole_heatmap.V2 <- function (object, signaling = NULL, 
                                                  pattern = c("outgoing", 
                                                              "incoming", "all"), 
                                                  slot.name = "netP", 
                                                  color.use = NULL, 
                                                  sources.use = NULL,
                                                  color.heatmap = "BuGn", title = NULL, width = 10, height = 8, 
                                                  font.size = 8, font.size.title = 10, cluster.rows = FALSE, 
                                                  cluster.cols = FALSE) 
{
  pattern <- match.arg(pattern)
  if (length(slot(object, slot.name)$centr) == 0) {
    stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
  }
  centr <- slot(object, slot.name)$centr
  if (!is.null(sources.use)){
    outgoing <- matrix(0, nrow = length(object@idents[sources.use]), ncol = length(centr))
    incoming <- matrix(0, nrow = length(object@idents[sources.use]), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents)[sources.use], names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
      outgoing[, i] <- centr[[i]]$outdeg[sources.use]
      incoming[, i] <- centr[[i]]$indeg[sources.use]
    }
  }else{
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
      outgoing[, i] <- centr[[i]]$outdeg
      incoming[, i] <- centr[[i]]$indeg
    }
  }

  if (pattern == "outgoing") {
    mat <- t(outgoing)
    legend.name <- "Outgoing"
  }
  else if (pattern == "incoming") {
    mat <- t(incoming)
    legend.name <- "Incoming"
  }
  else if (pattern == "all") {
    mat <- t(outgoing + incoming)
    legend.name <- "Overall"
  }
  if (is.null(title)) {
    title <- paste0(legend.name, " signaling patterns")
  }
  else {
    title <- paste0(paste0(legend.name, " signaling patterns"), 
                    " - ", title)
  }
  if (!is.null(signaling)) {
    mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
    mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
    idx <- match(rownames(mat1), signaling)
    mat[idx[!is.na(idx)], ] <- mat1
    dimnames(mat) <- list(signaling, colnames(mat1))
  }
  mat.ori <- mat
  mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
  mat[mat == 0] <- NA
 
  if (is.null(color.use)) {
    color.use <- scPalette(length(colnames(mat)))
  }
  color.heatmap.use = (grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, 
                                                                             name = color.heatmap))))(100)
  df <- data.frame(group = colnames(mat))
  rownames(df) <- colnames(mat)
  names(color.use) <- colnames(mat)
  col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                      which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                      simple_anno_size = grid::unit(0.2, "cm"))
  ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), 
                                                  border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                          show_annotation_name = FALSE)
  pSum <- rowSums(mat.ori)
  pSum.original <- pSum
  pSum <- -1/log(pSum)
  pSum[is.na(pSum)] <- 0
  idx1 <- which(is.infinite(pSum) | pSum < 0)
  if (length(idx1) > 0) {
    values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, 
                         length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(1:length(idx1), position)]
  }
  ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), 
                      show_annotation_name = FALSE)
  if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
    legend.break <- max(mat, na.rm = T)
  }
  else {
    legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
                      round(max(mat, na.rm = T), digits = 1))
  }
  ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
                name = "Relative strength", bottom_annotation = col_annotation, 
                top_annotation = ha2, right_annotation = ha1, cluster_rows = cluster.rows, 
                cluster_columns = cluster.rows, row_names_side = "left", 
                row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                column_names_gp = gpar(fontsize = font.size), width = unit(width, 
                                                                           "cm"), height = unit(height, "cm"), column_title = title, 
                column_title_gp = gpar(fontsize = font.size.title), 
                column_names_rot = 90, heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
                                                                                   fontface = "plain"), title_position = "leftcenter-rot", 
                                                                   border = NA, at = legend.break, legend_height = unit(20, 
                                                                                                                        "mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2, 
                                                                                                                                                                                 "mm")))
  return(ht1)
}

######### 15.Line plot 折线图
line_plot = function(data,x ="group" ,y = "Stress.associated",color = NA,
                     text.size = 8,
                     text.angle = 45,
                     text.hjust = 1,
                     legend.position = "none",
                     y.lab="Score",x.lab=NULL,title = NULL,my.test.method = "anova",label.size = 2.8,...){
  if(is.na(color)){
    color.set = c("#0072B5FF","#E18727FF","#BC3C29FF","#7876B1FF","#EE4C97FF","#FFDC91FF","#F39B7FB2","#8491B4B2","#6F99ADFF","#4DBBD5B2")
    color.set = color.set[1:length(unique(data[,x]))]
  }else{
    color.set = color
  }
  # data = data.frame(score = data[,x],group = data[,y])
  summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                        conf.interval=.95, .drop=TRUE,sd =TRUE,sem = FALSE ) {
    library(plyr)
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }
    
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                     c(N    = length2(xx[[col]], na.rm=na.rm),
                       mean = mean   (xx[[col]], na.rm=na.rm),
                       sd   = sd     (xx[[col]], na.rm=na.rm)
                     )
                   },
                   measurevar
    )
    
    # Rename the "mean" column    
    datac <- plyr::rename(datac, c("mean" = measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
  }
  if(T){
    mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                     axis.ticks = element_line(color = "black"),
                     axis.title = element_text(size = text.size,color ="black"), 
                     axis.text = element_text(size=text.size,color = "black"),
                     axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                     #axis.line = element_line(color = "black"),
                     #axis.ticks = element_line(color = "black"),
                     #panel.grid.minor.y = element_blank(),
                     #panel.grid.minor.x = element_blank(),
                     panel.grid=element_blank(), # 去网格线
                     legend.position = legend.position,
                     legend.text = element_text(size= text.size),
                     legend.title= element_text(size= text.size),
                     # panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
                     strip.background = element_rect(color="black",size= 1, linetype="solid") # fill="#FC4E07", 
    )
  }
  Data_summary <- summarySE(data, measurevar=y, groupvars=x)
  head(Data_summary)
  
  p1 = ggplot(data, aes(x = !!as.name(x), y = !!as.name(y), fill=!!as.name(x))) +
    geom_line(data = Data_summary,aes(group = 1)) +
    geom_errorbar(data = Data_summary,aes(ymin = !!as.name(y)-se, ymax= !!as.name(y)+se), 
                  width= 0.2, 
                  position= position_dodge(0.5), 
                  color="black",
                  alpha = 0.8) +
    theme_bw() +
    scale_fill_manual(values = color.set) + 
    geom_jitter( size =2.5,
                 alpha = 1,
                 shape = 21,
                 width = 0.2) +
    # geom_point(aes(x= group, y= !!as.name(y)),pch=19,
    #            position=position_dodge(0.5),size= 1,alpha = 1)+
    geom_point(data = Data_summary,aes(x= group, y= !!as.name(y)),pch=19,
               position=position_dodge(0.5),size= 1,alpha = 1)+ 
    labs(y= y.lab,x= x.lab,title = title) + 
    mytheme + 
    stat_compare_means(aes(group = group),
                       label = "p.format",
                       method = my.test.method,
                       size = label.size,
                       # label.y=label.y,
                       label.x = 1.5,
                       hide.ns = T)
  p1 
}

######### 16. 小提琴图
VlnPlot_3 = function(object,features,
                      group,
                      text.size = 8,
                      text.angle = 45,
                      text.hjust = 1,
                      legend.position = "right",
                      switch = "left",
                      fill.cols = NULL,
                      cols=NULL,
                      widths  = c(3,0.08),
                      heights = c(3,0.08),
                      legend.title = "Ave.exp",
                      x.lab=NULL,y.lab=NULL,title =NULL,...){
  library(tidyverse)
  library(grDevices)
  library(dplyr)
  if(T){
    mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                     axis.ticks = element_line(color = "black"),
                     axis.title = element_text(size = text.size,color ="black"), 
                     axis.text = element_text(size=text.size,color = "black"),
                     axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                     #axis.line = element_line(color = "black"),
                     #axis.ticks = element_line(color = "black"),
                     #panel.grid.minor.y = element_blank(),
                     #panel.grid.minor.x = element_blank(),
                     panel.grid=element_blank(), # 去网格线
                     legend.position = legend.position,
                     legend.text = element_text(size= text.size),
                     legend.title= element_text(size= text.size),
                     # panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
                     strip.background = element_rect(color="black",size= 1, linetype="solid") # fill="#FC4E07", 
    )
  }
  #从Seurat对象中提取细胞注释以及基因表达量
  features = features[features%in%row.names(object)]
  vln.dat=FetchData(object,c(features,group))
  colnames(vln.dat)[length(features)+1] = "celltype.sub"
  # 定义因子顺序，防止画图时对细胞注释进行重排
  # vln.dat$celltype.sub
  vln.dat=vln.dat[order(vln.dat$celltype.sub),]
  vln.dat.melt=vln.dat %>% 
    reshape2::melt(,features) %>%
    dplyr::rename("Gene"="variable") %>%
    dplyr::group_by(celltype.sub,Gene) %>%
    dplyr::mutate(fillcolor=mean(value))
  
  if(switch == "right"){
    switch = "x"
  }else{
    switch = "y"
  }
  if (is.null(fill.cols)) {
    # 小提琴图的填充颜色
    pal=colorRampPalette(RColorBrewer::brewer.pal(n = 9,name="YlOrRd"))(100)
    # 堆积小提琴图
    p1 = ggplot(vln.dat.melt,aes(x=celltype.sub,y=value,fill=fillcolor))+
      # 把小提琴图的外缘轮廓去除
      geom_violin(linetype="blank",scale = "width")+
      scale_fill_gradientn(colors=pal,name=legend.title)+
      facet_grid(Gene~.,switch = switch)+mytheme+
      theme(panel.grid = element_blank(),
            strip.text.y.left = element_text(angle = 0,hjust=1),
            strip.text.y = element_text(angle = 0,hjust=1),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            legend.position = legend.position
            
      )
    # p1
  }else{
    p1 = ggplot(vln.dat.melt,aes(x=celltype.sub,y=value,fill=celltype.sub))+
      # 把小提琴图的外缘轮廓去除
      geom_violin(linetype="blank",scale = "width",aes(fill=celltype.sub))+
      scale_fill_manual(values = fill.cols ,name=legend.title)+
      facet_grid(Gene~.,switch = switch)+mytheme+
      theme(panel.grid = element_blank(),
            strip.text.y.left = element_text(angle = 0,hjust=1),
            strip.text.y = element_text(angle = 0,hjust=1),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_blank(),
            legend.position = legend.position
            
      )
  }
  
  if(switch == "x"){
    p1 = p1 + theme(strip.text.y = element_text(angle = 0,hjust=0))
  }
  
  p1 = p1+labs(x=x.lab,y=y.lab,title =title)
  
  
  # 我们用geom_tile()完成细胞注释，当然也可以向上面一样用geom_segment或者geom_bar来完成
  p3=ggplot(vln.dat%>%
              dplyr::select(celltype.sub)%>%
              unique()%>%
              dplyr::mutate(value='A'),
            aes(x=celltype.sub,y=value,fill=celltype.sub))+
    geom_tile()+
    scale_x_discrete(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0))+
    # facet_grid(.~celltype.sub,scales = "free",space = "free",switch = 'x')+
    mytheme+
    theme(panel.background = element_blank(),
          strip.text = element_text(angle = 0,hjust = 0.5,vjust = 1),
          strip.background = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          # axis.text.x = element_blank(),
          legend.position = "none")
  p3
  if(!is.null(cols)){
    p3 = p3+scale_fill_manual(values = cols)
  }
  (p1+ theme(plot.margin = unit(c(0,30,0,0), "pt"))) + 
    (p3 + theme(plot.margin = unit(c(0,30,0,0), "pt")))+ 
    plot_layout(ncol = 1, widths  = widths,heights = heights)
  # ggsave(filename = "Outplot/Vis/Step3.celltype.EMT.markers.plot.pdf",
  #        height = 10,width = 10,units = "cm")
}

VlnPlot_2 <- function(object,features, cols = NULL,idents = NULL, pt.size = 0,
                      sort = FALSE, assay = NULL, group.by = NULL, split.by = NULL, 
                      adjust = 1, y.max = NULL, same.y.lims = FALSE, log = FALSE, 
                      ncol = NULL,  split.plot = FALSE, stack = FALSE, 
                      combine = TRUE, fill.by = "feature", flip = FALSE, add.noise = TRUE, 
                      raster = NULL,x.lab=NULL,y.lab="Expression level",
                      legend.position = "none",...
){
  library(ggplot2)
  library(ggpubr)
  VlnPlot(object = object, features = features, cols = cols, 
          pt.size = pt.size,idents = idents, 
          sort = sort, assay = assay, group.by = group.by, split.by = split.by, 
          adjust = adjust, y.max = y.max, same.y.lims = same.y.lims, log = log, 
          ncol = ncol, split.plot = split.plot, stack = stack, 
          combine = combine, fill.by = fill.by, flip = flip, add.noise = add.noise, raster = raster,...
  )&
    labs(x=x.lab,y=y.lab)&
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.9)) &
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,position = position_dodge(0.9))&
    mytheme&theme(legend.position = legend.position)
}

#### 17. harmony deBug
RunHarmony.Seurat <- function(
    object,
    group.by.vars,
    reduction = 'pca',
    dims.use = NULL,
    theta = NULL,
    lambda = NULL,
    sigma = 0.1,
    nclust = NULL,
    tau = 0,
    block.size = 0.05,
    max.iter.harmony = 10,
    max.iter.cluster = 20,
    epsilon.cluster = 1e-5,
    epsilon.harmony = 1e-4,
    plot_convergence = FALSE,
    verbose = TRUE,
    reference_values = NULL,
    reduction.save = "harmony",
    assay.use = NULL,
    project.dim = TRUE,
    ...
) {
  if (!requireNamespace('Seurat', quietly = TRUE)) {
    stop("Running Harmony on a Seurat object requires Seurat")
  }
  assay.use <- assay.use %||% Seurat::DefaultAssay(object)
  if (reduction == "pca" && !reduction %in% Seurat::Reductions(object = object)) {
    if (isTRUE(x = verbose)) {
      message("Harmony needs PCA. Trying to run PCA now.")
    }
    object <- tryCatch(
      expr = Seurat::RunPCA(
        object = object,
        assay = assay.use,
        verbose = verbose,
        reduction.name = reduction
      ),
      error = function(...) {
        stop("Harmony needs PCA. Tried to run PCA and failed.")
      }
    )
  }
  if (!reduction %in% Seurat::Reductions(object = object)) {
    stop("Requested dimension reduction is not present in the Seurat object")
  }
  embedding <- Seurat::Embeddings(object, reduction = reduction)
  if (is.null(dims.use)) {
    dims.use <- seq_len(ncol(embedding))
  }
  dims_avail <- seq_len(ncol(embedding))
  if (!all(dims.use %in% dims_avail)) {
    stop("trying to use more dimensions than computed. Rereun dimension reduction
         with more dimensions or run Harmony with fewer dimensions")
  }
  if (length(dims.use) == 1) {
    stop("only specified one dimension in dims.use")
  }
  metavars_df <- Seurat::FetchData(
    object,
    group.by.vars,
    cells = Seurat::Cells(x = object[[reduction]])
  )

  harmonyEmbed <- HarmonyMatrix(
    embedding,
    metavars_df,
    group.by.vars,
    FALSE,
    0,
    theta,
    lambda,
    sigma,
    nclust,
    tau,
    block.size,
    max.iter.harmony,
    max.iter.cluster,
    epsilon.cluster,
    epsilon.harmony,
    plot_convergence,
    FALSE,
    verbose,
    reference_values
  )

  # reduction.key <- Seurat::Key(reduction.save, quiet = TRUE)
  reduction.key = reduction.save
  rownames(harmonyEmbed) <- rownames(embedding)
  colnames(harmonyEmbed) <- paste0(reduction.key, seq_len(ncol(harmonyEmbed)))

  object[[reduction.save]] <- Seurat::CreateDimReducObject(
    embeddings = harmonyEmbed,
    stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)),
    assay = Seurat::DefaultAssay(object = object[[reduction]]),
    key = reduction.key
  )
  if (project.dim) {
    object <- Seurat::ProjectDim(
      object,
      reduction = reduction.save,
      overwrite = TRUE,
      verbose = FALSE
    )
  }
  return(object)
}

#### 18. 相关性分析
batch_cor <- function(gene, exprSet, 
                      cor_gene = NULL,method="sp", n.core=2){
  library(future.apply)
  plan("multisession", workers = n.core)
  plan()
  #设置可用的内存
  options(future.globals.maxSize = 10 * 1024^3)
  y <- as.numeric(exprSet[gene,])
  
  if(is.null(cor_gene)){
    cor_gene <-  row.names(exprSet)
    cor_gene <- cor_gene[!cor_gene%in%gene]
  }
  
  data =  do.call(rbind,future_lapply(cor_gene, function(x){
    dd  <- cor.test(as.numeric(exprSet[x,]),y,method=method)
    data.frame(gene=gene,mRNAs=x,cor=dd$estimate,p.value=dd$p.value )
  }))
  data$FDR = p.adjust(data$p.value,method = "BH")
  return(data)
}

############ CellphoneDB Vis
# mypvals <- read.table("./Step8.CellphoneDB/GSM5573466_Output/pvalues.txt",header = T,sep = "\t",stringsAsFactors = F)
# mymeans <- read.table("./Step8.CellphoneDB/GSM5573466_Output/means.txt",header = T,sep = "\t",stringsAsFactors = F) 

cellphoneDB_Dotplot <- function(pvals.data,
                                means.data,
                                target.cells_1,
                                target.cells_2 = NA,
                                gene_a = NA,
                                gene_b = NA,
                                p.cutoff = 0.05,
                                xlab=NULL,
                                ylab=NULL,
                                title=NULL,
                                text.size = 8,
                                text.angle = 45,
                                text.hjust = 1,
                                text.vjust = 0.5,
                                legend.position = "right",
                                filter.means = 0,
                                ...
){
  library(stringr)
  library(tidyverse)  
  library(dplyr)
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"), 
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust,vjust=text.vjust),
                   #axis.line = element_line(color = "black"),
                   #axis.ticks = element_line(color = "black"),
                   #panel.grid.minor.y = element_blank(),
                   #panel.grid.minor.x = element_blank(),
                   panel.grid=element_blank(), # 去网格线
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size),
                   panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
                   strip.background = element_rect(color="black",size= 1, linetype="solid") # fill="#FC4E07", 
  )
  colnames(pvals.data) = str_replace(string = colnames(pvals.data),
                                     pattern = "\\.",
                                     replacement = "\\_")
  
  colnames(means.data) = str_replace(string = colnames(means.data),
                                     pattern = "\\.",
                                     replacement = "\\_")
  
  # 调整配体-受体顺序，配体在前，受体在后
  order_sequence <- function(df){
    da<-data.frame()
    for(i in 1:length(df$gene_a)){
      sub_data <- df[i,]
      if(sub_data$receptor_b=="False"){
        if(sub_data$receptor_a=="True"){
          old_names<- colnames(sub_data)
          my_list<-strsplit(old_names[-c(1:11)],split="\\_")
          my_character <- paste(sapply(my_list,'[[',2L),
                                sapply(my_list,'[[',1L),sep='_')
          new_names<-c(names(sub_data)[1:4],"gene_b","gene_a","secreted",
                       "receptor_b","receptor_a","annotation_strategy",
                       "is_integrin",my_character)
          sub_data=dplyr::select(sub_data,new_names)
          names(sub_data)<-old_names
          da=rbind(da,sub_data)
        }
      }
      else
      {
        da=rbind(da,sub_data)
      }
    }
    return(da)
  }
  
  # df<-subset(pvals.data,receptor_a=="True"&receptor_b=="False"|receptor_a=="False"&receptor_b=="True")
  # df <- df %>% dplyr::mutate(na_count=rowSums(is.na(df)|df=="")) %>% subset(na_count==0) %>% dplyr::select(-na_count)
  # dim(df)
  # means.data <- order_sequence(df) %>% tidyr::unite(Pairs,gene_a,gene_b)
  # pvals.data <- order_sequence(pvals.data) %>% tidyr::unite(Pairs,gene_a,gene_b)
  
  if(is.na(target.cells_2[1])){
    kp = grepl(pattern = target.cells_1[1], colnames(pvals.data))
    for (i in target.cells_1[-1]) {
      tpm = grepl(pattern = i, colnames(pvals.data))
      kp =  kp|tpm
    }
  }else{
    kp_1 = grepl(pattern = target.cells_1[1], colnames(pvals.data))
    for (i in target.cells_1) {
      tpm = grepl(pattern = i, colnames(pvals.data))
      kp_1 =  kp_1|tpm
    }
    
    kp_2 = grepl(pattern = target.cells_2[1], colnames(pvals.data))
    for (i in target.cells_2) {
      tpm = grepl(pattern = i, colnames(pvals.data))
      kp_2 =  kp_2|tpm
    }
    kp = kp_1 & kp_2
  }
  
  pos = (1:ncol(pvals.data))[kp] 
  choose_pvalues <- pvals.data[,c(c(1,5,6,8,9),pos)]
  choose_means <- means.data[,c(c(1,5,6,8,9),pos)]
  
  if(!is.na(p.cutoff)){
    logi <- apply(choose_pvalues[,5:ncol(choose_pvalues)]<p.cutoff, 1, sum) 
    # 只保留具有细胞特异性的一些相互作用对
    choose_pvalues <- choose_pvalues[logi>=1,]
  }
  
  # 去掉空值
  logi1 <- choose_pvalues$gene_a != ""
  logi2 <- choose_pvalues$gene_b != ""
  logi <- logi1 & logi2
  choose_pvalues <- choose_pvalues[logi,]
  
  # 同样的条件保留choose_means
  choose_means <- choose_means[choose_means$id_cp_interaction %in% 
                                 choose_pvalues$id_cp_interaction,]
  
  # 将choose_pvalues和choose_means数据宽转长
  meansdf <- choose_means %>% reshape2::melt()
  meansdf <- data.frame(interacting_pair = paste0(meansdf$gene_a,"_",meansdf$gene_b),
                        CC = meansdf$variable,
                        means = meansdf$value)
  pvalsdf <- choose_pvalues %>% reshape2::melt()
  pvalsdf <- data.frame(interacting_pair = paste0(pvalsdf$gene_a,"_",pvalsdf$gene_b),
                        CC = pvalsdf$variable,
                        pvals = pvalsdf$value)
  
  # 合并p值和mean文件
  pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
  meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
  pldf <- merge(pvalsdf,meansdf,by = "joinlab")
  colnames(pldf) =   c("joinlab","interacting_pair","Celltypes",
                       "pvals", "interacting_pair.means", 
                       "CC.means", "means")
  
  if(!is.na(gene_a[1])){
    gene_a.df = str_split(pldf$interacting_pair,pattern = "\\_",simplify = T)[,1]
    pldf = pldf[gene_a.df %in%gene_a,]
  }
  
  if(!is.na(gene_b[1])){
    gene_b.df = str_split(pldf$interacting_pair,pattern = "\\_",simplify = T)[,2]
    pldf = pldf[gene_b.df %in%gene_b,]
  }
  # dotplot可视化
  head(pldf)
  
  if(!is.na(p.cutoff)){
    pcc =  pldf%>% filter(means > filter.means) %>% 
      ggplot(aes(Celltypes,interacting_pair) )+ 
      geom_point(aes(color=means,size=-log10(pvals+0.0001))) +
      scale_size_continuous(range = c(1,3),name="-log10(pvals)")+
      scale_color_gradient2(high="red",mid = "yellow",
                            low ="darkblue",
                            midpoint = 1,name="Means")+ 
      theme_bw()+ 
      labs(y= ylab,x= xlab,title = title) + 
      # scale_color_manual(values = rainbow(100))+
      mytheme
  }else{
    pldf$Sig = ifelse(pldf$pvals<0.05,"p<0.05","NS")
    pcc =  pldf %>% filter(means > filter.means) %>% 
      ggplot(aes(Celltypes,interacting_pair))+ 
      geom_point(shape=21,aes(fill=means,
                              size=-log10(pvals+0.0001),
                              color=Sig)) +
      scale_color_manual(values = c("white","black"))+
      scale_fill_gradient(
        low = "#0072B5FF",
        mid = "white",
        high = "#FF7F00",
        midpoint = 0,name="Means"
      )+ 
      scale_size(name = "-log10(FDR)",range = c(1,4))+
      theme_bw()+ 
      labs(y= ylab,x= xlab,title = title) + 
      # scale_color_manual(values = rainbow(100))+
      mytheme
  }
  return(pcc)
}


#### 19. pseudoBulk -------------------------
#### pseudoBulk
pseudoBulk = function(seurat.data,
                      group,
                      celltype = NA,...
){
  seurat.data.list = SplitObject(seurat.data,group)
  
  if(is.na(celltype)){
    res.pseudo = lapply(seurat.data.list, function(tem.data){
      sampleID = tem.data@meta.data[,group] %>% unique()
      print(paste("Start:",sampleID,sep = " "))
      
      count.data = tem.data@assays$RNA@counts %>% rowSums() %>% as.data.frame()
      colnames(count.data) = sampleID
      return(count.data)
    })
    res.pseudo = do.call("cbind",res.pseudo)
  }else{
    res.pseudo = lapply(seurat.data.list, function(tem.data){
      sampleID = tem.data@meta.data[,group] %>% unique()
      print(paste("Start:",sampleID,sep = " "))
      tem.data.list = SplitObject(tem.data,celltype)
      
      tem.pseudo = lapply(tem.data.list, function(cell_type){
        count.data = cell_type@assays$RNA@counts %>% rowSums() %>% as.data.frame()
        colnames(count.data) = unique(cell_type@meta.data[,celltype])
        return(count.data)
      })
      tem.pseudo = do.call("cbind",tem.pseudo)
      colnames(tem.pseudo) = colnames(tem.pseudo)
      return(tem.pseudo)
    })
    res.pseudo = do.call("cbind",res.pseudo)
  }
  print("Finish!")
  return(res.pseudo)
}


###### DESeq2
get.DESeq2 = function(input.data, #输入目标表达矩阵
                      group_list, #设置分组，因子数据
                      contrast = NA,
                      NA.remove = TRUE, #DESeq2会产生NA值， NA.remove = TRUE移除NA
                      ...){
  # 加载包
  library(DESeq2)
  if(is.na(contrast)){
    contrast = rev(levels(group_list))
  }
  # 第一步，构建DESeq2的DESeq对象
  colData <- data.frame(row.names=colnames(input.data),group_list=group_list)
  dds <- DESeqDataSetFromMatrix(countData = input.data,colData = colData,design = ~ group_list)
  
  # 第二步，进行差异表达分析
  dds2 <- DESeq(dds)
  
  # 提取差异分析结果，trt组对untrt组的差异分析结果
  table(group_list)
  tmp <- results(dds2,contrast=c("group_list",contrast))
  # tmp <- results(dds2)
  DEG_DESeq2 <- as.data.frame(tmp[order(tmp$padj),])
  head(DEG_DESeq2)
  
  if(NA.remove){
    # 去除差异分析结果中包含NA值的行
    DEG_DESeq2 = na.omit(DEG_DESeq2)
  }
  DEG_DESeq2$Gene = row.names(DEG_DESeq2)
  return(DEG_DESeq2)
}

###### MetaCells
MetaCells <- function(count.data,
                      meta.data, 
                      cell.type, 
                      N=20, seed=1, 
                      return.seurat = FALSE) {
  cell.types <- table(meta.data[[cell.type]])
  ## drop the cell states less than N cells
  cell.types <- cell.types[cell.types >= N]
  cell.types <- sort(cell.types, decreasing = TRUE)
  new.data <- lapply(names(cell.types), function(cc) {
    ## get cells in same cell types and tissue (cell state)
    new.meta.data <- subset(meta.data, get(cell.type) == cc)
    ## set the seed to make the shuffle function reproducible
    set.seed(seed)
    new.meta.data <- new.meta.data[permute::shuffle(rownames(new.meta.data)), ]
    ## generate PoolID for pooling
    new.meta.data$PoolID <- floor(0:(nrow(new.meta.data)-1) / N)
    ## pooling cells with the same PoolID
    pool.ids <- table(new.meta.data$PoolID)
    pool.ids <- pool.ids[pool.ids >= N]
    pool.ids <- sort(names(pool.ids))
    ## rowSums the expression profiles
    pseudoBulk.data <- lapply(pool.ids, function(x) {
      select.cells <- rownames(subset(new.meta.data, PoolID == x))
      new.data <- count.data[, select.cells]
      return(rowSums(new.data))
    })
    pseudoBulk.data <- do.call(cbind, pseudoBulk.data)
    colnames(pseudoBulk.data) <- paste(cc, pool.ids, sep = ".")
    return(pseudoBulk.data)
  })
  ## merge the expression profiles for each cell state
  new.data <- do.call(cbind, new.data)
  
  if(!return.seurat){
    return(new.data)
  }else{
    seurat.data = CreateSeuratObject(counts = new.data)
    return(seurat.data)
  }
}

#### 20.Scissor
## 出图函数
scissor_barplot <- function(seurat.data,
                            celltype,
                            scissor_name = "scissor",
                            pseudo_count = 0.001,
                            width = 0.8,
                            text.size = 8,
                            text.angle = 45,
                            text.hjust = 1,
                            legend.position = "right",
                            kepp_zero = TRUE, legend.key.size = 0.5,
                            color_fill = "0",
                            color_low = 'royalblue',color_high = 'indianred1',
                            x.lab = NULL, y.lab = "log2_ratio", title = "Scissor"){
  library(aplot)
  library(ggplotify)
  library(data.table)
  library(cowplot)
  
  if(T){
    mytheme <- theme_bw () + theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                                   axis.ticks = element_line(color = "black"),
                                   axis.title = element_text(size = text.size,color ="black"), 
                                   axis.text = element_text(size=text.size,color = "black"),
                                   # axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                                   #axis.line = element_line(color = "black"),
                                   #axis.ticks = element_line(color = "black"),
                                   #panel.grid.minor.y = element_blank(),
                                   #panel.grid.minor.x = element_blank(),
                                   panel.grid=element_blank(), # 去网格线
                                   legend.position = legend.position,
                                   legend.text = element_text(size= text.size),
                                   legend.title= element_text(size= text.size)
                                   # panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")
                                   # strip.background = element_rect(color="black",size= 1, linetype="solid") # fill="#FC4E07", 
    )
  }
  
  dat.plot = seurat.data@meta.data
  dat.plot$Celltype = dat.plot[,celltype]
  dat.plot$scissor = dat.plot[,scissor_name]
  dat.plot = data.table(dat.plot)
  dat.plot <- dat.plot[,.(N=.N),by=c("scissor","Celltype")]
  dat.plot
  
  sum.data <- dat.plot[,{.(NTotal = sum(.SD$N))},by=c("Celltype")]
  
  suppressMessages({dat.plot = left_join(dat.plot,sum.data)})
  # dat.plot
  
  dat.plot$freq = dat.plot$N / dat.plot$NTotal *100
  # dat.plot
  
  ##  1. "Scissor+","Scissor-"
  plot_Data = filter(dat.plot,scissor %in% c("Scissor+","Scissor-")) %>% 
    arrange(desc(freq))
  plot_Data$log2freq = log2(plot_Data$freq)
  head(plot_Data)
  
  library(tidyr)
  plot_Data = as.data.frame(plot_Data) %>% 
    pivot_wider(., 
                id_cols = c("Celltype"),         # 列标识符（保持不变的列）
                names_from = "scissor",      # 用于生成新列名的列
                values_from = "freq"     # 包含数据值的列
    )
  plot_Data[is.na(plot_Data)] = 0
  
  # Step 1: 计算对数比例
  plot_Data$log_ratio <- log2((plot_Data$`Scissor+`+pseudo_count) /( plot_Data$`Scissor-`+pseudo_count)) 

  if(!kepp_zero){
    (plot_Data = plot_Data %>%   arrange(desc(log_ratio)))
    plot_Data$Celltype = factor(plot_Data$Celltype,levels = rev(plot_Data$Celltype))
  }else{
    plot_Data2 = data.frame(Celltype = unique(dat.plot$Celltype))
    
    suppressMessages({plot_Data2 = left_join(plot_Data2,plot_Data)})
    plot_Data2[is.na(plot_Data2)] = 0
    
    (plot_Data2 = plot_Data2 %>% arrange(desc(log_ratio)))
    plot_Data2$Celltype = factor(plot_Data2$Celltype,levels = rev(plot_Data2$Celltype))
    plot_Data = plot_Data2
  }
  
  # 绘制柱状图
  p1 <- ggplot(plot_Data, aes(x = Celltype, y = log_ratio, fill = log_ratio)) +
    geom_bar(stat = "identity",width = width) +
    labs(x = x.lab, y = y.lab, title = title) +
    mytheme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = legend.position
    )  + theme(legend.key.size = unit(legend.key.size, "cm")) # 使x轴标签倾斜显示
    
  ## 2. Scissor+
  plot_Data = filter(dat.plot,scissor %in% c("Scissor+")) %>% 
    arrange(desc(freq)) %>% as.data.frame
  
  if(!kepp_zero){
    (plot_Data = plot_Data %>% arrange(desc(freq)))
    plot_Data$Celltype = factor(plot_Data$Celltype,levels = rev(plot_Data$Celltype))
  }else{
    plot_Data2 = data.frame(scissor = "Scissor+",
                            Celltype = unique(dat.plot$Celltype))
    suppressMessages({plot_Data2 = left_join(plot_Data2,plot_Data)})
    plot_Data2$freq[is.na(plot_Data2$freq)]  = 0
    (plot_Data2 = plot_Data2 %>% arrange(desc(freq)))
    plot_Data2$Celltype = factor(plot_Data2$Celltype,levels = rev(plot_Data2$Celltype))
    plot_Data = plot_Data2
    ct.leves = levels(plot_Data$Celltype)
  }
  
  # 绘制柱状图
  p2 <- ggplot(plot_Data, aes(x = Celltype, y = freq, fill = freq)) +
    geom_bar(stat = "identity",width = width) +
    labs(x = x.lab, y = "Cellular fraction (%)", title = "Scissor+ cells") +
    # scale_fill_gradient(low = color_low,high = color_high)+
    mytheme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = legend.position
    ) + theme(legend.key.size = unit(legend.key.size, "cm"))  # 使x轴标签倾斜显示
  
  ## 3. Scissor-
  plot_Data = filter(dat.plot,scissor %in% c("Scissor-")) %>% 
    arrange(desc(freq)) %>% as.data.frame
  
  plot_Data2 = data.frame(scissor = "Scissor-",
                          Celltype = unique(dat.plot$Celltype))
  
  if(!kepp_zero){
    (plot_Data = plot_Data %>% arrange(desc(freq)))
    plot_Data$Celltype = factor(plot_Data$Celltype,levels = rev(plot_Data$Celltype))
  }else{
    plot_Data2 = data.frame(scissor = "Scissor-",
                            Celltype = unique(dat.plot$Celltype))
    suppressMessages({plot_Data2 = left_join(plot_Data2,plot_Data)})
    plot_Data2$freq[is.na(plot_Data2$freq)]  = 0
    (plot_Data2 = plot_Data2 %>% arrange(desc(freq)))
    plot_Data2$Celltype = factor(plot_Data2$Celltype,levels = ct.leves)
    plot_Data = plot_Data2
  }
  
  # 绘制柱状图
  p3 <- ggplot(plot_Data, aes(x = Celltype, y = freq, fill = freq)) +
    geom_bar(stat = "identity",width = width) +
    labs(x = x.lab, y = "Cellular fraction (%)", title = "Scissor- cells") +
    mytheme +
    # scale_fill_gradient(low = color_low,high = color_high)+ 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = legend.position
    ) + theme(legend.key.size = unit(legend.key.size, "cm")) # 使x轴标签倾斜显示
  
  
  if(color_fill == "1"){
    temp <- RColorBrewer::brewer.pal(11, "Spectral")
    temp[6] <- "gold"
    rbPal <- grDevices::colorRampPalette(temp)
    
    p1 = p1 + scale_fill_gradientn(name = "Log2_ratio", 
                                   colours = rev(rbPal(50)), 
                                   guide = ggplot2::guide_colourbar(ticks.colour = "black", 
                                                                    ticks.linewidth = 0.5, frame.colour = "black"))
    
    p2 = p2 + scale_fill_gradientn(name = "Freq", 
                                   colours = rev(rbPal(50)), 
                                   guide = ggplot2::guide_colourbar(ticks.colour = "black", 
                                                                    ticks.linewidth = 0.5, frame.colour = "black"))
    
    p3 = p3 + scale_fill_gradientn(name = "Freq", 
                                   colours = rev(rbPal(50)), 
                                   guide = ggplot2::guide_colourbar(ticks.colour = "black", 
                                                                    ticks.linewidth = 0.5, frame.colour = "black"))
  }else if(color_fill == "0"){
    p1 = p1 + scale_fill_gradient(low = color_low,high = color_high)
    p2 = p2 + scale_fill_gradient(low = color_low,high = color_high)
    p3 = p3 + scale_fill_gradient(low = color_low,high = color_high)
  }else{
    p1 = p1 + scale_fill_gradientn(name = "Freq", 
                                   colours = viridis::viridis(7), 
                                   guide = ggplot2::guide_colourbar(ticks.colour = "black", 
                                                                    ticks.linewidth = 0.5, frame.colour = "black"))
    p2 = p2 + scale_fill_gradientn(name = "Freq", 
                                   colours = viridis::viridis(7), 
                                   guide = ggplot2::guide_colourbar(ticks.colour = "black", 
                                                                    ticks.linewidth = 0.5, frame.colour = "black"))
    p3 = p3 + scale_fill_gradientn(name = "Freq", 
                                   colours = viridis::viridis(7), 
                                   guide = ggplot2::guide_colourbar(ticks.colour = "black", 
                                                                    ticks.linewidth = 0.5, frame.colour = "black"))
  }
  
  return(list("Both" = p1,
              "Postive" = p2,
              "Negative" = p3))
}

#### 21. pySCENIC Rank plot
PlotRegulonRank <- function(rssMat, cell.type, topn=5, point.size = 1,front.size = 8) {
  data <- data.frame(
    Regulons = 1:nrow(rssMat),
    RSS = sort(rssMat[, cell.type], decreasing = T),
    label = sub("(+)", "", names(sort(rssMat[, cell.type], decreasing = T)), fixed = T)
  )
  
  data$pt.col <- ifelse(data$Regulons <= topn, "#007D9B", "#BECEE3")
  data <- head(data, n=200)
  data.label <- head(data, n=topn)
  
  ggplot(data, aes(Regulons, RSS)) + 
    geom_point(size=point.size, color=data$pt.col) + 
    ggrepel::geom_text_repel(inherit.aes = FALSE, data = data.label, aes(Regulons, RSS, label=label), size=2.8) + 
    ggtitle(cell.type) + ylab("Specificity score") + 
    theme_bw(base_size = front.size) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          axis.ticks = element_line(color="black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = .5)
    )
}

#### 22. 去批次      
scRemoveBatch = function(seurat.data, 
                         batchID = "orig.ident", 
                         methods = "harmony", n.pcs = 30){
  
  library(Seurat)
  
  ### 0. 不整合
  if(methods %in% c("none")){
    seurat_int <- seurat.data %>% NormalizeData(verbose = F) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>% 
      ScaleData(verbose = F) %>%
      RunPCA(npcs = n.pcs, verbose = F)%>% 
      RunUMAP(reduction = "pca", dims = 1:20, verbose = F)
  }
  
  ######多种算法整合去批次
  ### 1. harmony
  if(methods %in% c("harmony")){
    library(harmony)
    seurat.data <- seurat.data %>% NormalizeData(verbose = F) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>% 
      ScaleData(verbose = F) %>%
      RunPCA(npcs = n.pcs, verbose = F)
    
    seurat.data <- seurat.data %>% RunHarmony(batchID, plot_convergence = T)
    
    seurat_int <- seurat.data %>% 
      RunUMAP(reduction = "harmony", dims = 1:20, verbose = F)
  }
  
  ### 2. CCA
  if(methods %in% c("CCA")){
    seurat_list = SplitObject(seurat.data, split.by = batchID)
    rm(seurat.data)
    
    # normalize and identify variable features for each dataset independently
    seurat_list <- lapply(X = seurat_list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
      return(x)
    })
    
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = seurat_list)
    
    #找到整合的锚点anchors
    seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, 
                                             anchor.features = features, dims = 1:30)
    rm(seurat_list)
    
    # 利用anchors进行整合
    seurat_int <- IntegrateData(anchorset = seurat_anchors, dims = 1:n.pcs) %>% 
      ScaleData() %>% 
      RunPCA(npcs = n.pcs, verbose = FALSE) %>% 
      RunUMAP(reduction = "pca", dims = 1:n.pcs)
  }
  
  ### 3. RPCA
  if(methods %in% c("RPCA")){
    seurat_list = SplitObject(seurat.data, split.by = batchID)
    rm(seurat.data)
    
    # normalize and identify variable features for each dataset independently
    seurat_list <- lapply(X = seurat_list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
      return(x)
    })
    
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = seurat_list)
    seurat_list <- lapply(X = seurat_list, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE) %>% 
        RunPCA( features = features, verbose = FALSE)
      return(x)
    })
    
    #Perform integration
    seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, 
                                             anchor.features = features, 
                                             reduction = "rpca", 
                                             k.anchor = 20 # 5 by default. Increasing this parameter to 20 will assist in aligning these populations.
    )
    rm(seurat_list)
    
    seurat_int <- IntegrateData(anchorset = seurat_anchors, dims = 1:n.pcs) %>% 
      ScaleData() %>% 
      RunPCA(npcs = n.pcs, verbose = FALSE) %>% 
      RunUMAP(reduction = "pca", dims = 1:n.pcs)
  }
  
  ### 4. fastMNN
  if(methods %in% c("fastMNN")){
    library(SeuratWrappers)
    seurat_list = SplitObject(seurat.data, split.by = batchID)
    rm(seurat.data)
    
    seurat_list <- lapply(seurat_list, FUN = function(x){
      x = x %>% NormalizeData() %>% 
        FindVariableFeatures()
      return(x)
    })
    
    seurat_int <- RunFastMNN(object.list = seurat_list) %>% 
      RunUMAP(reduction = "mnn", dims = 1:n.pcs)
  }
  
  ### 5. LIGER
  if(methods %in% c("LIGER")){
    library(rliger)
    seurat_int <- seurat.data %>% NormalizeData(verbose = F) %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>% 
      ScaleData(verbose = F, do.center = F) %>%
      RunPCA(npcs = n.pcs, verbose = F)
    rm(seurat.data)
    
    #k是矩阵分解的因子数，一般取值20-40
    seurat_int <- RunOptimizeALS(seurat_int, k = 30, lambda = 5, split.by = batchID) %>% 
      RunQuantileNorm(split.by = batchID)
    
    #### 4.降维聚类：进行 UMAP 和 clustering:
    n.pcs = ncol(seurat_int[["iNMF"]])
    seurat_int <- seurat_int %>% 
      RunUMAP(reduction = "iNMF", dims = 1:n.pcs, verbose = F)
  }
  return(seurat_int)
}