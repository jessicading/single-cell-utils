library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
library(scales)
library(viridis)
library(varhandle)
library(Seurat)
library(VennDiagram)
library(RColorBrewer)
library(pheatmap)
library(grid)
library(ggpubr)
library(ggfittext)
library(stringr)
library(ggrepel)
library(ComplexHeatmap)
source("~/Desktop/Yang_Lab/Mergeomics.R")
source("~/Desktop/Yang_Lab/Resources/genesets/R-tomfunctions.R")

# meta_data_1 will be in the columns
cellTypeCount <- function(seuratObject, meta_data_1 = "Cell_type", meta_data_2 = "orig.ident", includePerc = TRUE){
  Idents(seuratObject) <- meta_data_1
  idents_table <- as.data.frame(table(Idents(seuratObject), seuratObject@meta.data[,meta_data_2]))
  idents_table$id <- rep(levels(Idents(seuratObject)),times = length(unique(seuratObject@meta.data[,meta_data_2])))
  idents_table <- dcast(data = idents_table, formula = id~Var2, fun.aggregate = sum,value.var = "Freq")
  if(includePerc){
    idents_table_prop <- as.data.frame(prop.table(table(Idents(seuratObject), seuratObject@meta.data[,meta_data_2]), margin=2))
    idents_table_prop$id <- rep(levels(Idents(seuratObject)),times = length(unique(seuratObject@meta.data[,meta_data_2])))
    idents_table_prop <- dcast(data = idents_table_prop, formula = id~Var2, fun.aggregate = sum,value.var = "Freq")
    
    total_idents = data.frame("Cell_Type" = idents_table$id, stringsAsFactors = FALSE)
    for(i in 2:(length(unique(seuratObject@meta.data[,meta_data_2]))+1)){
      total_idents[,colnames(idents_table)[i]] = idents_table[,i]
      total_idents[,paste0(colnames(idents_table_prop)[i], "_perc")] = idents_table_prop[,i]
    }
    total_idents = rbind(total_idents, c("Totals", colSums(total_idents[,2:ncol(total_idents)])))
    for(j in seq(from = 3, to = (length(unique(seuratObject@meta.data[,meta_data_2]))*2)+1, by = 2)){
      total_idents[,j] = percent(as.numeric(total_idents[,j]), accuracy = 0.01)
    }
    colnames(total_idents)[1] <- meta_data_1
    return(total_idents)
  }
  else{
    colnames(idents_table)[1] <- meta_data_1
    return(idents_table)
  }
}

# receive dataframe from cellTypeCount with just counts
# meta_data_2 is same as meta_data_2 from cellTypeCount but can be named anything
normalizedDistributionProportionPlot <- function(data, 
                                                 meta_data_2 = "Sample",
                                                 addCounts = TRUE, 
                                                 cts_include = NULL, 
                                                 custom_colors = NULL,
                                                 coord_flip = TRUE,
                                                 meta_data_2_levels = NULL){
  # melt data
  data <- readRDS("./Potential_Paper_Plots/HP_CellTypeCounts.rds")
  data <- readRDS("./Potential_Paper_Plots/HYP_CellTypeCounts.rds")
  # data <- as.data.frame(t(data))
  # Samples <- rownames(data)[2:length(rownames(data))]
  # celltypes <- as.character(data["Cell_typev2",])
  # data <- data[-1,]
  # colnames(data) <- celltypes
  # data$Sample <- Samples
  # data <- data[,c(ncol(data), 1:(ncol(data)-1))]
  melted <- melt(data, id.vars = colnames(data)[1])
  melted = melted[,c(2,1,3)]
  #colnames(melted) <- c(colnames(data)[1], meta_data_2,"Counts")
  colnames(melted) <- c("Sample", meta_data_2,"Counts")
  melted$Sample <- gsub("TGW2","5XFAD1", melted$Sample)
  melted$Sample <- gsub("TGW3","5XFAD2", melted$Sample)
  melted$Sample <- gsub("TGW5","5XFAD3", melted$Sample)
  melted$Sample <- gsub("TGW6","5XFAD4", melted$Sample)
  melted$Sample <- gsub("WTW","WT", melted$Sample)
  melted$Cell_typev2 <- gsub("\\<Endothelium\\>", "Endothelium1",melted$Cell_typev2)
  melted$Cell_typev2 <- gsub("Endothelium 2", "Endothelium2",melted$Cell_typev2)
  melted$Cell_typev2 <- gsub("\\<OPC\\>", "OPC1",melted$Cell_typev2)
  
  melted$Counts <- as.numeric(melted$Counts) 
  # add percentage as feature
  # percentages = melted$value[grepl("_perc", melted$variable)]
  # melted = melted[!grepl("_perc", melted$variable),]
  # melted$Percentage = percentages
  # melted = melted[!grepl("Totals", melted$Cell_Type),]
  percentages = c()
  # for(samp in unique(melted[,meta_data_2])){
  #   total = sum(melted$Counts[melted[,meta_data_2]==samp])
  #   for(ct in unique(melted[,colnames(melted)[1]])){
  #     percentages = append(percentages, melted$Counts[melted[,meta_data_2]==samp & melted[,colnames(melted)[1]]==ct]/total)
  #   }
  # }
  melted = melted[order(melted$Cell_typev2),]
  for(samp in unique(melted[,"Cell_typev2"])){
    total = sum(melted$Counts[melted[,"Cell_typev2"]==samp])
    for(ct in unique(melted[,"Sample"])){
      percentages = append(percentages, melted$Counts[melted[,"Cell_typev2"]==samp & melted[,"Sample"]==ct]/total)
    }
  }
  melted$Percentage = percentages*100
  
  # order by most abundant cell type
  #     add up cell type counts across all samples
  count_totals = c()
  for(ct in unique(melted[,colnames(data)[1]])){
    count_totals = append(count_totals, sum(as.numeric(melted$Counts[melted[,colnames(data)[1]]==ct])))
  }
  cell_type_order = unique(melted[,colnames(data)[1]])[order(count_totals, decreasing = TRUE)]
  
  melted[,colnames(data)[1]] <- factor(melted[,colnames(data)[1]], levels = rev(cell_type_order))
  
  if(addCounts){
    df <- melted %>%
      group_by_at(meta_data_2) %>%
      arrange(!!! rlang::syms(meta_data_2), desc(!!! rlang::syms(colnames(data)[1]))) %>%
      mutate(lab_ypos = cumsum(Percentage) - 0.5 * Percentage)
    # df <- melted %>%
    #   group_by(`Cell Type`) %>%
    #   arrange(`Cell Type`, desc(condition)) %>%
    #   mutate(lab_ypos = cumsum(Percentage) - 0.5 * Percentage)
    
    df$Counts[df$Percentage<1] <- ""
    
    if(!is.null(cts_include)){
      cts_lose_text = setdiff(df[,colnames(data)[1]], cts_include)
      for(iter in 1:nrow(df)){
        if(sum(df[,colnames(data)[1]][iter]==cts_lose_text)>0){
          df$Counts[iter] <- ""
        }
      }
    }
    
    if(!is.null(meta_data_2_levels)){
      df[,meta_data_2] <- factor(df, levels = meta_data_2_levels)
    }

    # pal = c("5XFAD"="#F8766D",
    #         "WT"="#00BFC4")
    
    g <- ggplot(df, aes_string(x=paste0("`",meta_data_2,"`"), y="Percentage")) + geom_col(aes_string(fill=colnames(data)[1])) +
      geom_text(aes(y=lab_ypos, label=Counts, group = colnames(data)[1]), color = "white") + coord_flip() + 
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), legend.position = "right",
            axis.text.y = element_text(margin = margin(r=5), size=14, face = "bold"),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            #legend.title = element_text(face = "bold", size=13)) +
            legend.title = element_blank()) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 102)) #+
    # scale_fill_manual(
    #   values = pal,
    #   limits = names(pal))
  } 
  else{
    # pal = c("5XFAD"="salmon",
    #         "WT"="cyan",
    #         "Hippocampus"="salmon",
    #         "Hypothalamus"="cyan")
    pal = c("5XFAD1"="#800026",
            "5XFAD2"="#BD0026",
            "5XFAD3"="#E31A1C",
            "5XFAD4"="#FC4E2A",
            "WT1"="#1A1A1A",
            "WT2"="#4D4D4D",
            "WT3"="#878787")
    celltypes <- c("Microglia","Astrocyte","Ependyma","OPC1","Oligodendrocyte","Endothelium1","Neuron","Pericyte")
    celltypes <- c("Microglia","Astrocyte","Ependyma","OPC1","Oligodendrocyte","Endothelium1","Neuron","Pericyte","Tanycyte")
    melted = melted[melted$Cell_typev2 %in% celltypes,]
    melted$Cell_typev2 <- factor(melted$Cell_typev2, levels=celltypes)
    # df = melted
    # df$Counts[df$Percentage<1] <- ""
    g <- ggplot(melted, aes_string(x=paste0("`",meta_data_2,"`"), y="Percentage")) + geom_col(aes_string(fill=colnames(melted)[1])) +
      #coord_flip() + 
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(),
            legend.position = "right",
            #axis.text.y = element_text(margin = margin(r=5), size=14, face = "bold"),
            axis.text.x = element_text(angle=60, vjust=1, hjust=1, size=12),
            #axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
            #legend.title = element_text(face = "bold", size=13)) +
            #legend.title = element_blank()) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 102)) +
      guides(fill=guide_legend(title="Sample")) +
      xlab("Cell Type") +
      scale_fill_manual(
        values = pal,
        limits = names(pal))
    pdf("./Potential_Paper_Plots/Figure1/Sample_Cell_Type_Proportions_HYP.pdf", height=4, width = 4.3)
    g
    dev.off()
  }
  if(coord_flip){
    g <- g + coord_flip()
  }
  results <- list(g, df)
  return(results)
}

# receive dataframe from cellTypeCount with just counts
# meta_data_2 is same as meta_data_2 from cellTypeCount but can be named anything
normalizedDistributionProportionPlotTwoTissues <- function(data, 
                                                           meta_data_2 = "Sample",
                                                           addCounts = TRUE, 
                                                           cts_include = NULL, 
                                                           custom_colors = NULL,
                                                           coord_flip = TRUE,
                                                           meta_data_2_levels = NULL){
  # melt data
  melted <- melt(data, id.vars = colnames(data)[1])
  colnames(melted) <- c(colnames(data)[1], meta_data_2,"Counts")
  melted$Counts <- as.numeric(melted$Counts) 
  # add percentage as feature
  # percentages = melted$value[grepl("_perc", melted$variable)]
  # melted = melted[!grepl("_perc", melted$variable),]
  # melted$Percentage = percentages
  # melted = melted[!grepl("Totals", melted$Cell_Type),]
  percentages = c()
  for(samp in unique(melted[,meta_data_2])){
    total = sum(melted$Counts[melted[,meta_data_2]==samp])
    for(ct in unique(melted[,colnames(data)[1]])){
      percentages = append(percentages, melted$Counts[melted[,meta_data_2]==samp & melted[,colnames(data)[1]]==ct]/total)
    }
  }
  melted$Percentage = percentages*100
  
  # order by most abundant cell type
  #     add up cell type counts across all samples
  count_totals = c()
  for(ct in unique(melted[,colnames(data)[1]])){
    count_totals = append(count_totals, sum(as.numeric(melted$Counts[melted[,colnames(data)[1]]==ct])))
  }
  cell_type_order = unique(melted[,colnames(data)[1]])[order(count_totals, decreasing = TRUE)]
  
  melted[,colnames(data)[1]] <- factor(melted[,colnames(data)[1]], levels = rev(cell_type_order))
  
  if(addCounts){
    df <- melted %>%
      group_by_at(meta_data_2) %>%
      arrange(!!! rlang::syms(meta_data_2), desc(!!! rlang::syms(colnames(data)[1]))) %>%
      mutate(lab_ypos = cumsum(Percentage) - 0.5 * Percentage)
    # df <- melted %>%
    #   group_by(`Cell Type`) %>%
    #   arrange(`Cell Type`, desc(condition)) %>%
    #   mutate(lab_ypos = cumsum(Percentage) - 0.5 * Percentage)
    
    df$Counts[df$Percentage<1] <- ""
    
    if(!is.null(cts_include)){
      cts_lose_text = setdiff(df[,colnames(data)[1]], cts_include)
      for(iter in 1:nrow(df)){
        if(sum(df[,colnames(data)[1]][iter]==cts_lose_text)>0){
          df$Counts[iter] <- ""
        }
      }
    }
    
    if(!is.null(meta_data_2_levels)){
      df[,meta_data_2] <- factor(df, levels = meta_data_2_levels)
    }
    
    g <- ggplot(df, aes_string(x=paste0("`",meta_data_2,"`"), y="Percentage")) + geom_col(aes_string(fill=colnames(data)[1])) +
      geom_text(aes(y=lab_ypos, label=Counts, group = colnames(data)[1]), color = "white") + coord_flip() + 
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), legend.position = "top",
            axis.text.y = element_text(margin = margin(r=5), size=14, face = "bold"),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            #legend.title = element_text(face = "bold", size=13)) +
            legend.title = element_blank()) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 102)) #+
    # scale_fill_manual(
    #   values = pal,
    #   limits = names(pal))
  }
  else{
    df = melted
    df$Counts[df$Percentage<1] <- ""
    g <- ggplot(df, aes_string(x=paste0("`",meta_data_2,"`"), y="Percentage")) + geom_col(aes_string(fill=colnames(data)[1])) +
      coord_flip() + 
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), legend.position = "top",
            axis.text.y = element_text(margin = margin(r=5), size=14, face = "bold"),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            #legend.title = element_text(face = "bold", size=13)) +
            legend.title = element_blank()) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 102))
  }
  if(coord_flip){
    g <- g + coord_flip()
  }
  results <- list(g, df)
  return(results)
}

# the point of this function is to look quickly at average expression for the individual samples to see if there is a 
# dominating sample that is driving differential gene expression, for example
# "ident_meta" can also be any metadata
drawAvgGeneExpHeatMap <- function(seuratObject, genes = seuratObject@assays$RNA@var.features[1:50], 
                                  cell_type_meta = "Cell_type", cell_types = NULL, 
                                  ident_meta = "orig.ident", idents = NULL){
  if(!is.null(cell_types)){
    Idents(seuratObject) <- cell_type_meta
    seuratObject = subset(seuratObject, idents = cell_types)
  }
  if(!is.null(idents)){
    Idents(seuratObject) <- ident_meta
    seuratObject = subset(seuratObject, idents = idents)
  }
  
  avg_expr_data = data.frame(stringsAsFactors = FALSE)
  if(length(unique(seuratObject@meta.data[,cell_type_meta]))>1 & length(unique(seuratObject@meta.data[,ident_meta]))>1){
    for(ct in unique(seuratObject@meta.data[,cell_type_meta])){
      for(ident in unique(seuratObject@meta.data[,ident_meta])){
        for(gene in genes){
          temp = data.frame(cell_type_meta = ct,
                            ident_meta = ident,
                            "GENE" = gene,
                            avg_exp = mean(seuratObject@assays$RNA@data[gene,][seuratObject@meta.data[,cell_type_meta]==ct &
                                                                                 seuratObject@meta.data[,ident_meta]==ident]))
          avg_expr_data = rbind(avg_expr_data, temp)
        }
      }
    }
    colnames(avg_expr_data) = c(cell_type_meta, ident_meta, "GENE","AvgExpr")
    gene_heat <- ggplot(avg_expr_data, aes(x = GENE, y = ident_meta, fill=Expression)) +
      geom_tile()+
      scale_fill_gradient2(low="blue",mid="white",high="red") +
      theme_bw() +
      theme(panel.border = element_blank(), 
            axis.title.x = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(size = 12,face = "bold", angle = 90, vjust = .5,hjust = 0),
            panel.grid = element_blank(),
            axis.title.y = element_blank(),
            legend.title = element_text(size=9),
            legend.position = "bottom") +
      scale_x_discrete(position = "top") +
      facet_grid(Cell_type + orig.ident~., scales = "free", space = "free") + #facet by group
      theme(strip.background = element_blank(), #remove background for facet labels
            panel.border = element_rect(colour = "white", fill = NA), #add black border
            panel.spacing.y = unit(0, "lines"), 
            strip.placement = "outside",
            strip.text = element_text(angle = 180, face = "bold", margin = margin(t=0,b=0,r=0,l=0),
                                      hjust = .5, size = 12)) 
  }
  else if(length(unique(seuratObject@meta.data[,cell_type_meta]))>1){ # one cell type, multiple identities
    gene_heat <- ggplot(avg_expr_data, aes_string(x = GENE, y = ident_meta, fill=Expression)) +
      geom_tile()+
      scale_fill_gradient2(low="blue",mid="white",high="red") +
      theme_bw() +
      theme(panel.border = element_blank(), 
            axis.title.x = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_text(size = 12, face = "bold"),
            axis.text.x = element_text(size = 12,face = "bold", angle = 90, vjust = .5,hjust = 0),
            panel.grid = element_blank(),
            axis.title.y = element_blank(),
            legend.title = element_text(size=9),
            legend.position = "bottom") +
      scale_x_discrete(position = "top")
  }
  
  
  # Idents(seuratObject) <- ident_meta
  # if(length(unique(seuratObject@meta.data[,cell_type_meta]))>1 & length(unique(seuratObject@meta.data[,ident_meta]))>1){
  #   avg_expr <- AverageExpression(object = seuratObject, features = genes, return.seurat = TRUE, add.ident = cell_type_meta)
  # }
  # else{
  #   avg_expr <- AverageExpression(object = seuratObject, features = genes, return.seurat = TRUE)
  # }
  # 
  # avg_expr_data <- avg_expr@assays$RNA@scale.data
  # avg_expr_data = setNames(melt(avg_expr_data), c("GENE","Condition","Expression"))
  # avg_expr_data$Condition <- as.character(avg_expr_data$Condition)
  
  # if(length(unique(seuratObject@meta.data[,cell_type_meta]))>1 & length(unique(seuratObject@meta.data[,ident_meta]))>1){
  #   origidents = c()
  #   cts = c()
  #   for(row in 1:nrow(avg_expr_data)){
  #     origidents[row] = unlist(strsplit(avg_expr_data$Condition[row], split = "_"))[1]
  #     cts[row] = unlist(strsplit(avg_expr_data$Condition[row], split = "_"))[2]
  #   }
  #   avg_expr_data$orig.ident = origidents
  #   avg_expr_data$Cell_type = cts
  # }

  # else{
  #   gene_heat <- ggplot(avg_expr_data, aes(x = GENE, y = Condition, fill=Expression)) +
  #     geom_tile()+
  #     scale_fill_gradient2(low="blue",mid="white",high="red") +
  #     theme_bw() +
  #     theme(panel.border = element_blank(), 
  #           axis.title.x = element_blank(),
  #           axis.ticks = element_blank(),
  #           axis.text = element_text(size = 12, face = "bold"),
  #           axis.text.x = element_text(size = 12,face = "bold", angle = 90, vjust = .5,hjust = 0),
  #           panel.grid = element_blank(),
  #           axis.title.y = element_blank(),
  #           legend.title = element_text(size=9),
  #           legend.position = "bottom") +
  #     scale_x_discrete(position = "top")
  # }
  
  return(gene_heat)
  
}

# must make column names the same
# if comparing mouse with human, put mouse FIRST (DEG_result.1)
# convertGenes is a vector with the same length to convert the mouse genes to human genes or vice versa
compareSecondStudy <- function(DEG_result.1, 
                               DEG_result.2, 
                               DEG_result.3 = NULL,
                               study_names, 
                               cell_type = NULL, 
                               convertGenes = FALSE, # boolean vector
                               FDR_threshold = 0.05,
                               logFC_threshold = .1,
                               file_output_name = "Comparison.png"){
  # make melted dataframe 
  DEG_result.1 = DEG_result.1[abs(DEG_result.1$avg_logFC)>.1 & DEG_result.1$p_val_adj<0.05,]
  DEG_result.2 = DEG_result.2[abs(DEG_result.2$avg_logFC)>.1 & DEG_result.2$p_val_adj<0.05,]
  
  if(convertHuman){
    # convertedToHuman <- convertMouseGeneList(DEG_result.1$GENE)
    # DEG_result.1$HUMAN <- convertedToHuman$HGNC.symbol[match(DEG_result.1$GENE, convertedToHuman$MGI.symbol)]
    # DEG_result.1$HUMAN[which(is.na(DEG_result.1$HUMAN))] <- toupper(DEG_result.1$GENE[which(is.na(DEG_result.1$HUMAN))])
    # DEG_result.1$GENE <- DEG_result.1$HUMAN
    # 
    # convertedToHuman <- convertMouseGeneList(DEG_result.3$GENE)
    # DEG_result.3$HUMAN <- convertedToHuman$HGNC.symbol[match(DEG_result.3$GENE, convertedToHuman$MGI.symbol)]
    # DEG_result.3$HUMAN[which(is.na(DEG_result.3$HUMAN))] <- toupper(DEG_result.3$GENE[which(is.na(DEG_result.3$HUMAN))])
    # DEG_result.3$GENE <- DEG_result.3$HUMAN
    
  }
  
  # for(set in c(DEG_result.1, DEG_result.2, DEG_result.3)){
  #   # add negative for negative logFC
  #   gene_info = c()
  #   for(iter in 1:length(set$GENE)){
  #     if(set$avg_logFC[iter]<0){
  #       gene_info[iter] = concatenate("-",set$GENE[iter])
  #     }
  #     else next
  #   }
  #   set$gene_info = gene_info
  # }
  
  gene_info = c()
  for(iter in 1:length(DEG_result.1$GENE)){
    if(DEG_result.1$avg_logFC[iter]<0){
      gene_info[iter] = concatenate(c("-",DEG_result.1$GENE[iter]), mysep = "")
    }
    else{
      gene_info[iter] = DEG_result.1$GENE[iter]
    }
  }
  DEG_result.1$gene_info = gene_info
  
  gene_info = c()
  for(iter in 1:length(DEG_result.2$GENE)){
    if(DEG_result.2$avg_logFC[iter]<0){
      gene_info[iter] = concatenate(c("-",DEG_result.2$GENE[iter]), mysep = "")
    }
    else{
      gene_info[iter] = DEG_result.2$GENE[iter]
    }
  }
  DEG_result.2$gene_info = gene_info
  
  gene_info = c()
  for(iter in 1:length(DEG_result.3$GENE)){
    if(DEG_result.3$avg_logFC[iter]<0){
      gene_info[iter] = concatenate(c("-",DEG_result.3$GENE[iter]), mysep = "")
    }
    else{
      gene_info[iter] = DEG_result.3$GENE[iter]
    }
  }
  DEG_result.3$gene_info = gene_info
  
  venn.diagram(
    x = list(DEG_result.1$gene_info, DEG_result.2$gene_info, DEG_result.3$gene_info),
    category.names = study_names,
    filename = file_output_name,
    output = TRUE,
    
    # Output features
    imagetype="png" ,
    height = 480 , 
    width = 600 , 
    resolution = 500,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = myCol,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
  )
}

# list must have names if genes are to be converted to human genes
# important that strings are not factors!!
compareStudiesVennDiagram <- function(study_list,
                               study_names, 
                               convertMouseToHumanGenes = NULL, # boolean vector, yes if mouse to human
                               FDR_threshold = 0.05,
                               logFC_threshold = .1,
                               addDirection = FALSE,
                               file_output_name = "Comparison.png"){
  # make melted dataframe
  for(iter in 1:length(study_list)){
    study_list[[iter]] = study_list[[iter]][abs(study_list[[iter]]$avg_logFC)>logFC_threshold & 
                                              study_list[[iter]]$p_val_adj<FDR_threshold,]
  }
  
  if(!is.null(convertMouseToHumanGenes)){
    studiesToConvert = names(study_list)[convertMouseToHumanGenes]
    for(study in studiesToConvert){
      study_list[[study]] <- convertDfGeneColumnMouseHuman(df = study_list[[study]], toSpecies = "human", forPathway = FALSE)
    }
  }

  if(addDirection){
    for(iter in 1:length(study_list)){
      gene_info = c()
      for(gene in 1:length(study_list[[iter]]$GENE)){
        if(study_list[[iter]]$avg_logFC[gene]<0){
          gene_info[gene] = paste0("-",study_list[[iter]]$GENE[gene])
        }
        else{
          gene_info[gene] = study_list[[iter]]$GENE[gene]
        }
      }
      study_list[[iter]]$gene_info = gene_info
    }
  }
  
  if(addDirection){
    final_list = list()
    for(iter in 1:length(study_list)){
      final_list[[iter]] = study_list[[iter]]$gene_info
    }
  }
  else{
    final_list = list()
    for(iter in 1:length(study_list)){
      final_list[[iter]] = study_list[[iter]]$GENE
    }
  }
  
  venn.diagram(
    x = final_list,
    #x = list(study_list[[1]]$gene_info, study_list[[2]]$gene_info),
    category.names = study_names,
    filename = file_output_name,
    output = TRUE,
    
    # Output features
    imagetype="png" ,
    height = 480 ,
    width = 480 ,
    resolution = 5000,
    compression = "lzw",

    # Circles
    #lwd = 1,
    lty = 'blank',
    fill = c("#87a2b5", "#94db82"),
    scaled = FALSE,

    # Numbers
    cex = .05,
    fontface = "bold",
    fontfamily = "sans",

    # Set names
    cat.cex = 0.02,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-15, 15),
    cat.dist = c(0.03, 0.03),
    cat.fontfamily = "sans"
    #rotation = 1
  )
  
}



# looks for each cell type
# use for only two comparisons
compareStudiesTableDEGs <- function(combined_DEGs, comparison = "Tissue", 
                                    FDR_threshold = 0.05, lfc_threshold = 0.1){
  combined_DEGs = combined_DEGs[combined_DEGs$p_val_adj<FDR_threshold & abs(combined_DEGs$avg_logFC)>lfc_threshold,]
  combined_DEGs = combined_DEGs[order(abs(combined_DEGs$avg_logFC), decreasing = TRUE),] # so that large effect ones show up first
  combined_DEGs$GENE_Direction = ifelse(combined_DEGs$avg_logFC>0, 
                                        paste0(combined_DEGs$GENE,"_UP"),
                                        paste0(combined_DEGs$GENE,"_DOWN"))
  result = data.frame()
  melted_result = data.frame()
  for(ct in unique(combined_DEGs$Cell_type)){
    group_genes = list()
    group_genes_direction = list()
    temp = data.frame("Cell_type" = ct, stringsAsFactors = FALSE)
    for(el in unique(combined_DEGs[,comparison])){
      group_genes[[el]] = combined_DEGs$GENE[combined_DEGs$Cell_type==ct & combined_DEGs[,comparison]==el]
      group_genes_direction[[el]] = combined_DEGs$GENE_Direction[combined_DEGs$Cell_type==ct & combined_DEGs[,comparison]==el]
      temp[,el] = paste0("(",length(group_genes[[el]]),") ", concatenate(group_genes[[el]], mysep = ", "))
    }
    specific1 = setdiff(group_genes[[1]], group_genes[[2]])
    temp[,paste0("Specific_", unique(combined_DEGs[,comparison])[1])] = paste0("(",length(specific1),") ",
                                                                               concatenate(specific1, mysep = ", "))
    if(length(specific1)>1){
      melt_temp = data.frame("Cell_type" = ct, 
                             "Type"=paste0("Specific_", unique(combined_DEGs[,comparison])[1]), 
                             "GENE"= specific1, stringsAsFactors = FALSE)
      melted_result = rbind(melted_result, melt_temp)
    }
    specific2 = setdiff(group_genes[[2]], group_genes[[1]])
    temp[,paste0("Specific_", unique(combined_DEGs[,comparison])[2])] = paste0("(",length(specific2),") ",
                                                                               concatenate(specific2, mysep = ", "))
    if(length(specific2)>1){
      melt_temp = data.frame("Cell_type" = ct, 
                             "Type"=paste0("Specific_", unique(combined_DEGs[,comparison])[2]), 
                             "GENE"= specific2, stringsAsFactors = FALSE)
      melted_result = rbind(melted_result, melt_temp)
    }
    shared = intersect(group_genes[[1]], group_genes[[2]])
    shared_direction = intersect(group_genes_direction[[1]], group_genes_direction[[2]])
    shared_direction_stripped <- gsub("_DOWN","",shared_direction)
    shared_direction_stripped <- gsub("_UP","",shared_direction_stripped)
    shared_discordant = setdiff(shared, shared_direction_stripped)
    temp[,paste0("Shared")] = paste0("(",length(setdiff(shared, shared_discordant)),") ",
                                     concatenate(setdiff(shared, shared_discordant), mysep = ", "))
    temp[,"Shared_Direction_Discordant"] = paste0("(",length(shared_discordant),") ",
                                        concatenate(shared_discordant, mysep = ", "))
    if(length(shared)>1){
      melt_temp = data.frame("Cell_type" = ct, 
                             "Type"= "Shared", 
                             "GENE"= setdiff(shared, shared_discordant), 
                             stringsAsFactors = FALSE)
      melted_result = rbind(melted_result, melt_temp)
    }
    if(length(shared_discordant)>1){
      melt_temp = data.frame("Cell_type" = ct, 
                             "Type"= "Shared_discordant", 
                             "GENE"= shared_discordant, 
                             stringsAsFactors = FALSE)
      melted_result = rbind(melted_result, melt_temp)
    }
    result = rbind(result, temp)
  }
  return(list(result, melted_result))
}

constructHierarchicalDEGPlot <- function(
  DEG_df, 
  genes=NULL, 
  ct_specific = FALSE, 
  numGenes = 30, 
  ct_list=NULL, 
  show_significance = FALSE,
  orderByPval = FALSE,
  annotation.bar = TRUE
) {
  if(!is.null(ct_list)){
    selected = data.frame()
    for(ct in ct_list){
      selected = rbind(selected, DEG_df[DEG_df$Cell_type==ct,])
    }
    DEG_df = selected
    rm(selected)
  }
  else{
    ct_list = unique(DEG_df$Cell_type)
  }
  if(is.null(genes)){
    # get top 20 DEGs and then 3 cell-type specific DEGs if any
    if(!orderByPval){
      DEG_df = DEG_df[order(abs(DEG_df$avg_logFC), decreasing = TRUE),]
    }
    else{
      DEG_df = DEG_df[order(DEG_df$p_val_adj, decreasing = TRUE),]
    }
    genes = c()
    for(g in DEG_df$GENE){
      genes = append(genes, g)
      genes = unique(genes)
      if(length(genes)==numGenes) break
    }
    if(ct_specific){
      specific_genes = list()
      for(ct in unique(DEG_df$Cell_type)){
        specific_genes[[ct]] = c("init")
      }
      sig_DEG_df <- DEG_df[DEG_df$p_val_adj<0.05 & abs(DEG_df$avg_logFC>0.25),]
      for(gene in unique(sig_DEG_df$GENE)){
        if(sum(sig_DEG_df$GENE==gene)==1){
          if(length(specific_genes[[sig_DEG_df$Cell_type[sig_DEG_df$GENE==gene]]])<4){
            specific_genes[[sig_DEG_df$Cell_type[sig_DEG_df$GENE==gene]]] = append(specific_genes[[sig_DEG_df$Cell_type[sig_DEG_df$GENE==gene]]],
                                                                                   gene)
          }
        }
      }
      rm(sig_DEG_df)
      
      show_specific_genes = c()
      for(ct in names(specific_genes)){
        show_specific_genes = c(show_specific_genes, specific_genes[[ct]])
      }
      show_specific_genes = show_specific_genes[show_specific_genes!="init"]
      genes = append(genes, show_specific_genes)
    }
    
    selected = data.frame()
    for(g in genes){
      selected = rbind(selected, DEG_df[DEG_df$GENE==g,])
    }
  }
  else{
    selected = data.frame()
    for(g in genes){
      selected = rbind(selected, DEG_df[DEG_df$GENE==g,])
    }
  }
  rm(DEG_df)
  
  selected = selected[,c("Cell_type","GENE","avg_logFC","p_val_adj")]
  
  # if empty DEG values, fill with zeroes
  if((length(unique(selected$Cell_type))*length(genes))!=nrow(selected)){
    for(gene in unique(selected$GENE)){
      if(sum(selected$GENE==gene)!=length(selected$Cell_type)){
        # get cts that are missing the value
        cts_missing = setdiff(unique(selected$Cell_type), selected$Cell_type[selected$GENE==gene])
        for(ct in cts_missing){
          selected = rbind(selected, data.frame("Cell_type"=ct, "GENE"=gene, "avg_logFC"=0, 
                                                "p_val_adj" = 1,
                                                stringsAsFactors = FALSE))
        }
      }
      else next
    }
  }
  
  # convert DEG_df to a matrix
  DEG_matrix <- matrix(nrow = length(genes), ncol = length(ct_list))
  colnames(DEG_matrix) = ct_list
  rownames(DEG_matrix) = genes
  for(ct in ct_list){
    for(gene in genes){
      DEG_matrix[gene,ct] = selected$avg_logFC[selected$Cell_type==ct & 
                                              selected$GENE==gene]
    }
  }
  
  # https://jcoliver.github.io/learn-r/008-ggplot-dendrograms-and-heatmaps.html
  DEG.dendro <- as.dendrogram(hclust(d = dist(DEG_matrix)))
  dendro.plot <- ggdendrogram(data = DEG.dendro, rotate = TRUE)
  
  selected$Cell_type <- factor(selected$Cell_type, levels = ct_list)
  
  heatmap <- ggplot(data = selected, aes(x = Cell_type, y = GENE, fill = avg_logFC)) +
    geom_tile() +
    scale_fill_gradient2(low="dodgerblue1",high="red", mid = "white") +
    theme_bw() +
    scale_x_discrete(position = "top")
  if(annotation.bar){
    pbuild <- ggplot_build(plot = heatmap)
    y.range = diff(x = pbuild$layout$panel_params[[1]]$y.range)
    y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range
    y.max <- y.pos + 0.02 * y.range
    cols = matrix(hue_pal(h=c(0,360))(length(ct_list)), nrow = 1)
    names(cols) = levels(selected$Cell_type)
    heatmap <- heatmap + annotation_raster(raster = (x = cols),
                                           xmin = -Inf,
                                           xmax = Inf, 
                                           ymin = y.pos, 
                                           ymax = y.max) +
      coord_cartesian(ylim = c(0, y.max), clip = 'off') +
      scale_color_manual(values = cols)
  }
  
  metadata <- data.frame(ct_list, row.names = ct_list)
  colnames(metadata) = "Cell Type"
  
  colfunc <- colorRampPalette(c("blue","white","red"))
  heat <- pheatmap(mat = DEG_matrix, 
                   annotation_col = metadata, 
                   cluster_cols = FALSE, 
                   color = colfunc(10))
    

}

cellTypeSpecificDEGs <- function(DEG_df){
  
}

# same number of genes for each plot
# Use GetAssayData - use two group bars
# marker_genes is a list for each of the cts_show
# Seurat object has the marker genes in the misc slot
markerGenePlot <- function(
  seuratObject, 
  cts_show, 
  subsample = FALSE, 
  subsample_size = 50, 
  marker_genes, 
  nMarkers = 10,
  metadata = c("Condition","Tissue", "Cell_type"),
  metadata_display_name = c("Condition","Tissue", "Cell Type")
){
  if(is.null(marker_genes)){
    markers <- seuratObject@misc$markers
    markers <- markers[markers$avg_logFC>0 & markers$p_val_adj<0.05,]
    # get top 10 from each
    marker_genes = list()
    all_marker = c()
    for(ct in cts_show){
      marker_genes[[ct]] = markers$gene[markers$cluster==ct][1:nMarkers]
      all_marker = append(all_marker, markers$gene[markers$cluster==ct][1:nMarkers])
    }
    dup_genes = all_marker[duplicated(all_marker)]
    if(length(dup_genes)>0){
      for(dgene in dup_genes){
        first_occurence = FALSE
        if(!first_occurence){
          for(ct in names(marker_genes)){
            if(!first_occurence){
              if(dgene %in% marker_genes[[ct]]){
                  first_occurence = TRUE
                }
              else next
              }
          }
        }
        else{
          if(dgene %in% marker_genes[[ct]]){
            marker_genes[[ct]] = marker_genes[[ct]][marker_genes[[ct]]!=dgene]
          }
          else next
        }
      }
    }
    for(ct in names(marker_genes)){
      marker_genes[[ct]] = marker_genes[[ct]][!is.na(marker_genes[[ct]])]
    }
    all_marker = unique(all_marker)
    all_marker = all_marker[!is.na(all_marker)]
  }
  else{
    all_marker = c()
    for(it in 1:length(marker_genes)){
      all_marker = append(all_marker, marker_genes[[it]])
    }
  }
  
  # extract scaled data 
  for(gene in all_marker){
    if(!(gene %in% rownames(seuratObject@assays$RNA@scale.data))){
      cat(gene, " not in scaled data\n")
    }
  }
  
  meta.data <- seuratObject@meta.data[,c("condition","Tissue","Cell_typev2")]

  data <- seuratObject@assays$RNA@scale.data[all_marker,]
  
  
  meta.data$Condition <- as.character(meta.data$condition)
  meta.data$Condition <- gsub("WTW", "WT",meta.data$Condition)
  meta.data$Condition <- gsub("TGW", "5XFAD",meta.data$Condition)
  
  meta.data$condition <- NULL
  meta.data$`Cell type` <- as.character(meta.data$Cell_typev2)
  meta.data$Cell_typev2 <- NULL
  ordered_data = data.frame(stringsAsFactors = FALSE)
  first = FALSE
  for(ct in cts_show){
    for(tissue in c("Hippocampus","Hypothalamus"))
      for(cond in c("WT","5XFAD")){
        # get cells that meeting position
        cells = rownames(meta.data)[meta.data$`Cell type`==ct & 
                                      meta.data$Condition==cond &
                                      meta.data$Tissue==tissue]
        if(!first){
          ordered_data <- rbind(ordered_data, data[,cells])
          first = TRUE
        }
        else{
          ordered_data <- cbind(ordered_data, data[,cells])
        }
      }
  }
  ordered_data <- as.matrix(ordered_data)
  rm(data)
  
  colfunc <- colorRampPalette(c("blue", "white","red"))
  heat <- pheatmap(mat = ordered_data, 
                   annotation_col = meta.data, 
                   cluster_cols = FALSE, 
                   color = magma(n = 50), 
                   show_colnames = FALSE, scale = "row")
  return(heat)
}

# use runCellTypeDEGsMain_mod to be inclusive of more DEGs (no logFC threshold)
# report "lost" DEGs
# report offending sample - sample with least effect
# two outputs - (1) - numbers table
#               (2) - for every gene that was taken out, reason and offending sample
compareGlobalAndSampleWiseDEGs <- function(global_DEGs, sample_wise_DEGs, 
                                           lfc_cutoffs = c(0.1,0.25), 
                                           pval_adj_cutoffs = c(0.05,0.01),
                                           pvalType = "minPadj",
                                           conditions = c("5XFAD", "WT"),
                                           nTotalSamples = 6){
  quantified = data.frame(stringsAsFactors = FALSE)
  # geneDetails = data.frame(stringsAsFactors = FALSE)
  sample_wise_DEGs_melt = data.frame(stringsAsFactors = FALSE)
  names(sample_wise_DEGs[[1]]) <- gsub(" ", "_", names(sample_wise_DEGs[[1]]))
  cts_leave_out = c()
  for(iter in 1:length(names(sample_wise_DEGs[[1]]))){
    temp = as.data.frame(sample_wise_DEGs[[1]][iter])
    if(length(colnames(temp))<(nTotalSamples*3+5)){
      cts_leave_out = append(cts_leave_out, iter)
      next
    }
    temp$Cell_type = names(sample_wise_DEGs[[1]])[iter]
    temp$GENE = rownames(temp)
    colnames(temp) <- gsub(paste0(names(sample_wise_DEGs[[1]])[iter],"."), "", colnames(temp))
    sample_wise_DEGs_melt = rbind(sample_wise_DEGs_melt, temp)
  }
  if(length(cts_leave_out)>0){
    final_celltypes = names(sample_wise_DEGs[[1]])[-cts_leave_out]
  }
  else{
    final_celltypes = names(sample_wise_DEGs[[1]])
  }
  for(ct in final_celltypes){
    for(lfc in lfc_cutoffs){
      for(fdr in pval_adj_cutoffs){
        GenesGlobal = global_DEGs$GENE[global_DEGs$Cell_type==gsub("_"," ",ct) & 
                                          abs(global_DEGs$avg_logFC)>lfc &
                                          global_DEGs$p_val_adj<fdr]
        GenesSampleWise = sample_wise_DEGs_melt$GENE[sample_wise_DEGs_melt$Cell_type==ct &
                                                       sample_wise_DEGs_melt[,pvalType]<fdr &
                                                       rowSums(abs(sample_wise_DEGs_melt[,grep("avg_logFC", 
                                                                                               colnames(sample_wise_DEGs_melt))])>lfc)==nTotalSamples]
        temp_q = data.frame("Cell Type" = ct, 
                            "lfc_threshold" = lfc,
                            "FDR_threshold" = fdr,
                            "nGeneGlobal" = length(GenesGlobal),
                            "nGeneSampleWise" = length(GenesSampleWise),
                            "ConservedGenes" = paste0("(",length(intersect(GenesSampleWise, GenesGlobal)), ") ", 
                                                      concatenate(intersect(GenesSampleWise, GenesGlobal), mysep = ", ")),
                            "LostGlobalGenes" = paste0("(",length(setdiff(GenesGlobal, GenesSampleWise)), ") ", 
                                                       concatenate(setdiff(GenesGlobal, GenesSampleWise), mysep = ", ")),
                            "OnlySampleWise" = paste0("(",length(setdiff(GenesSampleWise, GenesGlobal)), ") ", 
                                                      concatenate(setdiff(GenesSampleWise, GenesGlobal), mysep = ", "))
        )
        quantified <- rbind(quantified, temp_q)
      }
    }
  }
  colnames(quantified) = gsub("_", " ", colnames(quantified))
  colnames(quantified) = gsub(".", " ", colnames(quantified), fixed = TRUE)
  return(quantified)
}

# Due to constraints of the manipulations in this function, it will get rid of
# underscores in cell type names
# filters based on the meta p value type. Many DEGs may be significant by the meta p-value
# but not by the normal global DEG analysis. To filter by more "consistent" DEGs, use the
# sampleWiseFilteredDEGsConsistent function.
sampleWiseFilteredDEGs <- function(global_DEGs, sample_wise_DEGs,
                                   lfc_threshold = .1,FDR_threshold=0.05,
                                   pval_type = "ConsP", nTotalSamples = 6){
  # make reference for "actually" significant genes
  sample_wise_DEGs_melt = data.frame(stringsAsFactors = FALSE)
  names(sample_wise_DEGs[[1]]) <- gsub(" ", "_", names(sample_wise_DEGs[[1]]))
  for(iter in 1:length(names(sample_wise_DEGs[[1]]))){
    temp = as.data.frame(sample_wise_DEGs[[1]][iter])
    if(length(colnames(temp))<(nTotalSamples*3+5)){
      next # should not consider because not enough cells in a sample for this cell type
    }
    temp$Cell_type = names(sample_wise_DEGs[[1]])[iter]
    temp$GENE = rownames(temp)
    colnames(temp) <- gsub(paste0(names(sample_wise_DEGs[[1]])[iter],"."), "", colnames(temp))
    sample_wise_DEGs_melt = rbind(sample_wise_DEGs_melt, temp)
  }
  sample_wise_DEGs_melt = sample_wise_DEGs_melt[sample_wise_DEGs_melt[,pval_type]<FDR_threshold &
                                                  rowSums(abs(sample_wise_DEGs_melt[,grep("avg_logFC", 
                                                                                          colnames(sample_wise_DEGs_melt))])>lfc_threshold)==nTotalSamples,
                                                c("Cell_type","GENE", pval_type)]

  # filter global DEGs for only those included in the sample_wise_DEGs_melt
  
  sample_wise_DEGs_melt$Cell_type <- gsub("_", " ", sample_wise_DEGs_melt$Cell_type)
  filtered_DEGs <- merge(sample_wise_DEGs_melt, global_DEGs, by = c("Cell_type","GENE"))
  
  return(filtered_DEGs)
}

# not based on the meta p type. based on the p adj for each sample DEG run
sampleWiseFilteredDEGsConsistent <- function(global_DEGs, sample_wise_DEGs, 
                                             nTotalSamples=6, nSampleSignificance=6, 
                                             lfc_threshold=.1, pval_type = "ConsP"){
  sample_wise_DEGs_melt = data.frame(stringsAsFactors = FALSE)
  names(sample_wise_DEGs[[1]]) <- gsub(" ", "_", names(sample_wise_DEGs[[1]]))
  for(iter in 1:length(names(sample_wise_DEGs[[1]]))){
    temp = as.data.frame(sample_wise_DEGs[[1]][iter])
    if(length(colnames(temp))<(nTotalSamples*3+5)){
      next # should not consider because not enough cells in a sample for this cell type
    }
    temp$Cell_type = names(sample_wise_DEGs[[1]])[iter]
    temp$GENE = rownames(temp)
    colnames(temp) <- gsub(paste0(names(sample_wise_DEGs[[1]])[iter],"."), "", colnames(temp))
    sample_wise_DEGs_melt = rbind(sample_wise_DEGs_melt, temp)
  }
  sample_wise_DEGs_melt = sample_wise_DEGs_melt[rowSums(abs(sample_wise_DEGs_melt[,grep("avg_logFC", 
                                                                                        colnames(sample_wise_DEGs_melt))])>lfc_threshold)==nTotalSamples,]
  sample_wise_DEGs_melt = sample_wise_DEGs_melt[rowSums(sample_wise_DEGs_melt[,grep("p_val_adj", 
                                                                                      colnames(sample_wise_DEGs_melt))]<0.05)==nSampleSignificance,
                                                c("Cell_type","GENE", pval_type)]
  
  sample_wise_DEGs_melt$Cell_type <- gsub("_", " ", sample_wise_DEGs_melt$Cell_type)
  filtered_DEGs <- merge(sample_wise_DEGs_melt, global_DEGs, by = c("Cell_type","GENE"))
  return(filtered_DEGs)
}




# must have "GENE" column
# single-cell RNA-seq gene name convention applies 
#       -so mitochondrial genes are "mt-"/"MT-"
# if converting mouse genes for pathway enrichment, set forPathway = TRUE
#       This is because pathway databases will represent MT-ND1 as ND1, for example
convertDfGeneColumnMouseHuman <- function(df, toSpecies="human", forPathway=FALSE){
  MT_gene_names = c("MT-ND1"="ND1",  "MT-ND2"="ND2",  "MT-CO1"="COX1",  "MT-CO2"="COX2",  "MT-ATP8"="ATP8", 
                    "MT-ATP6"="ATP6", "MT-CO3"="COX3","MT-ND3" ="ND3", "MT-ND4L"="ND4L", "MT-ND4"="ND4",
                    "MT-ND5"="ND5",  "MT-ND6"="ND6",  "MT-CYB"="CYB")
  if(toSpecies=="mouse"){
    convertedtoMouse <- convertMouseGeneList(df$GENE)
    df$MOUSE = convertedtoMouse$MGI.symbol[match(df$GENE, convertedtoMouse$HGNC.symbol)]
    df$MOUSE[which(is.na(df$MOUSE))] <- tolower(df$GENE[which(is.na(df$MOUSE))])
    df$GENE <- df$MOUSE
    df$GENE <- paste0(toupper(x=substr(df$GENE, start = 1, stop = 1)),tolower(substring(df$GENE, first = 2)))
    # correct for mt- genes
    MT_mouse_gene_names = MT_gene_names
    new_names = c()
    for(name in names(MT_gene_names)){
      new_names = append(new_names, 
                         paste0(tolower(unlist(strsplit(name, split = "-"))[1]), "-",
                                paste0(toupper(x=substr(unlist(strsplit(name, split = "-"))[2], start = 1, stop = 1)),
                                       tolower(substring(unlist(strsplit(name, split = "-"))[2], first = 2)))))
    }
    names(MT_mouse_gene_names) = new_names
    
    for(gene in 1:nrow(df)){
      if(sum(MT_mouse_gene_names==df$gene[gene])>0){
        df$GENE[gene] = names(MT_mouse_gene_names)[MT_mouse_gene_names==df$gene[gene]]
      }
      else{
        next
      }
    }
    return(df)
  }
  else if(toSpecies=="human"){
    convertedToHuman <- convertMouseGeneList(df$GENE)
    df$HUMAN <- convertedToHuman$HGNC.symbol[match(df$GENE, convertedToHuman$MGI.symbol)]
    df$HUMAN[which(is.na(df$HUMAN))] <- toupper(df$GENE[which(is.na(df$HUMAN))])
    df$GENE <- df$HUMAN
    df$HUMAN <- NULL
    if(sum(is.na(df$GENE))>0) cat("Warning: NAs created.\n")
    if(forPathway){ # convert MT- genes
      for(gene in 1:nrow(df)){
        if(grepl("^MT-",df$GENE[gene])){
          df$GENE[gene] = MT_gene_names[df$GENE[gene]]
        }
        else{
          next
        }
      }
    }
    return(df)
  }
  else cat("Put a valid toSpecies value - 'human' or 'mouse'\n")
}

convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , 
                   mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  # Print the first 6 genes found to the screen
  print(head(genesV2))
  return(genesV2)
}

pathwayEnrichment <- function(modules_list, resources_path, output_Dir){
  # pVal this is our method of switching off between FDR and pval threshold
  list_of_pathway_databases = list.files(path = resources_path, pattern = "*.txt", full.names = TRUE)
  ifelse(!dir.exists(file.path(output_Dir)), dir.create(file.path(output_Dir),recursive = T), FALSE)
  result = list()
  all_result = list()
  description_file = data.frame(stringsAsFactors = FALSE)
  annotations = c()
  for(module in unique(modules_list$MODULE)){
    cat(module, "\n")
    deg_list = modules_list[modules_list$MODULE==module,] # why is this creating NAs?
    deg_list = deg_list[!is.na(deg_list$MODULE),]
    deg_list[["gene"]] <- deg_list$GENE
    deg_list[["module"]] <- deg_list$MODULE
    deg_list <- deg_list[,c("module","gene")]
    unique_module1 <- unique(deg_list$module)
    module1_len <- length(unique(deg_list$module))
    x <- list()
    for(z in 1:length(list_of_pathway_databases)){
      pathway_database <- list_of_pathway_databases[z]
      database_name = unlist(strsplit(pathway_database,"/"))
      database_name = database_name[length(database_name)]
      database_name = unlist(strsplit(database_name, ".txt"))[1]
      print(database_name)
      
      Module2 <- tool.read(pathway_database)
      Unique_module2 <- unique(Module2$module)
      
      Module2_len <- length(unique(Module2$module))
      
      data_matrix_for_enrichment <- data.frame()
      List_initial <- 1
      # go through the different databases
      for(k in 1:Module2_len){
        data_matrix_for_enrichment[List_initial,1] <- unique_module1
        data_matrix_for_enrichment[List_initial,2] <- length(deg_list$gene[which(match(deg_list$module, unique_module1)>0)])
        
        Overlapped_genes <- intersect(deg_list$gene[which(match(deg_list$module, unique_module1)>0)],
                                      Module2$gene[which(match(Module2$module, Unique_module2[k])>0)])
        data_matrix_for_enrichment[List_initial,3] <- length(Overlapped_genes)
        data_matrix_for_enrichment[List_initial,4] <- length(Module2$gene[which(match(Module2$module, Unique_module2[k])>0)])
        
        data_matrix_for_enrichment[List_initial, 5] <- 20000
        data_matrix_for_enrichment[List_initial,6] <- Unique_module2[k]
        
        if(length(Overlapped_genes)){
          data_matrix_for_enrichment[List_initial,7] <- paste(Overlapped_genes, collapse = ",")
        } else{
          data_matrix_for_enrichment[List_initial,7] <- c("NULL")
        }
        data_matrix_for_enrichment[List_initial,8] <- database_name
        
        List_initial=List_initial + 1
      }
      # write.table(data_matrix_for_enrichment, file = paste0("Csvs/pathway_enrichment/temp1.",module,".dat"), quote = FALSE, sep = "\t",
      #             row.names = FALSE, col.names = FALSE)
      write.table(data_matrix_for_enrichment, file = paste0(output_Dir,"/",module,".dat"), quote = FALSE, sep = "\t",
                  row.names = FALSE, col.names = FALSE)
      # record_mat <- read.table(paste0("Csvs/pathway_enrichment/temp1.",module,".dat"))
      record_mat <- read.table(paste0(output_Dir,"/",module,".dat"))
      record_length <- dim(record_mat)
      enrichment_score <- data.frame()
      
      for(i in 1:record_length[1]){
        enrichment_score[i,1] <- record_mat[i,1]
        enrichment_score[i,2] <- phyper(record_mat[i,3], record_mat[i,4], record_mat[i,5]-record_mat[i,4], record_mat[i,2], lower.tail = FALSE)
        enrichment_score[i,3] <- record_mat[i,3]/record_mat[i,2]*record_mat[i,5]/record_mat[i,4]
        enrichment_score[i,4] <- 0
        enrichment_score[i,5]<-record_mat[i,8]
        enrichment_score[i,6]<-record_mat[i,6]
        enrichment_score[i,7]<-record_mat[i,3]
        enrichment_score[i,8]<-record_mat[i,7]
      }
      enrichment_score[,4]<-p.adjust(enrichment_score[,2], 'bonferroni')
      colnames(enrichment_score) <- c("Module","Pval","Enrichment","FDR","PathwaySource","Pathway","nOverlap","Overlap")
      enrichment_score$ModuleGeneCount = length(deg_list$gene)
      enrichment_score <- enrichment_score[order(enrichment_score$FDR),]
      x[[database_name]] <- data.frame(enrichment_score)
      
      if(z==1){
        all_pathways_df <- enrichment_score
      }else{
        all_pathways_df <- rbind(all_pathways_df, enrichment_score)
      }
    }
    
    all_pathways_df <- all_pathways_df[order(all_pathways_df$FDR),]
    all_pathways_df <- all_pathways_df[which(all_pathways_df$FDR < 0.05),]
    all_pathways_df <- all_pathways_df[which(all_pathways_df$nOverlap > 3),]
    x[["Combined"]] <- data.frame(all_pathways_df)
    all_pathways_df <- data.frame(all_pathways_df)
    result[[module]] <- all_pathways_df[order(all_pathways_df$FDR),]
    # x[["Combined"]][["Disease_Model"]] <- rep(celltype_cluster, nrow(x[["Combined"]]))
    
    total = data.frame(stringsAsFactors = FALSE)
    for(iter in 1:(length(x)-1)){
      total = rbind(total, x[[iter]])
    }
    total = total[order(total$nOverlap, decreasing = TRUE),]
    write.table(total, paste0(output_Dir,"/",module, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
    all_result[[module]] = total[order(total$FDR),]
    
    if(sum(total$FDR<0.05)>0){
      if(sum(total$FDR<0.05)<10){
        annotations = append(annotations, concatenate(total$Pathway[total$FDR<0.05], mysep = ", "))
      }
      else{
        annotations = append(annotations, concatenate(total$Pathway[total$FDR<0.05][1:10], mysep = ", "))
      }
    }
    else{
      annotations = append(annotations, paste0("(ns) ",concatenate(total$Pathway[1:3], mysep = ", ")))
    }
  }
  # DESCR_file = data.frame("MODULE"=unique(modules_list$MODULE), 
  #                         "ANNOTATION" = annotations, 
  #                         stringsAsFactors = FALSE)
  # write.table(DESCR_file, paste0(output_Dir,"/","DESCR.txt"),row.names = FALSE, quote = FALSE, sep = "\t")

}


summatedPathwaylogFCTableFigure <- function(module_df, DEG_df, onlySig = FALSE, 
                                            pathways = NULL, nPathways = 20,
                                            FDR_threshold = 0.05, logFC_threshold = .1){
  
  # use only significant genes to do the pathway enrichment
  DEG_df <- DEG_df[DEG_df$p_val_adj<FDR_threshold & abs(DEG_df$avg_logFC)>logFC_threshold,]
  DEG_mod = DEG_df[,c("Cell_type","GENE")]
  
  # when calculating the total logFC, take all genes into account (unless onlySig is TRUE)
  
  
  
  # plot the top different pathways
  #     overlap cannot share more than 4 genes?
  
  
    
}

# without integration
# if already have isolated cell type, then set onlyModifyPlots = TRUE and put 
#     "SeuratObj_celltype.rds" in the Subclustering folder for that celltype in your current working directory
subclusterAnalysis <- function(seuratObject, celltype, 
                               cellTypeMeta = "Cell_typev2",
                               resolution = "RNA_snn_res.0.5",
                               proportionMetaData = c("condition","Tissue"),
                               orderByMeta = NULL, orderFirst = NULL, 
                               conditions = list(c("WTW","TGW"),c("Hippocampus,Hypothalamus")),
                               renameConditions = list(c("WT","5XFAD"),c("HP","HYP")),
                               nMarkers = 6,
                               onlyModifyPlots = FALSE,
                               markerGeneDims = c(5,15),
                               markerPathwayDims = c(5,18),
                               reduction = "umap",
                               combinePropPlots=TRUE, pathway_resource = "./", 
                               mouseStudy = FALSE,
                               markerGeneMinMaxExpression = c(-2,2)){
  # subset
  Idents(seuratObject) <- cellTypeMeta
  if(length(unique(seuratObject@active.ident))!=1){ # not already subsetted
    seuratObject <- subset(seuratObject, idents = celltype)
  }
  
  ifelse(!dir.exists(file.path(paste0("./Subclustering/",celltype))), dir.create(file.path(paste0("./Subclustering/",celltype)),recursive = T), FALSE)
  
  if(!onlyModifyPlots){
    # Run RNA Clustering
    DefaultAssay(seuratObject) <- "RNA"
    seuratObject <- NormalizeData(seuratObject)
    seuratObject <- FindVariableFeatures(seuratObject)
    seuratObject <- ScaleData(seuratObject, features = rownames(seuratObject))
    seuratObject <- ScaleData(seuratObject)
    seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject))
    seuratObject <- FindNeighbors(seuratObject, dims = 1:30, reduction="pca",k.param = 25)
    resUse=seq(0.5,4,by=0.5)
    seuratObject <- FindClusters(object = seuratObject, resolution = resUse, verbose = T, reduction = "pca")
    seuratObject <- RunUMAP(seuratObject, dims = 1:30, reduction = "pca")
    seuratObject <- RunTSNE(object = seuratObject, reduction = "pca", dims = 1:30)
    
    pdf(file=paste0("./Subclustering/",celltype,"/DimPlot_",resolution,".pdf"))
    print(DimPlot(object = seuratObject, reduction = reduction, group.by = resolution, label = TRUE))
    dev.off()
    
    for(groupby in proportionMetaData){
      pdf(file=paste0("./Subclustering/",celltype,"/DimPlot_",groupby,".pdf"))
      print(DimPlot(object = seuratObject, reduction = reduction, group.by = groupby))
      dev.off()
    }
    
    Idents(seuratObject) <- resolution
    seuratObject@misc$markers <- FindAllMarkers(seuratObject, logfc.threshold = .1)
    saveRDS(seuratObject, paste0("./Subclustering/",celltype,"/SeuratObj_",celltype,".rds"))
    
  }
  else{
    seuratObject <- readRDS(paste0("./Subclustering/",celltype,"/SeuratObj_",celltype,".rds"))
  }

  # Plot of marker genes for each sub cluster -------------------------------------------------
  markers <- markersForPathway
  
  markers <- seuratObject@misc$markers
  markers <- markers[markers$avg_logFC>.1,]
  markers <- markers[markers$p_val_adj<0.05,]
  
  markersForPathway <- seuratObject@misc$markers
  markersForPathway <- markersForPathway[markersForPathway$p_val_adj<0.05,] # need to be changed!!
  # get top marker genes from each
  marker_genes = list()
  all_marker = c()
  for(ct in unique(markers$cluster)){
    marker_genes[[ct]] = markers$gene[markers$cluster==ct][1:nMarkers]
    all_marker = append(all_marker, markers$gene[markers$cluster==ct][1:nMarkers])
  }
  dup_genes = all_marker[duplicated(all_marker)]
  
  # so that all marker genes are unique
  # add the next top gene
  if(length(dup_genes)>0){
    cat("The duplicated marker genes are \n")
    print(dup_genes)
    for(dgene in dup_genes){
      for(ct in 1:length(marker_genes)){
        if(dgene %in% marker_genes[[ct]]){
          # can "keep" this marker gene
          # but need to take out of other clusters
          for(clust in (ct+1):length(marker_genes)){
            if(dgene %in% marker_genes[[clust]]){
              marker_genes[[clust]] = marker_genes[[clust]][marker_genes[[clust]]!=dgene]
              for(gene in markers$gene[markers$cluster==names(marker_genes)[clust]]){
                if(!(gene %in% all_marker)){
                  marker_genes[[clust]] = append(marker_genes[[clust]], gene) # replace with next marker gene
                  # update all_marker
                  all_marker = c()
                  for(ct in unique(seuratObject@active.ident)){
                    all_marker = append(all_marker, marker_genes[[ct]])
                  }
                  break
                }
                else next
              }
            }
            else next
          }
        }
        else next
      }
    }
  }
  
  for(gene in all_marker){
    if(!(gene %in% rownames(seuratObject@assays$RNA@scale.data))){
      cat(gene, " not in scaled data\n. You must scale all genes.\n")
    }
  }
  
  # order by condition proportion
  if(!is.null(orderByMeta)){
    counts <- cellTypeCount(seuratObject = seuratObject, 
                            meta_data_1 = orderByMeta, meta_data_2 = resolution, 
                            includePerc = FALSE)
    melted <- melt(counts, id.vars = colnames(counts)[1])
    colnames(melted) <- c(colnames(counts)[1], resolution,"Counts")
    melted$Counts <- as.numeric(melted$Counts) 
    percentages = c()
    for(samp in unique(melted[,resolution])){
      total = sum(melted$Counts[melted[,resolution]==samp])
      for(ct in unique(melted[,colnames(counts)[1]])){
        percentages = append(percentages, melted$Counts[melted[,resolution]==samp & melted[,colnames(counts)[1]]==ct]/total)
      }
    }
    melted$Percentage = percentages*100
    ordered_clusters <- unique(melted[,resolution])[order(melted$Percentage[melted$condition==orderFirst],
                                                          decreasing = TRUE)]
    ordered_clusters = as.character(unfactor(ordered_clusters))
    # rename clusters
    iter = 0
    for(clust in ordered_clusters){
      seuratObject@meta.data[,resolution] = gsub(paste0("\\<",clust,"\\>"), paste0(iter,"_new"), 
                                                 seuratObject@meta.data[,resolution])
      iter = iter + 1
    }
    seuratObject@meta.data[,resolution] = gsub("_new", "", seuratObject@meta.data[,resolution])
    seuratObject@meta.data[,resolution] = factor(seuratObject@meta.data[,resolution], 
                                                 levels = as.character(0:(length(unique(seuratObject@meta.data[,resolution]))-1)))
    
    # new DimPlot
    pdf(file=paste0("./Subclustering/",celltype,"/ordered_DimPlot_",resolution,".pdf"))
    print(DimPlot(object = seuratObject, reduction = reduction, group.by = resolution, label=TRUE))
    dev.off()
  }
  
  meta.data <- data.frame("Cluster" = seuratObject@meta.data[,c(resolution)])
  rownames(meta.data) = rownames(seuratObject@meta.data)
  
  data <- seuratObject@assays$RNA@scale.data[all_marker,]
  
  ordered_data = data.frame(stringsAsFactors = FALSE)
  first = FALSE
  for(ct in levels(seuratObject@meta.data[,resolution])){
    cells = rownames(meta.data)[meta.data[,"Cluster"]==ct]
    if(!first){
      ordered_data <- rbind(ordered_data, data[,cells])
      first = TRUE
    }
    else{
      ordered_data <- cbind(ordered_data, data[,cells])
    }
  }

  ordered_data <- as.matrix(ordered_data)
  ordered_data <- MinMax(ordered_data, min = markerGeneMinMaxExpression[1], max = markerGeneMinMaxExpression[2])
  rm(data)
  
  pdf(file=paste0("./Subclustering/",celltype,"/MarkerGenes_",resolution,".pdf"), height = 12)
  pheatmap(mat = ordered_data, 
           annotation_col = meta.data, 
           cluster_cols = FALSE, 
           # color = magma(n = 50), 
           show_colnames = FALSE, ) 
  dev.off()
  
  png(file=paste0("./Subclustering/",celltype,"/MinMax5_MarkerGenes_",resolution,".png"), 
      width = 800, height = 1200, units = "px")
  pheatmap(mat = ordered_data, 
           annotation_col = meta.data, 
           cluster_cols = FALSE, 
           # color = magma(n = 50), 
           show_colnames = FALSE, scale = "column") 
  dev.off()
  
  # Proportion plots ---------
  combined_df = list()
  for(meta in proportionMetaData){
    counts <- cellTypeCount(seuratObject = seuratObject, 
                            meta_data_1 = meta, meta_data_2 = resolution, 
                            includePerc = FALSE)
    propPlot <- normalizedDistributionProportionPlot(data = counts, meta_data_2 = resolution)
    #pdf(file=paste0("./Subclustering/",celltype,"/Proportion_Plot_",meta,".pdf")) 
    pdf(paste0("microglia_subclusters_prop",meta,"_SampleIntegrateRemoveDoub.pdf"))
    print(propPlot[[1]])
    dev.off()
    combined_df[[meta]] <- as.data.frame(propPlot[[2]])
  }
  
  if(combinePropPlots){
    toPlot <- data.frame()
    for(meta in names(combined_df)){
      prop_df = combined_df[[meta]]
      colnames(prop_df)[1] = "Variable"
      prop_df$Meta = meta
      toPlot = rbind(toPlot, prop_df)
    }
    combined <- ggplot(toPlot, aes(x = as.numeric(interaction(Meta,integrated_snn_res.0.5)), y = Percentage, fill = Variable)) + 
      geom_col(stat="identity",color="white", position = "stack") +
      scale_x_continuous(breaks=((as.numeric(unique(toPlot[,resolution]))*2-.5)),
                         labels=as.character(unique(toPlot[,resolution])), expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      facet_grid(.~integrated_snn_res.0.5, scales = "free", space = "free") +
      xlab("Cluster") +
      theme_bw() +
      theme(panel.border = element_blank(), 
            axis.text.y = element_text(size=10),
            axis.text.x = element_text(size=10),
            legend.title = element_blank(), 
            strip.text = element_blank(),
            strip.background = element_blank())
    
    pdf(file="combined_prop_plot_microglia_integrate_animal.pdf", 
        height = 3)
    # pdf(file=paste0("./Subclustering/",celltype,"/Combined_Proportion_Plot.pdf"), 
    #     height = 3)
    print(combined)
    dev.off()
  }
  
  # Cdh5
  
  # Plot of marker pathways for each sub cluster -------------------------------------------------
  Idents(seuratObject) <- resolution
  #markersForPathway <- FindAllMarkers(seuratObject, logfc.threshold = .1)

  markersForPathway$GENE <- markersForPathway$gene
  markersForPathway <- convertDfGeneColumnMouseHuman(df = markersForPathway, toSpecies = "human", forPathway = TRUE)
  markersForPathway = markersForPathway[markersForPathway$p_val_adj<0.05,]
  
  modfile_pathway = data.frame("MODULE"=markers$cluster, 
                               "GENE"=markers$gene, 
                               "avg_logFC" = markers$avg_logFC, stringsAsFactors = FALSE)
  modfile_pathway$MOUSE = modfile_pathway$GENE
  if(mouseStudy){
    modfile_pathway <- convertDfGeneColumnMouseHuman(df = modfile_pathway, 
                                                     toSpecies = "human", forPathway = TRUE)
    for(ct in modfile_pathway$MODULE){
      if(sum(duplicated(modfile_pathway$GENE[modfile_pathway$MODULE==ct]))>0){
        cat("Cluster ", ct, "\n")
        cat("There were duplicate human mappings for different mouse genes:\n")
        print(modfile_pathway[duplicated(modfile_pathway$GENE[modfile_pathway$MODULE==ct]),])
      }
    }
  }
  pathway_dir = paste0("./Subclustering/",celltype,"/Marker_Pathway_Enrichment")
  pathwayEnrichment(modules_list = modfile_pathway[,c("MODULE","GENE")], 
                    resources_path = pathway_resource, output_Dir = pathway_dir)
  pathway_files = list.files(pathway_dir)[grep(".txt", list.files(pathway_dir))]
  pathway_df = data.frame(stringsAsFactors = FALSE)
  for(file in pathway_files){
    temp = read.delim(paste0(pathway_dir,"/",file), stringsAsFactors = FALSE, header = TRUE)
    pathway_df = rbind(pathway_df, temp)
  }
  pathway_df$`negLog10FDR` = -log10(pathway_df$FDR)
  
  # pick top 25 pathways 
  
  # attempting to pick non-redundant pathways - putting this on hold 
  # ordered_pathway_df <- pathway_df[order(pathway_df$nOverlap, decreasing = TRUE),]
  # topPathways = c()
  # first = TRUE
  # for(path in ordered_pathway_df$Pathway){
  #   if(first){
  #     topPathways = append(topPathways, path)
  #     first = FALSE
  #     next
  #   }
  #   else{
  #     onListGenes = unlist(strsplit(ordered_pathway_df$Overlap[ordered_pathway_df$Pathway==path][1]))
  #   }
  # }
  
  pathway_df <- pathway_df_add_logFC(pathway_df = pathway_df, logFC_df = markersForPathway, 
                                     cluster_meta = "cluster", 
                                     calculateOnlyOverlapLogFC = FALSE, compiled_pathways = compiled_pathways)
  
  topPathways <- unique(pathway_df$Pathway[order(abs(pathway_df$avg_logFC), decreasing = TRUE)])[1:50]
  
  pathway_heat <- pathway_df[pathway_df$Pathway %in% topPathways,c("Module","Pathway","avg_logFC")]
  pathway_heat_mat = matrix(ncol = length(levels(seuratObject@meta.data[,resolution])),
                            nrow = length(topPathways))
  colnames(pathway_heat_mat) <- levels(seuratObject@meta.data[,resolution])
  rownames(pathway_heat_mat) <- topPathways
  for(col in colnames(pathway_heat_mat)){
    for(row in rownames(pathway_heat_mat)){
      pathway_heat_mat[row,col] = pathway_heat$avg_logFC[pathway_heat$Module==col & pathway_heat$Pathway==row]
    }
  }
  
  rownames(pathway_heat_mat) <- make.names(rownames(pathway_heat_mat))
  rownames(pathway_heat_mat) <- makePathwayPretty(rownames(pathway_heat_mat))
  
  meta.data <- data.frame("Cluster" = levels(seuratObject@meta.data[,resolution]))
  rownames(meta.data) <- levels(seuratObject@meta.data[,resolution])
  
  pdf(file=paste0("./Subclustering/",celltype,"/Pathway_heat_",resolution,".pdf"), width = 15, height = 12)
  print(pheatmap(mat = pathway_heat_mat, 
           cluster_cols = FALSE, annotation_col = meta.data,
           scale = "row", border_color = "grey60", cellheight = 12, cellwidth = 12,
           annotation_names_col = TRUE))
  dev.off()
  
}

# if calculateOnlyOverlapLogFC = FALSE, better to input an "inclusive" logFC_df 
pathway_df_add_logFC <- function(pathway_df = pathway_df, 
                                 logFC_df, 
                                 orderBy = "Overlap", cluster_meta = "cluster",
                                 calculateOnlyOverlapLogFC = FALSE, 
                                 compiled_pathways = "compiled_pathways.txt",
                                 onlySig=FALSE){
  if(calculateOnlyOverlapLogFC){
    avg_logFC = c()
    for(row in 1:nrow(pathway_df)){
      if(pathway_df$nOverlap[row]>0){
        ct_logFC_all = logFC_df[logFC_df[,cluster_meta]==pathway_df$Module[row],]
        genes = unlist(strsplit(pathway_df$Overlap[row], split = ","))
        for(gene in 1:length(genes)){
          #cat(genes[gene], "\n")
          if(sum(MT_gene_names==genes[gene])>0){
            #genes[gene] = names(MT_gene_names)[grep(genes[gene], MT_gene_names, fixed = TRUE)]
            tryCatch({genes[gene] = names(MT_gene_names)[grep(paste0("\\b",genes[gene],"\\b"), MT_gene_names)]}, 
                     warning=function(w) print(row))
            if(length(names(MT_gene_names)[grep(genes[gene], MT_gene_names, fixed = TRUE)]>1)){
              cat(row, "\n")
            }
          }
          else{
            next
          }
        }
        logFCs = ct_logFC_all$avg_logFC[match(x = genes, table = ct_logFC_all$GENE)]
        avg_logFC[row] = sum(logFCs[!is.na(logFCs)])
        if(length(logFCs)!=length(genes)){
          cat("Did not find equal number of logFC values in row ",row,"\n")
        }
      }
      else{
        avg_logFC[row] = 0
      }
    }
    pathway_df$avg_logFC = avg_logFC
  }
  else if(!calculateOnlyOverlapLogFC & onlySig){
    avg_logFC = c()
    compiled_pathways_df = read.delim(compiled_pathways, header = TRUE, stringsAsFactors = FALSE)
    for(row in 1:nrow(pathway_df)){
      pathway_genes = compiled_pathways_df$GENE[compiled_pathways_df$MODULE==pathway_df$Pathway[row]]
      ct_logFC_all = logFC_df[logFC_df[,cluster_meta]==pathway_df$Module[row] & logFC_df$p_val_adj<0.05,]
      logFCs = ct_logFC_all$avg_logFC[match(x = pathway_genes, table = ct_logFC_all$GENE_human)]
      avg_logFC[row] = sum(logFCs[!is.na(logFCs)])
    }
    pathway_df$avg_logFC = avg_logFC
    return(pathway_df)
  }
  else{
    avg_logFC = c()
    compiled_pathways_df = read.delim(compiled_pathways, header = TRUE, stringsAsFactors = FALSE)
    for(row in 1:nrow(pathway_df)){
      pathway_genes = compiled_pathways_df$GENE[compiled_pathways_df$MODULE==pathway_df$Pathway[row]]
      ct_logFC_all = logFC_df[logFC_df[,cluster_meta]==pathway_df$Module[row],]
      logFCs = ct_logFC_all$avg_logFC[match(x = pathway_genes, table = ct_logFC_all$GENE_human)]
      avg_logFC[row] = sum(logFCs[!is.na(logFCs)])
    }
    pathway_df$avg_logFC = avg_logFC
    return(pathway_df)
  }
}

# no title for conditions - distinguish with colors
stackedVlnPlotsCompareConditions <- function(seuratObject, features, celltypes, cellTypeMeta = "Cell_type",
                                             conditionMeta = "condition",
                                             conditions = c("WT","5XFAD"), 
                                             conditionColors = c("deepskyblue2","brown1")){
  first = TRUE
  Idents(seuratObject) <- cellTypeMeta
  Idents(seuratObject) <- "Cell_type_tissuev2"
  toPlot_genes <- list()
  for(feat in features){
    toPlot_celltypes <- list()
    for(ct in celltypes){
      if(first){
        plt <- ggarrange(VlnPlot(subset(seuratObject, idents = paste0(ct,"_Hippocampus")),features = feat,
                                 group.by = conditionMeta, sort = conditions,
                                 cols = conditionColors, pt.size = .05) +
                            NoLegend() + ggtitle("Hippocampus") +
                            theme(axis.line = element_blank(), axis.title = element_blank(), 
                            axis.text = element_blank(), axis.ticks = element_blank()),
                         VlnPlot(subset(seuratObject, idents = paste0(ct,"_Hypothalamus")),features = feat,
                                 group.by = conditionMeta, sort = conditions,
                                 cols = conditionColors, pt.size = .05, ) +
                           NoLegend() + ggtitle("Hypothalamus") +
                           scale_color_manual(values=c("brown")) +
                           theme(axis.line = element_blank(), axis.title = element_blank(), 
                                 axis.text = element_blank(), axis.ticks = element_blank(),
                                 plot.title = element_text(colour = "brown"))) 
        toPlot_celltypes[[ct]] <- annotate_figure(plt, top = text_grob(ct, size = 30, face = "bold"))
      }
      else{
        plt <- ggarrange(VlnPlot(subset(seuratObject, idents = paste0(ct,"_Hippocampus")),features = feat,
                                 group.by = conditionMeta, sort = conditions,
                                 cols = conditionColors, pt.size = .05) +
                           NoLegend() + 
                           theme(axis.line = element_blank(), axis.title = element_blank(), 
                                 axis.text = element_blank(), axis.ticks = element_blank(),
                                 plot.title = element_blank()),
                         VlnPlot(subset(seuratObject, idents = paste0(ct,"_Hypothalamus")),features = feat,
                                 group.by = conditionMeta, sort = conditions,
                                 cols = conditionColors, pt.size = .05) +
                           NoLegend() +
                           theme(axis.line = element_blank(), axis.title = element_blank(), 
                                 axis.text = element_blank(), axis.ticks = element_blank(),
                                 plot.title = element_blank())) 
        toPlot_celltypes[[ct]] <- plt 
      }
    }
    first = FALSE
    ct_plots <- ggarrange(plotlist = toPlot_celltypes, nrow = 1)
    toPlot_genes[[feat]] <- annotate_figure(ct_plots, 
                                            left = text_grob(feat, size = 30, face = "bold", just="right", rot = 90))
  }
  grid <- plot_grid(plotlist = toPlot_genes, nrow = length(features), ncol = 1, align = "v")
  png(file="test.png", 
      width = 5000, height = 1500, units = "px")
  grid
  dev.off()
  
}


# stackedVlnPlotsCompareConditionsv1 <- function(seuratObject, features, celltypes, cellTypeMeta = "Cell_type",
#                                              conditionMeta = "condition",
#                                              conditions = c("WT","5XFAD"), 
#                                              conditionColors = c("deepskyblue2","brown1")){
#   first = TRUE
#   Idents(seuratObject) <- cellTypeMeta
#   Idents(seuratObject) <- "Cell_type_tissuev2"
#   toPlot_genes <- list()
#   for(feat in features){
#     toPlot_celltypes <- list()
#     for(ct in celltypes){
#       if(first){
#         plt <- ggarrange(VlnPlot(subset(seuratObject, idents = paste0(ct,"_Hippocampus")),features = feat,
#                                  group.by = conditionMeta, sort = conditions,
#                                  cols = conditionColors, pt.size = .05) +
#                            NoLegend() + ggtitle("Hippocampus") +
#                            theme(axis.line = element_blank(), axis.title = element_blank(), 
#                                  axis.text = element_blank(), axis.ticks = element_blank()),
#                          VlnPlot(subset(seuratObject, idents = paste0(ct,"_Hypothalamus")),features = feat,
#                                  group.by = conditionMeta, sort = conditions,
#                                  cols = conditionColors, pt.size = .05) +
#                            NoLegend() + ggtitle("Hypothalamus") +
#                            theme(axis.line = element_blank(), axis.title = element_blank(), 
#                                  axis.text = element_blank(), axis.ticks = element_blank()))
#       }
#       if(first){
#         toPlot_celltypes[[ct]] <- plt + ggtitle(ct) + theme(plot.title = element_text(size=10))
#       }
#       else{
#         toPlot_celltypes[[ct]] <- plt + theme(plot.title = element_blank())
#       }
#     }
#     first = FALSE
#     ct_plots <- ggarrange(plotlist = toPlot_celltypes, nrow = 1)
#     toPlot_genes[[feat]] <- annotate_figure(ct_plots, 
#                                             left = text_grob(feat, size = 30, face = "bold", just="right", rot = 90))
#   }
#   grid <- plot_grid(plotlist = toPlot_genes, nrow = length(features), ncol = 1, align = "v")
#   png(file="test.png", 
#       width = 5000, height = 1500, units = "px")
#   grid
#   dev.off()
#   
# }
# 
# stackedVlnPlotsCompareConditions <- function(seuratObject, features, celltypes, cellTypeMeta = "Cell_type",
#                                              conditionMeta = "condition",
#                                              conditions = c("WT","5XFAD"), 
#                                              conditionColors = c("deepskyblue2","brown1")){
#   first = TRUE
#   Idents(seuratObject) <- cellTypeMeta
#   Idents(seuratObject) <- "Cell_type_tissuev2"
#   first = TRUE
#   toPlot_genes <- list()
#   for(feat in features){
#     toPlot_celltypes <- list()
#     for(ct in celltypes){
#       if(first){
#         plt <- ggarrange(VlnPlot(subset(HP, idents = ct),features = feat,
#                                  group.by = conditionMeta, sort = conditions,
#                                  cols = conditionColors, pt.size = .05) +
#                            NoLegend() + ggtitle("Hippocampus") +
#                            theme(axis.line = element_blank(), axis.title = element_blank(), 
#                                  axis.text = element_blank(), axis.ticks = element_blank(),
#                                  plot.title = element_text(size=15)),
#                          VlnPlot(subset(HYP, idents = ct),features = feat,
#                                  group.by = conditionMeta, sort = conditions,
#                                  cols = conditionColors, pt.size = .05, ) +
#                            NoLegend() + ggtitle("Hypothalamus") +
#                            scale_color_manual(values=c("brown")) +
#                            theme(axis.line = element_blank(), axis.title = element_blank(), 
#                                  axis.text = element_blank(), axis.ticks = element_blank(),
#                                  plot.title = element_text(colour = "brown", size=15))) 
#         toPlot_celltypes[[ct]] <- annotate_figure(plt, top = text_grob(ct, size = 30, face = "bold"))
#       }
#       else{
#         plt <- ggarrange(VlnPlot(subset(HP, idents = ct),features = feat,
#                                  group.by = conditionMeta, sort = conditions,
#                                  cols = conditionColors, pt.size = .05) +
#                            NoLegend() + 
#                            theme(axis.line = element_blank(), axis.title = element_blank(), 
#                                  axis.text = element_blank(), axis.ticks = element_blank(),
#                                  plot.title = element_blank()),
#                          VlnPlot(subset(HYP, idents = ct),features = feat,
#                                  group.by = conditionMeta, sort = conditions,
#                                  cols = conditionColors, pt.size = .05) +
#                            NoLegend() +
#                            theme(axis.line = element_blank(), axis.title = element_blank(), 
#                                  axis.text = element_blank(), axis.ticks = element_blank(),
#                                  plot.title = element_blank())) 
#         toPlot_celltypes[[ct]] <- plt 
#       }
#     }
#     first = FALSE
#     ct_plots <- ggarrange(plotlist = toPlot_celltypes, nrow = 1)
#     toPlot_genes[[feat]] <- annotate_figure(ct_plots, 
#                                             left = text_grob(feat, size = 30, face = "bold", just = "center", rot = 90))
#   }
#   grid <- plot_grid(plotlist = toPlot_genes, nrow = length(features), ncol = 1, align = "v")
#   png(file="testagain23.png", 
#       width = 5000, height = 2400, units = "px")
#   grid
#   dev.off()
#   
# }
# 
# 
# stackedVlnPlotsCompareConditions <- function(seuratObject, features, celltypes, cellTypeMeta = "Cell_type",
#                                              conditionMeta = "condition",
#                                              conditions = c("WT","5XFAD"), 
#                                              conditionColors = c("deepskyblue2","brown1")){
#   first = TRUE
#   Idents(seuratObject) <- cellTypeMeta
#   Idents(seuratObject) <- "Cell_type_tissuev2"
#   first = TRUE
#   toPlot_genes <- list()
#   
#     toPlot_celltypes <- list()
#     for(ct in celltypes){
#       ct_HP <- subset(HP, idents = ct)
#       ct_HYP <- subset(HYP, idents = ct)
#       for(feat in features){
#       if(first){
#         
#         plt <- ggarrange(VlnPlot(ct_HP,features = feat,
#                                  group.by = conditionMeta, sort = conditions,
#                                  cols = conditionColors, pt.size = .05) +
#                            NoLegend() + ggtitle("Hippocampus") +
#                            theme(axis.line = element_blank(), axis.title = element_blank(), 
#                                  axis.text = element_blank(), axis.ticks = element_blank(),
#                                  plot.title = element_text(size=20)),
#                          VlnPlot(subset(ct_HYP, idents = ct),features = feat,
#                                  group.by = conditionMeta, sort = conditions,
#                                  cols = conditionColors, pt.size = .05, ) +
#                            NoLegend() + ggtitle("Hypothalamus") +
#                            scale_color_manual(values=c("brown")) +
#                            theme(axis.line = element_blank(), axis.title = element_blank(), 
#                                  axis.text = element_blank(), axis.ticks = element_blank(),
#                                  plot.title = element_text(colour = "brown", size=20)))
#         toPlot_celltypes[[feat]] <- annotate_figure(plt, left = text_grob(feat, size = 30, face = "bold", just = "center", rot = 90))
#       }
#       else{
#         plt <- ggarrange(VlnPlot(ct_HP,features = feat,
#                                  group.by = conditionMeta, sort = conditions,
#                                  cols = conditionColors, pt.size = .05) +
#                            NoLegend() + 
#                            theme(axis.line = element_blank(), axis.title = element_blank(), 
#                                  axis.text = element_blank(), axis.ticks = element_blank(),
#                                  plot.title = element_blank()),
#                          VlnPlot(ct_HYP,features = feat,
#                                  group.by = conditionMeta, sort = conditions,
#                                  cols = conditionColors, pt.size = .05) +
#                            NoLegend() +
#                            theme(axis.line = element_blank(), axis.title = element_blank(), 
#                                  axis.text = element_blank(), axis.ticks = element_blank(),
#                                  plot.title = element_blank())) 
#         toPlot_celltypes[[feat]] <- plt 
#       }
#         first = FALSE
#     }
#     
#     toPlot_genes[[ct]] <- annotate_figure(plt, top = text_grob(ct, size = 30, face = "bold"))
#   }
#   grid <- plot_grid(plotlist = toPlot_genes, nrow = 1, ncol = length(features), align = "v")
#   png(file="testagain2.png", 
#       width = 5000, height = 2400, units = "px")
#   grid
#   dev.off()
# }
# 
# 
# # make version that is faster and utilizes annotate_figure
stackedVlnPlotsCompareConditions <- function(seuratObject, features, celltypes, cellTypeMeta = "Cell_type",
                                             conditionMeta = "condition",
                                             conditions = c("WT","5XFAD"),
                                             conditionColors = c("brown1","deepskyblue2")){
  # Idents(seuratObject) <- cellTypeMeta
  # Idents(seuratObject) <- "Cell_type_tissuev2"
  first_gene = TRUE
  toPlot_genes <- list()
  toPlot_celltypes <- list()
  for(ct in celltypes){
    first_ct = TRUE
    ct_HP <- subset(HP, idents = ct)
    ct_HYP <- subset(HYP, idents = ct)
    for(feat in features){
      if(first_ct){ # label the gene once only
        plt <- VlnPlot(ct_HP,features = feat,
                                 group.by = conditionMeta, sort = conditions,
                                 cols = conditionColors, pt.size = .05) +
                           NoLegend() +
                           # theme(axis.line = element_blank(), axis.title = element_blank(),
                           #       axis.text = element_blank(), axis.ticks = element_blank(),
                           #       plot.title = element_blank()) +
          theme(axis.line.x = element_blank(), axis.title = element_blank(),
                axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                plot.title = element_blank())
        
        plt_HYP <- VlnPlot(subset(ct_HYP, idents = ct),features = feat,
                                 group.by = conditionMeta, sort = conditions,
                                 cols = conditionColors, pt.size = .05, ) +
                           NoLegend() +
                           # theme(axis.line = element_blank(), axis.title = element_blank(),
                           #       axis.text = element_blank(), axis.ticks = element_blank(),
                           #       plot.title = element_blank())
          theme(axis.line.x = element_blank(), axis.title = element_blank(),
                axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                plot.title = element_blank())
        
        # plt$layers[[2]]$inherit.aes <- TRUE
        # plt_HYP$layers[[2]]$inherit.aes <- TRUE
        # plt_HYP$layers[[2]]$geom$draw_key <- TRUE
        # plt_HYP$layers[[2]]$geom$default_aes$colour <- "coral4"
        # plt$layers[[2]]$geom$default_aes$colour <- "black"
        tissues <- plot_grid(annotate_figure(plt, top = text_grob("Hippocampus", size=20, face = "bold")),
                             annotate_figure(plt_HYP, top = text_grob("Hypothalamus", size=20, face = "bold", color = "coral4")))
        if(first_gene){
          tissues <- annotate_figure(tissues, left = text_grob(feat, size = 30, face = "bold", just = "center", rot = 90))
          toPlot_celltypes[[feat]] <- tissues
        }
        else{
          toPlot_celltypes[[feat]] <- tissues
        }
        first_ct = FALSE
      }
      else{
        plt <- VlnPlot(ct_HP,features = feat,
                       group.by = conditionMeta, sort = conditions,
                       cols = conditionColors, pt.size = .05) +
          NoLegend() +
          # theme(axis.line = element_blank(), axis.title = element_blank(),
          #       axis.text = element_blank(), axis.ticks = element_blank(),
          #       plot.title = element_blank())
          theme(axis.line.x = element_blank(), axis.title = element_blank(),
                axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                plot.title = element_blank())
        plt_HYP <- VlnPlot(subset(ct_HYP, idents = ct),features = feat,
                           group.by = conditionMeta, sort = conditions,
                           cols = conditionColors, pt.size = .05) +
          NoLegend() +
            theme(axis.line.x = element_blank(), axis.title = element_blank(),
                  axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                  plot.title = element_blank())
        #plt_HYP$layers[[2]]$geom$default_aes$colour <- "coral4"
          # theme(axis.line = element_blank(), axis.title = element_blank(),
          #       axis.text = element_blank(), axis.ticks = element_blank(),
          #       plot.title = element_blank())
          
        tissues <- ggarrange(plt, plt_HYP)
        if(first_gene){
          tissues <- annotate_figure(tissues, left = text_grob(feat, size = 30, face = "bold", just = "center", rot = 90))
          toPlot_celltypes[[feat]] <- tissues
        }
        else{
          toPlot_celltypes[[feat]] <- tissues
        }
      }
    }
    first_gene = FALSE
    toPlot_genes[[ct]] <- annotate_figure(plot_grid(plotlist = toPlot_celltypes, ncol = 1), top = text_grob(ct, size = 25, face = "bold"))
  }
  grid <- plot_grid(plotlist = toPlot_genes, nrow = 1, ncol = length(toPlot_genes), align = "v")
  grid <- plot_grid(plotlist = toPlot_genes, nrow = 1, ncol = 2, align = "v")
  grid <- ggarrange(plotlist = toPlot_genes, nrow = 1, ncol = length(features))
  png(file="Shared_genes.png",
      width = 5000, height = 5000, units = "px")
  grid
  dev.off()
}

# takes in melted data frame and converts to a matrix
# First variable (column) are the rows
# Second variable are the columns
# Third variable is the values
convertToMatrix <- function(df){
  DEG_mat <- matrix(nrow = length(unique(df[,1])),
                    ncol = length(unique(df[,2])))
  rownames(DEG_mat) <- unique(df[,1])
  colnames(DEG_mat) <- unique(df[,2])
  for(col in colnames(DEG_mat)){
    for(row in rownames(DEG_mat)){
      if(length(df[,3][df[,1]==row & df[,2]==col])==0){
        DEG_mat[row,col] = 0 
      }
      else DEG_mat[row,col] = df[,3][df[,1]==row & df[,2]==col]
    }
  }
  return(DEG_mat)
}

library(RColorBrewer)
myCol <- c(brewer.pal(12, "Set3"), brewer.pal(3, "Set1"))
stackedVlnPlot <- function(seuratObject, genes, cellTypes, cellTypeColors){
  for(i in 1:length(celltypes)){
    single_cells <- names(Idents(seuratObject))[which(Idents(seuratObject) == celltypes[i])]
    for(j in 1:length(genes)){
      value <- GetAssayData(object = seuratObject, slot = "data")[genes[j],single_cells] # try with normalized data?
      gene <- rep(genes[j],length(value))
      cell_type <- rep(celltypes[i],length(value))
      if(j==1){
        cell_type_matrix <- data.frame(value,gene,cell_type)
      }
      else{
        cell_type_matrix <- rbind(cell_type_matrix,data.frame(value,gene,cell_type))
      }
    }
    if(i==1){
      DGE_matrix <- cell_type_matrix
    } else{
      DGE_matrix <- rbind(DGE_matrix,cell_type_matrix)
    }
  }
  DGE_matrix$cell_type <- factor(DGE_matrix$cell_type, levels = celltypes)
  gg <- ggplot(DGE_matrix, aes(x=cell_type, y=value, fill=cell_type))
  gg <- gg + geom_violin(position=position_dodge(1.0),scale="width") 
  gg <- gg + scale_fill_manual(values=myCol)
  gg <- gg + facet_wrap(~interaction(gene),ncol = 1,strip.position="left",scales = "free_y")
  gg <- gg + labs(x="")
  gg <- gg + theme_bw()
  gg <- gg + theme(strip.background = element_blank(),
                   legend.position="right",
                   axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.text.x = element_custom(colour="black", fill=myCol),
                   panel.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.spacing.y = unit(0.6,"lines"),
                   panel.border = element_blank(), 
                   #axis.line.x = element_line(),
                   strip.text.y = element_text(angle=180,hjust=1,vjust = 0.5, size = 10),
                   plot.margin = margin(t = 140))
  gg <- gg + scale_x_discrete(position = "top")
  gg
  return(gg)
}


###### MODIFED from https://stackoverflow.com/questions/45956883/set-text-background-to-ggplot-axis-text ######
element_custom <- function(...) {
  structure(list(...), class = c("element_custom", "element_blank"))
}

element_grob.element_custom <- function(element, label, x, y, ...)  {
  tg <- textGrob(label, x=x, gp=gpar(col=element$colour, fontsize=15), 
                 rot = 90, just = "left", y=unit(.25,"line"))
  padding <- unit(1,"line")
  rg <- rectGrob(x=x,width=unit(1,"line"), height=grobHeight(tg) + padding, 
                 gp=gpar(fill = element$fill, col=NA, lineheight=3), hjust = .5,vjust=0)
  gTree(children=gList(rg, tg), height=grobHeight(tg) + padding, cl="custom_axis")
}

####### END OF MODIFIED #####



# gg <- ggplot(DGE_matrix, aes(x=gene, y=value, fill=cell_type)) 
# gg <- gg + geom_violin(position=position_dodge(.25),scale="width") 
# gg <- gg + scale_fill_manual(values=myCol)
# gg <- gg + facet_wrap(~cell_type,ncol = 1,strip.position = "left", scales = "free_y") 
# gg <- gg + labs(x="")
# gg <- gg + theme_bw()
# gg <- gg + theme(
#                  strip.text.y = element_text(angle = 0, hjust=1, size=10, face="bold"),
#                  legend.position="none",
#                  axis.title.y = element_blank(),
#                  axis.text.y = element_blank(),
#                  axis.ticks.y = element_blank(),
#                  axis.ticks.x = element_blank(),
#                  axis.text.x = element_text(angle=360,hjust=1,vjust = 0.5, size = 12),
#                  panel.background = element_blank(),
#                  panel.grid.major = element_blank(),
#                  panel.grid.minor = element_blank(),
#                  panel.spacing.y = unit(0.1,"lines"),
#                  panel.border = element_blank(), 
#                  axis.line.x = element_line(),
#                  plot.margin = margin(t = 20))
# 
# 
# g <- ggplot_gtable(ggplot_build(gg))
# strip_both <- which(grepl('strip-', g$layout$name))
# fills <- myCol
# cols <- c(rep("black",12), rep("white",3))
# k <- 1
# for (i in strip_both) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   l <- which(grepl('title', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   m <- which(grepl('text',g$grobs[[i]]$grobs[[1]]$children[[l]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[l]]$children[[m]]$gp$col <- cols[k]
#   k <- k+1
# }
# pdf("./Neuron_subcluster/Neuron_MarkerGene_ColoredLabels_flip2.pdf", width=7, height = 7)
# grid.draw(g)
# dev.off()

# https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html

simpleStackedVlnPlot <- function(seuratObject, genes, celltypes=NULL, cellTypeMeta=NULL){
  if(!is.null(cellTypeMeta)) Idents(seuratObject) <- cellTypeMeta
  if(is.null(celltypes)) celltypes = unfactor(unique(seuratObject@active.ident))
  
  for(i in 1:length(celltypes)){
    single_cells <- names(Idents(seuratObject))[which(Idents(seuratObject) == celltypes[i])]
    for(j in 1:length(genes)){
      value <- GetAssayData(object = seuratObject, slot = "data")[genes[j],single_cells] # try with normalized data?
      gene <- rep(genes[j],length(value))
      cell_type <- rep(celltypes[i],length(value))
      if(j==1){
        cell_type_matrix <- data.frame(value,gene,cell_type)
      }
      else{
        cell_type_matrix <- rbind(cell_type_matrix,data.frame(value,gene,cell_type))
      }
    }
    if(i==1){
      DGE_matrix <- cell_type_matrix
    } else{
      DGE_matrix <- rbind(DGE_matrix,cell_type_matrix)
    }
  }
  DGE_matrix$cell_type <- factor(DGE_matrix$cell_type, levels = celltypes)
  gg <- ggplot(DGE_matrix, aes(x=cell_type, y=value, fill=cell_type))
  gg <- gg + geom_violin(position=position_dodge(1.0),scale="width") 
  gg <- gg + scale_fill_manual(values=myCol)
  gg <- gg + facet_wrap(~interaction(gene),ncol = 1,strip.position="left",scales = "free_y")
  gg <- gg + labs(x="")
  gg <- gg + theme_bw()
  gg <- gg + theme(strip.background = element_blank(),
                   legend.position="right",
                   axis.title.y = element_blank(),
                   #axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.ticks.x = element_blank(),
                   #axis.text.x = element_custom(colour="black", fill=myCol),
                   panel.grid.major.x = element_blank(), 
                   panel.background = element_blank(),
                   #panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.spacing.y = unit(0.6,"lines"),
                   panel.border = element_blank(), strip.placement = "inside",
                   #axis.line.x = element_line(),
                   strip.text.y = element_text(angle=180,hjust=1,vjust = 0.5, size = 10, ),
                   plot.margin = margin(t = 140)) + scale_y_continuous(position = "right")
  gg <- gg + scale_x_discrete(position = "top")
  gg
  return(gg)
}


# only positive markers
findUniqueMarkers <- function(markers_df, lfc=.25, FDR=0.05){
  markers_df %>% 
    filter(avg_logFC>lfc, p_val_adj<FDR) -> markers_df
  unique_df = data.frame(stringsAsFactors = FALSE)
  for(gene in unique(markers_df$gene)){
    if(sum(markers_df$gene %in% gene)==1){
      temp = markers_df[markers_df$gene==gene,]
      unique_df = rbind(unique_df, temp)
    }
    else next
  }
  result <- list("Unique_Markers"=unique_df, 
                 "Melted_Unique"=unmelt(df = unique_df, id_column = "cluster", value_column = "gene"))
}

processSeurat <- function(seuratObject, combineLevel, name, skipIntegration = FALSE){
  if(!skipIntegration){
    DefaultAssay(seuratObject) <- "RNA"
    seuratObject <- runSampleCCA(seuratObject=seuratObject,combineLevel=combineLevel,numFeatures=2000,numDims=30,
                                 tissue=name,kParam=25,resUse=seq(0.5,4,by=0.5))
    seuratObject@reductions$IntegrateSampleUMAP <- seuratObject@reductions$umap
    seuratObject@reductions$IntegrateSampleTSNE <- seuratObject@reductions$tsne
  }
  DefaultAssay(seuratObject) <- "RNA"
  seuratObject <- NormalizeData(seuratObject)
  seuratObject <- FindVariableFeatures(seuratObject)
  seuratObject <- ScaleData(seuratObject)
  seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject))
  seuratObject <- FindNeighbors(seuratObject, dims = 1:30, reduction="pca",k.param = 25)
  resUse=seq(0.5,4,by=0.5)
  seuratObject <- FindClusters(object = seuratObject, resolution = resUse, verbose = T, reduction = "pca")
  seuratObject <- RunUMAP(seuratObject, dims = 1:30, reduction = "pca")
  seuratObject <- RunTSNE(object = seuratObject, reduction = "pca", dims = 1:30)
  return(seuratObject)
}


concatenate=function(myvect, mysep="")
{
  if(length(myvect)==0) return(myvect)
  if(length(myvect)==1) return(myvect)
  string = ""
  for(item in myvect){
    string = paste(string, item, sep = mysep)
  }
  string = substring(string, first=(nchar(mysep)+1))
  return(string)
}

describe_manipulation <- function(DEG_df, comparison = "Comparison",
                                  logFC_threshold = 0.25, p_threshold = 0.05){
  DEG_df = DEG_df[abs(DEG_df$avg_logFC)>logFC_threshold & DEG_df$p_val_adj<p_threshold,]
  cell_type_perturbed = list()
  # narrow down to "multiple perturbed subset"
  result = data.frame(stringsAsFactors = FALSE)
  for(ct in unique(DEG_df$Cell_type)){
    genes_perturbed = c()
    for(gene in unique(DEG_df$GENE)){
      if(sum(DEG_df$GENE[DEG_df$Cell_type==ct]==gene)>1){
        genes_perturbed = append(genes_perturbed, gene)
      }
      else{
        next
      }
    }
    if(length(genes_perturbed)==0) next
    cell_type_perturbed[[ct]] = genes_perturbed
    celltype_result = data.frame("Cell_type" = ct,
                                 "GENE"=genes_perturbed, stringsAsFactors = FALSE)
    Comparisons = unique(DEG_df[,comparison])
    for(comp in Comparisons){
      vect_FC = c()
      for(g in genes_perturbed){
        if(sum(DEG_df$GENE[DEG_df[,comparison]==comp & DEG_df$Cell_type==ct]==g)>0){
          vect_FC = append(vect_FC, DEG_df$avg_logFC[DEG_df[,comparison]==comp & DEG_df$GENE==g & DEG_df$Cell_type==ct])
        }
        else{
          vect_FC = append(vect_FC, "None") # may have fc but not significant
        }
      }
      celltype_result[,comp] = vect_FC
    }
    result = rbind(result, celltype_result)
  }
  return(result)
}

describe_effect <- function(df, comparisons, corresponding_sign){
  Effect_vect = c()
  for(n in 1:nrow(df)){
    effect = c()
    if(df[,comparisons[1]][n]=="None"){
      effect = append(effect, "NONE")
    }
    else if(df[,comparisons[1]][n]>0){
      effect = append(effect, "UP")
    }
    else{
      effect = append(effect, "DOWN")
    }
    effect = append(effect, corresponding_sign[1])
    if(df[,comparisons[2]][n]=="None"){
      effect = append(effect, "NONE")
    }
    else if(df[,comparisons[2]][n]>0){
      effect = append(effect, "UP")
    }
    else{
      effect = append(effect, "DOWN")
    }
    effect = append(effect, corresponding_sign[2])
    
    effect = concatenate(effect, mysep = "_")
    Effect_vect = append(Effect_vect, effect)
  }
  return(Effect_vect)
}


#https://github.com/SydneyBioX/scdney/blob/master/R/scDC_noClustering.R
scDC_noClustering <- function (cellTypes = NULL,
                               subject = NULL,
                               calCI = TRUE,
                               calCI_method = c("BCa", "multinom", "percentile"),
                               nboot = 10000,
                               conf_level = 0.95,
                               ncores = 1,
                               verbose = TRUE)
{
  
  x <- 1:length(cellTypes)
  n <- length(x)
  
  if(length(cellTypes)!=length(subject)){
    stop("the vector length of cell type info and subject info don't match!")
  }
  
  if(calCI){
    calCI_method <- match.arg(calCI_method, choices = c("BCa", "multinom",
                                                        "percentile"), several.ok = TRUE)
    if(is.null(nboot)){
      warnings("number of bootstrap is set as 10000")
      nboot = 10000
    }
    
  }
  
  if(!calCI&is.null(nboot)){
    warnings("number of bootstrap is set as 100")
    nboot = 100
  }
  
  cellTypes <- as.character(cellTypes)
  subject <- as.character(subject)
  
  tab <- table(cellTypes, subject)
  info <-  reshape2::melt(tab)
  
  if(verbose){
    print("Calculating sample proportion...")
  }
  df <- data.frame(cellTypes = cellTypes, subject = subject)
  # thetahat <- .calculateProp(x, df)
  calProp_hat <- .calculateProp(x, df)
  thetahat <- calProp_hat$prop
  nhat <- calProp_hat$count
  
  if(verbose){
    print("Calculating bootstrap proportion...")
  }
  bootsam <- do.call(rbind, .stratifiedBootstrap(1:length(cellTypes),
                                                 strata = subject, times = nboot))
  
  # thetastar <-  do.call(cbind, parallel::mclapply(1:nboot, function(i){
  #   .calculateProp(bootsam[i,], df)
  # }, mc.cores = ncores))
  calProp_star <- parallel::mclapply(1:nboot, function(i){
    .calculateProp(bootsam[i,], df)
  }, mc.cores = ncores)
  
  thetastar <-  do.call(cbind, lapply(calProp_star, "[[", "prop"))
  nstar <- do.call(cbind, lapply(calProp_star, "[[", "count"))
  
  res_multinom <- NULL
  res_BCa <- NULL
  res_percentile <- NULL
  
  if(calCI){
    alpha <-  c((1-conf_level)/2, 1-(1-conf_level)/2)
    if (!all(alpha < 1) || !all(alpha > 0))
      stop("All elements of alpha must be in (0,1)")
    alpha_sorted <- sort(alpha)
    if (nboot <= 1/min(alpha_sorted[1], 1 - alpha_sorted[length(alpha_sorted)]))
      warning("nboot is not large enough to estimate your chosen alpha.")
    
    if(verbose){
      print(paste("Calculating", calCI_method, "..."))
    }
    if("multinom" %in% calCI_method){
      
      
      ### Calculating the CI based on the multinomial
      confpoints = lapply(1:ncol(tab), function(i){
        multi_res <- DescTools::MultinomCI(tab[,i], conf.level = conf_level)
        multi_res <- multi_res[,-1]
        colnames(multi_res) <- c("conf_low", "conf_high")
        multi_res
      })
      res_multinom <- info[,1:2]
      res_multinom <- cbind(res_multinom, do.call(rbind, confpoints))
      rownames(res_multinom) <- 1:nrow(res_multinom)
      res_multinom$method <- "multinom"
      # return(list(results = res,
      #             confpoints = confpoints,
      #             thetastar = thetastar,
      #             thetahat = thetahat,
      #             nstar = nstar,
      #             nhat = nhat,
      #             info = info))
      
    }
    if("BCa" %in% calCI_method){
      
      if(verbose){
        print("Calculating z0 ...")
      }
      
      z0 <- sapply(1:nrow(thetastar), function(i) qnorm(sum(thetastar[i,] < thetahat[i])/nboot))
      if(verbose){
        print("Calculating acc ...")
      }
      
      u <- list()
      u <- parallel::mclapply(1:n, function(i){
        .calculateProp(x[-i], df)$prop
      }, mc.cores = ncores)
      
      u <- do.call(cbind, u)
      uu <- matrix(rep(rowMeans(u), ncol(u)), ncol = ncol(u)) - u
      acc <- apply(uu, 1, function(x) sum(x * x * x)/(6 * (sum(x * x))^1.5))
      
      zalpha <- qnorm(alpha)
      tt <- list()
      for(i in 1:length(z0)){
        tt[[i]] <- pnorm(z0[i] + (z0[i] + zalpha)/(1 - acc[i] * (z0[i] + zalpha)))
      }
      
      confpoints <- lapply(1:length(tt), function(i) quantile(x = thetastar[i,], probs = tt[[i]], type = 1))
      confpoints_original <- confpoints
      confpoints <- lapply(confpoints, function(x) {
        names(x) <- NULL
        x})
      confpoints <- lapply(confpoints, function(x) cbind(alpha, x))
      confpoints <- lapply(confpoints, function(x){
        dimnames(x)[[2]] <- c("alpha", "bca point")
        x
      })
      
      res_BCa <- info[,1:2]
      res_BCa$conf_low <- do.call(rbind, lapply(confpoints, function(x)x[1,2]))
      res_BCa$conf_high <- do.call(rbind, lapply(confpoints, function(x)x[2,2]))
      colnames(res_BCa) <- c("cellTypes", "subject", "conf_low", "conf_high")
      res_BCa$method <- "BCa"
      # return(list(result = res,
      #             confpoints = confpoints,
      #             z0 = z0,
      #             acc = acc,
      #             u = u,
      #             thetastar = thetastar,
      #             thetahat = thetahat,
      #             nstar = nstar,
      #             nhat = nhat,
      #             info = info,
      #             confpoints_original = confpoints_original,
      #             tt = tt))
    }
    
    
    if("percentile" %in% calCI_method){
      
      confpoints <- apply(thetastar, 1, function(x)quantile(x, probs = alpha))
      
      res_percentile <- info[,1:2]
      res_percentile$conf_low <- confpoints[1,]
      res_percentile$conf_high <- confpoints[2,]
      res_percentile$method <- "percentile"
      
      # return(list(results = res,
      #             confpoints = confpoints,
      #             thetastar = thetastar,
      #             thetahat = thetahat,
      #             nstar = nstar,
      #             nhat = nhat,
      #             info = info))
    }
  }
  
  res <- rbind(res_BCa, res_percentile, res_multinom)
  
  
  return(list(results = res,
              thetastar = thetastar,
              thetahat = thetahat,
              nstar = nstar,
              nhat = nhat,
              info = info))
  
  
}

.calculateProp <- function(x, xdata){
  tab <- table(xdata[x,]$cellTypes, xdata[x,]$subject)
  tab_prop <- c(tab/matrix(rep(colSums(tab), nrow(tab)), nrow = nrow(tab), byrow = T))
  tab_count <- c(tab)
  return(list(count = tab_count, prop = tab_prop))
}

.stratifiedBootstrap <- function(idx, strata, times = 1000){
  strata_list <- unique(strata)
  index_list <- lapply(1:times, function(x){
    unlist(lapply(strata_list, function(s){
      idx_strata = idx[strata == s]
      sample(idx_strata, length(idx_strata), replace = T)
    }))
  })
  return(index_list)
}

barplotCI <- function(res, condition){
  df_toPlot <- res$results
  df_toPlot$median <- apply(res$thetastar, 1, median)
  df_toPlot$cond <- condition
  df_toPlot$method <- factor(df_toPlot$method, levels = c("BCa", "percentile", "multinom"))
  n_method <- length(unique(df_toPlot$method))
  
  g_bar <- ggplot2::ggplot(df_toPlot, aes(x = subject, y = median, fill = cond)) +
    ggplot2::geom_bar(stat="identity", position = "dodge", alpha = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::ylab("Proportion") +
    ggplot2::geom_errorbar(aes(ymin=conf_low, ymax=conf_high, color = method), width=.3,lwd = 1,
                           position=position_dodge(width = 0.5)) +
    ggplot2::theme(axis.text.x = element_text(angle = 90), text = element_text(size = 12)) +
    ggplot2::scale_color_manual(values = .CIbarColor(n_method)) +
    ggplot2::facet_wrap(~cellTypes, ncol = 4) +
    ggplot2::coord_flip()+
    ggplot2::ylim(c(0,1))+
    NULL
  
  g_bar
  return(g_bar)
  
  
}

densityCI <- function(res, condition){
  df_toPlot <- reshape2::melt(res$thetastar)
  df_toPlot$cellTypes <- res$info[,"cellTypes"][df_toPlot$Var1]
  df_toPlot$subject <- res$info[,"subject"][df_toPlot$Var1]
  
  df_toPlot$cond <- condition
  
  conf_line <- res$results
  conf_line$cond <- condition
  
  conf_line$method <- factor(conf_line$method, levels = c("BCa", "percentile", "multinom"))
  n_method <- length(unique(conf_line$method))
  
  g_density <- ggplot2::ggplot(df_toPlot, aes(x = value, y = subject, fill = cond)) +
    stat_density_ridges(alpha = 0.5) +
    ggplot2::geom_segment(data = conf_line, aes(x = conf_low, xend = conf_low, y = as.numeric(subject),
                                                yend = as.numeric(subject) + .9,
                                                color = method, linetype = method),
                          lwd = 1, alpha = 0.8) +
    ggplot2::geom_segment(data = conf_line, aes(x = conf_high, xend = conf_high, y = as.numeric(subject),
                                                yend = as.numeric(subject) + .9,
                                                color = method, linetype = method),
                          lwd = 1, alpha = 0.8) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = element_text(size = 12)) +
    ggplot2::scale_color_manual(values = .CIbarColor(n_method)) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::xlab("Proportion") +
    ggplot2::facet_wrap(~cellTypes, ncol = 4, scales = "free_x")
  
  g_density
  return(g_density)
  
}

fitGLM <- function(res, condition, subject_effect = TRUE, pairwise = TRUE, fixed_only = FALSE, verbose = TRUE){
  
  fit_random <- list()
  fit_fixed <- list()
  for(i in 1:ncol(res$nstar)){
    # idx <- indexes_list[[i]]
    if(verbose){
      if(i%%10==0){
        print(paste("fitting GLM...", i))
      }
    }
    
    
    glm_df <-  cbind(res$info[,1:2], res$nstar[,i])
    
    # glm_df <- melt(glm_df)
    colnames(glm_df) <- c("cellTypes", "subject", "cell_count")
    glm_df$cond <- condition
    
    if(subject_effect){
      if(pairwise){
        
        fit_fixed[[i]] <- stats::glm(cell_count ~ cellTypes + cond +  cellTypes:cond + subject,
                                     data = glm_df,
                                     family = poisson(link=log))
        if(!fixed_only){
          fit_random[[i]] <- lme4::glmer(cell_count ~ cellTypes + cond +  cellTypes:cond + (1 | subject ),
                                         data = glm_df, family = poisson(link=log),
                                         control = glmerControl(nAGQ = 0L))
        }
      }else{
        fit_fixed[[i]] <- stats::glm(cell_count ~ cellTypes + cond +  cellTypes:cond + subject,
                                     data = glm_df, family = poisson(link=log))
        if(!fixed_only){
          
          fit_random[[i]] <- lme4::glmer(cell_count ~ cellTypes + cond +
                                           cellTypes:cond + (1 | subject ), data = glm_df,
                                         family = poisson(link=log),
                                         control = glmerControl(nAGQ = 0L))
        }
      }
    }else{
      fit_fixed[[i]] <- stats::glm(cell_count ~ cellTypes + cond +  cellTypes:cond, data = glm_df,
                                   family = poisson(link=log))
    }
    
  }
  
  if(!subject_effect){
    fixed_only = TRUE
  }
  
  if(!fixed_only){
    pool_res_random = mice::pool(fit_random)
    pool_res_fixed = mice::pool(fit_fixed)
    return(list(pool_res_random = pool_res_random,
                pool_res_fixed = pool_res_fixed,
                fit_random = fit_random,
                fit_fixed = fit_fixed))
  }else{
    pool_res_fixed = mice::pool(fit_fixed)
    return(list(pool_res_fixed = pool_res_fixed,
                fit_fixed = fit_fixed))
  }
  
}

.CIbarColor <- function(n){
  if(n==1){
    colour <- "black"
  }else{
    colour <- c("red", "blue", "purple")
  }
  return(colour)
}



compareStudiesTableDEGsDirection <- function(combined_DEGs, comparison = "Tissue", 
                                    FDR_threshold = 0.05, lfc_threshold = 0.1){
  combined_DEGs = combined_DEGs[combined_DEGs$p_val_adj<FDR_threshold & abs(combined_DEGs$avg_logFC)>lfc_threshold,]
  combined_DEGs = combined_DEGs[order(abs(combined_DEGs$avg_logFC), decreasing = TRUE),] # so that large effect ones show up first
  combined_DEGs$GENE_Direction = ifelse(combined_DEGs$avg_logFC>0, 
                                        paste0(combined_DEGs$GENE,"_UP"),
                                        paste0(combined_DEGs$GENE,"_DOWN"))
  result = data.frame()
  melted_result = data.frame()
  for(ct in unique(combined_DEGs$Cell_type)){
    group_genes = list()
    group_genes_direction = list()
    temp = data.frame("Cell_type" = ct, stringsAsFactors = FALSE)
    for(el in unique(combined_DEGs[,comparison])){
      group_genes[[el]] = combined_DEGs$GENE[combined_DEGs$Cell_type==ct & combined_DEGs[,comparison]==el]
      group_genes_direction[[el]] = combined_DEGs$GENE_Direction[combined_DEGs$Cell_type==ct & combined_DEGs[,comparison]==el]
      temp[,el] = paste0("(",length(group_genes[[el]]),") ", concatenate(group_genes[[el]], mysep = ", "))
    }
    specific1 = setdiff(group_genes[[1]], group_genes[[2]])
    temp[,paste0("Specific_", unique(combined_DEGs[,comparison])[1])] = paste0("(",length(specific1),") ",
                                                                               concatenate(specific1, mysep = ", "))
    if(length(specific1)>1){
      melt_temp = data.frame("Cell_type" = ct, 
                             "Type"=paste0("Specific_", unique(combined_DEGs[,comparison])[1]), 
                             "GENE"= specific1, stringsAsFactors = FALSE)
      melted_result = rbind(melted_result, melt_temp)
    }
    specific2 = setdiff(group_genes[[2]], group_genes[[1]])
    temp[,paste0("Specific_", unique(combined_DEGs[,comparison])[2])] = paste0("(",length(specific2),") ",
                                                                               concatenate(specific2, mysep = ", "))
    if(length(specific2)>1){
      melt_temp = data.frame("Cell_type" = ct, 
                             "Type"=paste0("Specific_", unique(combined_DEGs[,comparison])[2]), 
                             "GENE"= specific2, stringsAsFactors = FALSE)
      melted_result = rbind(melted_result, melt_temp)
    }
    shared = intersect(group_genes[[1]], group_genes[[2]])
    shared_direction = intersect(group_genes_direction[[1]], group_genes_direction[[2]])
    shared_direction_stripped <- gsub("_DOWN","",shared_direction)
    shared_direction_stripped <- gsub("_UP","",shared_direction_stripped)
    shared_discordant = setdiff(shared, shared_direction_stripped)
    new_shared = setdiff(shared, shared_discordant)
    new_shared = combined_DEGs[combined_DEGs$Tissue=="pfCtx" &
                                 combined_DEGs$Cell_type==ct,]$GENE_Direction[match(new_shared,
                                                                                    combined_DEGs[combined_DEGs$Tissue=="pfCtx"&
                                                                                                    combined_DEGs$Cell_type==ct,]$GENE)]
    temp[,paste0("Shared")] = paste0("(",length(new_shared),") ",
                                     concatenate(new_shared, mysep = ", "))
    shared_discordant = combined_DEGs[combined_DEGs$Tissue=="pfCtx"&
                                        combined_DEGs$Cell_type==ct,]$GENE_Direction[match(shared_discordant,
                                                                                           combined_DEGs[combined_DEGs$Tissue=="pfCtx"&
                                                                                                           combined_DEGs$Cell_type==ct,]$GENE)]
    temp[,"Shared_Direction_Discordant"] = paste0("(",length(shared_discordant),") ",
                                                  concatenate(shared_discordant, mysep = ", "))
    if(length(shared)>1){
      melt_temp = data.frame("Cell_type" = ct, 
                             "Type"= "Shared", 
                             "GENE"= new_shared, 
                             stringsAsFactors = FALSE)
      melted_result = rbind(melted_result, melt_temp)
    }
    if(length(shared_discordant)>1){
      melt_temp = data.frame("Cell_type" = ct, 
                             "Type"= "Shared_discordant", 
                             "GENE"= shared_discordant, 
                             stringsAsFactors = FALSE)
      melted_result = rbind(melted_result, melt_temp)
    }
    result = rbind(result, temp)
  }
  return(list(result, melted_result))
}

gene_dot_plot_custom <- function(
  DEG_df, 
  vars = c("GENE","Cell_type"), 
  var_list = NULL, 
  pdf_name = "DEG_dotPlot.pdf", 
  colors = c("blue","gray","red"),
  explore.show = "logFC"
){
  DEG_df$neglog10fdr = -log10(DEG_df$p_val_adj)
  DEG_df$neglog10fdr[DEG_df$neglog10fdr==Inf] <- 301
  
  if(!is.null(var_list)){
    selected = DEG_df[DEG_df[,vars[1]] %in% var_list[[1]],]
    if(length(var_list)>1){
      selected = selected[selected[,vars[2]] %in% var_list[[2]],]
    }
    for(v in 1:length(var_list)){
      if(length(unique(selected[,vars[v]]))!=length(var_list[[v]])){
        selected[,vars[v]] = factor(selected[,vars[v]], levels = unique(selected[,vars[v]]))
      }
      else{
        selected[,vars[v]] = factor(selected[,vars[v]], levels = var_list[[v]])
      }
    }
  }
  else if(length(unique(DEG_df$GENE))==1 & is.null(var_list)){
    selected = DEG_df
  }
  else{
    # select the top 30 most "regulated" genes, show all cell types
    if(explore.show == "logFC"){
      DEG_df = DEG_df[order(abs(DEG_df$avg_logFC), decreasing = TRUE),]
    }
    else{
      DEG_df = DEG_df[order(abs(DEG_df$neglog10fdr), decreasing = TRUE),]
    }
    genes = c()
    for(g in DEG_df$GENE){
      genes = append(genes, g)
      genes = unique(genes)
      if(length(genes)==30) break
    }
    for(g in genes){
      selected = rbind(selected, DEG_df[DEG_df$GENE==g,])
    }
  }
  
  selected$`p(Adj) < 0.05` = selected$p_val_adj < 0.05
  
  if(sum(selected$p_val_adj>0.05)==0){ 
    shapes = c(16)
  }
  else{
    shapes = c(18,16)
  }
  
  selected$GENE = factor(selected$GENE, levels = var_list[["GENE"]])
  
  D <- ggplot(selected, aes_string(x=vars[1], y=vars[2], size="neglog10fdr", color="avg_logFC")) + 
    geom_point(aes(shape=`p(Adj) < 0.05`)) + 
    #geom_point(data=selected[selected$Stringent=="Yes",],pch=21,colour="black", stroke=1) +
    theme_bw() +
    facet_grid(.~Tissue, scales = "free", space = "free") +
    theme(panel.border = element_blank(), #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, margin = margin(t=5), size=14, hjust = 1,vjust = 1),
          axis.text.y = element_text(size = 14, margin = margin(l=45)),
          #axis.text = element_text(face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_text(size=15),
          legend.text = element_text(size=10),
          plot.margin = margin(t=0),
          legend.position = "bottom") +
    theme(strip.background = element_blank(), #remove background for facet labels
          panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
          #panel.spacing = unit(0, "lines"),
          #strip.placement = "outside",
          #strip.text = element_blank())+
          strip.text = element_text(margin = margin(t=0,b=3,r=0,l=0),
                                    hjust = .5, size=14)) +
    labs(size=expression(-log[10]~p(Adj)), color=expression(logFC)) +
    scale_color_gradient2(low = colors[1],  mid=colors[2], high = colors[3],
                          breaks = c(-.5, 0, .5), labels = c(-.5, 0, .5)) +
    scale_size(range = c(3, 12), breaks = c(1.3,25,50,75,100)) +
    scale_shape_manual(values=shapes) +
    guides(shape = guide_legend(override.aes = list(size = 5)))
  
  return(D)
}

drawPathwayGenePlot <- function(pathway_df, 
                                pathway_gene_df, 
                                pathway, genes=NULL, 
                                min_max=c("min"=-8,"max"=8),
                                DEG_df){
  pathway_df <- pathway_df[pathway_df$Pathway==pathway,]
  if(is.null(genes)){
    pathway_gene_df <- pathway_gene_df[pathway_gene_df$Pathway==pathway,]
    pathway_gene_df$numAppearances <- sapply(pathway_gene_df$GENE, function(x){return(nrow(pathway_gene_df[pathway_gene_df$GENE==x,]))})
    pathway_gene_df <- pathway_gene_df[order(pathway_gene_df$numAppearances, decreasing = TRUE),]
    #pathway_gene_df <- pathway_gene_df[order(abs(pathway_gene_df$avg_logFC), decreasing = TRUE),]
    genes <- unique(pathway_gene_df$GENE)[1:20]
  }
  
  # cluster cell types based on DEGs
  DEG_df = DEG_df[DEG_df$GENE %in% genes,]
  DEG_mat <- convertToMatrix(df = DEG_df[DEG_df$Tissue=="Hippocampus",c("Cell_type","GENE","avg_logFC")])
  dendro <- as.dendrogram(hclust(d = dist(x = DEG_mat)))
  ct_order <- order.dendrogram(dendro)
  ct_levels <- rownames(DEG_mat)[ct_order]
  ct_levels <- ct_levels[ct_levels!="Choroid Plexus"]
  ct_levels <- c(ct_levels, "Choroid Plexus","Tanycyte")
  
  pathway_df$`Cell type` <- factor(pathway_df$`Cell type`, levels = ct_levels)
  DEG_df$Cell_type <- factor(DEG_df$Cell_type, levels = ct_levels)
  
  
  pathway_df$avg_logFC <- MinMax(pathway_df$avg_logFC, min = min_max["min"], max = min_max["max"])
  heat <- ggplot(pathway_df, aes(x=`Cell type`, y=Pathway, fill=avg_logFC, label=nOverlap)) +
    geom_tile(colour="grey", size=.7) +
    geom_text(size=4, aes(colour=abs(avg_logFC)>7), fontface="bold") +
    theme_bw() +
    facet_grid(.~Tissue, scales = "free", space = "free") +
    #scale_fill_gradient2(low="blue",high="red", mid = "white", breaks = c(0), labels=c("")) + 
    scale_fill_gradient2(low="blue",high="red", mid = "white") + 
    theme(strip.background = element_blank(), #remove background for facet labels
          panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
          #panel.spacing = unit(0, "lines"),
          #strip.placement = "outside",
          plot.margin = margin(t=0,b=0),
          strip.text = element_text(margin = margin(t=0,b=3,r=0,l=0),
                                    hjust = .5, size=14),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 11.5),
          axis.text.y = element_text(size = 15),
          panel.grid = element_blank(), axis.line = element_blank(), 
          legend.title = element_text(margin=margin(t=0,b=3,r=0,l=0)),
          legend.position = "none",
          axis.title = element_blank(),
          axis.text.x.bottom = element_blank(),
          axis.ticks = element_blank())+
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_color_manual(values = c("black", "white"))+
    labs(fill="Summed logFC")
  
  DEG_df$neglog10fdr = -log10(DEG_df$p_val_adj)
  DEG_df$neglog10fdr[DEG_df$neglog10fdr==Inf] <- 301
  
  #DEG_df$neglog10fdr <- MinMax(DEG_df$neglog10fdr, min = 0, max = 75)
  
  # cluster genes
  DEG_mat <- convertToMatrix(df = DEG_df[,c("GENE","Ct_tissue","avg_logFC")])
  dendro <- as.dendrogram(hclust(d = dist(x = DEG_mat)))
  gene_order <- order.dendrogram(dendro)
  gene_levels <- rownames(DEG_mat)[gene_order]
  DEG_df$GENE <- factor(DEG_df$GENE, levels = rev(gene_levels))
  
  
  DEG_df$`p(Adj) < 0.05` = DEG_df$p_val_adj < 0.05
  
  shapes = c(18,16)
  
  DEG_df$avg_logFC <- MinMax(DEG_df$avg_logFC, min = -1, max = 1)

  dot <- ggplot(DEG_df, aes_string(x="Cell_type", y="GENE", size="neglog10fdr", color="avg_logFC")) + 
    geom_point(aes(shape=`p(Adj) < 0.05`)) + 
    #geom_point(data=selected[selected$Stringent=="Yes",],pch=21,colour="black", stroke=1) +
    theme_bw() +
    facet_grid(.~Tissue, scales = "free", space = "free") +
    theme(panel.border = element_blank(), #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, margin = margin(t=0), size=14, hjust = 1,vjust = 1),
          axis.text.y = element_text(size = 14, margin = margin(l=80)),
          #axis.text = element_text(face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_text(size=15),
          legend.text = element_text(size=10),
          plot.margin = margin(t=0),
          legend.position = "bottom") +
    theme(strip.background = element_blank(), #remove background for facet labels
          panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
          #panel.spacing = unit(0, "lines"),
          #strip.placement = "outside",
          strip.text = element_blank())+
    # strip.text = element_text(margin = margin(t=0,b=3,r=0,l=0),
    #                           hjust = .5, size=14)) +
    labs(size=expression(-log[10]~p["Adj"]), color=expression(logFC), shape=expression(p["Adj"] < 0.05)) +
    scale_color_gradient2(low = colors[1],  mid=colors[2], high = colors[3],
                          breaks = c(-.5, 0, .5), labels = c(-.5, 0, .5)) +
    scale_size(range = c(3, 9), breaks = c(1.3,25,50,75,100)) +
    scale_shape_manual(values=shapes) +
    guides(shape = guide_legend(override.aes = list(size = 5)))
  
  gb1 <- ggplot_build(heat)
  gb2 <- ggplot_build(dot)
  n1 <- length(gb1$panel$ranges[[1]]$y.labels)
  n2 <- length(gb2$panel$ranges[[1]]$y.labels)
  gA <- ggplot_gtable(gb1)
  gB <- ggplot_gtable(gb2)
  #g <- gtable:::rbind_gtable(gA, gB, "last")
  g <- gtable:::rbind_gtable(gA, gB)
  
  return(g)
  
}

drawPathwayGenePlot_first <- function(pathway_df, 
                                pathway_gene_df, 
                                pathway, genes=NULL, 
                                min_max=c("min"=-8,"max"=8),
                                DEG_df){
  pathway_df <- pathway_df[pathway_df$Pathway==pathway,]
  if(is.null(genes)){
    pathway_gene_df <- pathway_gene_df[pathway_gene_df$Pathway==pathway,]
    pathway_gene_df$numAppearances <- sapply(pathway_gene_df$GENE, function(x){return(nrow(pathway_gene_df[pathway_gene_df$GENE==x,]))})
    pathway_gene_df <- pathway_gene_df[order(pathway_gene_df$numAppearances, decreasing = TRUE),]
    #pathway_gene_df <- pathway_gene_df[order(abs(pathway_gene_df$avg_logFC), decreasing = TRUE),]
    genes <- unique(pathway_gene_df$GENE)[1:20]
  }
  
  # cluster cell types based on DEGs
   DEG_df = DEG_df[DEG_df$GENE %in% genes,]
  # DEG_mat <- convertToMatrix(df = DEG_df[DEG_df$Tissue=="Hippocampus",c("Cell_type","GENE","avg_logFC")])
  # dendro <- as.dendrogram(hclust(d = dist(x = DEG_mat)))
  # ct_order <- order.dendrogram(dendro)
  # ct_levels <- rownames(DEG_mat)[ct_order]
  # ct_levels <- ct_levels[ct_levels!="Choroid Plexus"]
  # ct_levels <- c(ct_levels, "Choroid Plexus","Tanycyte")
  
  pathway_df$`Cell type` <- factor(pathway_df$`Cell type`, levels = ct_levels)
  DEG_df$Cell_type <- factor(DEG_df$Cell_type, levels = ct_levels)
  
  
  pathway_df$avg_logFC <- MinMax(pathway_df$avg_logFC, min = min_max["min"], max = min_max["max"])
  heat <- ggplot(pathway_df, aes(x=`Cell type`, y=Pathway, fill=avg_logFC, label=nOverlap)) +
    geom_tile(colour="grey", size=.7) +
    geom_text(size=4, aes(colour=abs(avg_logFC)>7), fontface="bold") +
    theme_bw() +
    facet_grid(.~Tissue, scales = "free", space = "free") +
    #scale_fill_gradient2(low="blue",high="red", mid = "white", breaks = c(0), labels=c("")) + 
    scale_fill_gradient2(low="blue",high="red", mid = "white") + 
    theme(strip.background = element_blank(), #remove background for facet labels
          panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
          #panel.spacing = unit(0, "lines"),
          #strip.placement = "outside",
          plot.margin = margin(t=0,b=0),
          #strip.text = element_blank(),
          strip.text = element_text(margin = margin(t=0,b=3,r=0,l=0),
                                    hjust = .5, size=14, face = "bold"),
          axis.text.x = element_blank(),
          #axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 11.5),
          axis.text.y = element_text(size = 15, face = "bold"),
          panel.grid = element_blank(), axis.line = element_blank(), 
          #legend.title = element_text(margin=margin(t=0,b=3,r=0,l=0)),
          legend.position = "none",
          axis.title = element_blank(),
          axis.text.x.bottom = element_blank(),
          axis.ticks = element_blank())+
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_color_manual(values = c("black", "white"))+
    labs(fill="Summed logFC")
  
  DEG_df$neglog10fdr = -log10(DEG_df$p_val_adj)
  DEG_df$neglog10fdr[DEG_df$neglog10fdr==Inf] <- 301
  
  #DEG_df$neglog10fdr <- MinMax(DEG_df$neglog10fdr, min = 0, max = 75)
  
  # cluster genes
  DEG_mat <- convertToMatrix(df = DEG_df[,c("GENE","Ct_tissue","avg_logFC")])
  dendro <- as.dendrogram(hclust(d = dist(x = DEG_mat)))
  gene_order <- order.dendrogram(dendro)
  gene_levels <- rownames(DEG_mat)[gene_order]
  DEG_df$GENE <- factor(DEG_df$GENE, levels = rev(gene_levels))
  
  
  DEG_df$`p(Adj) < 0.05` = DEG_df$p_val_adj < 0.05
  
  shapes = c(18,16)
  
  DEG_df$avg_logFC <- MinMax(DEG_df$avg_logFC, min = -.7, max = .7)
  
  dot <- ggplot(DEG_df, aes_string(x="Cell_type", y="GENE", size="neglog10fdr", color="avg_logFC")) + 
    geom_point(aes(shape=`p(Adj) < 0.05`)) + 
    #geom_point(data=selected[selected$Stringent=="Yes",],pch=21,colour="black", stroke=1) +
    theme_bw() +
    facet_grid(.~Tissue, scales = "free", space = "free") +
    theme(panel.border = element_blank(), #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_blank(),
          #axis.text.x = element_text(angle = 45, margin = margin(t=0), size=14, hjust = 1,vjust = 1),
          axis.text.y = element_text(size = 14, margin = margin(l=150)),
          #axis.text = element_text(face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_text(size=15),
          legend.text = element_text(size=10),
          plot.margin = margin(t=0),
          legend.position = "none",
          axis.ticks = element_blank()) +
    theme(strip.background = element_blank(), #remove background for facet labels
          panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
          #panel.spacing = unit(0, "lines"),
          #strip.placement = "outside",
          strip.text = element_blank())+
    # strip.text = element_text(margin = margin(t=0,b=3,r=0,l=0),
    #                           hjust = .5, size=14)) +
    labs(size=expression(-log[10]~p["Adj"]), color=expression(logFC), shape=expression(p["Adj"] < 0.05)) +
    scale_color_gradient2(low="blue",high="red", mid = "gray",
                          breaks = c(-.5, 0, .5), labels = c(-.5, 0, .5)) +
    scale_size(range = c(3, 9), breaks = c(1.3,25,50,75,100)) +
    scale_shape_manual(values=shapes) +
    guides(shape = guide_legend(override.aes = list(size = 5)))
  
  gb1 <- ggplot_build(heat)
  gb2 <- ggplot_build(dot)
  n1 <- length(gb1$panel$ranges[[1]]$y.labels)
  n2 <- length(gb2$panel$ranges[[1]]$y.labels)
  gA <- ggplot_gtable(gb1)
  gB <- ggplot_gtable(gb2)
  #g <- gtable:::rbind_gtable(gA, gB, "last")
  g <- gtable:::rbind_gtable(gA, gB)
  
  return(g)
  
}

drawPathwayGenePlot_middle <- function(pathway_df, 
                                pathway_gene_df, 
                                pathway, genes=NULL, 
                                min_max=c("min"=-8,"max"=8),
                                DEG_df){
  pathway_df <- pathway_df[pathway_df$Pathway==pathway,]
  if(is.null(genes)){
    pathway_gene_df <- pathway_gene_df[pathway_gene_df$Pathway==pathway,]
    pathway_gene_df$numAppearances <- sapply(pathway_gene_df$GENE, function(x){return(nrow(pathway_gene_df[pathway_gene_df$GENE==x,]))})
    pathway_gene_df <- pathway_gene_df[order(pathway_gene_df$numAppearances, decreasing = TRUE),]
    #pathway_gene_df <- pathway_gene_df[order(abs(pathway_gene_df$avg_logFC), decreasing = TRUE),]
    genes <- unique(pathway_gene_df$GENE)[1:20]
  }
  
  # cluster cell types based on DEGs
  DEG_df = DEG_df[DEG_df$GENE %in% genes,]
  # DEG_mat <- convertToMatrix(df = DEG_df[DEG_df$Tissue=="Hippocampus",c("Cell_type","GENE","avg_logFC")])
  # dendro <- as.dendrogram(hclust(d = dist(x = DEG_mat)))
  # ct_order <- order.dendrogram(dendro)
  # ct_levels <- rownames(DEG_mat)[ct_order]
  # ct_levels <- ct_levels[ct_levels!="Choroid Plexus"]
  # ct_levels <- c(ct_levels, "Choroid Plexus","Tanycyte")
  
  pathway_df$`Cell type` <- factor(pathway_df$`Cell type`, levels = ct_levels)
  DEG_df$Cell_type <- factor(DEG_df$Cell_type, levels = ct_levels)
  
  
  pathway_df$avg_logFC <- MinMax(pathway_df$avg_logFC, min = min_max["min"], max = min_max["max"])
  heat <- ggplot(pathway_df, aes(x=`Cell type`, y=Pathway, fill=avg_logFC, label=nOverlap)) +
    geom_tile(colour="grey", size=.7) +
    geom_text(size=4, aes(colour=abs(avg_logFC)>7), fontface="bold") +
    theme_bw() +
    facet_grid(.~Tissue, scales = "free", space = "free") +
    #scale_fill_gradient2(low="blue",high="red", mid = "white", breaks = c(0), labels=c("")) + 
    scale_fill_gradient2(low="blue",high="red", mid = "white") + 
    theme(strip.background = element_blank(), #remove background for facet labels
          panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
          #panel.spacing = unit(0, "lines"),
          #strip.placement = "outside",
          plot.margin = margin(t=0,b=0),
          strip.text = element_blank(),
          # strip.text = element_text(margin = margin(t=0,b=3,r=0,l=0),
          #                           hjust = .5, size=14),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 11.5),
          axis.text.y = element_text(size = 15, face = "bold"),
          panel.grid = element_blank(), axis.line = element_blank(), 
          legend.title = element_text(margin=margin(t=0,b=3,r=0,l=0)),
          legend.position = "none",
          axis.title = element_blank(),
          axis.text.x.bottom = element_blank(),
          axis.ticks = element_blank())+
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_color_manual(values = c("black", "white"))+
    labs(fill="Summed logFC")
  
  DEG_df$neglog10fdr = -log10(DEG_df$p_val_adj)
  DEG_df$neglog10fdr[DEG_df$neglog10fdr==Inf] <- 301
  
  #DEG_df$neglog10fdr <- MinMax(DEG_df$neglog10fdr, min = 0, max = 75)
  
  # cluster genes
  DEG_mat <- convertToMatrix(df = DEG_df[,c("GENE","Ct_tissue","avg_logFC")])
  dendro <- as.dendrogram(hclust(d = dist(x = DEG_mat)))
  gene_order <- order.dendrogram(dendro)
  gene_levels <- rownames(DEG_mat)[gene_order]
  DEG_df$GENE <- factor(DEG_df$GENE, levels = rev(gene_levels))
  
  
  DEG_df$`p(Adj) < 0.05` = DEG_df$p_val_adj < 0.05
  
  shapes = c(18,16)
  
  DEG_df$avg_logFC <- MinMax(DEG_df$avg_logFC, min = -.7, max = .7)
  
  dot <- ggplot(DEG_df, aes_string(x="Cell_type", y="GENE", size="neglog10fdr", color="avg_logFC")) + 
    geom_point(aes(shape=`p(Adj) < 0.05`)) + 
    #geom_point(data=selected[selected$Stringent=="Yes",],pch=21,colour="black", stroke=1) +
    theme_bw() +
    facet_grid(.~Tissue, scales = "free", space = "free") +
    theme(panel.border = element_blank(), #panel.grid.major = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          #axis.text.x = element_text(angle = 45, margin = margin(t=0), size=14, hjust = 1,vjust = 1),
          axis.text.y = element_text(size = 14, margin = margin(l=150)),
          #axis.text = element_text(face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_text(size=15),
          legend.text = element_text(size=10),
          plot.margin = margin(t=0),
          legend.position = "none",
          axis.ticks = element_blank()) +
    theme(strip.background = element_blank(), #remove background for facet labels
          panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
          #panel.spacing = unit(0, "lines"),
          #strip.placement = "outside",
          strip.text = element_blank())+
    # strip.text = element_text(margin = margin(t=0,b=3,r=0,l=0),
    #                           hjust = .5, size=14)) +
    labs(size=expression(-log[10]~p["Adj"]), color=expression(logFC), shape=expression(p["Adj"] < 0.05)) +
    scale_color_gradient2(low = colors[1],  mid=colors[2], high = colors[3],
                          breaks = c(-.5, 0, .5), labels = c(-.5, 0, .5)) +
    scale_size(range = c(3, 9), breaks = c(1.3,25,50,75,100)) +
    scale_shape_manual(values=shapes) +
    guides(shape = guide_legend(override.aes = list(size = 5)))
  
  gb1 <- ggplot_build(heat)
  gb2 <- ggplot_build(dot)
  n1 <- length(gb1$panel$ranges[[1]]$y.labels)
  n2 <- length(gb2$panel$ranges[[1]]$y.labels)
  gA <- ggplot_gtable(gb1)
  gB <- ggplot_gtable(gb2)
  #g <- gtable:::rbind_gtable(gA, gB, "last")
  g <- gtable:::rbind_gtable(gA, gB)
  
  return(g)
  
}

drawPathwayGenePlot_last <- function(pathway_df, 
                                       pathway_gene_df, 
                                       pathway, genes=NULL, 
                                       min_max=c("min"=-8,"max"=8),
                                       DEG_df){
  pathway_df <- pathway_df[pathway_df$Pathway==pathway,]
  if(is.null(genes)){
    pathway_gene_df <- pathway_gene_df[pathway_gene_df$Pathway==pathway,]
    pathway_gene_df$numAppearances <- sapply(pathway_gene_df$GENE, function(x){return(nrow(pathway_gene_df[pathway_gene_df$GENE==x,]))})
    pathway_gene_df <- pathway_gene_df[order(pathway_gene_df$numAppearances, decreasing = TRUE),]
    #pathway_gene_df <- pathway_gene_df[order(abs(pathway_gene_df$avg_logFC), decreasing = TRUE),]
    genes <- unique(pathway_gene_df$GENE)[1:20]
  }
  
  # cluster cell types based on DEGs
  DEG_df = DEG_df[DEG_df$GENE %in% genes,]
  # DEG_mat <- convertToMatrix(df = DEG_df[DEG_df$Tissue=="Hippocampus",c("Cell_type","GENE","avg_logFC")])
  # dendro <- as.dendrogram(hclust(d = dist(x = DEG_mat)))
  # ct_order <- order.dendrogram(dendro)
  # ct_levels <- rownames(DEG_mat)[ct_order]
  # ct_levels <- ct_levels[ct_levels!="Choroid Plexus"]
  # ct_levels <- c(ct_levels, "Choroid Plexus","Tanycyte")
  
  pathway_df$`Cell type` <- factor(pathway_df$`Cell type`, levels = ct_levels)
  DEG_df$Cell_type <- factor(DEG_df$Cell_type, levels = ct_levels)
  
  
  pathway_df$avg_logFC <- MinMax(pathway_df$avg_logFC, min = min_max["min"], max = min_max["max"])
  heat <- ggplot(pathway_df, aes(x=`Cell type`, y=Pathway, fill=avg_logFC, label=nOverlap)) +
    geom_tile(colour="grey", size=.7) +
    geom_text(size=4, aes(colour=abs(avg_logFC)>7), fontface="bold") +
    theme_bw() +
    facet_grid(.~Tissue, scales = "free", space = "free") +
    #scale_fill_gradient2(low="blue",high="red", mid = "white", breaks = c(0), labels=c("")) + 
    scale_fill_gradient2(low="blue",high="red", mid = "white") + 
    theme(strip.background = element_blank(), #remove background for facet labels
          panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
          #panel.spacing = unit(0, "lines"),
          #strip.placement = "outside",
          plot.margin = margin(t=0,b=0),
          strip.text = element_blank(),
          # strip.text = element_text(margin = margin(t=0,b=3,r=0,l=0),
          #                           hjust = .5, size=14),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 11.5),
          axis.text.y = element_text(size = 15, face = "bold"),
          panel.grid = element_blank(), axis.line = element_blank(), 
          legend.title = element_text(margin=margin(t=0,b=3,r=0,l=0)),
          legend.position = "none",
          axis.title = element_blank(),
          axis.text.x.bottom = element_blank(),
          axis.ticks = element_blank())+
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_color_manual(values = c("black", "white"))+
    labs(fill="Summed logFC")
  
  DEG_df$neglog10fdr = -log10(DEG_df$p_val_adj)
  DEG_df$neglog10fdr[DEG_df$neglog10fdr==Inf] <- 301
  
  #DEG_df$neglog10fdr <- MinMax(DEG_df$neglog10fdr, min = 0, max = 75)
  
  # cluster genes
  DEG_mat <- convertToMatrix(df = DEG_df[,c("GENE","Ct_tissue","avg_logFC")])
  dendro <- as.dendrogram(hclust(d = dist(x = DEG_mat)))
  gene_order <- order.dendrogram(dendro)
  gene_levels <- rownames(DEG_mat)[gene_order]
  DEG_df$GENE <- factor(DEG_df$GENE, levels = rev(gene_levels))
  
  
  DEG_df$`p(Adj) < 0.05` = DEG_df$p_val_adj < 0.05
  
  shapes = c(18,16)
  
  DEG_df$avg_logFC <- MinMax(DEG_df$avg_logFC, min = -.7, max = .7)
  
  dot <- ggplot(DEG_df, aes_string(x="Cell_type", y="GENE", size="neglog10fdr", color="avg_logFC")) + 
    geom_point(aes(shape=`p(Adj) < 0.05`)) + 
    #geom_point(data=selected[selected$Stringent=="Yes",],pch=21,colour="black", stroke=1) +
    theme_bw() +
    facet_grid(.~Tissue, scales = "free", space = "free") +
    theme(panel.border = element_blank(), #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, margin = margin(t=0), size=14, hjust = 1,vjust = 1),
          axis.text.y = element_text(size = 14, margin = margin(l=150)),
          #axis.text = element_text(face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_text(size=15),
          legend.text = element_text(size=10),
          plot.margin = margin(t=0),
          legend.position = "bottom") +
    theme(strip.background = element_blank(), #remove background for facet labels
          panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
          #panel.spacing = unit(0, "lines"),
          #strip.placement = "outside",
          strip.text = element_blank())+
    # strip.text = element_text(margin = margin(t=0,b=3,r=0,l=0),
    #                           hjust = .5, size=14)) +
    labs(size=expression(-log[10]~p["Adj"]), color=expression(logFC), shape=expression(p["Adj"] < 0.05)) +
    scale_color_gradient2(low = colors[1],  mid=colors[2], high = colors[3],
                          breaks = c(-.5, 0, .5), labels = c(-.5, 0, .5)) +
    scale_size(range = c(3, 9), breaks = c(1.3,25,50,75,100)) +
    scale_shape_manual(values=shapes) +
    guides(shape = guide_legend(override.aes = list(size = 5)))
  
  gb1 <- ggplot_build(heat)
  gb2 <- ggplot_build(dot)
  n1 <- length(gb1$panel$ranges[[1]]$y.labels)
  n2 <- length(gb2$panel$ranges[[1]]$y.labels)
  gA <- ggplot_gtable(gb1)
  gB <- ggplot_gtable(gb2)
  #g <- gtable:::rbind_gtable(gA, gB, "last")
  g <- gtable:::rbind_gtable(gA, gB)
  
  return(g)
  
}

gene_dot_plot_custom <- function(DEG_df, genes, min_max=c("min"=-1,"max"=1), 
                                 colors = c("blue","gray","red")){
  # cluster cell types based on DEGs
  DEG_df = DEG_df[DEG_df$GENE %in% genes,]
  DEG_mat <- convertToMatrix(df = DEG_df[DEG_df$Tissue=="Hippocampus",c("Cell_type","GENE","avg_logFC")])
  dendro <- as.dendrogram(hclust(d = dist(x = DEG_mat)))
  ct_order <- order.dendrogram(dendro)
  ct_levels <- rownames(DEG_mat)[ct_order]
  ct_levels <- ct_levels[ct_levels!="Choroid Plexus"]
  ct_levels <- c(ct_levels, "Choroid Plexus","Tanycyte")
  DEG_df$Cell_type <- factor(DEG_df$Cell_type, levels = ct_levels)
  
  # cluster genes
  DEG_mat <- convertToMatrix(df = DEG_df[,c("GENE","Ct_tissue","avg_logFC")])
  dendro <- as.dendrogram(hclust(d = dist(x = DEG_mat)))
  gene_order <- order.dendrogram(dendro)
  gene_levels <- rownames(DEG_mat)[gene_order]
  DEG_df$GENE <- factor(DEG_df$GENE, levels = rev(gene_levels))
  
  DEG_df$neglog10fdr = -log10(DEG_df$p_val_adj)
  DEG_df$neglog10fdr[DEG_df$neglog10fdr==Inf] <- 301
  DEG_df$`p(Adj) < 0.05` = DEG_df$p_val_adj < 0.05
  
  shapes = c(18,16)
  
  DEG_df$avg_logFC <-  MinMax(DEG_df$avg_logFC, min = min_max["min"], max = min_max["max"])
  
  dot <- ggplot(DEG_df, aes_string(x="Cell_type", y="GENE", size="neglog10fdr", color="avg_logFC")) + 
    geom_point(aes(shape=`p(Adj) < 0.05`)) + 
    #geom_point(data=selected[selected$Stringent=="Yes",],pch=21,colour="black", stroke=1) +
    theme_bw() +
    facet_grid(.~Tissue, scales = "free", space = "free") +
    theme(panel.border = element_blank(), #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, margin = margin(t=0), size=14, hjust = 1,vjust = 1),
          axis.text.y = element_text(size = 14, margin = margin(l=0)),
          #axis.text = element_text(face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_text(size=15),
          legend.text = element_text(size=10),
          plot.margin = margin(t=0),
          legend.position = "bottom") +
    theme(strip.background = element_blank(), #remove background for facet labels
          panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
          #panel.spacing = unit(0, "lines"),
          #strip.placement = "outside",
          #strip.text = element_blank())+
    strip.text = element_text(margin = margin(t=0,b=3,r=0,l=0),
                              hjust = .5, size=14)) +
    labs(size=expression(-log[10]~p["Adj"]), color=expression(logFC), shape=expression(p["Adj"] < 0.05)) +
    scale_color_gradient2(low = colors[1],  mid=colors[2], high = colors[3],
                          breaks = c(-.5, 0, .5), labels = c(-.5, 0, .5)) +
    scale_size(range = c(3, 9), breaks = c(1.3,25,50,75,100)) +
    scale_shape_manual(values=shapes) +
    guides(shape = guide_legend(override.aes = list(size = 5)))
  
  return(dot)
}



concordanceTable <- function(DEG1, DEG2, names = c("Study1","Study2")){
  
}

element_custom <- function(...) {
  structure(list(...), class = c("element_custom", "element_blank"))
}

# x axis labels
element_grob.element_custom <- function(element, label, x, y, ...)  {
  tg <- textGrob(label, x=x, gp=gpar(col=element$colour), 
                 rot = 45, just = "right", y=unit(-.25,"line"))
  padding <- unit(1,"line")
  rg <- rectGrob(x=x,width=unit(1,"line"), height=grobHeight(tg) + padding , 
                 gp=gpar(fill = element$fill, col=NA, alpha=.5, rot=45),
                 hjust = .5,vjust=1)
  gTree(children=gList(rg, tg), height=grobHeight(tg) + padding, cl="custom_axis")
}

# y axis labels
element_grob.element_custom <- function(element, label, x, y, ...)  {
  tg <- textGrob(label, y=y, gp=gpar(col=element$colour), 
                 rot = 0, just="right",x=unit(-.5,"line"))
  padding <- unit(1,"line")
  rg <- rectGrob(y=y,width=unit(1,"line"), height=grobHeight(tg) + padding , 
                 gp=gpar(fill = element$fill, col=NA, alpha=.5),
                 hjust = .5,vjust=1)
  gTree(children=gList(rg, tg), height=grobHeight(tg) + padding, cl="custom_axis")
}

element_grob.element_custom <- function(element, label, x, y, ...)  {
  tg <- textGrob(label, y=y, gp=gpar(col=element$colour), just = "right",
                 x=unit(-.1,"line"))
  padding <- unit(1,"line")
  rg <- rectGrob(y=y,width=grobWidth(tg), height=unit(1,"line"), 
                 gp=gpar(fill = element$fill, col=NA, alpha=0.4), just = "right",
                 x=unit(-.1,"line"))
  gTree(children=gList(rg, tg), width=grobWidth(tg) + padding, cl="custom_axis")
}
#widthDetails.custom_axis <- function(x) x$width #+ unit(2,"mm")

cellTypeFills = c("#33A02C","#1F78B4","#FF7F00","#984EA3",
                  "#E31A1C","#5E4FA2","#FB9A99","#FDBF6F")

gene_dot_plot_ct_specific <- function(DEG_df,genes, ct_levels, 
                                      min_max=c("min"=-1,"max"=1),
                                      colors = c("blue","gray","red"),
                                      cellTypeFills){
  # cluster cell types based on DEGs
  DEG_df = DEG_df[DEG_df$GENE %in% genes,]
  DEG_df = DEG_df[DEG_df$Cell_type %in% ct_levels,]
  
  DEG_df$neglog10fdr = -log10(DEG_df$p_val_adj)
  DEG_df$neglog10fdr[DEG_df$neglog10fdr==Inf] <- 301
  DEG_df$`p(Adj) < 0.05` = DEG_df$p_val_adj < 0.05
  
  shapes = c(18,16)
  
  DEG_df$avg_logFC <-  MinMax(DEG_df$avg_logFC, min = min_max["min"], max = min_max["max"])
  
  DEG_df$Tissue <- gsub("Hippocampus","HP", DEG_df$Tissue)
  DEG_df$Tissue <- gsub("Hypothalamus","HYP", DEG_df$Tissue)
  
  DEG_df$Cell_type <- factor(DEG_df$Cell_type, levels = ct_levels)
  DEG_df$GENE <- factor(DEG_df$GENE, levels = rev(genes))
  
  dot <- ggplot(DEG_df, aes_string(x="Tissue", y="GENE", size="neglog10fdr", color="avg_logFC")) + 
    geom_point(aes(shape=`p(Adj) < 0.05`)) + 
    #geom_point(data=selected[selected$Stringent=="Yes",],pch=21,colour="black", stroke=1) +
    theme_bw() +
    #facet_grid(Cell_type~., switch = "y") +
    facet_grid(.~Cell_type, scales = "free_y", space = "free") +
    theme(panel.border = element_blank(), #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 0, margin = margin(t=0), size=12, hjust = .5,vjust = 1),
          #axis.text.x = element_custom(fill=cellTypeFills),
          #axis.text.y = element_text(size = 13, margin = margin(l=0)),
          axis.text.y = element_custom(fill=cellTypeFills),
          #axis.text.y = element_blank(),
          #axis.text = element_text(face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_text(size=15),
          legend.text = element_text(size=10),
          plot.margin = margin(l=60),
          legend.position = "bottom") +
    theme(#strip.background = element_blank(), #remove background for facet labels
          panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
          #panel.spacing = unit(0, "lines"),
          #strip.placement = "outside",
          #strip.text = element_blank())+
          strip.placement = "outside",
          strip.text.x = element_text(angle = 90, hjust = 0, vjust = .5, size=13))+
          #strip.text.y.left=element_text(angle = 0, hjust=1, size=15)) +
    labs(size=expression(-log[10]~p["Adj"]), color=expression(logFC), shape=expression(p["Adj"] < 0.05)) +
    scale_color_gradient2(low = colors[1],  mid=colors[2], high = colors[3],
                          breaks = c(-.5, 0, .5), labels = c(-.5, 0, .5)) +
    scale_size(range = c(3, 9), breaks = c(1.3,25,50,75,100)) +
    scale_shape_manual(values=shapes) +
    guides(shape = guide_legend(override.aes = list(size = 5)))+
    scale_x_discrete(position = "top") 
  
  g <- ggplot_gtable(ggplot_build(dot))
  strip_both <- which(grepl('strip', g$layout$name))
  fills <- rev(unique(cellTypeFills))
  #cols <- c(rep("black",12), rep("white",3))
  k <- 1
  for (i in strip_both) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    #l <- which(grepl('title', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    g$grobs[[i]]$grobs[[1]]$children[[j]]$width <- unit(2,"lines")
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$alpha <- .4
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$lwd <- 0
    #m <- which(grepl('text',g$grobs[[i]]$grobs[[1]]$children[[l]]$childrenOrder))
    #g$grobs[[i]]$grobs[[1]]$children[[l]]$children[[m]]$gp$col <- cols[k]
    k <- k+1
  }
  

  
  #grid.draw(g)
  return(g)
}

gene_dot_plot_ct_tiss_specific <- function(DEG_df,genes, ct_levels, 
                                      min_max=c("min"=-1,"max"=1),
                                      colors = c("blue","gray","red"),
                                      cellTypeFills){
  
  DEG_df = DEG_df[DEG_df$GENE %in% genes,]
  DEG_df = DEG_df[DEG_df$Cell_type %in% ct_levels,]
  
  DEG_df$neglog10fdr = -log10(DEG_df$p_val_adj)
  DEG_df$neglog10fdr[DEG_df$neglog10fdr==Inf] <- 301
  DEG_df$`p(Adj) < 0.05` = DEG_df$p_val_adj < 0.05
  
  shapes = c(18,16)
  
  DEG_df$avg_logFC <-  MinMax(DEG_df$avg_logFC, min = min_max["min"], max = min_max["max"])
  
  DEG_df$Tissue <- gsub("Hippocampus","Hip", DEG_df$Tissue)
  DEG_df$Tissue <- gsub("Hypothalamus","Hyp", DEG_df$Tissue)
  
  DEG_df$Cell_type <- factor(DEG_df$Cell_type, levels = ct_levels)
  DEG_df$GENE <- factor(DEG_df$GENE, levels = rev(genes))
  
  dot <- ggplot(DEG_df, aes_string(x="Tissue", y="GENE", size="neglog10fdr", color="avg_logFC")) + 
    geom_point(aes(shape=`p(Adj) < 0.05`)) + 
    #geom_point(data=selected[selected$Stringent=="Yes",],pch=21,colour="black", stroke=1) +
    theme_bw() +
    #facet_grid(Cell_type~., switch = "y") +
    facet_grid(.~Cell_type, scales = "free_y", space = "free") +
    theme(panel.border = element_blank(), #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 0, margin = margin(t=0), size=12, hjust = .5,vjust = 1),
          #axis.text.x = element_custom(fill=cellTypeFills),
          axis.text.y = element_text(size = 13, margin = margin(l=0)),
          #axis.text.y = element_custom(fill=cellTypeFills),
          #axis.text.y = element_blank(),
          #axis.text = element_text(face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_text(size=15),
          legend.text = element_text(size=10),
          plot.margin = margin(l=5),
          legend.position = "bottom") +
    theme(strip.background = element_blank(), #remove background for facet labels
      panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
      #panel.spacing = unit(0, "lines"),
      #strip.placement = "outside",
      #strip.text = element_blank())+
      strip.placement = "outside",
      strip.text.x = element_text(angle = 90, hjust = 0, vjust = .5, size=13))+
    #strip.text.y.left=element_text(angle = 0, hjust=1, size=15)) +
    labs(size=expression(-log[10]~p["Adj"]), color=expression(logFC), shape=expression(p["Adj"] < 0.05)) +
    scale_color_gradient2(low = colors[1],  mid=colors[2], high = colors[3],
                          breaks = c(-.5, 0, .5), labels = c(-.5, 0, .5)) +
    scale_size(range = c(3, 9), breaks = c(1.3,25,50,75,100)) +
    scale_shape_manual(values=shapes) +
    guides(shape = guide_legend(override.aes = list(size = 5)))+
    scale_x_discrete(position = "top") 
  
  
  #grid.draw(g)
  return(dot)
}
logFC <- read.delim("./DEGs/HP_HYP_DEGs_new.txt")
DEG_df <- logFC



addInfoAndTrim <- function(DEG_df, 
                           FDR_cutoff=0.05, 
                           multiple_comparisons=TRUE, 
                           convertToHuman=FALSE, 
                           lfc_cutoff=0.1){
  DEG_df = DEG_df[DEG_df$p_val_adj<FDR_cutoff,]
  DEG_df = DEG_df[abs(DEG_df$avg_logFC)>lfc_cutoff,]
  module = c()
  module_ct = c()
  module_comp = c()
  for(row in 1:nrow(DEG_df)){
    if(DEG_df$avg_logFC[row]>0){
      if(multiple_comparisons){
        module[row] = paste(DEG_df$Comparison[row], DEG_df$Cell_type[row], "UP", sep = ".")
        module_ct[row] = paste(DEG_df$Cell_type[row], "UP", sep = ".")
        module_comp[row] = paste(DEG_df$Comparison[row], "UP", sep = ".")
      }
      else{
        module[row] = paste(DEG_df$Cell_type[row], "UP", sep = "_")
      }
    }
    else{
      if(multiple_comparisons){
        module[row] = paste(DEG_df$Comparison[row], DEG_df$Cell_type[row], "DOWN", sep = ".")
        module_ct[row] = paste(DEG_df$Cell_type[row], "DOWN", sep = ".")
        module_comp[row] = paste(DEG_df$Comparison[row], "DOWN", sep = ".")
      }
      else{
        module[row] = paste(DEG_df$Cell_type[row], "DOWN", sep = ".")
      }
    }
  }
  DEG_df$MODULE = module
  DEG_df$MODULE_ct = module_ct
  DEG_df$MODULE_comp = module_comp
  DEG_df$MODULE_no_direct = paste(DEG_df$Cell_type, DEG_df$Comparison)
  DEG_df = DEG_df[order(DEG_df$MODULE),]
  
  if(convertToHuman){
    DEG_human = DEG_df
    convertedToHuman <- convertMouseGeneList(DEG_human$GENE)
    DEG_human$GENE <- convertedToHuman$HGNC.symbol[match(DEG_human$GENE, convertedToHuman$MGI.symbol)]
    DEG_human$GENE[which(is.na(DEG_human$GENE))] <- toupper(DEG_df$GENE[which(is.na(DEG_human$GENE))])
    DEG_df = DEG_human
    rm(DEG_human)
    MT_gene_names = c("MT-ND1"="ND1",  "MT-ND2"="ND2",  "MT-CO1"="COX1",  "MT-CO2"="COX2",  "MT-ATP8"="ATP8", 
                      "MT-ATP6"="ATP6", "MT-CO3"="COX3","MT-ND3" ="ND3", "MT-ND4L"="ND4L", "MT-ND4"="ND4",
                      "MT-ND5"="ND5",  "MT-ND6"="ND6",  "MT-CYB"="CYB")
    for(gene in 1:nrow(DEG_df)){
      if(grepl("^MT-",DEG_df$GENE[gene])){
        DEG_df$GENE[gene] = MT_gene_names[DEG_df$GENE[gene]]
      }
      else{
        next
      }
    }
  }
  return(DEG_df)
}

drawTopDEGsDotPlot <- function(DEG_df, 
                               FDR_threshold=0.05, 
                               logFC_threshold=.25,
                               comparisons=NULL,
                               add_human_ortholog_column=FALSE,
                               add_expression_levels_column=FALSE,
                               expression_data,
                               addGWAS = FALSE,
                               GWAS_catalog_data){
  
}


# if specifying no genes, then DEG_df must have 'GENE' and 'Cell_type' columns
# and it will display only genes meeting 
gene_dot_plot <- function(
  DEG_df,
  vars = c("GENE","Cell_type"), 
  var_list = NULL, 
  num_genes_show = 30, # if var_list "GENE" not specified
  colors = c("blue","gray","red"),
  FDR_threshold = 0.05, # only used if no genes specified
  logFC_threshold = 0.25 # only used if no genes specified
  #logFC_min_max = c(-1,1)
){
  DEG_df$neglog10fdr = -log10(DEG_df$p_val_adj)
  DEG_df$neglog10fdr[DEG_df$neglog10fdr==Inf] <- 301
  
  selected = data.frame(stringsAsFactors = FALSE)
  if(!is.null(var_list)){ # user specifies specific variables to show
    temp = DEG_df[DEG_df[,names(var_list)[1]] %in% var_list[[1]],]
    if(length(var_list)>1){
      selected = temp[temp[,names(var_list)[2]] %in% var_list[[2]],]
    }
    else{
      selected = temp
    }
    for(v in 1:length(var_list)){
      if(length(unique(selected[,vars[v]]))!=length(var_list[v])){
        selected[,vars[v]] = factor(selected[,vars[v]], levels = unique(selected[,vars[v]]))
      }
      else{
        selected[,vars[v]] = factor(selected[,vars[v]], levels = var_list[[v]])
      }
    }
  }
  else if(length(unique(DEG_df$GENE))==1 & is.null(var_list)){ # single gene plot with multiple comparisons
    selected = DEG_df
  }
  else{
    # select the top 30 most "regulated" genes, show all cell types
    # subset DEG df to only significant DEGs
    DEG_df <- DEG_df[DEG_df$p_val_adj<FDR_threshold & abs(DEG_df$avg_logFC)>logFC_threshold,]
    DEG_df$numAppearances <- sapply(DEG_df$GENE, function(x){return(sum(DEG_df$GENE==x))})
    DEG_df <- DEG_df[order(DEG_df$numAppearances, decreasing = TRUE),]
    if(length(unique(DEG_df$GENE))<num_genes_show){
      cat("Number of significant genes less than specified num_genes_show.\n")
      cat("Proceeding to show ", length(unique(DEG_df$GENE)), " genes.\n")
      genes = unique(DEG_df$GENE)[1:length(unique(DEG_df$GENE))]
    }
    else{
      genes = unique(DEG_df$GENE)[1:num_genes_show]
    }
    selected = DEG_df[DEG_df$GENE %in% genes,]
    selected$GENE <- factor(selected$GENE, levels=genes)
  }
  
  selected$`p(Adj) < 0.05` = selected$p_val_adj < 0.05
  
  if(sum(selected$p_val_adj>0.05)==0){ 
    shapes = c(16)
  }
  else{
    shapes = c(18,16)
  }
  
  #selected$avg_logFC <- MinMax(selected$avg_logFC, min = logFC_min_max[1], max = logFC_min_max[2])
  
  D <- ggplot(selected, aes_string(x=vars[1], y=vars[2], size="neglog10fdr", color="avg_logFC")) + 
    geom_point(aes(shape=`p(Adj) < 0.05`)) + 
    #geom_point(data=selected[selected$Stringent=="Yes",],pch=21,colour="black", stroke=1) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, margin = margin(t=5), size=14, hjust = 1,vjust = 1),
          #axis.text.x = element_text(angle = 90, size=14, hjust = 0,vjust = 1),
          axis.text.y = element_text(size = 14),
          axis.text = element_text(face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_text(size=15),
          legend.text = element_text(size=10),
          legend.position = "bottom") +
    labs(size=expression(-log[10]~p["Adj"]), color=expression(logFC)) +
    scale_color_gradient2(low = colors[1],  mid=colors[2], high = colors[3]) +
                          #breaks = c(-.5, 0, .5), labels = c(-.5, 0, .5)) +
    scale_size(range = c(3, 12), breaks = c(1.3,25,50,75,100)) +
    scale_shape_manual(values=shapes) +
    guides(shape = guide_legend(override.aes = list(size = 5))) 
    #scale_x_discrete(position = "top")
  

  
  return(D)
}

reportCellTypeDiscordance <- function(DEG_df, ct_column="Cell_type",
                                      logFC_threshold=.25, FDR_threshold=0.05){
  DEG_df = DEG_df[abs(DEG_df$avg_logFC)>logFC_threshold & DEG_df$p_val_adj<FDR_threshold,]
  # narrow down to "multiple perturbed subset"
  result = data.frame(stringsAsFactors = FALSE)
  genes_perturbed = c()
  for(gene in unique(DEG_df$GENE)){
    if(sum(DEG_df$GENE==gene)>1){
      # to see if there are any varying directions
      logFCs =  DEG_df$avg_logFC[DEG_df$GENE==gene]
      if((sum(logFCs>0)>0 & sum(logFCs>0)!=length(logFCs)) | 
         (sum(logFCs<0)>0 & sum(logFCs<0)!=length(logFCs))){
        genes_perturbed = append(genes_perturbed, gene)
      }
      else next
    }
    else next
  }
  celltype_result = data.frame("GENE"=genes_perturbed, stringsAsFactors = FALSE)
  Comparisons = unique(DEG_df[,ct_column])
  for(comp in Comparisons){
    vect_FC = c()
    for(g in genes_perturbed){
      if(sum(DEG_df$GENE[DEG_df[,ct_column]==comp]==g)>0){
        vect_FC = append(vect_FC, DEG_df$avg_logFC[DEG_df[,ct_column]==comp & DEG_df$GENE==g])
      }
      else{
        vect_FC = append(vect_FC, "None") # may have fc but not significant
      }
    }
    celltype_result[,comp] = vect_FC
  }
  result = rbind(result, celltype_result)
  return(result)
}

# df is from describe_manulation() or reportCellTypeDiscordance()
addDescriptiveColumn <- function(df){
  # go through the entries in one row 
  descriptions <- c()
  for(row in 1:nrow(df)){
    col_descrips = c()
    for(col in 2:ncol(df)){
      if(df[row,col]!="None"){
        col_descrips = c(col_descrips, 
                         paste0(colnames(df)[col],"_",ifelse(df[row,col]>0, "UP","DOWN")))
      }
    }
    descriptions = c(descriptions, concatenate(col_descrips, mysep = ", "))
  }
  df$Description <- descriptions
  return(df)
}

prepareDEGsGmtFormat <- function(DEG_df, FDR_threshold=0.05, 
                                 lfc_threshold=0.1,
                                 convertToHuman=TRUE, name="na"){
  DEG_df  <- addInfoAndTrim(DEG_df = DEG_df, 
                               FDR_cutoff = FDR_threshold, 
                               multiple_comparisons = FALSE, 
                               convertToHuman = TRUE, 
                               lfc_cutoff = lfc_threshold)
  mat = matrix(ncol = max(table(DEG_DF$MODULE)) + 2,
               nrow = length(unique(DEG_DF$MODULE)))
  for(mod in 1:nrow(mat)){
    genes = DEG_DF$HUMAN[DEG_DF$MODULE==unique(DEG_DF$MODULE)[mod]]
    mat[mod,] = c(unique(DEG_DF$MODULE)[mod],name, genes,
                  rep("", (ncol(mat)-length(genes))-2))
  }
}

getHeirGeneOrder <- function(df){
  mat <- convertToMatrix(df)
  dendro <- as.dendrogram(hclust(d = dist(x = mat)))
  # just to get order
  ct_order <- order.dendrogram(dendro)
  ct_levels <- rownames(mat)[ct_order]
  return(ct_levels)
}

# color_by is column that contains meta data you want to color by
# if no color_by set, upregulated genes are red and downregulated genes are blue
# as written, cannot be more than 8 levels. You shouldn't visualize more than 8 anyway!
volcano_plot <- function(DEG_df, 
                         cell_type=NULL, 
                         color_by=NULL, 
                         colors = NULL, 
                         lbl_size = 3,
                         lfc = .25,
                         padjusted = 0.001){
  if(!is.null(cell_type)){
    DEG_df = DEG_df[DEG_df$Cell_type==cell_type,]
  }
  
  DEG_df$log10_pval_adj = -log10(DEG_df$p_val_adj)
  DEG_df$Sign <- (abs(DEG_df$avg_logFC) > lfc & DEG_df$log10_pval_adj > -log10(padjusted))
  DEG_df$GENE <- as.character(DEG_df$GENE)
  
  DEG_df = DEG_df[!grepl("^Rps", DEG_df$GENE),]
  DEG_df = DEG_df[!grepl("^Rpl", DEG_df$GENE),]
  
  if(!is.null(color_by)){
    if(is.null(colors)){
      if(length(unique(DEG_df[,color_by]))==2){
        colors = brewer.pal(3,"Set1")[1:2]
      } else{
        colors = brewer.pal(length(unique(DEG_df[,color_by])),"Set1")
      }
    }
    
    pal = colors
    names(pal) = unique(DEG_df[,color_by])
    
    v <- ggplot(DEG_df) + geom_point(aes_string(x="avg_logFC", y="log10_pval_adj", color=color_by)) +
      geom_hline(yintercept=-log10(padjusted), linetype="dashed", color = "gray") +
      geom_vline(xintercept = -lfc, color="gray",linetype = "dashed") +
      geom_vline(xintercept = lfc, color="gray", linetype = "dashed") +
      ylab(expression(-log[10]~FDR)) +
      xlab("logFC") + 
      geom_text_repel(aes(x = avg_logFC, y = log10_pval_adj, label = ifelse(Sign == T, GENE,"")),
                      point.padding = 0.5, segment.color = 'grey20',segment.size = 0.05,size = lbl_size) +
      theme_bw() +
      theme(legend.position = "right", panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.title = element_text(size=15), title = element_text(size=15)) +
      scale_color_manual(values = pal, limits = names(pal))
    return(v)
  } else{
    pal <- c(
      "UP" = "red",
      "DOWN" = "blue",
      "Not significant"="black"
    )
    direction = c()
    for(i in 1:length(DEG_df$p_val)){
      if(DEG_df$log10_pval_adj[i]>-log10(padjusted) & DEG_df$avg_logFC[i]>0 & abs(DEG_df$avg_logFC[i])>lfc) direction[i] = "UP"
      else if(DEG_df$log10_pval_adj[i]>-log10(padjusted) & DEG_df$avg_logFC[i]<0 & abs(DEG_df$avg_logFC[i])>lfc) direction[i] = "DOWN"
      else direction[i] = "Not significant"
    }
    DEG_df$direction = direction
    
    v <- ggplot(DEG_df) + geom_point(aes(x=avg_logFC, y=log10_pval_adj, color=direction)) +
      geom_hline(yintercept=-log10(padjusted), linetype="dashed", color = "gray") +
      geom_vline(xintercept = -lfc, color="gray",linetype = "dashed") +
      geom_vline(xintercept = lfc, color="gray", linetype = "dashed") +
      ylab(expression(-log[10]~FDR)) +
      xlab(expression(log[2]~FC)) + 
      geom_text_repel(aes(x = avg_logFC, y = log10_pval_adj, label = ifelse(Sign == T, GENE,"")),
                      point.padding = 0.5, segment.color = 'grey20',segment.size = 0.05,size = lbl_size) +
      theme_bw() +
      theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.title = element_text(size=15), title = element_text(size=15)) +
      scale_color_manual(values =pal, limits = names(pal))
    return(v)
  }
}


addStringency <- function(DEG_df, stringentDEG_df){
  # add Stringent column with "Yes"
  stringent = c()
  for(i in 1:nrow(DEG_df)){
    if(sum(rowSums(cbind(DEG_df$Cell_type[i]==stringentDEG_df$Cell_type,
                         DEG_df$GENE[i]==stringentDEG_df$GENE,
                         DEG_df$Comparison[i]==stringentDEG_df$Comparison))==3)>0){
      stringent[i] = "Yes"
    }
    else{
      stringent[i] = "No"
    }
  }
  DEG_df$Stringent = stringent
  return(DEG_df)
}

















