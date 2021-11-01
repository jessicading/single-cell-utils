# Sorry very messy 

library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
library(scales)
library(viridis)
library(varhandle)
library(Seurat)
#library(VennDiagram)
library(RColorBrewer)
#library(pheatmap)
library(grid)
library(ggpubr)
#library(ggfittext)
#library(stringr)
#library(ggrepel)
#library(ComplexHeatmap)
#library(WriteXLS)
#library(readxl)
library(metap)


# meta_data_1 will be in the rows, meta_data_2 will be in the columns
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
    colnames(total_idents) = meta_data_1
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
                                                 meta_data_1 = "Sample",
                                                 meta_data_2 = "Cell_type",
                                                 addCounts = FALSE, 
                                                 cts_include = NULL, 
                                                 custom_colors = NULL,
                                                 pal = NULL,
                                                 coord_flip = FALSE,
                                                 meta_data_2_levels = NULL, 
                                                 order_by_abundance=FALSE,
                                                 add_outline=FALSE){
  # melt data
  melted <- melt(data, id.vars = colnames(data)[1])
  colnames(melted) <- c(meta_data_1, meta_data_2, "Counts")
  
  melted$Counts <- as.numeric(melted$Counts) 
  # add percentage as feature
  percentages = c()
  melted = melted[order(melted[,meta_data_2]),]
  for(samp in unique(melted[,meta_data_2])){
    total = sum(melted$Counts[melted[,meta_data_2]==samp])
    for(ct in unique(melted[,meta_data_1])){
      percentages = append(percentages, melted$Counts[melted[,meta_data_2]==samp & melted[,meta_data_1]==ct]/total)
    }
  }
  melted$Percentage = percentages*100
  
  # order by most abundant cell type
  #     add up cell type counts across all samples
  if(order_by_abundance){
    count_totals = c()
    for(ct in unique(melted[,meta_data_1])){
      count_totals = append(count_totals, sum(as.numeric(melted$Counts[melted[,meta_data_1]==ct])))
    }
    cell_type_order = unique(melted[,meta_data_1])[order(count_totals, decreasing = TRUE)]
    
    melted[,meta_data_1] <- factor(melted[,meta_data_1], levels = rev(cell_type_order))
    
    # reorder pal if it exists
    pal <- pal[cell_type_order]
    
  }
  
  if(addCounts){
    df <- melted %>%
      group_by_at(meta_data_2) %>%
      arrange(!!! rlang::syms(meta_data_2), desc(!!! rlang::syms(colnames(data)[1]))) %>%
      mutate(lab_ypos = cumsum(Percentage) - 0.5 * Percentage)
    
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
      df <- as.data.frame(df)
      df[,meta_data_2] <- as.character(df[,meta_data_2])
      df[,meta_data_2] <- factor(df[,meta_data_2], levels = as.character(meta_data_2_levels))
    } 
    
    g <- ggplot(df, aes_string(x=paste0("`",meta_data_2,"`"), y="Percentage"))
    
    if(add_outline){
      g <- g + geom_col(aes_string(fill=colnames(data)[1]), colour="black")
    } else {
      g <- g + geom_col(aes_string(fill=colnames(data)[1]))
    }
    
    g <- g + geom_text(aes(y=lab_ypos, label=Counts, group = colnames(data)[1]), color = "black") + 
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            axis.line = element_line(),
            panel.grid.minor = element_blank(), legend.position = "right",
            #axis.text.y = element_text(margin = margin(r=5), size=14, face = "bold"),
            axis.text.x = element_text(angle=90, vjust=1, hjust=0.5, size=12),
            #axis.ticks = element_blank(),
            #axis.title.x = element_blank(),
            #axis.title.y = element_blank(),
            #legend.title = element_text(face = "bold", size=13)) +
            legend.title = element_blank()) +
      guides(fill=guide_legend(title=gsub("_"," ",meta_data_1))) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 102))
  } 
  else{
    g <- ggplot(melted, aes_string(x=paste0("`",meta_data_2,"`"), y="Percentage")) 
    if(add_outline){
      g <- g + geom_col(aes_string(fill=colnames(melted)[1]), colour="black")
    } else {
      g <- g + geom_col(aes_string(fill=colnames(melted)[1]))
    }
    g <- g + theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(),
            legend.position = "right",
            #axis.text.y = element_text(margin = margin(r=5), size=14, face = "bold"),
            axis.text.x = element_text(angle=90, vjust=.5, size=12)) +
            #axis.ticks = element_blank(),
            #axis.title.x = element_blank(),
            #axis.title.y = element_blank()) +
            #legend.title = element_text(face = "bold", size=13)) +
            #legend.title = element_blank()) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 102)) +
      guides(fill=guide_legend(title=gsub("_"," ",meta_data_1))) #+
      #xlab("Cell Type") 
    # pdf("./Potential_Paper_Plots/Figure1/Sample_Cell_Type_Proportions_HYP.pdf", height=4, width = 4.3)
    # g
    # dev.off()
  }
  
  if(!is.null(pal)){
    g <- g + scale_fill_manual(values = pal,limits = names(pal))
  }
  
  if(coord_flip){
    g <- g + coord_flip()
  }
  #results <- list(g, df)
  return(g)
}

# MODULE contains all of the meta data information without direction info and can be used for pathway enrichment
addInfoAndTrim <- function(DEG_df, FDR_cutoff=0.05, meta_data, convertToHuman=FALSE, lfc_cutoff=0.1){
  if(!("avg_logFC" %in% colnames(DEG_df))){
    fc_type <- "avg_log2FC"
  } else{
    fc_type <- "avg_logFC"
  }
  DEG_df = DEG_df[DEG_df$p_val_adj<FDR_cutoff,]
  DEG_df = DEG_df[abs(DEG_df[,fc_type])>lfc_cutoff,]
  DEG_df$Direction = ifelse(DEG_df[,fc_type]>0,"UP","DOWN")
  combinations = c()
  for(iter in 1:length(meta_data)){
    for(j in 1:ncol(combn(meta_data, m = iter))){
      combinations = c(combinations, concatenate(combn(meta_data, m = iter)[,j], mysep = "."))
    }
  }
  for(comb in setdiff(combinations, meta_data)){
    metas = unlist(strsplit(comb, split = ".", fixed = TRUE))
    for(iter in 1:length(metas)){
      if(iter==1) DEG_df[,comb] = DEG_df[,metas[iter]]
      else DEG_df[,comb] = paste(DEG_df[,comb], DEG_df[,metas[iter]], sep = ".")
    }
  }
  for(comb in combinations){
    DEG_df[,paste0(comb, ".Direction")] = paste(DEG_df[,comb], DEG_df$Direction, sep = ".")
  }
  
  # assuming that for pathway enrichment, the group to enrich will be the 'most descriptive' (containing all
  # info provided in the meta_data)
  longest = combinations[order(nchar(combinations), decreasing = TRUE)][1]
  DEG_df$MODULE = DEG_df[,longest] # basically making a duplicate of this column in the 'MODULE' column
  longest = paste0(longest, ".Direction")
  DEG_df$MODULE_direction = DEG_df[,longest]
  
  
  if(convertToHuman){
    temp <- convertDfGeneColumnMouseHuman(df = DEG_df, toSpecies = "human", forPathway = FALSE)
    DEG_df$HUMAN = temp$GENE
  }
  return(DEG_df)
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
                                           metapvalType = "minPadj",
                                           nTotalSamples = 6,
                                           output_name="Global_SampleWise_DEG_Comparison.txt"){
  if(!("avg_logFC" %in% colnames(global_DEGs))){
    fc_type <- "avg_log2FC"
  } else{
    fc_type <- "avg_logFC"
  }
  quantified = data.frame(stringsAsFactors = FALSE)
  all_sample_wise_DEGs = sample_wise_DEGs
  # geneDetails = data.frame(stringsAsFactors = FALSE)
  for(comp in unique(sample_wise_DEGs$Comparison)){
    sample_wise_DEGs <- all_sample_wise_DEGs[all_sample_wise_DEGs$Comparison==comp,]
    for(ct in unique(sample_wise_DEGs$Cell_type)){
      for(lfc in lfc_cutoffs){
        for(fdr in pval_adj_cutoffs){
          GenesGlobal = global_DEGs$GENE[global_DEGs$Cell_type==gsub("_"," ",ct) & 
                                           abs(global_DEGs[,fc_type])>lfc &
                                           global_DEGs$p_val_adj<fdr]
          GenesSampleWise = sample_wise_DEGs$GENE[sample_wise_DEGs$Cell_type==ct &
                                                    sample_wise_DEGs[,metapvalType]<fdr &
                                                    rowSums(abs(sample_wise_DEGs[,grep(fc_type, 
                                                                                       colnames(sample_wise_DEGs))])>lfc)>=nTotalSamples]
          temp_q = data.frame("Cell Type" = ct, 
                              "Comparison"=comp,
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
  }
  if(!is.null(output_name)){
    write.table(quantified, output_name, 
                row.names = FALSE,sep = "\t",quote = FALSE)
  }
  return(quantified)
}

# unlike compareGlobalAndSampleWiseDEGs(), does not go by the metap value because just filters based on
# sample wise p_val_adj
compareGlobalAndSampleWiseDEGsConsistent <- function(global_DEGs, sample_wise_DEGs, 
                                           lfc_cutoffs = c(0.1,0.25), 
                                           pval_adj_cutoffs = c(0.05,0.01),
                                           conditions = c("5XFAD", "WT"),
                                           nTotalSamples = 6,
                                           nSampleSignificance=6){
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
                                                       rowSums(sample_wise_DEGs_melt[,grep("p_val_adj", 
                                                                                           colnames(sample_wise_DEGs_melt))]<fdr)>=nSampleSignificance &
                                                       rowSums(abs(sample_wise_DEGs_melt[,grep("avg_logFC", 
                                                                                               colnames(sample_wise_DEGs_melt))])>lfc)>=nTotalSamples]
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
  
  if(!("avg_logFC" %in% colnames(global_DEGs))){
    fc_type <- "avg_log2FC"
  } else{
    fc_type <- "avg_logFC"
  }
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
                                                  rowSums(abs(sample_wise_DEGs_melt[,grep(fc_type, 
                                                                                          colnames(sample_wise_DEGs_melt))])>lfc_threshold)==nTotalSamples,
                                                c("Cell_type","GENE", pval_type)]
  
  # filter global DEGs for only those included in the sample_wise_DEGs_melt
  
  sample_wise_DEGs_melt$Cell_type <- gsub("_", " ", sample_wise_DEGs_melt$Cell_type)
  filtered_DEGs <- merge(sample_wise_DEGs_melt, global_DEGs, by = c("Cell_type","GENE"))
  
  return(filtered_DEGs)
}

# not based on the meta p type. based on the p adj for each sample DEG run
sampleWiseFilteredDEGsConsistent_old <- function(global_DEGs, sample_wise_DEGs, 
                                             nTotalSamples=6, nSampleSignificance=6, 
                                             fdr_threshold=0.1,
                                             lfc_threshold=.1, pval_type = "ConsP"){
  if(!("avg_logFC" %in% colnames(global_DEGs))){
    fc_type <- "avg_log2FC"
  } else{
    fc_type <- "avg_logFC"
  }
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
  sample_wise_DEGs_melt = sample_wise_DEGs_melt[rowSums(abs(sample_wise_DEGs_melt[,grep(fc_type, 
                                                                                        colnames(sample_wise_DEGs_melt))])>lfc_threshold)==nTotalSamples,]
  sample_wise_DEGs_melt = sample_wise_DEGs_melt[rowSums(sample_wise_DEGs_melt[,grep("p_val_adj", 
                                                                                      colnames(sample_wise_DEGs_melt))]<fdr_threshold)>=nSampleSignificance,
                                                c("Cell_type","GENE", pval_type)]
  
  sample_wise_DEGs_melt$Cell_type <- gsub("_", " ", sample_wise_DEGs_melt$Cell_type)
  filtered_DEGs <- merge(sample_wise_DEGs_melt, global_DEGs, by = c("Cell_type","GENE"))
  return(filtered_DEGs)
}

# some sample consistent DEGs are FDR=1 for normal deg analysis... idk why
sampleWiseFilteredDEGsConsistent <- function(global_DEGs, 
                                             sample_wise_DEGs, 
                                             nTotalSamples=6, 
                                             nSampleSignificance=6, 
                                             fdr_threshold=0.1,
                                             lfc_threshold=.1, 
                                             metapval_type = "ConsP"){
  
  if(!("avg_logFC" %in% colnames(global_DEGs))){
    fc_type <- "avg_log2FC"
  } else{
    fc_type <- "avg_logFC"
  }
  res <- data.frame()
  all_sample_wise_DEGs = sample_wise_DEGs
  all_global_DEGs <- global_DEGs
  for(comp in unique(sample_wise_DEGs$Comparison)){
    sample_wise_DEGs <- all_sample_wise_DEGs[all_sample_wise_DEGs$Comparison==comp,]
    global_DEGs <- all_global_DEGs[all_global_DEGs$Comparison==comp,]
    sample_wise_DEGs = sample_wise_DEGs[rowSums(abs(sample_wise_DEGs[,grep(fc_type, 
                                                                           colnames(sample_wise_DEGs))])>lfc_threshold)==nTotalSamples,]
    sample_wise_DEGs = sample_wise_DEGs[rowSums(sample_wise_DEGs[,grep("p_val_adj", 
                                                                       colnames(sample_wise_DEGs))]<fdr_threshold)>=nSampleSignificance,]
    
    filtered_DEGs <- merge(sample_wise_DEGs, global_DEGs, by = c("Cell_type","GENE","Comparison"))
    filtered_DEGs = filtered_DEGs[filtered_DEGs$p_val_adj<0.05,]
    res <- rbind(res, filtered_DEGs)
  }
  return(res)
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
                                 orderBy = "nOverlap", cluster_meta = "Cell_type",
                                 calculateOnlyOverlapLogFC = FALSE, 
                                 compiled_pathways,
                                 onlySig=FALSE){
  
  if(!("avg_logFC" %in% colnames(DEG_df))){
    fc_type <- "avg_log2FC"
  } else{
    fc_type <- "avg_logFC"
  }
  
  MT_gene_names = c("MT-ND1"="ND1",  "MT-ND2"="ND2",  "MT-CO1"="COX1",  "MT-CO2"="COX2",  "MT-ATP8"="ATP8", 
                    "MT-ATP6"="ATP6", "MT-CO3"="COX3","MT-ND3" ="ND3", "MT-ND4L"="ND4L", "MT-ND4"="ND4",
                    "MT-ND5"="ND5",  "MT-ND6"="ND6",  "MT-CYB"="CYB")
  
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
        logFCs = ct_logFC_all[,fc_type][match(x = genes, table = ct_logFC_all$GENE)]
        avg_logFC[row] = sum(logFCs[!is.na(logFCs)])
        if(length(logFCs)!=length(genes)){
          cat("Did not find equal number of logFC values in row ",row,"\n")
        }
      }
      else{
        avg_logFC[row] = 0
      }
    }
    pathway_df[,fc_type] = avg_logFC
    return(pathway_df)
  }
  else if(!calculateOnlyOverlapLogFC & onlySig){
    avg_logFC = c()
    for(row in 1:nrow(pathway_df)){
      pathway_genes = compiled_pathways$GENE[compiled_pathways$MODULE==pathway_df$Pathway[row]]
      ct_logFC_all = logFC_df[logFC_df[,cluster_meta]==pathway_df$Module[row] & logFC_df$p_val_adj<0.05,] # take only those that are significant
      logFCs = ct_logFC_all[,fc_type][match(x = pathway_genes, table = ct_logFC_all$GENE)]
      avg_logFC[row] = sum(logFCs[!is.na(logFCs)])
    }
    pathway_df[,fc_type] = avg_logFC
    return(pathway_df)
  }
  else{
    avg_logFC = c()
    for(row in 1:nrow(pathway_df)){
      pathway_genes = compiled_pathways$GENE[compiled_pathways$MODULE==pathway_df$Pathway[row]]
      ct_logFC_all = logFC_df[logFC_df[,cluster_meta]==pathway_df$Module[row],]
      logFCs = ct_logFC_all[,fc_type][match(x = pathway_genes, table = ct_logFC_all$GENE_human)]
      avg_logFC[row] = sum(logFCs[!is.na(logFCs)])
    }
    pathway_df[,fc_type] = avg_logFC
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

processSeurat <- function(seuratObject, 
                          combineLevel="orig.ident", 
                          name="single_cell_analysis", # should be tissue or descriptive name
                          skipIntegration = TRUE,
                          vars.to.regress = NULL,
                          resUse=seq(0.5,4,by=0.5),
                          nDims=30,
                          plot_meta=c("orig.ident"), # add conditions here (i.e., treatment, etc.)
                          plot_qcfeats = c("nCount_RNA","nFeature_RNA",
                                           "percent.ribo",
                                           "percent.mito", "percent.Hb"),
                          save_seurat=FALSE
                          ){
  
  obj_name = name
  if(!skipIntegration){
    DefaultAssay(seuratObject) <- "RNA"
    seuratObject <- runSampleCCA(seuratObject=seuratObject,
                                 combineLevel=combineLevel,
                                 numFeatures=2000,numDims=nDims,
                                 tissue=name,
                                 kParam=25,
                                 resUse=resUse,
                                 vars.to.regress = vars.to.regress)
    seuratObject@reductions$IntegrateSampleUMAP <- seuratObject@reductions$umap
    seuratObject@reductions$IntegrateSampleTSNE <- seuratObject@reductions$tsne
    obj_name = paste0("Integrated_",obj_name)
  }
  DefaultAssay(seuratObject) <- "RNA"
  seuratObject <- NormalizeData(seuratObject)
  seuratObject <- FindVariableFeatures(seuratObject)
  if(!is.null(vars.to.regress)){
    seuratObject <- ScaleData(seuratObject,vars.to.regress=vars.to.regress)
  } else {
    seuratObject <- ScaleData(seuratObject)
  }
  
  seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject))
  seuratObject <- FindNeighbors(seuratObject, dims = 1:nDims, reduction="pca",k.param = 25)
  seuratObject <- FindClusters(object = seuratObject, resolution = resUse, verbose = T, reduction = "pca")
  seuratObject <- RunUMAP(seuratObject, dims = 1:nDims, reduction = "pca")
  seuratObject <- RunTSNE(object = seuratObject, reduction = "pca", dims = 1:nDims)
  
  if(save_seurat){
    obj_output_name =  paste0("RNAClust_",obj_name,".rds")
    if(file.exists(obj_output_name)){
      iter = 1
      obj_output_name <- paste0(obj_output_name, iter)
      while(file.exists(obj_output_name)){
        obj_output_name <- paste0(substr(obj_output_name, 
                                         start = 1, 
                                         stop = nchar(obj_output_name)-1),
                                  iter)
        
        iter = iter + 1
      }
    }
    saveRDS(seuratObject, obj_output_name)
  }
  
  ifelse(!dir.exists(file.path(paste0("./Plots/rnaUMAP/resolution/",name,"/"))), dir.create(file.path(paste0("./Plots/rnaUMAP/resolution/",name,"/")),recursive = T), FALSE)
  ifelse(!dir.exists(file.path(paste0("./Plots/rnaTSNE/resolution/",name,"/"))), dir.create(file.path(paste0("./Plots/rnaTSNE/resolution/",name,"/")),recursive = T), FALSE)
  
  if(!(combineLevel %in% plot_meta)) plot_meta <- c(plot_meta, combineLevel)
  for(meta in plot_meta){
    pdf(file=paste0("./Plots/rnaTSNE/",name,"_",meta,".pdf"))
    print(DimPlot(object = seuratObject, reduction = "tsne", group.by = meta))
    dev.off()
  }
  for(meta in plot_meta){
    pdf(file=paste0("./Plots/rnaUMAP/",name,"_",meta,".pdf"))
    print(DimPlot(object = seuratObject, reduction = "umap", group.by = meta))
    dev.off()
  }
  
  # plot QC feature plots
  for(feat in plot_qcfeats){
    if(!(feat %in% colnames(seuratObject@meta.data))){
      plot_qcfeats <- plot_qcfeats[plot_qcfeats!=feat]
    }
  }
  wid = 7*length(plot_qcfeats)
  
  pdf(file=paste0("./Plots/rnaUMAP/QC_FeaturePlots_",name,"_umap.pdf"), width = wid)
  print(FeaturePlot(object = seuratObject, 
                    features = plot_qcfeats, 
                    reduction = "umap", 
                    ncol = length(plot_qcfeats)))
  dev.off()
  
  pdf(file=paste0("./Plots/rnaTSNE/QC_FeaturePlots_",name,"_tsne.pdf"), width = wid)
  print(FeaturePlot(object = seuratObject, 
                    features = plot_qcfeats, 
                    reduction = "tsne"))
  dev.off()

  resolutions <- as.character(resUse)
  for(reso in resolutions){
    pdf(file=paste0("./Plots/rnaTSNE/resolution/",name,"/PCA_tSNEres",reso,".pdf"))
    seuratObject <- SetIdent(seuratObject, value = paste0("RNA_snn_res.",reso))
    print(DimPlot(seuratObject, label = T, reduction = "tsne"))
    dev.off()
    pdf(file=paste0("./Plots/rnaUMAP/resolution/",name,"/PCA_UMAPres",reso,".pdf"))
    seuratObject <- SetIdent(seuratObject, value = paste0("RNA_snn_res.",reso))
    print(DimPlot(seuratObject, label = T, reduction = "umap"))
    dev.off()
  }
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
  if(!("avg_logFC" %in% colnames(DEG_df))){
    fc_type <- "avg_log2FC"
  } else{
    fc_type <- "avg_logFC"
  }
  
  DEG_df = DEG_df[abs(DEG_df[,fc_type])>logFC_threshold & DEG_df$p_val_adj<p_threshold,]
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
      subset <- DEG_df[DEG_df$Cell_type.Comparison==paste(ct,comp,sep = "."),]
      for(g in genes_perturbed){
        if(g %in% subset$GENE){
          vect_FC = c(vect_FC,subset[,fc_type][subset$GENE==g])
        } else {
          vect_FC = c(vect_FC, "None")
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
    new_shared = combined_DEGs[combined_DEGs$Tissue==unique(combined_DEGs$Tissue)[1] &
                                 combined_DEGs$Cell_type==ct,]$GENE_Direction[match(new_shared,
                                                                                    combined_DEGs[combined_DEGs$Tissue==unique(combined_DEGs$Tissue)[1] &
                                                                                                    combined_DEGs$Cell_type==ct,]$GENE)]
    temp[,paste0("Shared")] = paste0("(",length(new_shared),") ",
                                     concatenate(new_shared, mysep = ", "))
    shared_discordant = combined_DEGs[combined_DEGs$Tissue==unique(combined_DEGs$Tissue)[1] &
                                        combined_DEGs$Cell_type==ct,]$GENE_Direction[match(shared_discordant,
                                                                                           combined_DEGs[combined_DEGs$Tissue==unique(combined_DEGs$Tissue)[1] &
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
element_grob.element_custom_x <- function(element, label, x, y, ...)  {
  tg <- textGrob(label, x=x, gp=gpar(col=element$colour), 
                 rot = 45, just = "right", y=unit(-.25,"line"))
  padding <- unit(1,"line")
  rg <- rectGrob(x=x,width=unit(1,"line"), height=grobHeight(tg) + padding , 
                 gp=gpar(fill = element$fill, col=NA, alpha=.5, rot=45),
                 hjust = .5,vjust=1)
  gTree(children=gList(rg, tg), height=grobHeight(tg) + padding, cl="custom_axis")
}

# y axis labels
element_grob.element_custom_y <- function(element, label, x, y, ...)  {
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
gene_dot_plot <- function(
  DEG_df,
  vars = c("GENE","Cell_type"), 
  var_list = NULL, 
  num_genes_show = 30, # if var_list "GENE" not specified
  colors = c("blue","gray","red"),
  FDR_threshold = 0.05, # only used if no genes specified
  logFC_threshold = 0.25, # only used if no genes specified
  #logFC_min_max = c(-1,1)
  facet_by = NULL,
  cluster = NULL,
  cluster_by = NULL
){
  if(!("avg_log2FC" %in% colnames(DEG_df))){
    logFCtype = "avg_logFC"
  } else {
    logFCtype = "avg_log2FC"
  }
  
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
      if(length(unique(selected[,vars[v]]))!=length(var_list[[v]])){
        selected[,vars[v]] = factor(selected[,vars[v]], levels = unique(selected[,vars[v]]))
      }
      else{
        selected[,names(var_list)[v]] = factor(selected[,names(var_list)[v]], levels = var_list[[v]])
      }
    }
  }
  else if(length(unique(DEG_df$GENE))==1 & is.null(var_list)){ # single gene plot with multiple comparisons
    selected = DEG_df
  }
  else{
    # select the top 30 most "regulated" genes, show all cell types
    # subset DEG df to only significant DEGs
    DEG_df <- DEG_df[DEG_df$p_val_adj<FDR_threshold & abs(DEG_df$logFCtype)>logFC_threshold,]
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
  
  if(!is.null(cluster)){
    for(meta in cluster){
      if(!is.null(cluster_by)){
        meta_levels <- getHeirOrder(df = selected[,c(meta, cluster_by, logFCtype)])
      }
      else{
        meta_levels <- getHeirOrder(df = selected[,c(meta, setdiff(vars, meta), logFCtype)])
      }
      selected[,meta] = factor(selected[,meta], levels = meta_levels)
    }
  }
  
  if(is.null(facet_by)){
    D <- ggplot(selected, aes_string(x=vars[1], y=vars[2], size="neglog10fdr", color=logFCtype)) + 
      geom_point(aes(shape=`p(Adj) < 0.05`)) + 
      #geom_point(data=selected[selected$Stringent=="Yes",],pch=21,colour="black", stroke=1) +
      theme_bw() +
      theme(panel.border = element_blank(), #panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 45, size=14, hjust = 1,vjust = 1),
            #axis.text.x = element_text(angle = 90, size=14, hjust = 0,vjust = 1),
            axis.text.y = element_text(size = 14),
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
  }
  else{
    if(length(unique(selected[,facet_by]))==1){
      strip_x_angle = 0
    } else {
      strip_x_angle = 90
    }
    D <- ggplot(data = selected, aes_string(x=vars[1], y=vars[2], size="neglog10fdr", color=logFCtype)) +
      geom_point(aes(shape=`p(Adj) < 0.05`)) + 
      theme_bw() +
      facet_grid(reformulate(facet_by,"Tissue"), scales = "free_y", space = "free") + # hard coding Comparison! Switch in between this and Tissue...
      theme(panel.border = element_blank(), #panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 45, margin = margin(t=2), size=14, hjust = 1,vjust = 1),
            #axis.text.x = element_text(angle = 90, size=14, hjust = 0,vjust = 1),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.title = element_text(size=15),
            legend.text = element_text(size=10),
            legend.position = "bottom") +
      theme(strip.background = element_blank(), #remove background for facet labels
            panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
            #panel.spacing = unit(0, "lines"),
            #strip.placement = "outside",
            #strip.text = element_blank())+
            strip.placement = "outside",
            strip.text.x = element_text(angle = as.numeric(strip_x_angle), hjust = 0, size=18),
            strip.text.y = element_text(angle = 270, hjust = .5, size=16))+
      labs(size=expression(-log[10]~p["Adj"]), color=expression(logFC)) +
      scale_color_gradient2(low = colors[1],  mid=colors[2], high = colors[3]) +
      #breaks = c(-.5, 0, .5), labels = c(-.5, 0, .5)) +
      scale_size(range = c(3, 12), breaks = c(1.3,25,50,75,100)) +
      scale_shape_manual(values=shapes) +
      guides(shape = guide_legend(override.aes = list(size = 5))) 
  }
  return(D)
}

gene_dot_plot_log2FC <- function(
  DEG_df,
  vars = c("GENE","Cell_type"), 
  var_list = NULL, 
  num_genes_show = 30, # if var_list "GENE" not specified
  colors = c("blue","gray","red"),
  FDR_threshold = 0.05, # only used if no genes specified
  logFC_threshold = 0.25, # only used if no genes specified
  #logFC_min_max = c(-1,1)
  facet_by = NULL,
  cluster = NULL,
  cluster_by = NULL
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
      if(length(unique(selected[,vars[v]]))!=length(var_list[[v]])){
        selected[,vars[v]] = factor(selected[,vars[v]], levels = unique(selected[,vars[v]]))
      }
      else{
        selected[,names(var_list)[v]] = factor(selected[,names(var_list)[v]], levels = var_list[[v]])
      }
    }
  }
  else if(length(unique(DEG_df$GENE))==1 & is.null(var_list)){ # single gene plot with multiple comparisons
    selected = DEG_df
  }
  else{
    # select the top 30 most "regulated" genes, show all cell types
    # subset DEG df to only significant DEGs
    DEG_df <- DEG_df[DEG_df$p_val_adj<FDR_threshold & abs(DEG_df$avg_log2FC)>logFC_threshold,]
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
  
  if(!is.null(cluster)){
    for(meta in cluster){
      if(!is.null(cluster_by)){
        meta_levels <- getHeirOrder(df = selected[,c(meta, cluster_by, "avg_log2FC")])
      }
      else{
        meta_levels <- getHeirOrder(df = selected[,c(meta, setdiff(vars, meta), "avg_log2FC")])
      }
      selected[,meta] = factor(selected[,meta], levels = meta_levels)
    }
  }
  
  if(is.null(facet_by)){
    D <- ggplot(selected, aes_string(x=vars[1], y=vars[2], size="neglog10fdr", color="avg_log2FC")) + 
      geom_point(aes(shape=`p(Adj) < 0.05`)) + 
      #geom_point(data=selected[selected$Stringent=="Yes",],pch=21,colour="black", stroke=1) +
      theme_bw() +
      theme(panel.border = element_blank(), #panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 45, size=14, hjust = 1,vjust = 1),
            #axis.text.x = element_text(angle = 90, size=14, hjust = 0,vjust = 1),
            axis.text.y = element_text(size = 14),
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
  }
  else{
    D <- ggplot(data = selected, aes_string(x=vars[1], y=vars[2], size="neglog10fdr", color="avg_log2FC")) +
      geom_point(aes(shape=`p(Adj) < 0.05`)) + 
      theme_bw() +
      facet_grid(reformulate(facet_by,"Tissue"), scales = "free_y", space = "free") + # hard coding Comparison! Switch in between this and Tissue...
      theme(panel.border = element_blank(), #panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 45, margin = margin(t=5), size=14, hjust = 1,vjust = 1),
            #axis.text.x = element_text(angle = 90, size=14, hjust = 0,vjust = 1),
            axis.text.y = element_text(size = 14),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.title = element_text(size=15),
            legend.text = element_text(size=10),
            legend.position = "bottom") +
      theme(strip.background = element_blank(), #remove background for facet labels
            panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
            #panel.spacing = unit(0, "lines"),
            #strip.placement = "outside",
            #strip.text = element_blank())+
            strip.placement = "outside",
            strip.text.x = element_text(angle = 90, hjust = 0, size=18),
            strip.text.y = element_text(angle = 270, hjust = .5, size=16))+
      labs(size=expression(-log[10]~p["Adj"]), color=expression(logFC)) +
      scale_color_gradient2(low = colors[1],  mid=colors[2], high = colors[3]) +
      #breaks = c(-.5, 0, .5), labels = c(-.5, 0, .5)) +
      scale_size(range = c(3, 12), breaks = c(1.3,25,50,75,100)) +
      scale_shape_manual(values=shapes) +
      guides(shape = guide_legend(override.aes = list(size = 5))) 
  }
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
  return(mat)
}

getHeirOrder <- function(df){
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
                         padjusted = 0.001, 
                         indicate_stringency = FALSE,
                         only_show_stringent = FALSE){
  if(!is.null(cell_type)){
    DEG_df = DEG_df[DEG_df$Cell_type==cell_type,]
  }
  
  if(!("avg_logFC" %in% colnames(DEG_df))){
    colnames(DEG_df)[colnames(DEG_df)=="avg_log2FC"] <- "avg_logFC"
    label=expression(log[2]~FDR)
  } else{
    label="logFC"
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
      xlab(label) + 
      theme_bw() +
      theme(legend.position = "right", panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.title = element_text(size=15), title = element_text(size=15)) +
      scale_color_manual(values = pal, limits = names(pal))
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
    
    v <- ggplot(DEG_df) + geom_point(aes(x=avg_logFC, y=log10_pval_adj, fill=direction), cex=2, pch = 21) +
      geom_hline(yintercept=-log10(padjusted), linetype="dashed", color = "gray") +
      geom_vline(xintercept = -lfc, color="gray",linetype = "dashed") +
      geom_vline(xintercept = lfc, color="gray", linetype = "dashed") +
      ylab(expression(-log[10]~FDR)) +
      xlab(label) + 
      theme_bw() +
      theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            axis.title = element_text(size=15), title = element_text(size=15)) +
      scale_fill_manual(values = pal,  limits = names(pal))
      #scale_color_manual(values =pal, limits = names(pal))
  }
  if(indicate_stringency){
    if(only_show_stringent){
      v <- v + geom_text_repel(data = DEG_df[DEG_df$Stringent==T & DEG_df$Sign==T,], aes(x = avg_logFC, y = log10_pval_adj,
                                                                   label = GENE),
                               point.padding = 0.5, segment.color = 'gray14',segment.size = 0.1,size = lbl_size, max.overlaps = Inf) 
    }
    else{
      v <- v + geom_text_repel(data = DEG_df[DEG_df$Sign==T,], aes(x = avg_logFC, y = log10_pval_adj,
                                                                   label = GENE, 
                                                                   color = ifelse(DEG_df[DEG_df$Sign==T,]$Stringent==T,"black","gray54")),
                               point.padding = 0.5, segment.color = 'gray14',segment.size = 0.1,size = lbl_size, max.overlaps = Inf) +
        scale_colour_manual(values = c("black","gray54"), limits=c("black","gray54"))
    }
  }
  else{
    v <- v + geom_text_repel(aes(x = avg_logFC, y = log10_pval_adj, label = ifelse(Sign == T, GENE,"")),
                             point.padding = 0.5, segment.color = 'gray14',segment.size = 0.1,size = lbl_size,
                             max.overlaps = Inf)
  }
  return(v)
}

# # slow! should figure out how to speed up
# addStringency <- function(DEG_df, stringentDEG_df){
#   # add Stringent column with "Yes"
#   stringent = c()
#   for(i in 1:nrow(DEG_df)){
#     if(sum(rowSums(cbind(DEG_df$Cell_type[i]==stringentDEG_df$Cell_type,
#                          DEG_df$GENE[i]==stringentDEG_df$GENE,
#                          DEG_df$Comparison[i]==stringentDEG_df$Comparison,
#                          DEG_df$Tissue[i]==stringentDEG_df$Tissue))==4)>0){
#       stringent[i] = TRUE
#     }
#     else{
#       stringent[i] = FALSE
#     }
#   }
#   DEG_df$Stringent = stringent
#   return(DEG_df)
# }

addSampleWiseConsistency <- function(DEG_df, stringentDEG_df,match_columns=c("Cell_type")){
  # add Stringent column with "Yes"
  match_columns <- c(match_columns, "GENE")
  match <- c()
  for(col in match_columns){
    match <- paste(match, DEG_df[,col], sep = "_")
  }
  DEG_df$match <- match
  match <- c()
  for(col in match_columns){
    match <- paste(match, stringentDEG_df[,col], sep = "_")
  }
  stringentDEG_df$match <- match
  stringent = c()
  DEG_df$SampleWiseConsistent = ifelse(DEG_df$match %in% stringentDEG_df$match,"Yes","No")
  DEG_df$match <- NULL
  return(DEG_df)
}

makeMeanNormExprDf <- function(seuratObject, meta_data = c("Cell_type","condition"), genes, nonZero=FALSE){
  # make new meta data containing all labels
  if(length(meta_data)>1){
    for(iter in 1:length(meta_data)){
      if(iter==1) seuratObject$Combined_label = seuratObject@meta.data[,meta_data[iter]]
      else seuratObject$Combined_label = paste(seuratObject$Combined_label, seuratObject@meta.data[,meta_data[iter]], sep = ".")
    }
  }
  else{
    seuratObject$Combined_label = seuratObject@meta.data[,meta_data]
  }
  seuratObject$Combined_label <- as.character(seuratObject$Combined_label)
  
  all_norm_data = data.frame("GENE" = genes, stringsAsFactors = FALSE)
  for(meta in unique(seuratObject$Combined_label)){
    if(nonZero){
      all_norm_data[,paste0(meta,"_prct_expr")] = apply(X = seuratObject@assays$RNA@data[genes,colnames(seuratObject)[seuratObject$Combined_label==meta]], 
                                                        MARGIN = 1, 
                                                        FUN = function(x){
                                                          return(sum(x > 0) / length(x = x))
                                                        })
      all_norm_data[,paste0(meta,"_expression_level")] <- sapply(all_norm_data[,paste0(meta,"_prct_expr")], function(x){
        if(x<0.05) return("low")
        else if(x>=0.05 & x<0.1) return("med")
        else return("high")
      })
      all_norm_data[,paste0(meta,"_mean_in_expr_cells")] = apply(X = seuratObject@assays$RNA@data[genes,colnames(seuratObject)[seuratObject$Combined_label==meta]], 
                                                   MARGIN = 1, 
                                                   FUN = function(x){
                                                    return(mean(x[x>0]))
                                                   })
      all_norm_data[is.na(all_norm_data)] <- 0
      all_norm_data[,paste0(meta,"_prct_expr")] <- round(all_norm_data[,paste0(meta,"_prct_expr")], digits = 4)
      all_norm_data[,paste0(meta,"_mean_in_expr_cells")] <- round(all_norm_data[,paste0(meta,"_mean_in_expr_cells")], digits = 4)
      
      # all_norm_data[,paste0(meta,"_expr.score")] <- all_norm_data[,paste0(meta,"_mean")]*all_norm_data[,paste0(meta,"_pct.expr")]
      # quants <- quantile(all_norm_data[,paste0(meta,"_expr.score")], probs=c(.33,.66))
      # pdf(paste0(meta,"_expr.score.pdf"), width = 12)
      # hist(all_norm_data[,paste0(meta,"_expr.score")], breaks = 200)
      # abline(v=quants[1], col="blue")
      # abline(v=quants[2], col="red")
      # dev.off()
      # all_norm_data[,paste0(meta,"_discr.expr")] <- sapply(all_norm_data[,paste0(meta,"_expr.score")], function(x){
      #   if(x<=quants[1]) return("low")
      #   else if(x>quants[1] & x<=quants[2]) return("med")
      #   else return("high")
      # })
      # all_norm_data[,paste0(meta,"_mean")] <- round(all_norm_data[,paste0(meta,"_mean")], digits = 4)
      # all_norm_data[,paste0(meta,"_pct.expr")] <- round(all_norm_data[,paste0(meta,"_pct.expr")], digits = 4)
      # all_norm_data[,paste0(meta,"_expr.score")] <- round(all_norm_data[,paste0(meta,"_expr.score")], digits = 4)
      
      
    }
    else{
      all_norm_data[,meta] = round(rowMeans(as.matrix(seuratObject@assays$RNA@data[genes,
                                                                                   colnames(seuratObject)[seuratObject$Combined_label==meta]])), digits = 3)
      # removed any NaNs - this might be because a certain celltype is not found in one of the conditions, for ex.
      remove = c()
      for(col in 1:length(colnames(all_norm_data))){
        if(sum(is.nan(all_norm_data[,col]))>0){
          remove = c(remove, col)
        }
      }
      if(length(remove)>0) all_norm_data = all_norm_data[,-remove]
    }
  }
  # order in alphabetical order
  #all_norm_data <- all_norm_data[,c(1, (order(colnames(all_norm_data)[2:length(colnames(all_norm_data))])+1))]
  
  return(all_norm_data)
}

runDEGs <- function(seuratObject, 
                    comparisons,
                    logfc_threshold=0.1,
                    Cell_type_column="Cell_type",
                    Condition_column="condition",
                    skip="", 
                    save_progress=FALSE,
                    output_name="DEGs.txt"){
  
  logFC = data.frame(stringsAsFactors = FALSE)
  DefaultAssay(seuratObject) <- "RNA"
  seuratObject$Cell_type_condition <- paste(seuratObject@meta.data[,Cell_type_column],
                                            seuratObject@meta.data[,Condition_column], 
                                            sep = "_")
  Idents(seuratObject) <- "Cell_type_condition"
  if(save_progress){
    ifelse(!dir.exists(file.path("DEGs/")), dir.create(file.path("DEGs/"),recursive = T), FALSE)
  }
  for(i in 1:length(comparisons)){
    cat("Comparison: ", concatenate(gsub("_","",comparisons[[i]]), mysep = ".v."), "\n")
    for(celltype in setdiff(unique(seuratObject@meta.data[,Cell_type_column]), skip)){
      cat("Cell Type: ",celltype, "\n")
      DEGs <- FindMarkers(seuratObject, ident.1 = paste(celltype,comparisons[[i]][1], sep = "_"), 
                          ident.2 = paste(celltype,comparisons[[i]][2], sep="_"), 
                          logfc.threshold = logfc_threshold)
      DEGs$Cell_type = celltype
      DEGs$GENE <- rownames(DEGs)
      DEGs$Comparison <- concatenate(gsub("_","",comparisons[[i]]), mysep = ".v.")
      logFC  = rbind(logFC, DEGs)
      if(save_progress){
        write.table(logFC, paste0("DEGs/",output_name), row.names = FALSE, sep = "\t",quote = FALSE)
      }
    }
  }
  return(logFC)
}

makeGeneNumberTable <- function(pathway_df, comparisons, corresponding_names){
  table = data.frame(stringsAsFactors = FALSE)
  for(path in unique(pathway_df$Pathway)){
    for(ct in unique(pathway_df$Cell_type)){
      if(sum(pathway_df$Cell_type[pathway_df$Pathway==path]==ct)==0) next
      temp = data.frame("Pathway"=path, "Cell_type"=ct, stringsAsFactors = FALSE)
      for(comp in 1:length(comparisons)){
        temp[,paste0(corresponding_names[comp],"_UP")] = as.numeric(pathway_df$nOverlap[pathway_df$Pathway==path &
                                                                                          pathway_df$Cell_type==ct &
                                                                                          pathway_df$Comparison_direction==paste0(comparisons[comp],"_UP")])
        temp[,paste0(corresponding_names[comp],"_UP_genes")] = pathway_df$Overlap[pathway_df$Pathway==path &
                                                                                    pathway_df$Cell_type==ct &
                                                                                    pathway_df$Comparison_direction==paste0(comparisons[comp],"_UP")]
        temp[,paste0(corresponding_names[comp],"_DOWN")] = as.numeric(pathway_df$nOverlap[pathway_df$Pathway==path &
                                                                                            pathway_df$Cell_type==ct &
                                                                                            pathway_df$Comparison_direction==paste0(comparisons[comp],"_DOWN")])
        temp[,paste0(corresponding_names[comp],"_DOWN_genes")] = pathway_df$Overlap[pathway_df$Pathway==path &
                                                                                      pathway_df$Cell_type==ct &
                                                                                      pathway_df$Comparison_direction==paste0(comparisons[comp],"_DOWN")]
        temp[is.na(temp)] <- 0
        temp[,paste0(corresponding_names[comp], "_Percent_UP")] = (temp[,paste0(corresponding_names[comp],"_UP")]/(temp[,paste0(corresponding_names[comp],"_UP")]+temp[,paste0(corresponding_names[comp],"_DOWN")]))
      }
      if(length(comparisons)!=1){
        end_positions = returnNum(length(comparisons))
      }
      beginning_positions = setdiff(1:ncol(temp), end_positions)
      temp = temp[,c(beginning_positions, end_positions)]
      table = rbind(table, temp)
    }
  }
  return(table)
}

unmelt <- function(df, id_column, value_column){
  first = TRUE
  for(i in unique(df[,id_column])){
    column = df[,value_column][df[,id_column]==i]
    if(first){
      unmelt_df = data.frame(i=column, stringsAsFactors = FALSE)
      names(unmelt_df) <- i
      first= FALSE
    }
    else{
      if(nrow(unmelt_df)>length(column)){
        column = append(column, rep(" ", nrow(unmelt_df)-length(column)))
      }
      else if(nrow(unmelt_df)==length(column)){
        column = column
      }
      else{
        for(t in 1:(length(column)-nrow(unmelt_df))){
          unmelt_df = rbind(unmelt_df, " ")
        }
      }
      unmelt_df[,i] = column
    }
  }
  return(unmelt_df)
}

addGWASstudy <- function(df, GWAS_info){
  GWAS = read.delim(GWAS_info, stringsAsFactors = FALSE)
  GWAS_genes = data.frame()
  
  #get genes
  genes = unique(GWAS$MAPPED_GENE)
  genes = c(genes, unique(GWAS$REPORTED.GENE.S.))
  genes_dash = genes[grep(" - ", genes)]
  genes = genes[!grepl(" - ", genes)]
  genes_comma = genes[grep(",", genes)]
  genes = genes[!grepl(",", genes)]
  genes_semicolon = genes[grepl(";", genes)]
  genes = genes[!grepl(";", genes)]
  genes_x = genes[grepl(" x ", genes)]
  genes = genes[!grepl(" x ", genes)]
  
  
  for(gene in genes_dash){
    genes = append(genes, unlist(strsplit(gene, split = " - ")))
  }
  for(gene in genes_comma){
    genes = append(genes, unlist(strsplit(gene, split = ", ")))
  }
  for(gene in genes_semicolon){
    genes = append(genes, unlist(strsplit(gene, split = "; ")))
  }
  for(gene in genes_x){
    genes = append(genes, unlist(strsplit(gene, split = " x ")))
  }
  
  genes = unique(genes)
  genes = genes[genes!=""]
  genes = genes[genes!="Intergenic"]
  genes = genes[genes!="intergenic"]
  
  studies = c()
  studies2 = c()
  for(gene in genes){ # might introduce 'errors' with grep() with similar gene names
    mapped_genes = concatenate(unique(GWAS$PUBMEDID[grep(gene, GWAS$MAPPED_GENE)]), 
                               mysep = ", ")
    reported_genes = concatenate(unique(GWAS$PUBMEDID[grep(gene, GWAS$REPORTED.GENE.S.)]), 
                                 mysep = ", ")
    studies = c(studies, ifelse(length(mapped_genes)==0, "", mapped_genes))
    studies2 = c(studies2, ifelse(length(reported_genes)==0,"", reported_genes))
  }
  
  for(iter in 1:length(studies)){
    if(studies[iter]!=studies2[iter]){
      if(studies[iter]=="") studies[iter] = studies2[iter]
      else studies[iter] = paste(studies[iter], studies2[iter], sep = ", ")
    }
    
    # get rid of duplicate studies
    studies_sep = unlist(strsplit(studies[iter], split = ", "))
    studies[iter] = concatenate(unique(studies_sep), mysep = ", ")
  }
  
  GWAS_genes = data.frame("GENE"=genes, "PUBMED_ID" = studies, stringsAsFactors = FALSE)
  GWAS_genes = GWAS_genes[order(nchar(GWAS_genes$PUBMED_ID), decreasing = TRUE),]
  
  # make GWAS hit info vector
  GWAS_hit_info = c()
  first = TRUE
  for(gene in df$HUMAN){
    if(sum(GWAS_genes$GENE==gene)>0){
      GWAS_hit_info = append(GWAS_hit_info, GWAS_genes$PUBMED_ID[GWAS_genes$GENE==gene])
    }
    else{
      GWAS_hit_info = append(GWAS_hit_info, "")
    }
  }
  return(GWAS_hit_info)
}


annotate_gene <- function(df, protein_names="10090.protein.info.v11.0.txt"){
  protein_names = read.delim(protein_names, stringsAsFactors = FALSE, quote = "", header = TRUE)
  info = c()
  for(m in df$GENE){
    if(length(protein_names$annotation[protein_names$preferred_name==m])==0){
      if(length(protein_names$annotation[protein_names$preferred_name==gsub("[[:digit:]]","",m)])==0){
        info = append(info, "Not annotated")
        #cat(m, "\n")
      }
      else{
        info = append(info, paste0("Possible annotation: ",
                                   concatenate(protein_names$annotation[protein_names$preferred_name==gsub("[[:digit:]]","",m)],
                                               mysep = ",")))
      }
    }
    else if(length(protein_names$annotation[protein_names$preferred_name==m])>1){
      info = append(info, 
                    concatenate(protein_names$annotation[protein_names$preferred_name==m],
                                mysep = ","))
      cat("This protein had more than one annotation:", m, "\n")
    }
    else{
      info = append(info, protein_names$annotation[protein_names$preferred_name==m])
    }
  }
  return(info)
}


makePathwayPretty <- function(vector){
  vector <- gsub("HALLMARK_","", vector)
  vector <- gsub("KEGG_","", vector)
  vector <- gsub("REACTOME_","", vector)
  vector <- gsub("BIOCARTA_","", vector)
  vector <- gsub("MHC_CLASS_II_ANTIGEN_PRESENTATION","MHC II antigens", vector)
  vector <- gsub("ADIPOCYTOKINE_SIGNALING_PATHWAY","Adipocytokine signaling", vector)
  vector <- gsub("XENOBIOTIC_METABOLISM","Xenobiotic metabolism", vector)
  vector <- gsub("INTERFERON_ALPHA_RESPONSE","Interferon alpha response", vector)
  vector <- gsub("MITOCHONDRIAL_PROTEIN_IMPORT","Mitochondrial protein import", vector)
  vector <- gsub("TRANS_GOLGI_NETWORK_VESICLE_BUDDING","Golgi vesicle budding", vector)
  vector <- gsub("AMYLOIDS","Amyloids", vector)
  vector <- gsub("COMPLEMENT","Complement", vector)
  vector <- gsub("EPITHELIAL_MESENCHYMAL_TRANSITION","Epithelial mesenchymal transition", vector)
  vector <- gsub("INFLAMMATORY_RESPONSE","Inflammatory Response", vector)
  vector <- gsub("IL6_JAK_STAT3_SIGNALING","IL6 JAK STAT3 Signaling", vector)
  vector <- gsub("P53HYPOXIA_PATHWAY","P53 Hypoxia pathway", vector)
  vector <- gsub("HYPOXIA","Hypoxia", vector)
  vector <- gsub("TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT","TCA cycle and ETC", vector)
  vector <- gsub("METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES","Amino acid metabolism", vector)
  vector <- gsub("GLYCOSPHINGOLIPID_BIOSYNTHESIS_GLOBO_SERIES","Glycosphingolipid biosynthesis", vector)
  vector <- gsub("GLUCONEOGENESIS","Gluconeogenesis", vector)
  vector <- gsub("NDKDYNAMIN_PATHWAY","NDK dyamin pathway", vector)
  vector <- gsub("LYSOSOME","Lysosome", vector)
  vector <- gsub("IRON_UPTAKE_AND_TRANSPORT","Iron uptake and transport", vector)
  vector <- gsub("APOPTOSIS","Apoptosis", vector)
  vector <- gsub("PROTEIN_FOLDING","Protein folding", vector)
  vector <- gsub("NITROGEN_METABOLISM","Nitrogen Metabolism", vector)
  vector <- gsub("REACTIVE_OXIGEN_SPECIES_PATHWAY","Reactive oxygen species", vector)
  vector <- gsub("ESTROGEN_RESPONSE_LATE","Late estrogen response", vector)
  vector <- gsub("INTERFERON_GAMMA_RESPONSE","Interferon gamma response", vector)
  vector <- gsub("ADIPOGENESIS","Adipogenesis", vector)
  vector <- gsub("TNFA_SIGNALING_VIA_NFKB","TNFa signaling\nvia NFkB", vector)
  vector <- gsub("CHOLESTEROL_HOMEOSTASIS","Cholesterol\nhomeostasis", vector)
  vector <- gsub("GLYCOLYSIS","Glycolysis", vector)
  vector <- gsub("AXON_GUIDANCE","Axon guidance", vector)
  vector <- gsub("ION_CHANNEL_TRANSPORT","Ion channel transport", vector)
  vector <- gsub("ESTROGEN_RESPONSE_EARLY","Early estrogen response", vector)
  vector <- gsub("PPAR_SIGNALING_PATHWAY","PPAR signaling", vector)
  vector <- gsub("STEROID_BIOSYNTHESIS","Steroid biosynthesis", vector)
  vector <- gsub("OXIDATIVE_PHOSPHORYLATION","Oxidative\nphosphorylation", vector)
  vector <- gsub("METABOLISM_OF_LIPIDS_AND_LIPOPROTEINS","Lipid and lipoprotein metabolism", vector)
  vector <- gsub("TRANSMEMBRANE_TRANSPORT_OF_SMALL_MOLECULES","Small molecule transport", vector)
  vector <- gsub("NEURONAL_SYSTEM","Neuronal system", vector)
  vector <- gsub("HEMOSTASIS","Hemostasis", vector)
  vector <- gsub("TIGHT_JUNCTION","Tight junction", vector)
  vector <- gsub("REGULATION_OF_ORNITHINE_DECARBOXYLASE_ODC","Regulation of ODC", vector)
  vector <- gsub("PROTEASOME","Proteasome", vector)
  vector <- gsub("MRNA_SPLICING","mRNA splicing", vector)
  vector <- gsub("SPLICEOSOME","Spliceosome", vector)
  vector <- gsub("SIGNALING_BY_WNT","Wnt signaling", vector)
  vector <- gsub("IMMUNE_SYSTEM","Immune System", vector)
  vector <- gsub("PGC1A_PATHWAY","PGC-1a pathway", vector)
  vector <- gsub("BRANCHED_CHAIN_AMINO_ACID_CATABOLISM","Branched chain amino acid catabolism", vector)
  vector <- gsub("ZINC_TRANSPORTERS","Zinc transporters", vector)
  vector <- gsub("ABC_TRANSPORTERS","ABC transporters", vector)
  vector <- gsub("GLYCOSPHINGOLIPID_METABOLISM","Glycosphingolipid metabolism", vector)
  vector <- gsub("FORMATION_OF_TUBULIN_FOLDING_INTERMEDIATES_BY_CCT_TRIC","Formation of tubulin folding intermediates", vector)
  vector <- gsub("MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION","Mitochondrial fatty acid beta oxidation", vector)
  vector <- gsub("CHOLESTEROL_BIOSYNTHESIS","Cholesterol biosynthesis", vector)
  vector <- gsub("UNFOLDED_PROTEIN_RESPONSE","Unfolded Protein Response", vector)
  vector <- gsub("PEROXISOME","Peroxisome", vector)
  vector <- gsub("MEMBRANE_TRAFFICKING","Membrane Trafficking", vector)
  vector <- gsub("TRANSCRIPTION","Transcription", vector)
  vector <- gsub("SYNTHESIS_SECRETION_AND_DEACYLATION_OF_GHRELIN","Ghrelin secretion", vector)
  vector <- gsub("DNA_REPAIR","DNA repair", vector)
  vector <- gsub("GAP_JUNCTION","Gap junction", vector)
  vector <- gsub("CELL_CYCLE","Cell cycle", vector)
  vector <- gsub("INSULIN_SIGNALING_PATHWAY","Insulin signaling", vector)
  vector <- gsub("SMOOTH_MUSCLE_CONTRACTION","Smooth Muscle Contraction", vector)
  vector <- gsub("CELL_CYCLE","Cell cycle", vector)
  vector <- gsub("CHEMOKINE_SIGNALING_PATHWAY","Chemokine signaling", vector)
  vector <- gsub("Complement_AND_COAGULATION_CASCADES","Complement and coagulation", vector)
  vector <- gsub("ANTIGEN_PROCESSING_AND_PRESENTATION","Antigen processing", vector)
  vector <- gsub("MET_PATHWAY","Met pathway", vector)
  vector <- gsub("PATHOGENIC_ESCHERICHIA_COLI_INFECTION","Infection", vector)
  vector <- gsub("APICAL_JUNCTION","Apical Junction", vector)
  vector <- gsub("REGULATION_OF_ACTIN_CYTOSKELETON","Actin cytoskeleton regulation", vector)
  vector <- gsub("P53_PATHWAY","p53 Pathway", vector)
  vector <- gsub("MAPK_SIGNALING_PATHWAY","MAPK Signaling", vector)
  vector <- gsub("IL2_STAT5_SIGNALING","IL2 STAT5 Signaling", vector)
  vector <- gsub("TRANSLATION","Translation", vector)
  vector <- gsub("PARKINSONS_DISEASE","Parkinson's disease", vector)
  vector <- gsub("ALZHEIMERS_DISEASE","Alzheimer's disease", vector)
  vector <- gsub("ANDROGEN_RESPONSE","Androgen response", vector)
  vector <- gsub("FATTY_ACID_TRIACYLGLYCEROL_AND_KETONE_BODY_METABOLISM","Fatty acid metabolism", vector)
  vector <- gsub("MTORC1_SIGNALING","mTORC1 signaling", vector)
  vector <- gsub("KRAS_SIGNALING_UP","K-Ras signaling", vector)
  vector <- gsub("INSULIN_PATHWAY","Insulin pathway", vector)
  
  return(vector)
}


getInfo <- function(gene){
  protein_names <- read.delim("../10090.protein.info.v11.0.txt", header = TRUE, 
                              stringsAsFactors = FALSE, quote = "") # from String database
  if(length(gene)==1){
    return(protein_names$annotation[protein_names$preferred_name==gene])
  }
  else{
    genes = list()
    for(g in gene){
      genes[[g]] = protein_names$annotation[protein_names$preferred_name==g]
    }
    return(genes)
  }
}

getInfoHuman <- function(gene){
  protein_names <- read.delim("./9606.protein.info.v11.0.txt", header = TRUE, 
                              stringsAsFactors = FALSE, quote = "") # from String database
  if(length(gene)==1){
    return(protein_names$annotation[protein_names$preferred_name==gene])
  }
  else{
    genes = list()
    for(g in gene){
      genes[[g]] = protein_names$annotation[protein_names$preferred_name==g]
    }
    return(genes)
  }
}

MakeSpecificDEGsLogFCsDf <- function(specificDEGs, 
                                     DEG_df,
                                     column="Cell_type", 
                                     output_name = "Specific_DEG_summary.xls",
                                     add_very_specific=TRUE,
                                     very_specific_logFC_include=TRUE,
                                     very_specific_logFC_cutoff=0.1,
                                     very_specific_fdr_cutoff=0.05){
  
  if(!("avg_logFC" %in% colnames(DEG_df))){
    fc_type <- "avg_log2FC"
  } else{
    fc_type <- "avg_logFC"
  }
  
  specificDEGsSummary = list()
  for(ct in unique(DEG_df$Cell_type)){
    specific_subset = specificDEGs[grep(ct, specificDEGs$Cell_types),]
    if(length(specific_subset$GENE)==0) next
    specific_subset = specific_subset[order(specific_subset$nCellTypes, decreasing = TRUE),]
    first_cts = unlist(strsplit(specific_subset$Cell_types[1], split = ", "))
    first_cts = gsub(".DOWN", "", first_cts)
    first_cts = gsub(".UP","",first_cts)
    specific_df = data.frame(stringsAsFactors = FALSE)
    for(gene in unique(specific_subset$GENE)){
      temp = data.frame("GENE"=gene, stringsAsFactors = FALSE)
      celltypes = c(first_cts, setdiff(unique(DEG_df[,column]), first_cts))
      for(ct_tiss in celltypes){
        if(length(DEG_df[,fc_type][DEG_df[,column]==ct_tiss & DEG_df$GENE==gene])==0){
          temp[,paste0(ct_tiss,"_lfc")] = 0
          temp[,paste0(ct_tiss,"_fdr")] = 1
        }
        else{
          temp[,paste0(ct_tiss,"_lfc")] = DEG_df[,fc_type][DEG_df[,column]==ct_tiss & DEG_df$GENE==gene]
          temp[,paste0(ct_tiss,"_fdr")] = DEG_df$p_val_adj[DEG_df[,column]==ct_tiss & DEG_df$GENE==gene]
        }
      }
      specific_df = rbind(specific_df, temp)
    }
    specific_df = specific_df[order(abs(specific_df[,2]), decreasing = TRUE),]
    specificDEGsSummary[[ct]] <- specific_df
    if(add_very_specific){
      very_specific_df <- specific_df
      other_lfcs <- colnames(specific_df)[grepl("_lfc",colnames(specific_df))]
      other_lfcs <- other_lfcs[other_lfcs!=paste0(ct, "_lfc")]
      other_fdrs <- colnames(specific_df)[grepl("_fdr",colnames(specific_df))]
      other_fdrs <- other_fdrs[other_fdrs!=paste0(ct, "_fdr")]
      very_specific_df <- very_specific_df[apply(very_specific_df[,other_fdrs], 1,function(x){
        return(sum(x<very_specific_fdr_cutoff)==0)
      }),]
      specificDEGsSummary[[paste0(ct,"_fdr_filtered")]] <- very_specific_df
      if(very_specific_logFC_include){
        very_specific_df <- very_specific_df[apply(very_specific_df[,other_lfcs], 1,function(x){
          return(sum(abs(x)>very_specific_logFC_cutoff)==0)
        }),]
        specificDEGsSummary[[paste0(ct,"_fdr_lfc_filtered")]] <- very_specific_df
      }
    }
  }
  names(specificDEGsSummary) <- substr(names(specificDEGsSummary), start = 1, stop = 31)
  WriteXLS(specificDEGsSummary, output_name, names(specificDEGsSummary))
  return(specificDEGsSummary)
}


tool.read <- function(file, vars=NULL) {
  if(is.null(file)) return(data.frame())
  if(file == "") return(data.frame())
  dat <- read.delim(file=file, header=TRUE,
                    na.strings=c("NA", "NULL", "null", ""),
                    colClasses="character", comment.char="",
                    stringsAsFactors=FALSE)
  if(is.null(vars) == FALSE) dat <- dat[,vars]
  dat <- na.omit(dat)
  return(dat)
}

# @param DEG_df - can be either 'MODULE' 'GENE' file or DEG output from Seurat with genes in 'GENE' column
# @param MODULE_column - if 'Cell_type' is not the desired "module" to run the analysis, set the desired module column here
# @param resources_path - path to resources directory with txt files of resources in 'module' 'gene' format
# outputs the .txt reports for each "module"
# outputs an example heatmap of the top 50 consistent pathways
# returns a dataframe detailing each pathway enrichment result
# a note about logFC threshold:
#   i've seen opposite directions of pathway logFC when using 0.25 and 0.1 logFC thresholds, expectedly. something to consider when choosing a threshold
#   i.e. whether to be inclusive or only focus on DEGs with over 0.40546 logFC (this is a 1.5 fold change) (which you may claim are "valid" changes)
#   i'd prefer to be inclusive when discussing things on a pathway level
# if log2FC, correct labeling accordingly!
makePathwayEnrichmentDf <- function(DEG_df, 
                                    resources_path, 
                                    output_Dir="./PathwayEnrichmentResults", 
                                    convertToHuman=TRUE, 
                                    addlogFC=FALSE, 
                                    MODULE_column=NULL, 
                                    FDR_threshold=0.05, 
                                    logFC_threshold=0.1, 
                                    min_max=NULL, 
                                    makePlot=TRUE,
                                    plot_wid=15,
                                    plot_hgt=20,
                                    permute=FALSE,
                                    nperm=10,
                                    genepool){
  
  if("avg_log2FC" %in% colnames(DEG_df)){
    colnames(DEG_df)[colnames(DEG_df)=="avg_log2FC"] <- "avg_logFC"
  }
   # trim DEG output if necessary 
  if(length(colnames(DEG_df))>2){
    if(addlogFC){
      DEG_df_logFC_info = DEG_df # this is to add the logFC for the entire pathway
    }
    DEG_df = DEG_df[abs(DEG_df$avg_logFC)>logFC_threshold &
                      DEG_df$p_val_adj<FDR_threshold,]
    if(!is.null(MODULE_column)){
      DEG_df <- DEG_df[,c(MODULE_column,"GENE")]
    }
    else{
      DEG_df <- DEG_df[,c("Cell_type","GENE")]
    }
    colnames(DEG_df) <- c("MODULE","GENE")
  }else if(length(colnames(DEG_df))==2){
    colnames(DEG_df) <- c("MODULE","GENE")
  }
  else{
    cat("At least two columns are required, for example the geneset name in the 'MODULE' column and genes
        in the 'GENE' column.\n")
  }
  
  if(convertToHuman & addlogFC){
    temp <- convertDfGeneColumnMouseHuman(df = DEG_df_logFC_info, toSpecies = "human", forPathway = TRUE)
    DEG_df_logFC_info$HUMAN <- temp$GENE
    rm(temp)
    reference <- DEG_df_logFC_info[!duplicated(DEG_df_logFC_info[,c("GENE","HUMAN")]),]
    DEG_df$GENE <- DEG_df_logFC_info$HUMAN[match(DEG_df$GENE, DEG_df_logFC_info$GENE)]
    DEG_df_logFC_info$GENE <- DEG_df_logFC_info$HUMAN
  }
  if(convertToHuman & !addlogFC){
    DEG_df <- convertDfGeneColumnMouseHuman(df = DEG_df, toSpecies = "human", forPathway = TRUE)
  }
  # if(convertToHuman){
  #   DEG_df <- convertDfGeneColumnMouseHuman(df = DEG_df, toSpecies = "human", forPathway = TRUE)
  # }
  
  # enrichment doesn't work if you have " " in the module name
  space_modules = unique(DEG_df$MODULE)[grep(" ", unique(DEG_df$MODULE))]
  DEG_df$MODULE <- gsub(" ","_", DEG_df$MODULE)
  underscore_modules = gsub(" ","_", space_modules)
  
  # run pathway enrichment ------------------------------------------------------------------------------------------
  list_of_pathway_databases = list.files(path = resources_path, pattern = "*.txt", full.names = TRUE)
  ifelse(!dir.exists(file.path(output_Dir)), dir.create(file.path(output_Dir),recursive = T), FALSE)
  result = list()
  all_result = list()
  annotations = c()
  if(permute){
    genepool <- unique(c(genepool,DEG_df$GENE)) # genes from current module
  }
  for(module in unique(DEG_df$MODULE)){
    cat(module, "\n")
    deg_list = DEG_df[DEG_df$MODULE==module,] # why is this creating NAs?
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
      for(k in 1:Module2_len){ # go through each pathway
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
        
        if(permute & length(Overlapped_genes)>0){
          # get all genes from all pathways and add on module genes
          # randomly label the same amount of real overlapped genes as intersected
          nmodulegenes <- length(deg_list$gene[which(match(deg_list$module, unique_module1)>0)])
          # sample from gene pool ngenes times, find overlap, and rerun hypergeometric test, repeat 10000 times
          permutedenrichment <- c()
          for(i in 1:nperm){
            indices = sample(x = 1:length(genepool),size = nmodulegenes)
            permutedgenes <- genepool[indices]
            permutedgeneoverlap <- intersect(permutedgenes,
                                             Module2$gene[which(match(Module2$module, Unique_module2[k])>0)])
            
            enrichment <- length(permutedgeneoverlap)/nmodulegenes*25000/length(Module2$gene[which(match(Module2$module, Unique_module2[k])>0)])
            permutedenrichment <- c(permutedenrichment,enrichment)
          }
          noverlap <- length(Overlapped_genes)
          ngenepathway <- length(Module2$gene[which(match(Module2$module, Unique_module2[k])>0)])
          real_enrichment <- noverlap/nmodulegenes*25000/ngenepathway
          perm_p <- sum(permutedenrichment>real_enrichment)/length(permutedenrichment)
          data_matrix_for_enrichment[List_initial,9] <- perm_p
        } else if(permute){
          data_matrix_for_enrichment[List_initial,9] <- 1
        } else{
          
        }
        
        List_initial=List_initial + 1
      }
      
      
      
      write.table(data_matrix_for_enrichment, file = paste0(output_Dir,"/",module,".dat"), quote = FALSE, sep = "\t",
                  row.names = FALSE, col.names = FALSE)
      # record_mat <- read.table(paste0("Csvs/pathway_enrichment/temp1.",module,".dat"))
      record_mat <- read.table(paste0(output_Dir,"/",module,".dat"), sep = "\t")
      record_length <- dim(record_mat)
      enrichment_score <- data.frame()
      
      for(i in 1:record_length[1]){
        enrichment_score[i,1] <- record_mat[i,1]
        enrichment_score[i,2] <- phyper(record_mat[i,3], record_mat[i,4], record_mat[i,5]-record_mat[i,4], record_mat[i,2], lower.tail = FALSE)
        enrichment_score[i,3] <- record_mat[i,3]/record_mat[i,2]*record_mat[i,5]/record_mat[i,4]
        enrichment_score[i,4] <- 0
        enrichment_score[i,5]<-record_mat[i,8]
        enrichment_score[i,6]<-record_mat[i,6]
        enrichment_score[i,7]<-record_mat[i,4]
        enrichment_score[i,8]<-record_mat[i,3]
        enrichment_score[i,9]<-record_mat[i,7]
        if(permute){
          enrichment_score[i,10]<-record_mat[i,9]
        }
      }
      enrichment_score[,4]<-p.adjust(enrichment_score[,2], 'bonferroni')
      if(permute){
        colnames(enrichment_score) <- c("Module","Pval","Enrichment","FDR","PathwaySource","Pathway","PathwayGeneCount","nOverlap","Overlap","perm_p")
      } else {
        colnames(enrichment_score) <- c("Module","Pval","Enrichment","FDR","PathwaySource","Pathway","PathwayGeneCount","nOverlap","Overlap")
      }
      
      enrichment_score$ModuleGeneCount = length(deg_list$gene)
      enrichment_score <- enrichment_score[order(enrichment_score$FDR),]
      x[[database_name]] <- data.frame(enrichment_score)
      
      if(z==1){
        all_pathways_df <- enrichment_score
      }else{
        all_pathways_df <- rbind(all_pathways_df, enrichment_score)
      }
    }
    
    total = data.frame(stringsAsFactors = FALSE)
    for(iter in 1:length(x)){
      total = rbind(total, x[[iter]])
    }
    total = total[order(total$nOverlap, decreasing = TRUE),]
    write.table(total, paste0(output_Dir,"/",module, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
    
  }
  
  pathway_files = list.files(output_Dir)[grep(".txt", list.files(output_Dir))]
  pathway_df = data.frame(stringsAsFactors = FALSE)
  for(file in pathway_files){
    temp = read.delim(paste0(output_Dir,"/",file), stringsAsFactors = FALSE, header = TRUE)
    pathway_df = rbind(pathway_df, temp)
  }
  pathway_df = pathway_df[order(pathway_df$nOverlap, decreasing = TRUE),]
  pathway_df$`negLog10FDR` = -log10(pathway_df$FDR)
  if(!makePlot & !addlogFC){
    return(pathway_df)
  }
  
  # add summed logFC ------------------------------------------------------------------------------------------
  # can be very slow
  if(addlogFC){
    if(length(space_modules)>0){
      for(i in 1:length(space_modules)){
        pathway_df$Module <- gsub(underscore_modules[i], space_modules[i], pathway_df$Module)
      }
    }
    list_of_pathway_databases = list.files(path = resources_path, pattern = "*.txt", full.names = TRUE)
    compiled_pathways_df = data.frame(stringsAsFactors = FALSE)
    for(file in list_of_pathway_databases){
      compiled_pathways_df <- rbind(compiled_pathways_df, read.delim(file, stringsAsFactors = FALSE))
    }
    colnames(compiled_pathways_df) <- c("MODULE","GENE")
    
    if(!is.null(MODULE_column)){
      pathway_df <- pathway_df_add_logFC(pathway_df = pathway_df, 
                                         logFC_df = DEG_df_logFC_info, 
                                         cluster_meta = MODULE_column, 
                                         calculateOnlyOverlapLogFC = FALSE, 
                                         onlySig = TRUE, 
                                         compiled_pathways = compiled_pathways_df)
    }
    else{
      pathway_df <- pathway_df_add_logFC(pathway_df = pathway_df, 
                                         logFC_df = DEG_df_logFC_info,
                                         cluster_meta = "Cell_type",
                                         calculateOnlyOverlapLogFC = FALSE, onlySig = TRUE, 
                                         compiled_pathways = compiled_pathways_df)
    }
  }
  
  # plot heat map
  if(makePlot){
    # make simple heatmap --------------------------------------------------------------------------------------------------
    # pick top consistent significant pathways
    top_pathway_df <- pathway_df[pathway_df$FDR<0.05 & pathway_df$nOverlap>3,]
    top_pathway_df$numAppearances <- sapply(top_pathway_df$Pathway, function(x){
      return(nrow(top_pathway_df[top_pathway_df$Pathway==x,]))
    })
    top_pathway_df = top_pathway_df[order(top_pathway_df$numAppearances, decreasing = TRUE),]
    
    if(nrow(top_pathway_df)>50){
      topPathways <- unique(top_pathway_df$Pathway)[1:50]
    } else if(nrow(top_pathway_df)==0) {
      cat("No significant pathways to plot.\n")
      return(pathway_df)
    } else{
      topPathways <- unique(top_pathway_df$Pathway)
    }
    
    top_pathway_df = pathway_df[pathway_df$Pathway %in% topPathways,]
    top_pathway_df$Module <- paste0(top_pathway_df$Module, " (",top_pathway_df$ModuleGeneCount,")")
    
    # change pathway names to look better
    top_pathway_df$Pathway <- makePathwayPretty(top_pathway_df$Pathway)
    top_pathway_df <- top_pathway_df[!duplicated(top_pathway_df[,c("Module","Pathway")]),] # might be losing a "better" pathway here
    
    if(addlogFC){
      # hclust pathways by logFC
      pathway_heat_mat <- convertToMatrix(top_pathway_df[,c("Pathway","Module","avg_logFC")])
      pathway.dendro <- as.dendrogram(hclust(d = dist(x = pathway_heat_mat)))
      # get order
      pathway_order <- order.dendrogram(pathway.dendro)
      pathway_levels <- rownames(pathway_heat_mat)[pathway_order]
      
      # hclust pathways by Module
      if(length(unique(top_pathway_df$Module))>=2){
        pathway_heat_mat <- convertToMatrix(top_pathway_df[,c("Module","Pathway","nOverlap")])
        pathway.dendro <- as.dendrogram(hclust(d = dist(x = pathway_heat_mat)))
        ct_order <- order.dendrogram(pathway.dendro)
        ct_levels <- rownames(pathway_heat_mat)[ct_order]
        top_pathway_df$Module <- factor(top_pathway_df$Module, levels = ct_levels)
      }
    } else {
      # hclust pathways by nOverlap
      pathway_heat_mat <- convertToMatrix(top_pathway_df[,c("Pathway","Module","nOverlap")]) # can change to Enrichment or FDR
      pathway.dendro <- as.dendrogram(hclust(d = dist(x = pathway_heat_mat)))
      # get order
      pathway_order <- order.dendrogram(pathway.dendro)
      pathway_levels <- rownames(pathway_heat_mat)[pathway_order]
      
      # hclust pathways by Module
      if(length(unique(top_pathway_df$Module))>=2){
        pathway_heat_mat <- convertToMatrix(top_pathway_df[,c("Module","Pathway","nOverlap")])
        pathway.dendro <- as.dendrogram(hclust(d = dist(x = pathway_heat_mat)))
        ct_order <- order.dendrogram(pathway.dendro)
        ct_levels <- rownames(pathway_heat_mat)[ct_order]
        top_pathway_df$Module <- factor(top_pathway_df$Module, levels = ct_levels)
      }
    }
    
    top_pathway_df$Pathway <- factor(top_pathway_df$Pathway, levels = pathway_levels)

    if(addlogFC){
      if(!is.null(min_max)){
        top_pathway_df$avg_logFC <- MinMax(top_pathway_df$avg_logFC, min = min_max[1], max = min_max[2]) # need Seurat
      }
      heat <- ggplot(top_pathway_df, aes(x=Module, y=Pathway, fill=avg_logFC)) +#, label=nOverlap, )) +
        #geom_tile(aes(colour="gray"), size=.7) +
        geom_tile(colour="grey", size=.7) +
        #geom_text(size=2.7, aes(colour=FDR<0.05)) +
        theme_bw() +
        #facet_grid(.~Tissue, scales = "free", space = "free") +
        #scale_fill_gradient2(low="blue",high="red", mid = "white", breaks = c(0), labels=c("")) + 
        scale_fill_gradient2(low="blue",high="red", mid = "white") + 
        theme(#strip.background = element_blank(), #remove background for facet labels
          panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
          #panel.spacing = unit(0, "lines"),
          #strip.placement = "outside",
          # strip.text = element_text(margin = margin(t=0,b=3,r=0,l=0),
          #                           hjust = .5, size=14),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 11.5),
          axis.text.y = element_text(size = 11),
          panel.grid = element_blank(), axis.line = element_blank(), 
          legend.title = element_text(margin=margin(t=0,b=3,r=0,l=0)))+
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        #scale_color_manual(values = c("black", "darkgray"))+
        labs(fill="logFC")
      pdf(paste0(output_Dir,"/heatmap.pdf"), height = plot_hgt, width = plot_wid)
      print(heat)
      dev.off()
    } else {
      if(!is.null(min_max)){
        top_pathway_df$nOverlap <- MinMax(top_pathway_df$nOverlap, min = min_max[1], max = min_max[2]) # need Seurat
      }
      heat <- ggplot(top_pathway_df, aes(x=Module, y=Pathway, fill=nOverlap)) +#, label=nOverlap, )) +
        #geom_tile(aes(colour="gray"), size=.7) +
        geom_tile(colour="grey", size=.7) +
        #geom_text(size=2.7, aes(colour=FDR<0.05)) +
        theme_bw() +
        #facet_grid(.~Tissue, scales = "free", space = "free") +
        #scale_fill_gradient2(low="blue",high="red", mid = "white", breaks = c(0), labels=c("")) + 
        scale_fill_gradient2(low="blue",high="red", mid = "white") + 
        theme(#strip.background = element_blank(), #remove background for facet labels
          panel.border = element_rect(colour = "black", fill = NA, size = 1), #add black border
          #panel.spacing = unit(0, "lines"),
          #strip.placement = "outside",
          # strip.text = element_text(margin = margin(t=0,b=3,r=0,l=0),
          #                           hjust = .5, size=14),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 11.5),
          axis.text.y = element_text(size = 11),
          panel.grid = element_blank(), axis.line = element_blank(), 
          legend.title = element_text(margin=margin(t=0,b=3,r=0,l=0)))+
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        #scale_color_manual(values = c("black", "darkgray"))+
        labs(fill="nOverlap")
      pdf(paste0(output_Dir,"/heatmap.pdf"), height = plot_hgt, width = plot_wid)
      print(heat)
      dev.off()
    }
  }
  
  write.table(pathway_df,paste0(output_Dir,"/pathway_df.txt"),
              row.names = FALSE, sep = "\t", quote = FALSE)
  
  return(pathway_df)
  
}

summarizeTopDegs <- function(DEG_df, 
                             MODULE = "MODULE_direction",
                             comparisons = NULL, # useful to set for multiple comparisons
                             celltypes = NULL, 
                             output_name = "topDEGs.xls",
                             meanNormExprDf = NULL,
                             conditions = NULL, 
                             convertToHuman=FALSE,
                             gene_info = NULL,
                             GWAS_info = NULL){
  
  dir = do.call("paste",c(unlist(strsplit(output_name, "/"))[1:length(unlist(strsplit(output_name, "/")))-1],
                          list("sep"="/")))
  ifelse(!dir.exists(file.path(dir)), dir.create(file.path(dir),recursive = T), FALSE)
  
  topDEGs = list()
  if(convertToHuman){
    mouse_genes = data.frame("GENE" = unique(DEG_df$GENE), stringsAsFactors = FALSE)
    human <- convertDfGeneColumnMouseHuman(df = mouse_genes, toSpecies = "human", forPathway = TRUE)
    mouse_human_reference <- data.frame("MOUSE"=mouse_genes$GENE, "HUMAN" = human$GENE, stringsAsFactors = FALSE)
    mouse_human_reference$GENE = mouse_human_reference$MOUSE
    rm(human)
    rm(mouse_genes)
  }
  # by comparison
  # useful for multiple comparisons. can still be set if just one comparison specified
  if(is.null(comparisons)){
    if(!is.null(DEG_df$Comparison)){
      comparisons = unique(as.character(DEG_df$Comparison))
      names(comparisons) = unique(as.character(DEG_df$Comparison))
    }
    else{
      cat("Assuming only one comparison...\n")
      comparisons = c("DEG_analysis")
      names(comparisons) = "DEG_analysis"
    }
  }
  
  # can still have one comparison and set the comparisons
  for(comp in names(comparisons)){
    cat("**Now summarizing ", comp, "\n")
    if(length(comparisons)>1) subset = DEG_df[DEG_df$Comparison==comp,]
    else subset = DEG_df
    celltype_group = c()
    logFCs = c()
    nCelltypes = c()
    for(DEG in unique(subset$GENE)){
      if(sum(subset$GENE==DEG)>0){
        celltype_group = append(celltype_group,
                                concatenate(gsub(paste0(comp,"_"),"",subset[,MODULE][subset$GENE==DEG]),
                                            mysep = ", "))
        logFCs = append(logFCs, concatenate(gsub(paste0(comp,"_"),"",subset$avg_logFC[subset$GENE==DEG]), 
                                            mysep = ", "))
        nCelltypes = append(nCelltypes, sum(subset$GENE==DEG))
      }
    }
    top_DEGs = data.frame("GENE" = unique(subset$GENE),
                          "Cell_types"=celltype_group,
                          "logFCs" = logFCs,
                          "nCellTypes"=nCelltypes, stringsAsFactors = FALSE)
    if(!is.null(meanNormExprDf)){
      for(cond in  conditions[[comp]]){
        all_avgs = c()
        for(DEG in unique(subset$GENE)){
          # need to get modules that do not have the comparison in the name
          # get modules that match those in the meanNormExprDf
          modules = gsub(paste0(".",comp),"",paste0(subset[,MODULE][subset$GENE==DEG], ".", cond)) 
          modules = gsub(".UP","",gsub(".DOWN","", modules))
          avgs = c()
          for(mod in modules){
            avgs = c(avgs, meanNormExprDf[meanNormExprDf$GENE==DEG,mod])
          }
          if(length(avgs)==0) cat(DEG,"\n")
          all_avgs = c(all_avgs, concatenate(avgs, mysep = ", "))
        }
        top_DEGs[,paste0(cond, "_means")] = all_avgs
      }
    }
    top_DEGs = top_DEGs[order(top_DEGs$nCellTypes, decreasing = TRUE),]
    
    if(convertToHuman){
      top_DEGs$HUMAN <- mouse_human_reference$HUMAN[match(x = top_DEGs$GENE, table = mouse_human_reference$MOUSE)]
    }
    
    if(!is.null(GWAS_info)){
      # add GWAS study accession
      top_DEGs$`AD GWAS Pubmed IDs` = addGWASstudy(df = top_DEGs, GWAS_info = GWAS_info)
    }
    
    if(!is.null(gene_info)){
      # add info
      top_DEGs$INFO = annotate_gene(df = top_DEGs, protein_names = gene_info)
    }
    
    # since sheet is already named by comparison, redundant to have comparison in the Cell Types column
    top_DEGs$Cell_types = gsub(paste0(comp,"."),"", top_DEGs$Cell_types)
    
    topDEGs[[comparisons[comp]]] = top_DEGs
  }
  
  if(is.null(celltypes)){
    celltypes = unique(DEG_df$Cell_type)
  }
  
  if(length(comparisons)>1){
    # by cell type effect
    # rearrange data frame
    ordered = data.frame(stringsAsFactors = FALSE)
    for(comp in names(comparisons)){
      ordered = rbind(ordered, DEG_df[DEG_df$Comparison==comp,])
    }
    DEG_df <- ordered
  }
  
  # adding the periods before and after because of multiple cell type namings - Endothelium and Endothelium 2, for ex.
  # for(celltype in paste0(celltypes,".")){
  for(celltype in celltypes){
    cat("**Now summarizing ", celltype, "\n")
    # subset = DEG_df[grep(celltype, DEG_df[,MODULE], fixed = TRUE),]
    subset = DEG_df[DEG_df$Cell_type==celltype,]
    if(nrow(subset)==0) next
    celltype_group = c()
    nCelltypes = c()
    logFCs = c()
    for(DEG in unique(subset$GENE)){
      if(sum(subset$GENE==DEG)>0){
        celltype_group = append(celltype_group, 
                                concatenate(gsub(paste0(celltype,"_"),"",subset[,MODULE][subset$GENE==DEG]), 
                                            mysep = ", "))
        logFCs = append(logFCs, concatenate(gsub(paste0(celltype,"_"),"",subset$avg_logFC[subset$GENE==DEG]), 
                                            mysep = ", "))
        nCelltypes = append(nCelltypes, sum(subset$GENE==DEG))
      }
    }
    top_DEGs = data.frame("GENE" = unique(subset$GENE),
                          # called Comparisons because this was useful for multiple comparisons
                          # (seeing consistency of genes across comparisons)
                          "Comparisons"=celltype_group, 
                          "logFCs" = logFCs,
                          "nComparisons"=nCelltypes, stringsAsFactors = FALSE)
    if(!is.null(meanNormExprDf)){
      if(length(comparisons)>1){
        rownames(meanNormExprDf) <- meanNormExprDf$GENE
        toAdd <- meanNormExprDf[unique(subset$GENE), colnames(meanNormExprDf)[grep(celltype, colnames(meanNormExprDf))]]
        colnames(toAdd) <- gsub(celltype, "", colnames(toAdd)) 
        conds <- rev(unlist(rev(conditions)))
        conds <- unique(conds)
        toAdd <- toAdd[,conds]
        top_DEGs <- cbind(top_DEGs, toAdd)
      }
      else{
        for(cond in conditions[[comp]]){
          all_avgs = c()
          for(DEG in unique(subset$GENE)){
            modules = gsub(paste0(".",comp),"",paste0(subset[,MODULE][subset$GENE==DEG], ".", cond)) 
            modules = gsub(".UP","",gsub(".DOWN","", modules))
            avgs = c()
            for(mod in modules){
              avgs = c(avgs, meanNormExprDf[meanNormExprDf$GENE==DEG,mod])
            }
            all_avgs = c(all_avgs, concatenate(avgs, mysep = ", "))
          }
          top_DEGs[,paste0(cond, "_means")] = all_avgs
        }
      }
    }

    
    top_DEGs = top_DEGs[order(top_DEGs$nComparisons, decreasing = TRUE),]
    
    if(convertToHuman){
      top_DEGs$HUMAN <- mouse_human_reference$HUMAN[match(x = top_DEGs$GENE, table = mouse_human_reference$MOUSE)]
    }
    
    if(!is.null(GWAS_info)){
      # add GWAS study accession
      top_DEGs$`AD GWAS Pubmed IDs` = addGWASstudy(df = top_DEGs, GWAS_info = GWAS_info)
    }
    
    if(!is.null(gene_info)){
      # add info
      top_DEGs$INFO = annotate_gene(df = top_DEGs, protein_names = gene_info)
    }
    celltype = gsub(".", "",celltype, fixed = TRUE)
    
    if(length(comparisons)==1){
      top_DEGs$Comparisons = gsub(paste0(names(comparisons),"."), "", top_DEGs$Comparisons)
    }
    
    # since sheet is already named by celltype, redundant to have celltype in the Comparisons column
    top_DEGs$Comparisons = gsub(paste0(celltype, "."),"", top_DEGs$Comparisons)
    topDEGs[[celltype]] = top_DEGs
    
  }
  
  WriteXLS(topDEGs, output_name, gsub("\\.","",names(topDEGs)))
  return(topDEGs)
  
}

strippedDimPlot <- function(seuratObject, cols, reduction="tsne", group.by="Cell_type",
                            pdf_output_name = "DimPlot.pdf", height=7, width=7, pt.size=1.5){
  # new_order = names(cols)[sort(order(names(cols))[levels(seuratObject@meta.data[,group.by])])]
  # cols = cols[new_order]
  pdf(pdf_output_name, width = width, height = height)
  print(DimPlot(seuratObject, reduction = reduction, group.by = group.by,
                cols = cols, pt.size = pt.size) +
        theme(plot.title = element_text(size = 25, hjust = .5),
              axis.title = element_blank(), axis.text = element_blank(), 
              axis.ticks = element_blank(), axis.line = element_blank()) + NoLegend())
  dev.off()
  pdf(paste0(gsub(".pdf","",pdf_output_name),"_Legend.pdf"), width = width, height = height)
  print(DimPlot(seuratObject, reduction = reduction, group.by = group.by,
              cols = cols, pt.size = pt.size) +
        theme(plot.title = element_text(size = 25, hjust = .5),
              axis.title = element_blank(), axis.text = element_blank(), 
              axis.ticks = element_blank(), axis.line = element_blank()))
  dev.off()
}

# returns ggplot object, can add theme elements
# not directional (e.g., neg/pos)
# outputs pdf of plot
drawTilePlotWithText <- function(df,
                                 x_aes,
                                 y_aes,
                                 text_clm,
                                 text_size=4,
                                 fill_clm=NULL,
                                 fill_col="firebrick4",
                                 output_name="TilePlotWithText.pdf",
                                 save=TRUE,
                                 height_width=c(7,7)
){
  if(is.null(fill_clm)){
    fill_clm = text_clm
  }
  plt <- ggplot(df, aes_string(x=x_aes, y=y_aes, fill=fill_clm, label=text_clm)) +
    #geom_tile(aes(colour="gray"), size=.7) +
    geom_tile(colour="grey", size=.7) +
    geom_text(size=text_size, fontface="bold") +
    theme_bw() +
    scale_fill_gradient2(low="white",high=fill_col) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, size = 12),
          axis.text.y = element_text(size = 11),
          axis.line = element_blank(), 
          legend.title = element_text(margin=margin(t=0,b=3,r=0,l=0)),
          #legend.position = "none",
          axis.text.x.bottom = element_text(vjust = 0.5, hjust = 1))+
    scale_x_discrete(expand = c(0, 0), position = "bottom") +
    scale_y_discrete(expand = c(0, 0), position = "left") +
    scale_color_manual(values = c("black", "white"))
  #labs(fill="logFC")
  if(save){
    if(grepl("/",output_name)){
      # get dir
      dir <- do.call("paste",c(unlist(strsplit(output_name, split = "/"))[1:length(unlist(strsplit(output_name, split = "/")))-1],
                               list("sep"="/")))
      ifelse(!dir.exists(file.path(dir)), dir.create(file.path(dir),recursive = T), FALSE)
    }
    
    pdf(output_name, height=height_width[1], width = height_width[2])
    print(plt)
    dev.off()
  }
  
  return(plt)
}



theme_cleanDimPlot <- function(){ 
    theme(
      plot.title = element_text(size = 25), 
      axis.line = element_blank(), 
      axis.ticks = element_blank(), 
      axis.title = element_blank(),
      axis.text = element_blank(), 
      legend.position = "none"
    )
}

theme_cleanDimPlotLegend <- function(){ 
  theme(
    plot.title = element_text(size = 25), 
    axis.line = element_blank(), 
    axis.ticks = element_blank(), 
    axis.title = element_blank(),
    axis.text = element_blank()
  )
}

manyFeaturesPlot <- function(seuratObject, 
                             features, 
                             reduction = "tsne", 
                             ncol, 
                             nrow, 
                             h_w=NULL,
                             output_name="MultipleFeatures.png",
                             pt.size=1,
                             plot.title.size=15,
                             labels=NULL,
                             highlight_clusters = NULL,
                             cluster_meta = NULL,
                             png_resolution_factor=1,
                             display_plot=FALSE){
  require(ggpubr)
  feature_plots <- list()
  
  if(grepl("ntegrate",reduction)){ # to correct for Seurat bug
    seuratObject@reductions$tsne <- NULL
    seuratObject@reductions$umap <- NULL
  }
  
  for(iter in 1:length(features)){
    if(!is.null(labels)){
      name = paste(labels[iter], features[iter], sep = "_")
    } else {
      name = features[iter]
    }
    feature_plots[[name]] <- FeaturePlot(seuratObject, 
                                         reduction = reduction, 
                                         features = features[iter], 
                                         pt.size = pt.size) +
      theme_cleanDimPlot() +
      theme(plot.title = element_text(size = plot.title.size))
  }
  if(is.null(h_w)){
    if(ncol>nrow){
      wid = 1280
      hght = (1280*nrow)/ncol
      wid = wid*png_resolution_factor
      hght = hght*png_resolution_factor
    } else{
      hght = 1280
      wid = (1280*ncol)/nrow
    }
  } else {
    hght = h_w
  }
  
  if(!is.null(highlight_clusters)){
    Idents(seuratObject) <- cluster_meta
    for(c in highlight_clusters){
      feature_plots[[paste0("Cluster_",c)]] <- DimPlot(seuratObject, reduction = reduction, 
                                                       cells.highlight = colnames(seuratObject)[seuratObject@active.ident==c],
                                                       sizes.highlight = 0.3) + 
        ggtitle(paste0("Cluster_",c)) +
        theme_cleanDimPlot() +
        theme(plot.title = element_text(size = plot.title.size))
    }
    final_order <- list()
    highlight_plots <- names(feature_plots)[grep("Cluster",names(feature_plots))]
    start = 1
    end = ncol
    for(iter in 1:length(highlight_plots)){
      final_order[[highlight_plots[iter]]] <-  feature_plots[[highlight_plots[[iter]]]]
      for(add in start:end){
        final_order[[names(feature_plots)[add]]] <- feature_plots[[names(feature_plots)[add]]]
      }
      start = end + 1
      end = start + ncol - 1
    }
    ncol = ncol + 1
    feature_plots <- final_order
    rm(final_order)
  }
  
  if(!is.null(output_name)){
    png(output_name, width = wid, height = hght, units = "px")
    print(ggarrange(plotlist = feature_plots, ncol = ncol, nrow = nrow))
    dev.off()
  }
  
  if(display_plot){
    print(ggarrange(plotlist = feature_plots, ncol = ncol, nrow = nrow))
  }
}

ClusterAnalysis <- function(seuratObject, 
                            condition = c("condition","Sample"), 
                            sample_colors = c("5XFAD1"="#800026",
                                              "5XFAD2"="#BD0026",
                                              "5XFAD3"="#E31A1C",
                                              "5XFAD4"="#FC4E2A",
                                              "WT1"="#1A1A1A",
                                              "WT2"="#4D4D4D",
                                              "WT3"="#878787"),
                            show_clusters = c("integrated_snn_res.0.5",
                                              "RNA_snn_res.0.5"),
                            combineLevel = "orig.ident",
                            name,
                            saves_dir="Saves/",
                            plots_dir="Figures/",
                            markers_dir="Markers/",
                            process_seurat=FALSE,
                            nMarkers=10,
                            nClustersOneImage=5,
                            sample_proportion_plot_width=6,
                            runDEGsClusters=NULL,
                            runDEGsCondition="condition",
                            DEGComparisons=list("AD"=c("5XFAD","WT")),
                            save_seurat=FALSE,
                            png_resolution_factor=1){
  
  ifelse(!dir.exists(file.path(saves_dir)), dir.create(file.path(saves_dir),recursive = T), FALSE)
  ifelse(!dir.exists(file.path(plots_dir)), dir.create(file.path(plots_dir),recursive = T), FALSE)
  ifelse(!dir.exists(file.path(markers_dir)), dir.create(file.path(markers_dir),recursive = T), FALSE)
  plot_markers_dir <- paste0(plots_dir,"/Markers/")
  ifelse(!dir.exists(file.path(plot_markers_dir)), dir.create(file.path(plot_markers_dir),recursive = T), FALSE)
  
  if(sum(grepl("integrated", show_clusters))!=0){
    skipIntegration = FALSE
  } else {
    skipIntegration = TRUE
  }
  if(process_seurat){
    seuratObject <- processSeurat(seuratObject = seuratObject, 
                                  combineLevel = combineLevel, 
                                  name = name, 
                                  skipIntegration = skipIntegration)
    save_seurat = TRUE
  }
  for(clust in show_clusters){
    if(grepl("integrated",clust)){
      red = c("IntegrateSampleUMAP","IntegrateSampleTSNE")
    } else {
      red = c("umap","tsne")
    }
    for(r in red){
      plots <- list()
      plots[["clusters"]] <- DimPlot(seuratObject, reduction = r, group.by = clust, label = TRUE, label.size = 10) + theme_cleanDimPlot()
      for(c in condition){
        plots[[c]] <- DimPlot(seuratObject, reduction = r, group.by = c) + theme_cleanDimPlotLegend()
      }
      wid = 467*length(plots)
      png(paste0(plots_dir,name,"_",
                 gsub("_snn","",clust),"_",
                 gsub("IntegrateSample","",r),
                 ".png"), width = wid, height = 425)
      print(ggarrange(plotlist = plots, ncol = length(plots), nrow = 1))
      dev.off()
    }
    
    for(meta in condition){
      counts <- cellTypeCount(seuratObject = seuratObject, 
                              meta_data_1 = meta,
                              meta_data_2 = clust, includePerc = FALSE)
      if(grepl(".",meta,fixed = TRUE)){
        meta <- gsub(".","_",meta, fixed = TRUE)
        colnames(counts)[1] <- meta
      }
      pdf(paste0(plots_dir,name,"_",
                 gsub("_snn","",clust),"_",meta,
                 "_SampleProportions.pdf"), height = 4, width = sample_proportion_plot_width)
      print(normalizedDistributionProportionPlot(data = counts, 
                                           meta_data_1 = meta, 
                                           meta_data_2 = "Cluster",
                                           pal = sample_colors, addCounts = FALSE, 
                                           meta_data_2_levels = 1:length(unique(seuratObject@meta.data[,clust]))) + 
              theme(axis.text.x = element_text(angle = 0, size = 15)))
      dev.off()
    }
    
    Idents(seuratObject) <- clust
    markers <- FindAllMarkers(seuratObject, only.pos = TRUE)
    markers$res <- clust
    if(is.null(seuratObject@misc$markers)){
      seuratObject@misc$markers <- markers
    } else {
      seuratObject@misc$markers <- rbind(seuratObject@misc$markers, markers)
    }
    
  }
  # save before feature plots because of FeaturePlot reduction bug
  if(save_seurat) saveRDS(seuratObject, paste0(saves_dir, name, "_SeuratObject.rds"))
  write.table(seuratObject@misc$markers, paste0(markers_dir, name, "_markers.txt"), 
              row.names = FALSE, quote = FALSE, sep = "\t")
  
  # make marker plots
  # put integrated last
  integrated <- show_clusters[grep("integrated",show_clusters)]
  if(length(integrated)!=length(show_clusters)){
    nonintegrated <- setdiff(show_clusters, integrated)
    show_clusters <- c(nonintegrated, integrated)
  }
  for(clust in show_clusters){
    if(grepl("integrated",clust)){
      red = c("IntegrateSampleUMAP","IntegrateSampleTSNE")
      rna_umap_red <- seuratObject@reductions$umap
      rna_tsne_red <- seuratObject@reductions$tsne
      seuratObject@reductions$umap <- NULL
      seuratObject@reductions$tsne <- NULL
    } else {
      red = c("umap","tsne")
    }
    
    markers <- seuratObject@misc$markers[seuratObject@misc$markers$res==clust,]
    count = 0
    set = 1
    features <- c()
    clusters <- c()
    for(c in unique(markers$cluster)[order(unique(markers$cluster))]){
      count = count + 1
      clusters <- c(clusters, c)
      features <- c(features, markers$gene[markers$cluster==c][1:nMarkers])
      if(count==nClustersOneImage || c==unique(markers$cluster)[length(unique(markers$cluster))]){
        for(r in red){
          manyFeaturesPlot(seuratObject = seuratObject, features = features, 
                           ncol = nMarkers, nrow = nClustersOneImage, 
                           reduction = r, pt.size = 0.1, 
                           labels = 1:(nMarkers*nClustersOneImage),
                           highlight_clusters = clusters,
                           cluster_meta = clust,
                           output_name = paste0(plot_markers_dir, name, "_",
                                                gsub("_snn","",clust),"_",
                                                gsub("IntegrateSample","",r),
                                                "_Markers_Set",set,".png"),
                           png_resolution_factor=1)
        }
        set = set + 1
        features <- c()
        clusters <- c()
        count = 0
      }
    }
  }
  
  if(!is.null(runDEGsClusters)){
    # run DEGs
    seuratObject$Cell_type_condition <- paste(seuratObject@meta.data[,runDEGsClusters],
                                              seuratObject@meta.data[,runDEGsCondition], 
                                              sep = "_")
    seuratObject$Cell_type <- seuratObject@meta.data[,runDEGsClusters]
    Idents(seuratObject) <- "Cell_type_condition"
    DEGs <- runDEGs(seuratObject = seuratObject, comparisons = DEGComparisons)
    DEGs_sig <- DEGs[DEGs$p_val_adj<0.05,]
    if(length(unique(DEGs_sig$GENE))<80){
      pdf(paste0("Figures/",runDEGsClusters,"_",runDEGsCondition,"_DEGs.pdf"), 
          width = 16, height = 4)
      gene_dot_plot_log2FC(DEG_df = DEGs, var_list = list("GENE"=unique(DEGs_sig$GENE)[65:128]))
      dev.off()
    }
  }
  
  seuratObject@reductions$umap <- rna_umap_red
  seuratObject@reductions$tsne <- rna_tsne_red
  return(seuratObject)
}

# Calculate and plot marker genes (feature plots) (or just plot)

# @param markers: result output from FindAllMarkers()
# @param calculate_clusters: if markers is not set, calculates these (can be
#     more than 1) cluster markers
# @clusters_show: a list of clusters to show for each cluster assignment if
# only specific clusters are to be shown
#     ex. list("integrated_snn_res.0.5"=c(12,4,5,6,9))
makeMarkerFeaturePlots <- function(seuratObject,
                                   markers=NULL, 
                                   calculate_clusters=NULL,
                                   clusters_show=NULL,
                                   nMarkers=10, 
                                   nClustersOneImage=5,
                                   output_name="Figures/Markers/markers.png",
                                   red=c("IntegrateSampleUMAP"),
                                   png_resolution_factor=1,
                                   display_plot=FALSE){
  
  dir = do.call("paste",c(unlist(strsplit(output_name, "/"))[1:length(unlist(strsplit(output_name, "/")))-1],
                          list("sep"="/")))
  ifelse(!dir.exists(file.path(dir)), dir.create(file.path(dir),recursive = T), FALSE)
  
  if(is.null(markers)){
    markers <- data.frame(stringsAsFactors = FALSE)
    for(clust in calculate_clusters){
      Idents(seuratObject) <- clust
      temp <- FindAllMarkers(seuratObject, only.pos = TRUE)
      temp$res <- clust
      markers <- rbind(temp, markers)
    }
    if(!is.null(seuratObject@misc$markers)){
      seuratObject@misc$markers <- rbind(seuratObject@misc$markers,markers)
      seuratObject@misc$markers <- seuratObject@misc$markers[!duplicated(seuratObject@misc$markers),]
    } else {
      seuratObject@misc$markers <- markers
    }
    ifelse(!dir.exists(file.path("Markers/")), dir.create(file.path("Markers/"),recursive = T), FALSE)
    table_output_name <- gsub("Plots/","Markers/",gsub(".png",".txt",output_name))
    if(file.exists(table_output_name)){
      table_output_name <- paste0(gsub(".txt","",paste0(table_output_name,"_additional")),".txt")
      iter = 1
      while(file.exists(table_output_name)){
        table_output_name <- gsub(".txt","",table_output_name)
        if(iter==1){
          table_output_name <- paste0(table_output_name, iter)
        } else{
          table_output_name <- paste0(substr(table_output_name, 
                                             start = 1, 
                                             stop = nchar(table_output_name)-1),
                                      iter)
        }
        table_output_name <- paste0(table_output_name, ".txt")
        iter = iter + 1
      }
    }
    write.table(seuratObject@misc$markers,table_output_name,
                row.names = FALSE, sep = "\t", quote = FALSE)
  }
  assignments <- unique(markers$res)
  for(ass in assignments){
    markersPlot <- markers[markers$res==ass,]
    if(!is.null(clusters_show)){
      clusters <- clusters_show[[ass]]
    } else {
      clusters <- unique(markersPlot$cluster)
    }
    count = 0
    set = 1
    features <- c()
    clusters_plot <- c()
    for(c in clusters){
      count = count + 1
      clusters_plot <- c(clusters_plot, c)
      features <- c(features, markersPlot$gene[markersPlot$cluster==c][1:nMarkers])
      if(count==nClustersOneImage || c==unique(markersPlot$cluster)[length(unique(markersPlot$cluster))]){
        for(r in red){
          manyFeaturesPlot(seuratObject = seuratObject, features = features, 
                           ncol = nMarkers, nrow = nClustersOneImage, 
                           reduction = r, pt.size = 0.1, 
                           labels = 1:(nMarkers*nClustersOneImage),
                           highlight_clusters = clusters_plot,
                           cluster_meta = ass,
                           output_name = paste0(gsub(".png","",output_name),
                                                "_",ass,"_","Set",set,".png"),
                           png_resolution_factor=png_resolution_factor,
                           display_plot = display_plot)
        }
        set = set + 1
        features <- c()
        clusters_plot <- c()
        count = 0
      }
    }
  }
  if(!is.null(calculate_clusters)){
    return(seuratObject)
  }
}

# clusters out of order sometimes
# tried order(as.numeric(unique(seuratObject@active.ident))) - doesn't work
# levels() should work but doesn't
ClusterDistribitionAcrossMethods <- function(seuratObject,
                                             clustering1="RNA_snn_res.0.5", 
                                             clustering2="integrated_snn_res.0.5", 
                                             output_name="Plots/ClusterDistribitionAcrossMethods.png",
                                             reduction1="IntegrateSampleUMAP",
                                             reduction2="umap",
                                             plot_titles_size=15){
  
  dir = do.call("paste",c(unlist(strsplit(output_name, "/"))[1:length(unlist(strsplit(output_name, "/")))-1],
                          list("sep"="/")))
  ifelse(!dir.exists(file.path(dir)), dir.create(file.path(dir),recursive = T), FALSE)
  
  Idents(seuratObject) <- clustering1
  plots <- list()
  for(c in levels(seuratObject@active.ident)[order(as.numeric(levels(seuratObject@active.ident)))]){ # also sometimes out of order, have to do it two different ways
    plots[[paste0("Clustering1_",c)]] <- DimPlot(seuratObject, reduction = reduction1, 
                                             cells.highlight = colnames(seuratObject)[seuratObject@active.ident==c],
                                             sizes.highlight = 0.3) + 
      ggtitle(paste0("Cluster_",c)) +
      theme_cleanDimPlot() +
      theme(plot.title = element_text(size = plot_titles_size))
  }
  
  nClusters1 = length(unique(seuratObject@active.ident))
  nrow = 2
  
  Idents(seuratObject) <- clustering2
  for(c in levels(seuratObject@active.ident)[order(as.numeric(levels(seuratObject@active.ident)))]){ # still out of order?
    plots[[paste0("Clustering2_",c)]] <- DimPlot(seuratObject, reduction = reduction2, 
                                             cells.highlight = colnames(seuratObject)[seuratObject@active.ident==c],
                                             sizes.highlight = 0.3) + 
      ggtitle(paste0("Cluster_",c)) +
      theme_cleanDimPlot() +
      theme(plot.title = element_text(size = plot_titles_size))
  }
  
  nClusters2 = length(unique(seuratObject@active.ident))
  
  ncol = max(nClusters1, nClusters2)
  
  if(nClusters2>nClusters1){ # longer list of clusters goes first
    reordered <- list()
    new_order <- c((nClusters1+1):length(plots), 1:nClusters1)
    for(iter in 1:length(new_order)){
      reordered[[names(plots)[new_order[iter]]]] <- plots[[new_order[iter]]]
    }
  } else {
    reordered <- plots
    rm(plots)
  }
  
  wid = 1280
  hght = (1280*nrow)/(ncol-1)
  
  if(ncol>15){
    wid = wid*5
    hght = hght*5
  }
  
  png(output_name, width = wid, height = hght, units = "px")
  print(ggarrange(plotlist = reordered, ncol = ncol, nrow = nrow))
  dev.off()
  
}

markerGeneOverlap <- function(query_markers, 
                              reference_markers,
                              fdr_threshold=0.05,
                              logfc_threshold=0.25,
                              overlap_plot_wid=20,
                              overlap_plot_hgt=15){
  
  # warning: log2FC and logFC are different - adjust threshold accordingly
  if(!("avg_log2FC" %in% colnames(query_markers))){
    logFCtype = "avg_logFC"
  } else {
    logFCtype = "avg_log2FC"
  }
  query_markers <- query_markers[query_markers$p_val_adj<fdr_threshold &
                                   query_markers[,logFCtype]>logfc_threshold,]
  query_markers <- query_markers[,c("cluster","gene")]
  colnames(query_markers) <- c("MODULE","GENE")
  ref_markers <- read.delim(reference_markers, stringsAsFactors = FALSE)
  
  ifelse(!dir.exists(file.path("markers_temp")), 
         dir.create(file.path("markers_temp"),recursive = T), FALSE)
  if(length(colnames(ref_markers))==2){
    colnames(ref_markers) <- c("module","gene")
  } else{
    if((!("module" %in% colnames(ref_markers)))){
      colnames(ref_markers)[colnames(ref_markers)=="cluster"] <- "module"
    }
    if(ncol(ref_markers)>2){
      if(!("avg_log2FC" %in% colnames(ref_markers))){
        logFCtype = "avg_logFC"
      } else {
        logFCtype = "avg_log2FC"
      }
      ref_markers <- ref_markers[ref_markers$p_val_adj<fdr_threshold &
                                   ref_markers[,logFCtype]>logfc_threshold,]
    }
  }
  write.table(ref_markers, paste0("markers_temp/reference_markers.txt"),
              row.names = FALSE, quote = FALSE, sep = "\t")
  
  output_dir = "InitialCellTypeNaming"
  if(dir.exists(output_dir)){
    iter = 1
    output_dir <- paste0(output_dir, iter)
    while(dir.exists(output_dir)){
        output_dir <- paste0(substr(output_dir, 
                                    start = 1, 
                                    stop = nchar(output_dir)-1),
                                    iter)
      iter = iter + 1
    }
  }
  
  overlap <- makePathwayEnrichmentDf(DEG_df=query_markers, 
                                     resources_path="./markers_temp/",
                                     output_Dir=output_dir, 
                                     convertToHuman=FALSE,
                                     plot_wid = overlap_plot_wid, 
                                     plot_hgt = overlap_plot_hgt)
  predictions <- data.frame()
  for(mod in unique(overlap$Module)[order(unique(overlap$Module))]){
    subset <- overlap[overlap$Module==mod,]
    byOverlap <- subset$Pathway[subset$nOverlap==max(subset$nOverlap)]
    byOverlap <- do.call("paste",c(byOverlap, list("sep"=", ")))
    byEnrichment <- subset$Pathway[subset$Enrichment==max(subset$Enrichment)]
    byEnrichment <- do.call("paste",c(byEnrichment, list("sep"=", ")))
    byFDR <- subset$Pathway[subset$FDR==min(subset$FDR)]
    byFDR <- do.call("paste",c(byFDR, list("sep"=", ")))
    temp <- data.frame("Cluster"=mod, 
                       "ByOverlapPrediction"=byOverlap, 
                       "ByEnrichmentPrediction"=byEnrichment, 
                       "ByFDRPrediction"=byFDR)
    predictions <- rbind(predictions, temp)
  }
  
  return(list("predictions"=predictions,
              "overlap_results"=overlap))
  
}


runSampleWiseDEGs <- function(seuratObject,
                              comparisons = list("F effect"=c("ADF","AD")),
                              Cell_type_column="Cell_type",
                              Condition_column="condition",
                              sampleNameMetaEntry="orig.ident",
                              save_progress=FALSE,
                              name="SampleWiseDEGs",
                              logfc_threshold=0.05,
                              filter_fcdirection=TRUE){
  
  if(grepl("/",name)){
    dir = do.call("paste",c(unlist(strsplit(name, "/"))[1:length(unlist(strsplit(name, "/")))-1],
                            list("sep"="/")))
    ifelse(!dir.exists(file.path(dir)), dir.create(file.path(dir),recursive = T), FALSE)
  }
  
  res <- list()
  for(iter in 1:length(comparisons)){
    SampleWise_DEGs <- runCellTypeDEGsMain(seuratObject = seuratObject, 
                                           groupsCompare = comparisons[[iter]], 
                                           sampleNameMetaEntry = sampleNameMetaEntry, # can be orig.ident 
                                           compareMetaDataLabel = Condition_column, 
                                           cellMetaDataLabel = Cell_type_column,
                                           logfc_threshold = logfc_threshold,
                                           filter_fcdirection=filter_fcdirection)
    
    comp <- paste(comparisons[[iter]][1], comparisons[[iter]][2],sep=".v.")
    saveRDS(SampleWise_DEGs,paste0(name,"_Stringent_compiled_",comp,".rds"))
    
    temp <- data.frame()
    max_ncol = max(unlist(lapply(SampleWise_DEGs[[1]],ncol)))
    celltypes <- names(SampleWise_DEGs[[1]])[unlist(lapply(SampleWise_DEGs[[1]],ncol))==max_ncol]
    leftout <- setdiff(names(SampleWise_DEGs[[1]]),celltypes)
    cat("\nThese cell types are being left out in the combined data frame due to lack of samples/cells:\n",
        do.call("paste",c(leftout,list("sep"=", "))),"\n")
    
    for(ct in setdiff(celltypes, leftout)){
      degs <- SampleWise_DEGs[[1]][[ct]]
      degs$GENE <- rownames(degs)
      rownames(degs) <- NULL
      degs$Cell_type <- ct
      temp <- rbind(temp,degs)
    }
    temp$Comparison <- comp
    write.table(temp, paste0(name,"_Stringent_compiled_",comp,".txt"),
                row.names = FALSE, quote = FALSE, sep = "\t")
    
    for(ct in leftout){
      degs <- SampleWise_DEGs[[1]][[ct]]
      degs$GENE <- rownames(degs)
      rownames(degs) <- NULL
      degs$Cell_type <- ct
      degs$Comparison <- comp
      write.table(degs, paste0(name,"_Stringent_",gsub(" ","_",ct),"_",comp,".txt"),
                  row.names = FALSE, quote = FALSE, sep = "\t")
    }
    
    res[[comp]] <- temp
    
    # tryCatch({res <- rbind(res, temp)}, error=function(e){
    #   message("Cannot combine results from multiple comparisons. Samples/sample number may not be compatible.")
    #   message(e)
    #   return(NA)
    # })
    
    if(save_progress){
      saveRDS(res, paste0(name,"_Stringent_Results_in_progress.rds"))
    }
  }
  
  return(res)
}

CreateSeuratObjectFrom10X <- function(sample_dirs,
                                      meta_data){
  for(sample in sample_dirs){
    
    count.data <- Read10X(data.dir = sample)
    count.data@Dimnames[[2]] = paste0(sample,"_",count.data@Dimnames[[2]])
    sample.seurat <- CreateSeuratObject(counts = count.data, project = sample)
    rm(count.data)
    
    # meta data
    descriptions = unlist(strsplit(sample,split = "-"))
    if(length(descriptions)!=length(meta_data)){
      stop("meta_data length doesn't match sample directory descriptors")
    }
    for(iter in 1:length(meta_data)){
      sample.seurat@meta.data[,meta_data[iter]] <- descriptions[iter]
    }
    
    # merge
    if(!exists("seuratObject")){
      seuratObject <- sample.seurat
      firstSampleName = sample
      firstSample = TRUE
    } 
    else{
      if(firstSample==TRUE){
        seuratObject <- merge(x = seuratObject, y = sample.seurat)
        firstSample = FALSE
      } 
      else{
        seuratObject <- merge(x = seuratObject, y = sample.seurat)
      }
    }
    
    #clean up environment
    rm(sample.seurat)
    
  }
  return(seuratObject)
}

makePlot <- function(string_vector, title, labels){
  string_vector <- gsub("F","F on AD",string_vector)
  string_vector <- gsub("DHA","DHA on ADF",string_vector)
  string_vector <- gsub("NR","NR on ADF",string_vector)
  res <- data.frame(stringsAsFactors = FALSE)
  for(string in string_vector){
    var1 <- concatenate(unlist(strsplit(string,split = "_"))[1:2], mysep = "_")
    var2 <- concatenate(unlist(strsplit(string,split = "_"))[3:4], mysep = "_")
    temp <- data.frame("var1"=var1, "var2"=var2,"count"=unlist(strsplit(string,split = "_"))[5], stringsAsFactors = FALSE)
    res <- rbind(res, temp)
  }
  # res$var1 <- gsub("DOWN AD", expression(""%down%"AD"),res$var1)
  # res$var2 <- gsub("DOWN ", "Down in",res$var2)
  # res$var1 <- gsub("UP", "Up in",res$var1)
  # res$var2 <- gsub("UP", "Up in",res$var2)
  
  res$var1 <- factor(res$var1, levels = c(unique(res$var1[grep("UP",res$var1)]),
                                          unique(res$var1[grep("DOWN",res$var1)])))
  res$var2 <- factor(res$var2, levels = rev(c(unique(res$var2[grep("UP",res$var2)]),
                                              unique(res$var2[grep("DOWN",res$var2)]))))
  
  res$fill <- c(0,0,0,0)
  
  res$color = ifelse(grepl("Down",res$var1), "royalblue3","firebrick3")
  
  g <- ggplot(data = res, mapping = aes(x=var1, y=var2, label=count, fill=fill)) + 
    geom_tile(colour="black") +
    #geom_text(size=2.8) + 
    geom_text(size=5) + 
    #geom_tile() +
    theme_bw() +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, vjust = .5, colour = c("royalblue3","firebrick3"), face="bold"),
          #axis.text.x = element_custom(fills=c("red","blue")),
          axis.text.x = element_text(colour = c("firebrick3","royalblue3"), face="bold"),
          axis.text = element_text(margin=margin(t=0,b=0,r=0,l=0), size=10),
          axis.ticks = element_blank(),
          axis.title = element_blank(), legend.position = "none",
          plot.margin = margin(t=5,b=0,r=1,l=15))+
    #plot.margin = margin(t=20, l=30))+
    #text = element_text(lineheight = .2, size = 7))+
    scale_x_discrete(expand = c(0, 0), position = "top",
                     labels=c(bquote(bold(""%up%.(unlist(strsplit(as.character(res$var1[1]), split = "_"))[2]))),
                              bquote(bold(""%down%.(unlist(strsplit(as.character(res$var1[1]), split = "_"))[2]))))) +
    scale_y_discrete(expand = c(0, 0), 
                     labels=rev(c(bquote(bold(""%up%.(unlist(strsplit(as.character(res$var2[1]), split = "_"))[2]))),
                                  bquote(bold(""%down%.(unlist(strsplit(as.character(res$var2[1]), split = "_"))[2])))))) +
    scale_fill_gradient2(mid="white") #+ ggtitle(title)
  
  #g <- g + theme(axis.text.y = element_custom_y(fills=c("red","blue")))
  g
  return(g)
}

# text is the text inside the tile
# if size is set, pathway result will be a dot plot (size denotes aes to determine dot size)
# select_by can be "numAppearances" (how commonly enriched), "Pval", or "Enrichment"
# fill suggestions: Enrichment or nOverlap
# hclust_by suggestion: Enrichment or nOverlap
pathwayEnrichmentPlot <- function(pathway_df, npathways_plot=50, FDR_cutoff=0.05, 
                                  nOverlap_cutoff = 3, show_pathway_size=FALSE,
                                  fill="Enrichment", 
                                  text=NULL, text_size=4,
                                  size=NULL, 
                                  hclust_by="Pval", select_by="numAppearances",
                                  output_dir="./", 
                                  output_name="pathway_heatmap.pdf",
                                  max_fill_value=NULL,
                                  plot_wid=10, plot_hgt=13,
                                  nchar_pathway_name_cutoff=NULL){
  
  pathway_df$neglog10FDR <- -log10(pathway_df$FDR)
  pathway_df$neglog10FDR[pathway_df$neglog10FDR==Inf] <- 100
  pathway_df$`FDR < 0.05` = pathway_df$FDR < 0.05

  top_pathway_df <- pathway_df[pathway_df$FDR<FDR_cutoff & pathway_df$nOverlap>nOverlap_cutoff,]
  
  top_pathway_df$numAppearances <- sapply(top_pathway_df$Pathway, function(x){
    return(nrow(top_pathway_df[top_pathway_df$Pathway==x,]))
  })
  
  if(select_by=="Pval"){
    top_pathway_df = top_pathway_df[order(top_pathway_df[,select_by], decreasing = FALSE),]
  } else{
    top_pathway_df = top_pathway_df[order(top_pathway_df[,select_by], decreasing = TRUE),]
  }
  
  if(nrow(top_pathway_df)>npathways_plot){
    topPathways <- unique(top_pathway_df$Pathway)[1:npathways_plot]
  } else if(nrow(top_pathway_df)==0) {
    cat("No significant pathways to plot.\n")
  } else{
    topPathways <- unique(top_pathway_df$Pathway)
  }
  
  top_pathway_df = pathway_df[pathway_df$Pathway %in% topPathways,]
  top_pathway_df$Module <- paste0(top_pathway_df$Module, " (",top_pathway_df$ModuleGeneCount,")")
  
  # change pathway names to look better
  top_pathway_df$Pathway <- makePathwayPretty(top_pathway_df$Pathway)
  temp <- top_pathway_df[,c("PathwaySource","Pathway")]
  temp <- temp[!duplicated(temp),]
  dupped_pathways <- temp$Pathway[duplicated(temp$Pathway)]
  if(length(dupped_pathways)>0){
    indices <- which(top_pathway_df$Pathway %in% dupped_pathways)
    for(iter in indices){
      top_pathway_df$Pathway[iter] <- paste(top_pathway_df$PathwaySource[iter], top_pathway_df$Pathway[iter], sep = "_")
    }
  }
  
  if(!is.null(nchar_pathway_name_cutoff)){
    long_pathway_names <- which(nchar(top_pathway_df$Pathway)>nchar_pathway_name_cutoff)
    top_pathway_df$Pathway <- substr(top_pathway_df$Pathway, 1, 50)
    for(iter in long_pathway_names){
      top_pathway_df$Pathway[iter] <- paste0(top_pathway_df$Pathway[iter], "...")
    }
  }
  
  #top_pathway_df <- top_pathway_df[!duplicated(top_pathway_df[,c("Module","Pathway")]),]
  # hclust pathways by nOverlap
  pathway_heat_mat <- convertToMatrix(top_pathway_df[,c("Pathway","Module",hclust_by)]) # can change to Enrichment or FDR
  pathway.dendro <- as.dendrogram(hclust(d = dist(x = pathway_heat_mat)))
  # get order
  pathway_order <- order.dendrogram(pathway.dendro)
  pathway_levels <- rownames(pathway_heat_mat)[pathway_order]
  
  # hclust pathways by Module
  if(length(unique(top_pathway_df$Module))>=2){
    pathway_heat_mat <- convertToMatrix(top_pathway_df[,c("Module","Pathway",hclust_by)])
    pathway.dendro <- as.dendrogram(hclust(d = dist(x = pathway_heat_mat)))
    ct_order <- order.dendrogram(pathway.dendro)
    ct_levels <- rownames(pathway_heat_mat)[ct_order]
    top_pathway_df$Module <- factor(top_pathway_df$Module, levels = ct_levels)
  }
  
  top_pathway_df$Pathway <- factor(top_pathway_df$Pathway, levels = pathway_levels)
  
  if(show_pathway_size){
    top_pathway_df$Pathway <- paste0(top_pathway_df$Pathway, " (",top_pathway_df$PathwayGeneCount,")")
  }
  
  if(!is.null(max_fill_value)){
    top_pathway_df[,fill] <- MinMax(top_pathway_df[,fill], min = min(top_pathway_df[,fill]), max = max_fill_value)
  }
  
  if(!is.null(text)){
    heat <- ggplot(top_pathway_df, aes_string(x="Module", y="Pathway", fill=fill, label=text)) +
      geom_text(size=text_size, fontface="bold")
    
  }
  # } else if (!is.null(size)){
  #   heat <- ggplot(top_pathway_df, aes_string(x="Module", y="Pathway", size=size, fill=fill))
  #   heat <- heat + geom_point(aes(shape=`p(Adj) < 0.05`))
  # } 
    else{
    heat <- ggplot(top_pathway_df, aes_string(x="Module", y="Pathway", fill=fill))
  }
  
  if(!is.null(size)){
    heat <- heat + geom_point(aes(size=size, shape=`FDR < 0.05`))
  } else {
    heat <- heat + geom_tile(colour="grey", size=.7)
  }
  
  heat <- heat + theme_bw() +
    scale_fill_gradient2(low="blue",high="red", mid = "white") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 11.5),
          axis.text.y = element_text(size = 11),
          axis.line = element_blank(), 
          legend.title = element_text(margin=margin(t=0,b=3,r=0,l=0)))+
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(fill=fill)
  pdf(paste0(output_dir,"/",output_name,".pdf"), height = plot_hgt, width = plot_wid)
  print(heat)
  dev.off()
}


#### Doug's Functions ####
# Some modified, see original in HelperFunctions.R


#Get mito, ribo, and predicted % and add to seurat object meta data
getCellQuality <- function(seuratObject=dropEST.combined, 
                           feature_patterns=list("percent.mito"="^mt-",
                                                 "percent.ribo"=c("^Rps","^Rpl"),
                                                 "percent.pred"=c("^Gm1","^Gm2","^Gm3","^Gm4","^Gm5","^Gm6","^Gm7","^Gm8","^Gm9"),
                                                 "percent.Hb"=c("^Hba-","^Hbb-","^Hbq"))){ 
  for(iter in 1:length(feature_patterns)){
    pat <- do.call("paste", c(feature_patterns[[iter]], list("sep"="|")))
    feats <- grep(pattern = pat, x = rownames(x = seuratObject), value = TRUE)
    percentage <- Matrix::colSums(x = GetAssayData(object = seuratObject, slot = "counts")[feats,])/Matrix::colSums(x = GetAssayData(object = seuratObject, slot = "counts"))
    seuratObject[[names(feature_patterns)[iter]]] <- percentage
  }
  
  return(seuratObject)
}


cellQualityPlot <- function(seuratObject=dropEST.combined,fileName=NULL,H=9,W=40,
                            featuresPlot=c("nFeature_RNA","nCount_RNA","percent.mito",
                                           "percent.ribo","percent.pred","percent.Hb"),
                            identPlot="orig.ident",pointSize=0){
  for(feat in featuresPlot){
    if(!(feat %in% colnames(seuratObject@meta.data))){
      featuresPlot <- featuresPlot[featuresPlot!=feat]
    }
  }
  #Make directory
  filePath = unlist(strsplit(x=fileName,split="/"))
  filePath = paste0(filePath[1:(length(filePath)-1)],collapse = "/")
  ifelse(!dir.exists(file.path(filePath)), dir.create(file.path(filePath),recursive = T), FALSE)
  #Generate QCPlots
  seuratObject = SetIdent(seuratObject, value = identPlot)
  pdf(file=fileName,height=H,width=W)
  print(VlnPlot(object = seuratObject, features = featuresPlot ,pt.size = pointSize))
  dev.off()
}


runSampleCCA <- function(seuratObject=dropEST.combined.filtered,combineLevel="timpoint.condition",numFeatures=2000,numDims=30,
                         tissue=NULL,kParam=25,resUse=seq(0.5,4,by=0.5), vars.to.regress=NULL){
  seuratObject = SetIdent(seuratObject,value=combineLevel)
  seurat.list = list()
  for(i in 1:length(levels(seuratObject@active.ident))){
    seurat.list[[i]] = subset(seuratObject,idents=levels(seuratObject@active.ident)[i])
  }
  for (i in 1:length(x = seurat.list)) {
    seurat.list[[i]] <- NormalizeData(object = seurat.list[[i]], verbose = FALSE)
    seurat.list[[i]] <- FindVariableFeatures(object = seurat.list[[i]], selection.method = "vst", nfeatures = numFeatures, verbose = FALSE)
  }
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:numDims) #prob
  seuratIntegrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:numDims)
  library(ggplot2)
  library(cowplot)
  DefaultAssay(object = seuratIntegrated) <- "integrated"
  # Run the standard workflow for visualization and clustering
  seuratIntegrated <- ScaleData(object = seuratIntegrated, verbose = FALSE, vars.to.regress=vars.to.regress)
  seuratIntegrated <- RunPCA(object = seuratIntegrated, npcs = numDims, verbose = FALSE)
  seuratIntegrated <- RunUMAP(object = seuratIntegrated, reduction = "pca", dims = 1:numDims)
  seuratIntegrated <- RunTSNE(object = seuratIntegrated, reduction = "pca", dims = 1:numDims)
  seuratIntegrated <- FindNeighbors(object = seuratIntegrated, reduction = "pca", dims = 1:numDims, k.param = kParam)
  seuratIntegrated <- FindClusters(object = seuratIntegrated, resolution = resUse, verbose = T, reduction = "pca")
  ifelse(!dir.exists(file.path(paste0("./Plots/ccaTSNEcondition"))), dir.create(file.path(paste0("./Plots/ccaTSNEcondition")),recursive = T), FALSE)
  ifelse(!dir.exists(file.path(paste0("./Plots/ccaUMAPcondition"))), dir.create(file.path(paste0("./Plots/ccaUMAPcondition")),recursive = T), FALSE)
  ifelse(!dir.exists(file.path(paste0("./Plots/ccaTSNEresolution/",tissue))), dir.create(file.path(paste0("./Plots/ccaTSNEresolution/",tissue)),recursive = T), FALSE)
  ifelse(!dir.exists(file.path(paste0("./Plots/ccaUMAPresolution/",tissue))), dir.create(file.path(paste0("./Plots/ccaUMAPresolution/",tissue)),recursive = T), FALSE)
  pdf(file=paste0("./Plots/ccaUMAPcondition/",tissue,".pdf"))
  print(DimPlot(object = seuratIntegrated, reduction = "umap", group.by = combineLevel))
  dev.off()
  pdf(file=paste0("./Plots/ccaTSNEcondition/",tissue,".pdf"))
  print(DimPlot(object = seuratIntegrated, reduction = "tsne", group.by = combineLevel))
  dev.off()
  resolutions <- as.character(resUse)
  for(reso in resolutions){
    pdf(file=paste0("./Plots/ccaTSNEresolution/",tissue,"/PCA_tSNEres",reso,".pdf"))
    seuratIntegrated <- SetIdent(seuratIntegrated, value = paste0("integrated_snn_res.",reso))
    print(DimPlot(seuratIntegrated, label = T, reduction = "tsne"))
    dev.off()
    pdf(file=paste0("./Plots/ccaUMAPresolution/",tissue,"/PCA_UMAPres",reso,".pdf"))
    seuratIntegrated <- SetIdent(seuratIntegrated, value = paste0("integrated_snn_res.",reso))
    print(DimPlot(seuratIntegrated, label = T, reduction = "umap"))
    dev.off()
  }
  return(seuratIntegrated)
}


# Main DEG function -- can wrap it to include additional variable levels (e.g. timepoint and condition)
runCellTypeDEGsMain <- function(seuratObject=metaObject1Subset,metaData1=NULL,
                                cellMetaDataLabel="CellTypeRefined",metaObject1=NULL,
                                groupsCompare=NULL,sampleNameMetaEntry=NULL,compareMetaDataLabel=NULL,minCounts=5,
                                cellTypeDEGsNovelMethod=NULL,cellTypeDEGsStandard=NULL,logfc_threshold=0.1,
                                filter_fcdirection=TRUE){
  
  Idents(seuratObject) <- cellMetaDataLabel
  
  metaObject1Subset = seuratObject
  
  #iterate through each cell type
  for(cellType in levels(metaObject1Subset@active.ident)){
    print(cellType)
    #generate a cell type subset
    cellTypeSubset = subset(metaObject1Subset,idents = cellType)
    #get Group1 sample names
    group1Samples = unique(cellTypeSubset@meta.data[which(cellTypeSubset@meta.data[,compareMetaDataLabel]==groupsCompare[1]),sampleNameMetaEntry])
    #get Group2 sample names
    group2Samples = unique(cellTypeSubset@meta.data[which(cellTypeSubset@meta.data[,compareMetaDataLabel]==groupsCompare[2]),sampleNameMetaEntry])
    
    #set identity to sample level
    cellTypeSubset = SetIdent(cellTypeSubset, value = sampleNameMetaEntry)
    
    #remove samples which have less than 5 cells
    sampleCounts = table(cellTypeSubset@active.ident)
    samplesRemove = names(sampleCounts)[which(sampleCounts<minCounts)]
    
    group1Overlap = which(group1Samples %in% samplesRemove)
    if(length(group1Overlap)>0){
      group1Samples = group1Samples[-group1Overlap]
    }
    
    group2Overlap = which(group2Samples %in% samplesRemove)
    if(length(group2Overlap)>0){
      group2Samples = group2Samples[-group2Overlap]
    }
    
    if(length(group1Samples) > 0 & length(group2Samples) > 0){
      
      #list for overlapping genes
      allGenes = list()
      
      #get DEGs between each group1 sample and all group2 samples using Wilcoxon rank sum test
      group1DEGs = list()
      for(i in 1:length(group1Samples)){
        group1DEGs[[group1Samples[i]]] = FindMarkers(cellTypeSubset, ident.1 = group1Samples[i], ident.2 = group2Samples, logfc.threshold=logfc_threshold)
        colnames(group1DEGs[[group1Samples[i]]]) = paste0(group1Samples[i],colnames(group1DEGs[[group1Samples[i]]]))
        allGenes[[group1Samples[i]]] = rownames(group1DEGs[[group1Samples[i]]])
      }
      
      #get DEGs between each group2 sample and all group1 samples using Wilcoxon rank sum test
      group2DEGs = list()
      for(i in 1:length(group2Samples)){
        group2DEGs[[group2Samples[i]]] = FindMarkers(cellTypeSubset, ident.1 = group2Samples[i], ident.2 = group1Samples, logfc.threshold=logfc_threshold)
        colnames(group2DEGs[[group2Samples[i]]]) = paste0(group2Samples[i],colnames(group2DEGs[[group2Samples[i]]]))
        allGenes[[group2Samples[i]]] = rownames(group2DEGs[[group2Samples[i]]])
      }
      
      #get overlapping genes from all the statistical tests at the sample level
      overlappingGenes = Reduce(intersect, allGenes)
      
      #generate a matrix with the pvals, FC, and adj p vals
      for(i in 1:length(group1Samples)){
        if(i==1){
          combinedMatrix = group1DEGs[[group1Samples[i]]][overlappingGenes,c(1,2,5)]
        } else{
          combinedMatrix = cbind(combinedMatrix,group1DEGs[[group1Samples[i]]][overlappingGenes,c(1,2,5)])
        }
      }
      for(i in 1:length(group2Samples)){
        combinedMatrix = cbind(combinedMatrix,group2DEGs[[group2Samples[i]]][overlappingGenes,c(1,2,5)])
      }
      
      colsUseGroup1 = seq(from = 1, to = length(group1Samples)*3, by = 3)
      colsUseGroup2 = seq(from = length(group1Samples)*3+1, to = length(group1Samples)*3+1+(length(group2Samples)-1)*3, by = 3)
      
      minP = apply(combinedMatrix, 1, function(x) minimump(x[c(colsUseGroup1,colsUseGroup2)])$p)
      maxP = apply(combinedMatrix, 1, function(x) maximump(x[c(colsUseGroup1,colsUseGroup2)])$p)
      
      minPadj = p.adjust(minP,method = 'bonferroni',n = nrow(cellTypeSubset@assays$RNA@data))
      maxPadj = p.adjust(maxP,method = 'bonferroni',n = nrow(cellTypeSubset@assays$RNA@data))
      #add to the combined matrix
      combinedMatrix$minP = minP
      combinedMatrix$maxP = maxP
      combinedMatrix$minPadj = minPadj
      combinedMatrix$maxPadj = maxPadj
      
      #get the conversative p-val (max of minP and maxP)
      ConsP = apply(combinedMatrix, 1, function(x) max(x[c('minP','maxP')]))
      combinedMatrix$ConsP = ConsP
      
      if(filter_fcdirection){
        #only keep genes which show the FC in the same direction
        colsUseGroup1 = seq(from = 2, to = length(group1Samples)*3, by = 3)
        colsUseGroup2 = seq(from = length(group1Samples)*3+2, to = length(group1Samples)*3+2+(length(group2Samples)-1)*3, by = 3)
        
        if(length(colsUseGroup1)==1){
          goodGroup1 = rownames(combinedMatrix)
        } else{
          goodGroup1 = names(which(rowSums(combinedMatrix[,colsUseGroup1] > 0)==length(group1Samples) | rowSums(combinedMatrix[,colsUseGroup1] < 0)==length(group1Samples)))
        }
        if(length(colsUseGroup2)==1){
          goodGroup2 = rownames(combinedMatrix)
        } else{
          goodGroup2 = names(which(rowSums(combinedMatrix[,colsUseGroup2] > 0)==length(group2Samples) | rowSums(combinedMatrix[,colsUseGroup2] < 0)==length(group2Samples)))
        }
        correctDirection = intersect(goodGroup1,goodGroup2)
        combinedMatrix = combinedMatrix[correctDirection,]
      }
      
      
      #sort the matrix by the max pvalue
      combinedMatrix = combinedMatrix[order(combinedMatrix$ConsP),]
      
      if(!(is.null(metaData1))){
        cellTypeDEGsNovelMethod[[metaObject1]][[cellType]] = combinedMatrix
      } else{
        cellTypeDEGsNovelMethod[[cellType]] = combinedMatrix
      }
    }
    # #do traditional Seurat DEG identification
    # cellTypeSubset = SetIdent(cellTypeSubset, value = compareMetaDataLabel)
    # #only proceed if both have at least 3 cells
    # statusCounts = table(cellTypeSubset@active.ident)
    # if(length(which(statusCounts>=3))==2){
    #   StandardDEGs = FindMarkers(cellTypeSubset, ident.1 = groupsCompare[1], ident.2 = groupsCompare[2])
    #   #add cell type DEGs to matrices
    #   if(!(is.null(metaData1))){
    #     cellTypeDEGsStandard[[metaObject1]][[cellType]] = StandardDEGs
    #   } else{
    #     cellTypeDEGsStandard[[cellType]] = StandardDEGs
    #   }
    # }
  }
  #return(list(cellTypeDEGsNovelMethod,cellTypeDEGsStandard))
  return(list(cellTypeDEGsNovelMethod))
}


#Wrapper for the main DEG function -- allows us to add an additional meta data level (e.g. timepoint + condition rather than just condition)
runCellTypeDEGsNovelMethod <- function(seuratObject=dropEST.combined.filtered,metaData1=NULL,cellMetaDataLabel="CellTypeRefined",
                                       groupsCompare=NULL,sampleNameMetaEntry=NULL,compareMetaDataLabel=NULL,minCounts=5,tissue=NULL){
  
  #initialize a list to hold marker genes for all cell types
  cellTypeDEGsStandard = list()
  cellTypeDEGsNovelMethod = list()
  
  #iterate through the first supplied level of meta data
  if(!(is.null(metaData1))){
    seuratObject = SetIdent(seuratObject, value = metaData1)
    for(metaObject1 in levels(seuratObject@active.ident)){
      print(metaObject1)
      metaObject1Subset = subset(seuratObject,idents = metaObject1)
      metaObject1Subset = SetIdent(metaObject1Subset, value = cellMetaDataLabel)
      outputList = runCellTypeDEGsMain(seuratObject=metaObject1Subset,metaData1=metaData1,cellMetaDataLabel=cellMetaDataLabel,metaObject1=metaObject1,
                                       groupsCompare=groupsCompare,sampleNameMetaEntry=sampleNameMetaEntry,compareMetaDataLabel=compareMetaDataLabel,
                                       minCounts=minCounts,cellTypeDEGsNovelMethod=cellTypeDEGsNovelMethod,cellTypeDEGsStandard=cellTypeDEGsStandard)
      cellTypeDEGsNovelMethod = outputList[[1]]
      cellTypeDEGsStandard = outputList[[2]]
    }
  } else{
    metaObject1Subset = SetIdent(seuratObject, value = cellMetaDataLabel)
    outputList = runCellTypeDEGsMain(seuratObject=metaObject1Subset,metaData1=metaData1,cellMetaDataLabel=cellMetaDataLabel,metaObject1=NULL,
                                     groupsCompare=groupsCompare,sampleNameMetaEntry=sampleNameMetaEntry,compareMetaDataLabel=compareMetaDataLabel,
                                     minCounts=minCounts,cellTypeDEGsNovelMethod=cellTypeDEGsNovelMethod,cellTypeDEGsStandard=cellTypeDEGsStandard)
    cellTypeDEGsNovelMethod = outputList[[1]]
    cellTypeDEGsStandard = outputList[[2]]
  }
  #save DEGs
  saveRDS(cellTypeDEGsNovelMethod, file = paste0("./Saves/",tissue,".DEGsNovelMethod.rds"))
  saveRDS(cellTypeDEGsStandard, file = paste0("./Saves/",tissue,".DEGsStandardMethod.rds"))
}



#### Dan's Functions ####

export2SingleCellPortal<-function(seurat.obj,
                                  species, species_ontology_label,
                                  disease, disease_ontology_label,
                                  organ,organ.label,
                                  library_preparation_protocol,library_preparation_protocol_ontology_label,
                                  sex, donor_id,
                                  outdir="./",
                                  generateCountMatrix=F
){
  numb.rep<-length(colnames(seurat.obj))
  meta.scp<-seurat.obj@meta.data
  
  
  name<-c("TYPE",colnames(seurat.obj))
  biosample_id<-c("group",seurat.obj$orig.ident)
  nCount_RNA<-c("numeric",seurat.obj@meta.data$nCount_RNA)
  nFeature_RNA<-c("numeric",seurat.obj@meta.data$nFeature_RNA)
  percent_mt<-c("numeric",seurat.obj@meta.data$percent.mt)
  
  seurat_clusters<-c("group",seurat.obj$seurat_clusters)
  if("condition" %in% names(seurat.oec@meta.data)){
    condition<-c("group",seurat.obj$condition)
  }
  if("Cell_type" %in% names(seurat.oec@meta.data)){
    cell_type<-c("group",seurat.obj$Cell_type)
  }
  donor_id<-c("group",colnames(seurat.obj))
  species<-c("group",rep(species,numb.rep))
  species_ontology_label<-c("group",rep(species_ontology_label,numb.rep))
  disease<-c("group",rep(disease,numb.rep))
  disease_ontology_label<-c("group",rep(disease_ontology_label,numb.rep))
  organ<-c("group",rep(organ,numb.rep))
  organ.label<-c("group",rep(organ.label,numb.rep))
  library_preparation_protocol<-c("group",rep(library_preparation_protocol,numb.rep))
  library_preparation_protocol__ontology_label<-c("group",rep(library_preparation_protocol_ontology_label,numb.rep))
  sex<-c("group",rep(sex,numb.rep))
  
  
  meta.scp <- data.frame(NAME=name,biosample_id=biosample_id,
                         nCount_RNA=nCount_RNA,nFeature_RNA=nFeature_RNA,
                         percent_mt=percent_mt,seurat_clusters=seurat_clusters,condition=condition,
                         Cell_type=cell_type,donor_id=donor_id,
                         species=species,species__ontology_label=species_ontology_label,
                         disease=disease,disease__ontology_label=disease_ontology_label,
                         organ=organ,organ__ontology_label=organ.label,
                         library_preparation_protocol=library_preparation_protocol,
                         library_preparation_protocol__ontology_label=library_preparation_protocol__ontology_label,
                         sex=sex)
  DefaultAssay(seurat.obj)<-"RNA"
  data.table<- as.data.frame(GetAssayData(seurat.obj, slot = "data"))
  
  data.table<-cbind(GENE=rownames(data.table),data.table)
  
  write.table(data.table,file=paste0(outdir,"/","1.exp.matrix.tsv"),sep="\t",quote = F,row.names=F)
  
  
  if(generateCountMatrix){
    data.table<- as.data.frame(GetAssayData(seurat.obj, slot = "count"))
    data.table<-cbind(GENE=rownames(data.table),data.table)
    write.table(data.table,file=paste0(outdir,"/","1.1.count.matrix.tsv"),sep="\t",quote = F,row.names=F)
  }
  
  x<-c("numeric",seurat.obj@reductions$tsne@cell.embeddings[,1])
  y<-c("numeric",seurat.obj@reductions$tsne@cell.embeddings[,2])
  tsne<-data.frame(NAME=name, X=x,Y=y)
  write.table(tsne,
              file=paste0(outdir,"/","3.tsne.txt"),sep="\t",quote = F,row.names=F)
  
  x<-c("numeric",seurat.obj@reductions$umap@cell.embeddings[,1])
  y<-c("numeric",seurat.obj@reductions$umap@cell.embeddings[,2])
  umap<-data.frame(NAME=name, X=x,Y=y)
  write.table(umap,
              file=paste0(outdir,"/","3.umap.txt"),sep="\t",quote = F,row.names=F)
  write.table(meta.scp,
              file=paste0(outdir,"/","2.metadata.txt"),sep="\t",quote = F,row.names =F )
  print("File saved at :")
  print(paste0(outdir,"/","1.exp.tsv"))
  print(paste0(outdir,"/","2.metadata.txt"))
  print(paste0(outdir,"/","3.tsne.txt"))
  print(paste0(outdir,"/","3.umap.txt"))
}



#### Mergeomics Integration ####

# Marker Dependency Filtering
# 
# Filters markers based on dependency (such as linkage disequilibrium) file given. Prepares
# files for MSEA (mapping and loci files)
#
# @param MARFILE: association file with 'LOCUS' and 'VALUE' headers. Value reflects strength of association.
#                 (for example, can use -log10 pvalues)
# @param GENFILE: file maping loci to genes - 'GENE' 'LOCUS'
# @param LNKFILE: marker dependency file 'LOCUSa' 'LOCUSb' 'WEIGHT'
# @param trait_name: optional, becomes part of the output directory name with mapping_name
# @param mapping_name: optional, becomes part of the output directory name with trait_name
# @param nmax: percent of associations to consider
# @param ld_threshold: corresponds to LD threshold used - must be changed if not 50%, 
#                     used to name output files
# @param ldprune: path to ldprune program
#
# @return creates directory with name that combines trait and mapping information and that 
#         contains the mapping (-.g.txt) and loci (-.l.txt) file with percent associations
#         and ld threshold information appended (ex. 50.50.g.txt)
#
# @examples
# runMDF(MARFILE = "./GWAS/Kunkle_AD.txt",
#        GENFILE = "./mapping/Brain_Hippocampus.eQTL.txt", 
#        LNKFILE = "./linkage/LD50.1000G.CEU.txt", 
#        output_dir = "./MSEA/Data/",
#        ldprune = "./MDPRUNE/ldprune")
#
runMDF <-function(MARFILE,
                  GENFILE,
                  LNKFILE,
                  trait_name=NULL, 
                  mapping_name=NULL, 
                  output_dir="./MSEA/Data/",
                  nmax=0.5,
                  ld_threshold=50,
                  editFiles=FALSE,
                  mdprune="mdprune"){
  
  # change to appropriate headers
  if(editFiles){
    system(paste0("sed -i \"1s/.*/GENE\\tMARKER/\" ",GENFILE))
    system(paste0("sed -i \"1s/.*/MARKERa\\tMARKERb\\tWEIGHT/\" ",LNKFILE))
  }
  
  if(is.null(trait_name) & is.null(mapping_name)){
    trait_name = unlist(strsplit(MARFILE,"/"))[length(unlist(strsplit(MARFILE,"/")))]
    trait_name = gsub(".txt","",trait_name)
    
    mapping_name = unlist(strsplit(GENFILE,"/"))[length(unlist(strsplit(GENFILE,"/")))]
    mapping_name = gsub(".txt","",mapping_name)
  }
  
  label=paste(output_dir,trait_name,'.',mapping_name,sep="")
  ifelse(!dir.exists(label), dir.create(label, recursive = TRUE),FALSE)
  
  bash_file <- file(paste0(label,".bash"))
  writeLines(c(paste('MARFILE="',MARFILE,'"',sep=''),
               paste('GENFILE="',GENFILE,'"',sep=''),
               paste('LNKFILE="', LNKFILE,'"',sep=""),
               paste('OUTPATH="',output_dir,trait_name,'.',mapping_name,'/"',sep=""),
               paste('NTOP=',nmax,sep=""),
               paste("echo -e \"MARKER\\tVALUE\" > /tmp/header.txt"),
               paste("nice sort -r -g -k 2 $MARFILE > /tmp/sorted.txt"),
               paste("NMARKER=$(wc -l < /tmp/sorted.txt)"),
               paste("NMAX=$(echo \"($NTOP*$NMARKER)/1\" | bc)"),
               paste("nice head -n $NMAX /tmp/sorted.txt > /tmp/top.txt"),
               paste("cat /tmp/header.txt /tmp/top.txt > /tmp/subset.txt"),
               paste('nice ',mdprune," /tmp/subset.txt $GENFILE $LNKFILE $OUTPATH",sep="")),
             bash_file)
  close(bash_file)
  
  ifelse(!dir.exists(output_dir), dir.create(output_dir, recursive = TRUE),FALSE)
  
  name=paste(as.character(nmax*100), ld_threshold, sep=".")
  system(paste0("bash ", paste0(label,".bash")))
  system(paste("mv ",label,"/genes.txt ",label,"/",name,".g.txt",sep=""))
  system(paste("mv ",label,"/marker.txt ",label,"/",name,".m.txt",sep=""))
}

GWAS_MSEA_screen <- function(GWAS_dir="/u/project/xyang123/xyang123-NOBACKUP/icestrik/public_gwas/MSEA_Gtex_v7/",
                             SNP_map_method = "Adipose_Subcutaneous",
                             loci_pattern = ".l.txt",
                             gene_pattern = ".g.txt",
                             geneset="", 
                             info = NULL,
                             ntop = 50,
                             ld_threshold = 50,
                             output_dir = "./GWAS_Screen_Results/",
                             output_name = "Results",
                             max_module_genes = 5000,
                             trim = .005,
                             nperm = 10000,
                             permtype = "gene",
                             nGWASshow = 20){
  
  traits = list.files(GWAS_dir)[grep(SNP_map_method, list.files(GWAS_dir))]
  trim_traits = c()
  for(trait in traits){
    trim_traits = append(trim_traits, unlist(strsplit(trait, split = ".", fixed = TRUE))[1])
  }
  traits = trim_traits
  rm(trim_traits)
  para = paste0(ld_threshold,".",ntop)
  mapping = SNP_map_method
  
  ifelse(!dir.exists(output_dir), dir.create(output_dir), FALSE)
  setwd(output_dir)
  for(trait.item in traits){
    for (item.1 in mapping){
      for (item.2 in para){
        job.ssea <- list()
        job.ssea$label <- paste(trait.item,item.1,item.2,sep=".")
        job.ssea$folder <- output_dir
        job.ssea$genfile <- paste(GWAS_dir,trait.item,".",item.1,"/",item.2,gene_pattern,sep="")
        job.ssea$marfile <- paste(GWAS_dir,trait.item,".",item.1,"/",item.2,loci_pattern,sep="")		
        job.ssea$modfile <- geneset
        if(!is.null(info)){
          job.ssea$inffile <- info
        }
        job.ssea$permtype <- permtype	#optional
        job.ssea$maxgenes <- max_module_genes	#optional
        job.ssea$nperm <- nperm	#optional
        job.ssea <- ssea.start(job.ssea)
        job.ssea <- ssea.prepare(job.ssea)
        job.ssea <- ssea.control(job.ssea)
        job.ssea <- ssea.analyze(job.ssea,trim_start=trim,trim_end=(1-trim))
        job.ssea <- ssea.finish(job.ssea)
      }
    }
  }
  
  # make plot of results
  result_files <- list.files(paste0(output_dir, "msea/"))[grep("results", list.files(paste0(output_dir, "msea/")))]
  results_df = data.frame(stringsAsFactors = FALSE)
  for(result in result_files){
    msea <- read.delim(paste0(output_dir,"msea/", result), stringsAsFactors = FALSE)
    study = unlist(strsplit(result, split = ".", fixed = TRUE))[1]
    temp = data.frame("Study" = study, 
                      "Module" = msea$MODULE,
                      "FDR" = msea$FDR)
    results_df = rbind(results_df, temp)
  }
  
  results_df = results_df[results_df$Module!="_ctrlA",]
  results_df = results_df[results_df$Module!="_ctrlB",]
  
  results_df$log10FDR <- -log10(results_df$FDR)
  results_df <- results_df[order(results_df$log10FDR, decreasing = TRUE),]
  # get top GWAS studies 
  if(length(unique(results_df$Study))<nGWASshow){
    show_GWAS = unique(results_df$Study)
  } else {
    show_GWAS = unique(results_df$Study)[1:nGWASshow]
  }
  
  toPlot = data.frame(stringsAsFactors = FALSE)
  for(gwas in show_GWAS){
    toPlot = rbind(toPlot, results_df[results_df$Study==gwas,])
  }
  
  toPlot$FDR = toPlot$FDR/nGWASshow
  toPlot$log10FDR <- -log10(toPlot$FDR)
  
  heat <- ggplot(toPlot,aes(x=Module,y = Study, fill=log10FDR)) + geom_tile(colour="gray") + 
    scale_fill_gradientn(colors = c("white", "red")) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 15, vjust = 0, hjust = 0, 
                                     margin = margin(t=0,b=0,r=0,l=0)), 
          axis.text.y =  element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), panel.border = element_blank(),
          axis.title = element_text(size=15, face = "bold")) + 
    labs(fill=expression(-log[10]~FDR)) +
    scale_x_discrete(position = "top") 
  
  pdf(paste0(output_name,"_heatmap.pdf"), width = 10, height = 7)
  print(heat)
  dev.off()
  
  write.table(results_df,paste0(output_name,"_data.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
}


# nodes_file can either be .results.txt file from MSEA or a 'MODULE' 'NODE' file
runKDA <- function(MSEA_results_dir = NULL,
                   nodes_file, # can either be .results.txt file from MSEA or a 'MODULE' 'NODE' file
                   marker_set=NULL,
                   network, 
                   name = NULL,
                   trim = c(".results.txt",".50.20",".50.50",".txt",".results"), 
                   merge = FALSE, 
                   rcutoff=.33, 
                   info=NULL,
                   edgefactor = 0,
                   depth = 1,
                   direction = 0,
                   nperm = 10000,
                   nKDs_show = NULL,
                   return_job = FALSE
){
  
  job.kda <- list()
  job.kda$label<-"wKDA"
  
  if(!is.null(MSEA_results_dir)){
    results = list.files(MSEA_results_dir, 
                         full.names = TRUE)[grep(".results.txt",
                                                 list.files(MSEA_results_dir, 
                                                            full.names = TRUE))]
    if(length(results)>1) cat("Error: More than one MSEA results files detected.\n")
  } else{
    results = nodes_file
  }
  
  
  result <- read.delim(results, stringsAsFactors = FALSE)
  if(ncol(result)>4){ # result file from MSEA
    sig_mods = result$MODULE[result$FDR<0.05]
    sig_modfile = data.frame(stringsAsFactors = FALSE)
    
    modfile <- read.delim(marker_set, stringsAsFactors = FALSE)
    
    for(mod in sig_mods){
      sig_modfile = rbind(sig_modfile, modfile[modfile$MODULE==mod,])
    }
    colnames(sig_modfile) <- c("MODULE","NODE")
    write.table(sig_modfile, "Significant_module_nodes.txt", row.names = FALSE, quote = FALSE, sep = "\t")
    job.kda$modfile <- "./Significant_module_nodes.txt"
    forMerging = data.frame("MODULE"=unique(sig_modfile$MODULE), stringsAsFactors = FALSE)
  }
  else if(ncol(result)==4){ # output from merge modules
    result <- result[,c("MODULE","NODE")]
    write.table(result, "nodes_file_forKDA.txt", row.names = FALSE, quote = FALSE, sep = "\t")
    job.kda$modfile <- "nodes_file_forKDA.txt"
  }
  else{
    if(length(intersect(colnames(result),c("MODULE","NODE")))!=2){
      colnames(result) <- c("MODULE","NODE")
      write.table(result, "nodes_file_forKDA.txt", row.names = FALSE, quote = FALSE, sep = "\t")
      job.kda$modfile <- "nodes_file_forKDA.txt"
    }else{
      job.kda$modfile <- nodes_file
    }
    
    forMerging <- read.delim(nodes_file, stringsAsFactors = FALSE)
    forMerging = data.frame("MODULE"=unique(forMerging$MODULE), stringsAsFactors = FALSE)
  }
  
  if(merge){
    if(ncol(result)==2){
      modfile <- result
      colnames(modfile) <- c("MODULE","GENE")
      write.table(modfile,"modfile_forMerging.txt",
                  row.names = FALSE, quote = FALSE, sep = "\t")
      marker_set = "modfile_forMerging.txt"
    }else{
      modfile <- read.delim(marker_set, stringsAsFactors = FALSE)
    }
    
    if(is.null(info)){
      infofile = data.frame("MODULE"=unique(modfile$MODULE), "SOURCE"=unique(modfile$MODULE), "DESCR"=unique(modfile$MODULE))
      write.table(infofile,"infofile_forMerging.txt", row.names = FALSE, quote = FALSE, sep = "\t")
      merge_modules(name = "nodes", modules_df = forMerging, rcutoff = rcutoff, output_Dir = "./", 
                    modfile_path = marker_set, infofile_path = "./infofile_forMerging.txt")
    }
    else{
      merge_modules(name = "nodes", modules_df = forMerging, rcutoff = rcutoff, 
                    output_Dir = "./", modfile_path = marker_set, infofile_path = info)
    }
    job.kda$modfile <- "merged_nodes.mod.txt"
  }
  
  nodes_name = nodes_file
  nodes_name = unlist(strsplit(nodes_name, split = "/"))[length(unlist(strsplit(nodes_name, split = "/")))]
  for(i in trim){
    nodes_name = gsub(i,"",nodes_name)
  }
  
  network_name = gsub(".txt","",unlist(strsplit(network, split = "/"))[length(unlist(strsplit(network, split = "/")))])
  
  cat("\nNow analyzing:",nodes_name, "with", network_name, "\n")
  
  if(is.null(name)){
    name = paste(nodes_name, network_name,sep = "_")
  }
  
  job.kda$folder<- name ## parent folder for results
  
  ## Input a network
  ## columns: TAIL HEAD WEIGHT
  job.kda$netfile <- network
  
  
  ## "0" means we do not consider edge weights while 1 is opposite.
  job.kda$edgefactor <- edgefactor
  ## The searching depth for the KDA
  job.kda$depth <- depth
  ## 0 means we do not consider the directions of the regulatory interactions
  ## while 1 is opposite.
  job.kda$direction <- direction
  job.kda$nperm <- nperm
  
  moddata <- tool.read(job.kda$modfile)
  mod.names <- unique(moddata$MODULE)
  moddata <- moddata[which(!is.na(match(moddata$MODULE, mod.names))),]
  ## save this to a temporary file and set its path as new job.kda$modfile:
  tool.save(moddata, "subsetof.supersets.txt")
  job.kda$modfile <- "subsetof.supersets.txt"
  
  ## Run KDA
  job.kda <- kda.configure(job.kda)
  job.kda <- kda.start(job.kda)
  job.kda <- kda.prepare(job.kda)
  job.kda <- kda.analyze(job.kda)
  job.kda <- kda.finish(job.kda)
  
  if(!is.null(nKDs_show)){
    job.kda <- kda2cytoscape(job.kda, ndrivers = nKDs_show)
  }
  else{
    job.kda <- kda2cytoscape(job.kda)
  }
  
  if(return_job){
    saveRDS(job.kda, file = paste0(name, "_kda_job.rds"))
  }
}