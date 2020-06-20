################ FUNCTIONS TO SUMMARIZE DEGS ###########
# Written by Jessica Ding, updated August 2019 - not directly usable yet for others
# for all functions, input is a dataframe output from DEG analysis (FindMarkers)
#     Column Names: p_val   avg_logFC   pct.1   pct.2   p_val_adj   Cell_type   GENE    Comparison
# Comparison is an optional column if multiple comparisons are done

library(ggplot2)
library(varhandle)
library(WriteXLS)
library(reshape2)

# requires biomaRt
# human genes that map to different mouse genes: "ACAT2"  "SPCS3"  "FTL"    "HSPA1A" "BCL2A1" "UBA52"  "LILRB3" "CCL15"  "CD33"   "GBP6"   "IFITM1" "HBA2"
# the resulting dataframe can be used for pathway enrichment
addInfoAndTrim <- function(DEG_df, FDR_cutoff=0.05, multiple_comparisons=TRUE, convertToHuman=TRUE, lfc_cutoff=0.1){
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
        module[row] = paste(DEG_df$Cell_type[row], "DOWN", sep = "_")
      }
    }
  }
  DEG_df$MODULE = module
  DEG_df$MODULE_ct = module_ct
  DEG_df$MODULE_comp = module_comp
  DEG_df$MODULE_no_direct = paste(DEG_df$Cell_type, DEG_df$Comparison)
  DEG_df = DEG_df[order(DEG_df$MODULE),]
  
  if(convertToHuman){
    temp <- convertDfGeneColumnMouseHuman(df = DEG_df, toSpecies = "human", forPathway = FALSE)
    DEG_df$HUMAN = temp$GENE
  }
  return(DEG_df)
}

# MODULE contains all of the meta data information without direction info and can be used for pathway enrichment
addInfoAndTrim <- function(DEG_df, FDR_cutoff=0.05, meta_data, convertToHuman=TRUE, lfc_cutoff=0.1){
  DEG_df = DEG_df[DEG_df$p_val_adj<FDR_cutoff,]
  DEG_df = DEG_df[abs(DEG_df$avg_logFC)>lfc_cutoff,]
  DEG_df$Direction = ifelse(DEG_df$avg_logFC>0,"UP","DOWN")
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
  
  if(convertToHuman){
    temp <- convertDfGeneColumnMouseHuman(df = DEG_df, toSpecies = "human", forPathway = FALSE)
    DEG_df$HUMAN = temp$GENE
  }
  return(DEG_df)
}

# wrapper for pathway_enrichment
# identifer - string - if running multiple versions, add identifer to differentiate
run_pathway_enrichment <- function(DEG_df, FDR_cutoff=NULL, lfc_cutoff=NULL, module="MODULE", identifier){
  iter = 1
  deg_list = DEG_df
  if(!is.null(lfc_cutoff)){
    deg_list = deg_list[abs(deg_list$avg_logFC)>lfc_cutoff,]
  }
  deg_list$MODULE <- deg_list[,module]
  for(MOD in unique(deg_list$MODULE)){
    cat(iter,": ",MOD, "\n")
    pathway_enrichment(deg_list = deg_list, 
                       FDR_threshold = FDR_cutoff, 
                       pval_threshold = 0.05, 
                       celltype_cluster = MOD, 
                       pVal = FALSE,
                       identifier = identifier)
    iter = iter + 1
  }
}

# pathway_enrichment 
# input:
#     deg_list
#     FDR_threshold
#     pval_threshold
# JD modifications
# 1) deg_list is compiled list (all cell types and comparisons)
# 2) converted input deg list from mouse to human genes using biomaRt
# 3) deg_list has 4 additional columns: Cell Type, Gene, Comparison, MODULE
#        MODULE combines cell type, comparison, and direction of DEG
# finds the enriched pathways based on DEGs using hypergeometric approach
# note list_of_pathway_databases may need to be modified
# goal is to find the enriched pathways between treated and untreated looking at many combinations
# Output of this in compiled format still has enrichment that only has 1 overlap
pathway_enrichment <- function(deg_list, FDR_threshold=NULL, pval_threshold, celltype_cluster,identifier, pVal=FALSE){
  
  # pVal this is our method of switching off between FDR and pval threshold
  list_of_pathway_databases = list.files(path = "../Resources", pattern = "*.txt", full.names = TRUE)
  
  # only take specific module
  deg_list = deg_list[deg_list$MODULE==celltype_cluster,]
  
  #susbsets off of pval threshold or the adjusted pval
  # if(pVal == TRUE){
  #   deg_list <- deg_list[which(deg_list$p_val < pval_threshold),]
  # }else{
  #   if(!is.null(FDR_threshold)){
  #     deg_list <- deg_list[which(deg_list$p_val_adj < FDR_threshold),]
  #   }
  # }
  
  #remove ribosomal genes
  # rpsLoc = grepl("^RPS", deg_list$GENE)
  # deg_list <- deg_list[!rpsLoc,]
  # rplLoc = grepl("^RPL", deg_list$GENE)
  # deg_list <- deg_list[!rplLoc,]
  # if(length(deg_list[,1]) == 0){
  #   return(NULL)
  # }
  
  deg_list[["gene"]] <- deg_list$GENE
  deg_list[["gene"]] <- toupper(deg_list[["gene"]]) #changes the gene names to upper case for compatibility with database (redundant)
  deg_list[["module"]] <- rep(celltype_cluster, nrow(deg_list))
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
    ifelse(!dir.exists(file.path(paste0("./Csvs","_",identifier,"/pathway_enrichment"))), dir.create(file.path(paste0("./Csvs","_",identifier,"/pathway_enrichment")),recursive = T), FALSE)
    write.table(data_matrix_for_enrichment, file = paste0("./Csvs","_",identifier,"/pathway_enrichment/temp1.",celltype_cluster,".dat"), quote = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE)
    
    record_mat <- read.table(paste0("./Csvs","_",identifier,"/pathway_enrichment/temp1.",celltype_cluster,".dat"))
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
  x[["Combined"]][["Disease_Model"]] <- rep(celltype_cluster, nrow(x[["Combined"]]))
  
  total = data.frame(stringsAsFactors = FALSE)
  total = rbind(x[[1]], x[[2]])
  total = rbind(total, x[[3]])
  total = rbind(total, x[[4]])
  total$ModuleSize = nrow(deg_list)
  ifelse(!dir.exists(file.path(paste0("./consolidated_pathways","_",identifier))), dir.create(file.path(paste0("./consolidated_pathways","_",identifier)),recursive = T), FALSE)
  write.table(total, paste0("./consolidated_pathways","_",identifier,"/", celltype_cluster, ".txt"), row.names = FALSE, quote = FALSE, sep = "\t")
  
  if(pVal==TRUE){
    WriteXLS(x, paste('Csvs/pathway_enrichment/PVAL/',celltype_cluster,'_PVAL.',pval_threshold,'.xls',sep=""),names(x))
  }else{
    WriteXLS(x, paste("./Csvs","_",identifier,"/pathway_enrichment/",celltype_cluster,'_FDR.',FDR_threshold,'.xls',sep=""),names(x))
  }
  
}


# Generates a dot plot showing logFC and FDR of DEGs
# @param DEG_df: output from Seurat FindMarkers()
#          -must have columns named with those inputting into this function as vars
# @param vars: vector describing dimensions of the dot plot (corresponding to column names) 
#          first element is the dimension that is chosen (usually GENE)
#          (first element is on the x axis, second is on the y axis)
#
#          -Examples: c("GENE", "Cell_type") - across different genes on the x axis for different cell types    
#                     c("Cell_type","Comparison") - visualize single gene changing across different cell 
#                                                   types and comparisons
#                     c("GENE", "Comparison") - for a single cell type, visualize genes across comparisons
# @param var_list: [optional] list containing desired variables to appear in the specific orders. 
#                  list names must correspond with the column name
#                  -Example: list("GENE"=c("Apoe","Clu","Trem2"), "Cell_type"=c("Neuron","Astrocyte","Microglia"))
# @param pdf_name: name of output pdf
# @param pdf_height: height of output pdf
# @param pdf_width: width of output pdf
# @param colors: colors to use for logFC gradient (must be three - low, mid, high color)
# 
# (exploratory option) if var_list is not specified, will show top 30 regulated genes (exploratory) 
# (only works if "GENE" is a column name) (assumes no comparison dimension)
#
# Output: ggplot object and pdf of plot 
#
gene_dot_plot <- function(
  DEG_df, 
  vars = c("GENE","Cell_type"), 
  var_list = NULL, 
  pdf_name = "DEG_dotPlot.pdf", 
  colors = c("blue","gray","red"),
  explore.show = "logFC"
){
  DEG_df$neglog10fdr = -log10(DEG_df$p_val_adj)
  DEG_df$neglog10fdr[DEG_df$neglog10fdr==Inf] <- 301
  
  selected = data.frame(stringsAsFactors = FALSE)
  if(!is.null(var_list)){
    temp = data.frame(stringsAsFactors = FALSE)
    for(item in var_list[[1]]){
      temp = rbind(temp, DEG_df[DEG_df[,vars[1]]==item,])
    }
    if(length(var_list)>1){
      for(item2 in var_list[[2]]){
        selected = rbind(selected, temp[temp[,vars[2]]==item2,])
      }
    }
    else{
      selected = temp
    }
    for(v in 1:length(var_list)){
      if(length(unique(selected[,vars[v]]))!=length(var_list[v])){
        selected[,vars[v]] = factor(selected[,vars[v]], levels = unique(selected[,vars[v]]))
      }
      else{
        selected[,vars[v]] = factor(selected[,vars[v]], levels = var_list[v])
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
    geom_point(data=selected[selected$Stringent=="Yes",],pch=21,colour="black", stroke=1) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          #axis.text.x = element_text(angle = 45, margin = margin(t=5), size=14, hjust = 1,vjust = 1),
          axis.text.x = element_text(angle = 90, size=14, hjust = 0,vjust = 1),
          axis.text.y = element_text(size = 14),
          axis.text = element_text(face="bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_text(size=15),
          legend.text = element_text(size=10),
          legend.position = "bottom") +
    labs(size=expression(-log[10]~p(Adj)), color=expression(logFC)) +
    scale_color_gradient2(low = colors[1],  mid=colors[2], high = colors[3],
                          breaks = c(-.5, 0, .5), labels = c(-.5, 0, .5)) +
    scale_size(range = c(3, 12), breaks = c(1.3,25,50,75,100)) +
    scale_shape_manual(values=shapes) +
    guides(shape = guide_legend(override.aes = list(size = 5))) +
    scale_x_discrete(position = "top")
  
  return(D)
}



# this function will set active ident to last element in idents vector
makeNumbersTable <- function(seuratObject, idents){

}

# Example comparisons vector:
# comparisons = c("TGW.v.WTW"="AD effect",
#                 "TGF.v.TGW"="Fructose effect on AD",
#                 "TGDHAF.v.TGF"="DHA effect on AD Fructose",
#                 "TGDHAF.v.TGW"="DHA and Fructose effect on AD",
#                 "TGNRF.v.TGF"="NR effect on AD Fructose",
#                 "TGNRF.v.TGW"="NR and Fructose effect on AD")
#
# In the output Excel file, add navigation links using macros code:
# Sub CreateLinksToAllSheets()
# Dim sh As Worksheet
# Dim cell As Range
# For Each sh In ActiveWorkbook.Worksheets
# If ActiveSheet.Name <> sh.Name Then
# ActiveCell.Hyperlinks.Add Anchor:=Selection, Address:="", SubAddress:= _
# "'" & sh.Name & "'" & "!A1", TextToDisplay:=sh.Name
# ActiveCell.Offset(1, 0).Select
# End If
# Next sh
# End Sub
# had to change "_" to "." - make sure it is always _!!

summarizePathways <- function(path, comparisons=NULL, celltypes, add_direction, name, DEG_df = NULL){
  files = list.files(path)[!grepl("temp1", list.files(path))]
  pathway.enrichment = list()
  if(!is.null(comparisons)){
    for(comp in names(comparisons)){
      set = files[grepl(comp,files)]
      all = data.frame(stringsAsFactors = FALSE)
      for(s in set){
        pathways = read.delim(paste0(path,s), header = TRUE, stringsAsFactors = FALSE)
        # pathways = pathways[order(pathways$nOverlap, decreasing = TRUE),]
        # pathways = pathways[1:5,]
        pathways = pathways[pathways$nOverlap>0,]
        if(add_direction){
          new_overlap = c()
          for(g in pathways$Overlap){
            ct = gsub("_"," ",gsub(paste0("_",comp),"",gsub(".txt","",s)))
            genes = unlist(strsplit(g, split = ","))
            MT_gene_names = c("MT-ND1"="ND1",  "MT-ND2"="ND2",  "MT-CO1"="COX1",  "MT-CO2"="COX2",  "MT-ATP8"="ATP8", 
                              "MT-ATP6"="ATP6", "MT-CO3"="COX3","MT-ND3" ="ND3", "MT-ND4L"="ND4L", "MT-ND4"="ND4",
                              "MT-ND5"="ND5",  "MT-ND6"="ND6",  "MT-CYB"="CYB")
            # change to proper MT- name in DEG_df
            for(gene in 1:length(genes)){
              if(sum(MT_gene_names==genes[gene])>0){
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
            UP_genes = c("UP:")
            DOWN_genes = c("DOWN:")
            for(gene in genes){
              FC = DEG_df$avg_logFC[DEG_df$Comparison==comp &
                                                DEG_df$Cell_type==ct &
                                                DEG_df$GENE==gene]
              FC = FC[!is.na(FC)]
              if(length(FC)>1){
                dup_genes = append(dup_genes, gene)
                cat("dup\n")
              }
              if(FC>0){
                UP_genes = append(UP_genes, paste0(gene,","))
              }
              else{
                DOWN_genes = append(DOWN_genes, paste0(gene,","))
              }
            }
            UP_genes = append(UP_genes, DOWN_genes)
            detailed_overlap = concatenate(UP_genes, mysep = "")
            new_overlap = append(new_overlap, detailed_overlap)
          }
          pathways$Overlap = new_overlap
        }
        pathways$Module = gsub(paste0(comp,"_"), "", pathways$Module)
        all = rbind(all, pathways)
      }
      pathway.enrichment[[comparisons[comp]]] = all
    }
  }
  # get cell types
  
  # if(!add_direction){
  #   cellType_effect = c()
  #   top_pathways = c()
  #   for(c in celltypes){
  #     cellType_effect = append(cellType_effect, paste0(c,".UP"))
  #     cellType_effect = append(cellType_effect, paste0(c,".DOWN"))
  #   }
  #   celltypes = cellType_effect
  # }
  # for(effect in celltypes){
  #   set = files[grepl(effect,files)]
  #   all = data.frame(stringsAsFactors = FALSE)
  #   for(s in names(comparisons)){
  #     toSkip = c()
  #     pathway2 = paste0(path,set[grep(s,set)])
  #     if(pathway2==path) next
  #     skipOk = FALSE
  #     for(t in toSkip){
  #       if(sum(pathway2==t)>0) skipOk = TRUE
  #     }
  #     if(skipOk) next
  #     else if(length(pathway2)>1){
  #       for(p in pathway2){
  #         pathways = read.delim(p, header = TRUE, stringsAsFactors = FALSE)
  #         pathways = pathways[order(pathways$FDR),]
  #         # pathways = pathways[1:5,]
  #         top_pathways = append(top_pathways, pathways$Pathway)
  #         pathways = pathways[pathways$FDR<0.05,]
  #         pathways$Module = gsub(paste0("_",effect), "", pathways$Module)
  #         all = rbind(all, pathways)
  #       }
  #       toSkip = append(toSkip, pathway2)
  #     }
  #     else{
  #       pathways = read.delim(paste0(path,set[grep(s,set)]), header = TRUE, stringsAsFactors = FALSE)
  #       pathways = pathways[order(pathways$FDR),]
  #       # pathways = pathways[1:5,]
  #       top_pathways = append(top_pathways, pathways$Pathway)
  #       pathways = pathways[pathways$FDR<0.05,]
  #       pathways$Module = gsub(paste0("_",effect), "", pathways$Module)
  #       all = rbind(all, pathways)
  #     }
  #   }
  #   pathway.enrichment[[effect]] = all
  # }
  return(pathway.enrichment)
  WriteXLS(pathway.enrichment, name, names(pathway.enrichment))
}

# if one comparison, then do by cell type in long format. if there are many comparisons, then do by comparison
# c("TGW.v.WTW", "TGF.v.TGW", "TGDHAF.v.TGF","TGDHAF.v.TGW","TGNRF.v.TGF","TGNRF.v.TGW")
summarizeDEGslongFormat <- function(DEG_df, comparisons = NULL){
  
  if(!is.null(comparisons)){
    melted_DEGs = list()
    for(module in unique(DEG_df$MODULE_ct)){
      all_sig_DEGs_subset = unmelt(df = DEG_df[grepl(module,DEG_df$MODULE),], id_column = "MODULE", value_column = "GENE")
      colnames(all_sig_DEGs_subset) <- gsub(paste0("_",module), "", colnames(all_sig_DEGs_subset))
      all_sig_DEGs_subset = all_sig_DEGs_subset[,comparisons]
      melted_DEGs[[module]] = all_sig_DEGs_subset
    }
    WriteXLS(melted_DEGs, "all_melted_DEGs.xls", names(melted_DEGs))
  }
  else{
    for(module in unique(DEG_df$Cell_type)){
      all_sig_DEGs_subset = unmelt(df = DEG_df[grepl(module,DEG_df$MODULE),], id_column = "MODULE", value_column = "GENE")
      write.table(all_sig_DEGs_subset, row.names = FALSE, quote = FALSE, sep = "\t")
    }
  }
}

# if only one comparison, don't have to set comparisons if not adding mean normalized data
# if adding mean normalized data, then must have both comparisons and conditions
summarizeTopDegs <- function(DEG_df, comparisons = NULL, 
                             celltypes =unique(DEG_df$Cell_type), 
                             output_name = "topDEGs.xls",
                             meanNormDataDf = NULL,
                             conditions = list("TGW.v.WTW"=c("TGW","WTW")), convertToHuman=TRUE,
                             GWAS_info = NULL){
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
  if(!is.null(comparisons)){
    for(comp in names(comparisons)){
      subset = DEG_df[DEG_df$Comparison==comp,]
      celltype_group = c()
      logFCs = c()
      nCelltypes = c()
      for(DEG in unique(subset$GENE)){
        if(sum(subset$GENE==DEG)>0){
          #CHANGED TO MODULE_ct_tissue FOR CURRENT PURPOSES
          celltype_group = append(celltype_group,
                                  concatenate(gsub(paste0(comp,"_"),"",subset$MODULE_ct_tissue[subset$GENE==DEG]),
                                              mysep = ", "))
          logFCs = append(logFCs, concatenate(gsub(paste0(comp,"_"),"",subset$avg_logFC[subset$GENE==DEG]), 
                                              mysep = ", "))
          nCelltypes = append(nCelltypes, sum(subset$GENE==DEG))
        }
      }
      top_DEGs = data.frame("GENE" = unique(subset$GENE),
                            "Cell Types"=celltype_group,
                            "logFCs" = logFCs,
                            "nCellTypes"=nCelltypes, stringsAsFactors = FALSE)
      if(!is.null(meanNormDataDf)){
        for(cond in  conditions[[comp]]){
          all_avgs = c()
          for(DEG in unique(subset$GENE)){
            modules = paste0(subset$Ct_tissue[subset$GENE==DEG], ".", cond)
            avgs = c()
            for(mod in modules){
              avgs = append(avgs, meanNormDataDf[meanNormDataDf$GENE==DEG,mod])
            }
            all_avgs = append(all_avgs, concatenate(avgs, mysep = ", "))
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
        top_DEGs$`AD GWAS Pubmed IDs` = addGWASstudy(top_DEGs)
      }
      
      # add info
      top_DEGs$INFO = annotate_gene(top_DEGs)
      colnames(top_DEGs) <- gsub("TGW", "5XFAD", colnames(top_DEGs))
      colnames(top_DEGs) <- gsub("WTW", "WT", colnames(top_DEGs))
      topDEGs[[comparisons[comp]]] = top_DEGs
    }
    
    # by cell type effect
    # rearrange data frame
    ordered = data.frame(stringsAsFactors = FALSE)
    for(comp in names(comparisons)){
      ordered = rbind(ordered, DEG_df[DEG_df$Comparison==comp,])
    }
    
    # adding the periods before and after because of multiple cell type namings - Endothelium and Endothelium 2, for ex.
    for(celltype in paste0(".",celltypes,".")){
      cat("*********Now summarizing ", celltype, "\n")
      subset = ordered[grep(celltype, ordered$MODULE, fixed = TRUE),]
      if(nrow(subset)==0) next
      celltype_group = c()
      nCelltypes = c()
      logFCs = c()
      for(DEG in unique(subset$GENE)){
        if(sum(subset$GENE==DEG)>0){
          celltype_group = append(celltype_group, 
                                  concatenate(gsub(paste0(celltype,"_"),"",subset$MODULE[subset$GENE==DEG]), 
                                              mysep = ", "))
          logFCs = append(logFCs, concatenate(gsub(paste0(celltype,"_"),"",subset$avg_logFC[subset$GENE==DEG]), 
                                              mysep = ", "))
          nCelltypes = append(nCelltypes, sum(subset$GENE==DEG))
        }
      }
      top_DEGs = data.frame("GENE" = unique(subset$GENE),
                            "Comparisons"=celltype_group,
                            "logFCs" = logFCs,
                            "nComparisons"=nCelltypes, stringsAsFactors = FALSE)
      if(!is.null(meanNormDataDf)){
        for(cond in conditions[[comp]]){
          all_avgs = c()
          for(DEG in unique(subset$GENE)){
            modules = paste0(subset$Ct_tissue[subset$GENE==DEG], ".", cond)
            avgs = meanNormDataDf[meanNormDataDf$GENE==DEG,colnames(meanNormDataDf) %in% modules]
            all_avgs = append(all_avgs, concatenate(avgs, mysep = ", "))
          }
          top_DEGs[,paste0(cond, "_means")] = all_avgs
        }
      }
      
      top_DEGs = top_DEGs[order(top_DEGs$nComparisons, decreasing = TRUE),]
      
      if(convertToHuman){
        top_DEGs$HUMAN <- mouse_human_reference$HUMAN[match(x = top_DEGs$GENE, table = mouse_human_reference$MOUSE)]
      }
      
      if(!is.null(GWAS_info)){
        # add GWAS study accession
        top_DEGs$`AD GWAS Pubmed IDs` = addGWASstudy(top_DEGs)
      }
      
      # add info
      top_DEGs$INFO = annotate_gene(top_DEGs)
      celltype = gsub(".", "",celltype, fixed = TRUE)
      
      if(length(comparisons)==1){
        top_DEGs$Comparisons = gsub(paste0(comp,"."), "", top_DEGs$Comparisons)
        top_DEGs$Comparisons = gsub(paste0(celltype,"."), "", top_DEGs$Comparisons)
      }
      
      colnames(top_DEGs) <- gsub("TGW", "5XFAD", colnames(top_DEGs))
      colnames(top_DEGs) <- gsub("WTW", "WT", colnames(top_DEGs))
      topDEGs[[celltype]] = top_DEGs
    }
    # return(topDEGs)
    WriteXLS(topDEGs, output_name, gsub("\\.","",names(topDEGs)))
    return(topDEGs)
  }
  else{
    for(DEG in unique(DEG_df$GENE)){
      if(sum(DEG_df$GENE==DEG)>0){
        celltype_group = append(celltype_group, 
                                concatenate(gsub(paste0(celltype,"_"),"",DEG_df$MODULE[DEG_df$GENE==DEG]), 
                                            mysep = ", "))
        nCelltypes = append(nCelltypes, sum(DEG_df$GENE==DEG))
      }
    }
    top_DEGs = data.frame("GENE" = unique(DEG_df$GENE),
                          "Cell Types"=celltype_group,
                          "nCellTypes"=nCelltypes, stringsAsFactors = FALSE)
    top_DEGs = top_DEGs[order(-top_DEGs$nCellTypes),]
    # add human gene
    top_DEGs_human <- convertDfGeneColumnMouseHuman(df = top_DEGs, toSpecies = "human",forPathway = "TRUE")
    top_DEGs$HUMAN <- top_DEGs_human$GENE
    
    # add GWAS study accession
    top_DEGs$`GWAS Catalog Study Accession` = addGWASstudy(top_DEGs)
    
    # add info
    top_DEGs$INFO = annotate_gene(top_DEGs)
    
    write.table(top_DEGs, "top_DEGs.txt", row.names = FALSE, quote = FALSE, sep = "\t")
  }
}

summarizeTopDegs <- function(DEG_df, comparisons = NULL, 
                             celltypes =unique(DEG_df$Cell_type), 
                             output_name = "topDEGs.xls",
                             meanNormDataDf = NULL,
                             conditions = list("TGW.v.WTW"=c("TGW","WTW")), convertToHuman=TRUE,
                             gene_info = NULL,
                             GWAS_info = NULL){
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
  if(!is.null(comparisons)){
    for(comp in names(comparisons)){
      subset = DEG_df[DEG_df$Comparison==comp,]
      celltype_group = c()
      logFCs = c()
      nCelltypes = c()
      for(DEG in unique(subset$GENE)){
        if(sum(subset$GENE==DEG)>0){
          #CHANGED TO MODULE_ct_tissue FOR CURRENT PURPOSES
          celltype_group = append(celltype_group,
                                  concatenate(gsub(paste0(comp,"_"),"",subset$MODULE_ct_tissue[subset$GENE==DEG]),
                                              mysep = ", "))
          logFCs = append(logFCs, concatenate(gsub(paste0(comp,"_"),"",subset$avg_logFC[subset$GENE==DEG]), 
                                              mysep = ", "))
          nCelltypes = append(nCelltypes, sum(subset$GENE==DEG))
        }
      }
      top_DEGs = data.frame("GENE" = unique(subset$GENE),
                            "Cell Types"=celltype_group,
                            "logFCs" = logFCs,
                            "nCellTypes"=nCelltypes, stringsAsFactors = FALSE)
      if(!is.null(meanNormDataDf)){
        for(cond in  conditions[[comp]]){
          all_avgs = c()
          for(DEG in unique(subset$GENE)){
            modules = paste0(subset$Ct_tissue[subset$GENE==DEG], ".", cond)
            avgs = c()
            for(mod in modules){
              avgs = append(avgs, meanNormDataDf[meanNormDataDf$GENE==DEG,mod])
            }
            all_avgs = append(all_avgs, concatenate(avgs, mysep = ", "))
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
        top_DEGs$`AD GWAS Pubmed IDs` = addGWASstudy(top_DEGs)
      }
      
      if(!is.null(gene_info)){
        # add info
        top_DEGs$INFO = annotate_gene(top_DEGs)
      }
      colnames(top_DEGs) <- gsub("TGW", "5XFAD", colnames(top_DEGs))
      colnames(top_DEGs) <- gsub("WTW", "WT", colnames(top_DEGs))
      topDEGs[[comparisons[comp]]] = top_DEGs
    }
    
    # by cell type effect
    # rearrange data frame
    ordered = data.frame(stringsAsFactors = FALSE)
    for(comp in names(comparisons)){
      ordered = rbind(ordered, DEG_df[DEG_df$Comparison==comp,])
    }
    
    
      
      if(length(comparisons)==1){
        top_DEGs$Comparisons = gsub(paste0(comp,"."), "", top_DEGs$Comparisons)
        top_DEGs$Comparisons = gsub(paste0(celltype,"."), "", top_DEGs$Comparisons)
      }
  }
  
  # adding the periods before and after because of multiple cell type namings - Endothelium and Endothelium 2, for ex.
  for(celltype in paste0(".",celltypes,".")){
    cat("*********Now summarizing ", celltype, "\n")
    subset = ordered[grep(celltype, ordered$MODULE, fixed = TRUE),]
    if(nrow(subset)==0) next
    celltype_group = c()
    nCelltypes = c()
    logFCs = c()
    for(DEG in unique(subset$GENE)){
      if(sum(subset$GENE==DEG)>0){
        celltype_group = append(celltype_group, 
                                concatenate(gsub(paste0(celltype,"_"),"",subset$MODULE[subset$GENE==DEG]), 
                                            mysep = ", "))
        logFCs = append(logFCs, concatenate(gsub(paste0(celltype,"_"),"",subset$avg_logFC[subset$GENE==DEG]), 
                                            mysep = ", "))
        nCelltypes = append(nCelltypes, sum(subset$GENE==DEG))
      }
    }
    top_DEGs = data.frame("GENE" = unique(subset$GENE),
                          "Comparisons"=celltype_group,
                          "logFCs" = logFCs,
                          "nComparisons"=nCelltypes, stringsAsFactors = FALSE)
    if(!is.null(meanNormDataDf)){
      for(cond in conditions[[comp]]){
        all_avgs = c()
        for(DEG in unique(subset$GENE)){
          modules = paste0(subset$Ct_tissue[subset$GENE==DEG], ".", cond)
          avgs = meanNormDataDf[meanNormDataDf$GENE==DEG,colnames(meanNormDataDf) %in% modules]
          all_avgs = append(all_avgs, concatenate(avgs, mysep = ", "))
        }
        top_DEGs[,paste0(cond, "_means")] = all_avgs
      }
    }
    
    top_DEGs = top_DEGs[order(top_DEGs$nComparisons, decreasing = TRUE),]
    
    if(convertToHuman){
      top_DEGs$HUMAN <- mouse_human_reference$HUMAN[match(x = top_DEGs$GENE, table = mouse_human_reference$MOUSE)]
    }
    
    if(!is.null(GWAS_info)){
      # add GWAS study accession
      top_DEGs$`AD GWAS Pubmed IDs` = addGWASstudy(top_DEGs)
    }
    
    if(!is.null(gene_info)){
      # add info
      top_DEGs$INFO = annotate_gene(top_DEGs)
    }
    celltype = gsub(".", "",celltype, fixed = TRUE)
    
    topDEGs[[celltype]] = top_DEGs
  }
  
  WriteXLS(topDEGs, output_name, gsub("\\.","",names(topDEGs)))
  return(topDEGs)
  
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

addGWASstudy <- function(vector, GWAS_info){
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
  for(gene in vector$HUMAN){
    if(sum(GWAS_genes$GENE==gene)>0){
      GWAS_hit_info = append(GWAS_hit_info, GWAS_genes$PUBMED_ID[GWAS_genes$GENE==gene])
    }
    else{
      GWAS_hit_info = append(GWAS_hit_info, "")
    }
  }
  return(GWAS_hit_info)
}


annotate_gene <- function(vector, protein_names="10090.protein.info.v11.0.txt"){
  protein_names = read.delim(protein_names, stringsAsFactors = FALSE)
  info = c()
  for(m in vector$GENE){
    if(length(protein_names$annotation[protein_names$preferred_name==m])==0){
      if(length(protein_names$annotation[protein_names$preferred_name==gsub("[[:digit:]]","",m)])==0){
        info = append(info, "Not annotated")
        cat(m, "\n")
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

# pathway_dir: output of pathways in COMPILED format (from pathway enrichment function, "consolidated_pathways" directory)
# direction:
#   If FALSE, pathway modules do not have directionality information, intensity of tile color correlate to FDR significance
#   If TRUE, pathway modules do have directionality information (in _UP and _DOWN format) and tile color will be red/blue accordingly
#   and intensity of color will correlate to FDR significance
# cellTypes to display, by default does all cellTypes.
# topPathways - if TRUE, then uses top most significant pathways of the cell type (user defined) and not the pathway_list

pathway_heatmap <- function(pathway_dir, 
                            direction = TRUE,
                            trim=FALSE, 
                            cellTypes = NULL, 
                            pathway_list, 
                            topPathways = FALSE, topPathways_num = 30,
                            comparison_levels = NULL,
                            celltype_levels = NULL,
                            Min_Max = NULL,
                            pdf_width = 25,
                            pdf_height = 25,
                            name = "pathway_heatmap"){
  # compile all pathways into one df
  pathway_df = data.frame(stringsAsFactors = FALSE)
  for(file in list.files(pathway_dir)){
    temp = read.delim(paste0(pathway_dir,file), stringsAsFactors = FALSE, header = TRUE)
    pathway_df = rbind(pathway_df, temp)
  }
  pathway_df$`negLog10FDR` = -log10(pathway_df$FDR)
  
  if(!is.null(cellTypes)){
    subset = data.frame(stringsAsFactors = FALSE)
    for(c in cellTypes){
      temp = pathway_df[grep(c,pathway_df$Module),]
      subset = rbind(subset,temp)
    }
    pathway_df = subset
    rm(subset)
  }
  
  # select subset of pathways
  if(!topPathways){
    selected_pathways = data.frame(stringsAsFactors = FALSE)
    for(p in pathway_list){
      temp = pathway_df[pathway_df$Pathway==p,]
      selected_pathways = rbind(selected_pathways, temp)
    }
  }
  else{
    selected_pathways = data.frame(stringsAsFactors = FALSE)
    pathway_df = pathway_df[order(pathway_df$FDR),]
    top_pathways = pathway_df$Pathway[1:topPathways_num]
    for(p in top_pathways){
      temp = pathway_df[pathway_df$Pathway==p,]
      selected_pathways = rbind(selected_pathways, temp)
    }
  }
  rm(pathway_df)
  
  
  if(direction){
    # for trimming
    selected_pathways$MODULEv2 = gsub("_DOWN","",gsub("_UP","",selected_pathways$Module))
    
    # add pseudo directionality - change sign of the FDR
    for(mod in 1:length(selected_pathways$Module)){
      if(grepl("_DOWN",selected_pathways$Module[mod])){
        selected_pathways$negLog10FDR[mod] = -selected_pathways$negLog10FDR[mod]
      }
      else{
        next
      }
    }
  }
  
  # desired ordering that comparisons appear if both cell type and comparison order desired
  if(!is.null(comparison_levels) & !is.null(celltype_levels)){
    order = c()
    for(c in celltypes){
      for(comp in (comparison_levels)){
        module = unique(selected_pathways$Module[grep(c, selected_pathways$Module)])
        module = module[grep(comp, module)]
        order = append(order, module)
      }
    }
  }
  # ordering by only comparison
  if(!is.null(comparison_levels) & is.null(celltype_levels)){
    order = c()
    for(comp in (comparison_levels)){
      module = unique(selected_pathways$Module[grep(c, selected_pathways$Module)])
      module = module[grep(comp, module)]
      order = append(order, module)
    }
  }
  selected_pathways = selected_pathways[!duplicated(selected_pathways),]
  if(trim){
    to_trim = c()
    for(path in unique(selected_pathways$Pathway)){
      for(mod in unique(selected_pathways$MODULEv2)){
        positions = which(selected_pathways$Pathway==path & 
                            selected_pathways$MODULEv2==mod)
        position_UP = positions[grepl("UP",selected_pathways$Module[selected_pathways$MODULEv2==mod & selected_pathways$Pathway==path])]
        if(length(position_UP)>1) cat(path, mod)
                              
        position_DOWN = positions[grepl("DOWN",selected_pathways$Module[selected_pathways$MODULEv2==mod & selected_pathways$Pathway==path])]
        if(selected_pathways$FDR[position_UP]>selected_pathways$FDR[position_DOWN]){
          to_trim = append(to_trim, position_UP)
        }
        else if(selected_pathways$FDR[position_UP]<selected_pathways$FDR[position_DOWN]){
          to_trim = append(to_trim, position_DOWN)
        }
        else{
          to_trim = append(to_trim, position_DOWN)
        }
      }
    }
    selected_pathways = selected_pathways[-to_trim,]
    order = gsub("_DOWN","",order)
    order = gsub("_UP","",order)
    order = order[!duplicated(order)]
  }
  
  if(topPathways){
    selected_pathways$Pathway = factor(selected_pathways$Pathway, levels = rev(unique(selected_pathways$Pathway)))
  }
  if(!is.null(comparison_levels) | !is.null(celltype_levels)){
    if(!trim){
      selected_pathways$Module = factor(selected_pathways$Module, levels = order)
    }
    else{
      selected_pathways$MODULEv2 = factor(selected_pathways$MODULEv2, levels = order)
    }
  }
  if(!is.null(Min_Max)){
    selected_pathways$negLog10FDR = MinMax(selected_pathways$negLog10FDR, min = -Min_Max, max = Min_Max)
  }
  
  if(trim & direction){
    heat <- ggplot(selected_pathways,aes(x=MODULEv2,y = Pathway, fill=negLog10FDR)) + geom_tile() + 
      scale_fill_gradientn(colors = c("blue", "white", "red")) + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, size = 20, vjust = 0, hjust = 0, 
                                       margin = margin(t=0,b=0,r=0,l=0)), 
            axis.text.y = element_text(size=15),
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.title.x = element_blank()) + 
      labs(fill=expression(-log[10]~FDR)) +
      scale_x_discrete(position = "top") 
  }
  if(!trim & direction){
    cat("Using Module")
    heat <- ggplot(selected_pathways,aes(x=Module,y = Pathway, fill=negLog10FDR,  label=paste0(nOverlap,"/",ModuleSize))) +
      geom_tile() + 
      geom_fit_text() +
      scale_fill_gradientn(colors = c("blue", "white", "red")) + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, size = 20, vjust = 0, hjust = 0, 
                                       margin = margin(t=0,b=0,r=0,l=0)), 
            axis.text.y = element_text(size=15),
            panel.border = element_blank(), 
            axis.title.x = element_blank()) + 
      labs(fill=expression(-log[10]~FDR)) +
      scale_x_discrete(position = "top")
  }
  if(!direction){
    selected_pathways$Module = gsub(paste0(cellTypes,"_"), "",selected_pathways$Module)
    selected_pathways$Module = factor(selected_pathways$Module, levels = comparison_levels)
    heat <- ggplot(selected_pathways,aes(x=Module,y = Pathway, fill=negLog10FDR)) + geom_tile() + 
      scale_fill_gradientn(colors = c("white", "red")) + theme_bw() +
      theme(axis.text.x = element_text(angle = 60, size = 10, vjust = 0, hjust = 0,
                                       margin = margin(t=0,b=0,r=0,l=0)), 
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.title.x = element_blank()) + 
      labs(fill=expression(-log[10]~FDR)) +
      scale_x_discrete(position = "top")
  }
  
  pdf(paste0(name,".pdf"), width = pdf_width, height = pdf_height)
  print(heat)
  dev.off()
}


pathway_heatmapv2 <- function(pathway_dir, 
                            direction = TRUE,
                            trim=FALSE, 
                            cellTypes = NULL, 
                            pathway_list, 
                            topPathways = FALSE, topPathways_num = 30,
                            comparison_levels = NULL,
                            celltype_levels = NULL,
                            Min_Max = NULL,
                            pdf_width = 25,
                            pdf_height = 25,
                            name = "pathway_heatmap"){
  # compile all pathways into one df
  pathway_df = data.frame(stringsAsFactors = FALSE)
  for(file in list.files(pathway_dir)){
    temp = read.delim(paste0(pathway_dir,file), stringsAsFactors = FALSE, header = TRUE)
    pathway_df = rbind(pathway_df, temp)
  }
  pathway_df$`negLog10FDR` = -log10(pathway_df$FDR)
  
  if(!is.null(cellTypes)){
    subset = data.frame(stringsAsFactors = FALSE)
    for(c in cellTypes){
      temp = pathway_df[grep(c,pathway_df$Module, fixed = TRUE),]
      subset = rbind(subset,temp)
    }
    pathway_df = subset
    rm(subset)
  }
  
  # select subset of pathways
  if(!topPathways){
    selected_pathways = data.frame(stringsAsFactors = FALSE)
    for(p in pathway_list){
      temp = pathway_df[pathway_df$Pathway==p,]
      selected_pathways = rbind(selected_pathways, temp)
    }
  }
  else{
    selected_pathways = data.frame(stringsAsFactors = FALSE)
    pathway_df = pathway_df[order(pathway_df$nOverlap, decreasing = TRUE),]
    top_pathways = unique(pathway_df$Pathway)[1:topPathways_num]
    for(p in top_pathways){
      temp = pathway_df[pathway_df$Pathway==p,]
      selected_pathways = rbind(selected_pathways, temp)
    }
  }
  rm(pathway_df)
  
  
  if(direction){
    # for trimming
    selected_pathways$MODULEv2 = gsub("DOWN","",gsub("UP","",selected_pathways$Module))
    
    # add pseudo directionality - change sign of the FDR
    for(mod in 1:length(selected_pathways$Module)){
      if(grepl("DOWN",selected_pathways$Module[mod])){
        selected_pathways$negLog10FDR[mod] = -selected_pathways$negLog10FDR[mod]
        selected_pathways$nOverlap[mod] = -selected_pathways$nOverlap[mod]
      }
      else{
        next
      }
    }
  }
  
  # desired ordering that comparisons appear if both cell type and comparison order desired
  if(!is.null(comparison_levels) & !is.null(celltype_levels)){
    order = c()
    for(c in celltypes){
      for(comp in (comparison_levels)){
        module = unique(selected_pathways$Module[grep(c, selected_pathways$Module)])
        module = module[grep(comp, module)]
        order = append(order, module)
      }
    }
  }
  # ordering by only comparison
  if(!is.null(comparison_levels) & is.null(celltype_levels)){
    order = c()
    for(comp in (comparison_levels)){
      module = unique(selected_pathways$Module[grep(c, selected_pathways$Module)])
      module = module[grep(comp, module)]
      order = append(order, module)
    }
  }
  selected_pathways = selected_pathways[!duplicated(selected_pathways),]
  if(trim){
    to_trim = c()
    for(path in unique(selected_pathways$Pathway)){
      for(mod in unique(selected_pathways$MODULEv2)){
        positions = which(selected_pathways$Pathway==path & 
                            selected_pathways$MODULEv2==mod)
        position_UP = positions[grepl("UP",selected_pathways$Module[selected_pathways$MODULEv2==mod & selected_pathways$Pathway==path])]
        if(length(position_UP)>1) cat(path, mod)
        
        position_DOWN = positions[grepl("DOWN",selected_pathways$Module[selected_pathways$MODULEv2==mod & selected_pathways$Pathway==path])]
        if(selected_pathways$nOverlap[position_UP]<selected_pathways$nOverlap[position_DOWN]){
          to_trim = append(to_trim, position_UP)
        }
        else if(selected_pathways$nOverlap[position_UP]>selected_pathways$nOverlap[position_DOWN]){
          to_trim = append(to_trim, position_DOWN)
        }
        else{
          to_trim = append(to_trim, position_DOWN)
        }
      }
    }
    selected_pathways = selected_pathways[-to_trim,]
    order = gsub("_DOWN","",order)
    order = gsub("_UP","",order)
    order = order[!duplicated(order)]
  }
  
  if(topPathways){
    selected_pathways$Pathway = factor(selected_pathways$Pathway, levels = rev(unique(selected_pathways$Pathway)))
  }
  if(!is.null(comparison_levels) | !is.null(celltype_levels)){
    if(!trim){
      selected_pathways$Module = factor(selected_pathways$Module, levels = order)
    }
    else{
      selected_pathways$MODULEv2 = factor(selected_pathways$MODULEv2, levels = order)
    }
  }
  if(!is.null(Min_Max)){
    selected_pathways$negLog10FDR = MinMax(selected_pathways$negLog10FDR, min = -Min_Max, max = Min_Max)
  }
  
  if(trim & direction){
    heat <- ggplot(selected_pathways,aes(x=MODULEv2,y = Pathway, fill=negLog10FDR)) + geom_tile() + 
      scale_fill_gradientn(colors = c("blue", "white", "red")) + theme_bw() +
      theme(axis.text.x = element_text(angle = 90, size = 20, vjust = 0, hjust = 0, 
                                       margin = margin(t=0,b=0,r=0,l=0)), 
            axis.text.y = element_text(size=15),
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.title.x = element_blank()) + 
      labs(fill=expression(-log[10]~FDR)) +
      scale_x_discrete(position = "top") 
  }
  if(!trim & direction){
    write.table(selected_pathways, paste0(name,".txt"), row.names = FALSE, sep = "\t", quote = FALSE)
    cat("Using Module")
    #selected_pathways$Pathway <- makePathwayPretty(selected_pathways$Pathway)
    selected_pathways$Module <-gsub(cellTypes,"",selected_pathways$Module)
    selected_pathways$Module <-gsub("MCI.v.CON","MCI\n",selected_pathways$Module)
    selected_pathways$Module <-gsub("Yes.v.No","Supplement\n",selected_pathways$Module)
    heat <- ggplot(selected_pathways,aes(x=Module,y = Pathway, fill=nOverlap, label=paste0(nOverlap,"/",ModuleSize))) + 
      geom_tile() + 
      geom_fit_text() +
      scale_fill_gradient2(low="blue",high="red") + theme_bw() +
      theme(axis.text.x = element_text(angle = 0, size = 15, vjust = 0, hjust = .5, 
                                       margin = margin(t=0,b=0,r=0,l=0)), 
            axis.text.y = element_text(size=15),
            panel.border = element_blank(), 
            axis.title.x = element_blank()) + 
      labs(fill=expression(-log[10]~FDR)) +
      scale_x_discrete(position = "top") +
      ggtitle(gsub(".","",cellTypes)) +
      theme(legend.position="none")
    
    
  }
  if(!direction){
    selected_pathways$Module = gsub(paste0(cellTypes,"_"), "",selected_pathways$Module)
    selected_pathways$Module = factor(selected_pathways$Module, levels = comparison_levels)
    heat <- ggplot(selected_pathways,aes(x=Module,y = Pathway, fill=negLog10FDR, label=paste0(nOverlap,"/",ModuleSize))) + geom_tile() + 
      geom_fit_text() +
      scale_fill_gradientn(colors = c("white", "red")) + theme_bw() +
      theme(axis.text.x = element_text(angle = 60, size = 10, vjust = 0, hjust = 0,
                                       margin = margin(t=0,b=0,r=0,l=0)), 
            panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.title.x = element_blank()) + 
      labs(fill=expression(-log[10]~FDR)) +
      scale_x_discrete(position = "top") 
  }
  
  pdf(paste0(name,".pdf"), width = pdf_width, height = pdf_height)
  print(heat)
  dev.off()
}


# after getting pathways, add logFC information
# genes have to be compatible (must convert mouse to human if needed)
# makes many assumptions on string structure of file names, data frames (of cell types, etc.)
add_logFC_info_pathways <- function(pathway_dir, DEG_df, output_Dir){
  pathway_df = data.frame(stringsAsFactors = FALSE)
  for(file in list.files(pathway_dir)){
    pathway_df = read.delim(paste0(pathway_dir,file), stringsAsFactors = FALSE, header = TRUE)
    celltype = c()
    comparison = c()
    for(iter in 1:nrow(pathway_df)){
      module_vect = unlist(strsplit(pathway_df$Module[iter],split = "_"))
      celltype[iter] = concatenate(module_vect[1:(length(module_vect)-1)], mysep = "_")
      comparison[iter] = module_vect[length(module_vect)]
    }
    pathway_df$Cell_type = celltype
    pathway_df$Comparison = comparison
    
    # add up avg_LogFC
    
    pathway_df$Cell_type <- gsub("_", " ", pathway_df$Cell_type)
    MT_gene_names = c("MT-ND1"="ND1",  "MT-ND2"="ND2",  "MT-CO1"="COX1",  "MT-CO2"="COX2",  "MT-ATP8"="ATP8", 
                      "MT-ATP6"="ATP6", "MT-CO3"="COX3","MT-ND3" ="ND3", "MT-ND4L"="ND4L", "MT-ND4"="ND4",
                      "MT-ND5"="ND5",  "MT-ND6"="ND6",  "MT-CYB"="CYB")
    
    avg_logFC = c()
    for(row in 1:nrow(pathway_df)){
      if(pathway_df$nOverlap[row]>0){
        ct_logFC_all = DEG_df[DEG_df$Cell_type==pathway_df$Cell_type[row] & DEG_df$Comparison==pathway_df$Comparison[row],]
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
    pathway_df$Comparison <- gsub("Yes.v.No","MCI_Sup.v.MCI_No_Sup", pathway_df$Comparison)
    pathway_df$Module <- gsub("Yes.v.No","MCI_Sup.v.MCI_No_Sup", pathway_df$Module)
    pathway_df$CombinedLogFC = avg_logFC
    pathway_df = pathway_df[order(abs(pathway_df$CombinedLogFC), decreasing = TRUE),]
    file <- gsub("Yes.v.No","MCI_Sup.v.MCI_No_Sup", file)
    write.table(pathway_df, paste0(output_Dir,file), quote = FALSE, sep = "\t", row.names = FALSE)
  }
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
  vector <- gsub("TNFA_SIGNALING_VIA_NFKB","TNFa signaling via NFkB", vector)
  vector <- gsub("CHOLESTEROL_HOMEOSTASIS","Cholesterol homeostasis", vector)
  vector <- gsub("GLYCOLYSIS","Glycolysis", vector)
  vector <- gsub("AXON_GUIDANCE","Axon guidance", vector)
  vector <- gsub("ION_CHANNEL_TRANSPORT","Ion channel transport", vector)
  vector <- gsub("ESTROGEN_RESPONSE_EARLY","Early estrogen response", vector)
  vector <- gsub("PPAR_SIGNALING_PATHWAY","PPAR signaling", vector)
  vector <- gsub("STEROID_BIOSYNTHESIS","Steroid biosynthesis", vector)
  vector <- gsub("OXIDATIVE_PHOSPHORYLATION","Oxidative phosphorylation", vector)
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

getInfoHumam <- function(gene){
  protein_names <- read.delim("../9606.protein.info.v11.0.txt", header = TRUE, 
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
# gsl
# BiocManager::install(c("princurve", "NMF", 
#                        "phylobase", "howmany", "locfdr", "kernlab", "copula", "glmnet", "softImpute"))

