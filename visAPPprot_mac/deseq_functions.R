library(DESeq2)
library(tidyverse)
library(ggplot2)
library(VennDiagram)
library(ggfortify)
library(fgsea)
library(data.table)
library(Rcpp)
library(apeglm)
library(rjson)



# from pheatmap source code
cluster_mat = function(mat, distance, method){
    if(!(method %in% c("ward.D", "ward.D2", "ward", "single", "complete", "average", "mcquitty", "median", "centroid"))){
        stop("clustering method has to one form the list: 'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.")
    }
    if(!(distance[1] %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) & class(distance) != "dist"){
        stop("distance has to be a dissimilarity structure as produced by dist or one measure  form the list: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'")
    }
    if(distance[1] == "correlation"){
        d = as.dist(1 - cor(t(mat)))
    }
    else{
        if(class(distance) == "dist"){
            d = distance
        }
        else{
            d = dist(mat, method = distance)
        }
    }
    
    return(hclust(d, method = method))
}


rawdata <- function(data_dir, cur_dataset, exp_mat, levels, lfcshrink_type, cluster_rows, cluster_cols, folder_name, heatmap_fig_count) {
    
    cur_dataset_load <- paste(data_dir, cur_dataset, sep="")
    PatientCharacter <- read.csv(file=cur_dataset_load, header=TRUE, sep=",")
    
    expmat_load <- paste(data_dir, exp_mat, sep="")
    IgMexpmat <- read.csv(expmat_load, header=TRUE, row.names=1)
    
    
    Disease=factor(PatientCharacter$Disease)
    y = as.numeric(Disease == levels[2]) + 1
    
    sele=(!is.na(PatientCharacter$Disease))
    Disease=Disease[sele]
    
    IgMexpmat <- IgMexpmat %>% mutate(across(where(is.numeric), ~ .x+1))
    array_vst <- IgMexpmat
    
    # Hierarchical clustering using Complete Linkage
    if (cluster_rows & cluster_cols) 
    {
        d_proteins <- dist(array_vst, method = "euclidean")
        hc_proteins <- hclust(d_proteins, method = "complete" )
        
        d_patients <- dist(t(array_vst), method = "euclidean")
        hc_patients <- hclust(d_patients, method = "complete" )
        
        array_vst_cluster <- array_vst[unlist(hc_proteins["order"]), unlist(hc_patients["order"])]
        IgM_rownames <- row.names(array_vst_cluster)
    } else if (cluster_rows) 
    {
        hc_proteins = cluster_mat(array_vst, distance = "euclidean", method = "complete")
        
        array_vst_cluster <- array_vst[unlist(hc_proteins["order"]),]
        IgM_rownames <- row.names(array_vst_cluster)
    } else if (cluster_cols) 
    {
        hc_patients = cluster_mat(t(array_vst), distance = "euclidean", method = "complete")
        
        array_vst_cluster <- array_vst[, unlist(hc_patients["order"])]
        IgM_rownames <- row.names(array_vst_cluster)
    }
   
    res <- list(array_vst_cluster, IgM_rownames, array_vst)

}



deseq <- function(data_dir, cur_dataset, exp_mat, levels, lfcshrink_type, cluster_rows, cluster_cols, folder_name, heatmap_fig_count) {
    
    cur_dataset_load <- paste(data_dir, cur_dataset, sep="")
    PatientCharacter <- read.csv(file=cur_dataset_load, header=TRUE, sep=",")
    
    expmat_load <- paste(data_dir, exp_mat, sep="")
    IgMexpmat <- read.csv(expmat_load, header=TRUE, row.names=1)
    
    
    # Read in  data and Stromal Score
    Disease=factor(PatientCharacter$Disease)
    y = as.numeric(Disease == levels[2]) + 1

    
    sele=(!is.na(PatientCharacter$Disease))
    
    print(sele)
    
    # What columns are actually in the data?
    print(colnames(Disease))
    
    # Which ones are missing?
    setdiff(sele, colnames(Disease))
    
    Disease=Disease[sele]

    
    IgMexpmat <- IgMexpmat %>% mutate(across(where(is.numeric), ~ .x+1))
    count_input=IgMexpmat[,sele]
    
    x=PatientCharacter[sele,]
    output=c()
    for(i in 1:nrow(x)){
        if(x$Disease[i]==levels[2]){output=c(output,levels[2])}
        if(x$Disease[i]==levels[1]){output=c(output,levels[1])}
    }
    
    Disease=factor(output, levels = levels)
    
    
    # # Differential analysis
    samples <- data.frame(Disease)
    rownames(samples) <- PatientCharacter$TypeID[sele]
    dds <- DESeqDataSetFromMatrix(countData = round(count_input),
                                  colData = samples,
                                  design= ~ Disease)
    dds <- DESeq(dds)
    vst_dds <- varianceStabilizingTransformation(dds, blind = FALSE)
    array_vst <- assay(vst_dds)
    
    array_vst <- as.data.frame(array_vst)
    IgM_rownames <- row.names(array_vst)
    
    res1=results(dds,independentFiltering = FALSE)
    cat("Down-regulated:",sum(res1$log2FoldChange<(-1.0)&res1$padj<0.05,na.rm = TRUE),"\n")
    cat("Up-regulated:",sum(res1$log2FoldChange>1.0&res1$padj<0.05,na.rm = TRUE),"\n")
    coef = paste("Disease_", levels[2], '_vs_', levels[1], sep='')
    
    res2 = NULL
    
    if (lfcshrink_type == "GLM") {
        res2=res1
    } else {
        res2=lfcShrink(dds,coef = coef,type = lfcshrink_type)
    }
    
    # res2=lfcShrink(dds,coef = coef,type = lfcshrink_type)
    cat("Down-regulated:",sum(res2$log2FoldChange<(-1.0)&res2$padj<0.05,na.rm = TRUE),"\n")

    cat("Up-regulated:",sum(res2$log2FoldChange>1.0&res2$padj<0.05,na.rm = TRUE),"\n")

    
    res1name <- paste("./differential_expression_tables/", folder_name, sep = "", collapse = NULL)
    res1name <- paste(res1name, "/res1", sep = "", collapse = NULL)
    res1name <- paste(res1name, folder_name, sep = "_", collapse = NULL)
    res1name <- paste(res1name, "heatmap", sep = "_", collapse = NULL)
    res1name <- paste(res1name, lfcshrink_type, sep = "_", collapse = NULL)
    res1name <- paste(res1name, (heatmap_fig_count+1), sep = "_", collapse = NULL)
    res1name <- paste(res1name, "csv", sep = ".", collapse = NULL)
    
    write.csv(res1,file=res1name)
    
    # Hierarchical clustering using Complete Linkage
    if (cluster_rows & cluster_cols) {
        d_proteins <- dist(array_vst, method = "euclidean")
        hc_proteins <- hclust(d_proteins, method = "complete" )
        
        d_patients <- dist(t(array_vst), method = "euclidean")
        hc_patients <- hclust(d_patients, method = "complete" )
        
        array_vst_cluster <- array_vst[unlist(hc_proteins["order"]), unlist(hc_patients["order"])]
        IgM_rownames <- row.names(array_vst_cluster)
    }
    else if (cluster_rows) {
        hc_proteins = cluster_mat(array_vst, distance = "euclidean", method = "complete")
        
        array_vst_cluster <- array_vst[unlist(hc_proteins["order"]),]
        IgM_rownames <- row.names(array_vst_cluster)
    }
    else if (cluster_cols) {
        hc_patients = cluster_mat(t(array_vst), distance = "euclidean", method = "complete")
        
        array_vst_cluster <- array_vst[, unlist(hc_patients["order"])]
        IgM_rownames <- row.names(array_vst_cluster)
    }
    
    res <- list(array_vst_cluster, IgM_rownames, array_vst)
    
}

