library(data.table)
library(circlize)
library(tidyverse)
library(DESeq2)
library(VennDiagram)
library(fgsea)
library(Rcpp)
library(apeglm)
library(WGCNA)
library(cluster)
library(pheatmap)
library(jsonlite)
library(limma)
library(ashr)



wgcna_compute <- function(col_names, row_names, y, array_vst, levels, lfcshrink_type, cur_dataset, exp_mat, data_dir) {
        
    tryCatch({
        
        cur_dataset_load <- paste(data_dir, cur_dataset, sep="")
        PatientCharacter <- read.csv(file=cur_dataset_load, header=TRUE, sep=",")
        
        # Read in  data and Stromal Score
        Disease=factor(PatientCharacter$Disease)
    
        y <- unlist(y)
        array_vst <- t(array_vst)
        IgMexpmat_sub <- array_vst[unlist(row_names),unlist(col_names)] # col_names comes in from flask as a list
    
        # Langfelder, P. and Horvath, S., 2008. WGCNA: an R package for weighted correlation network analysis.
        # tutorial code
        
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        sft = pickSoftThreshold(IgMexpmat_sub, powerVector = powers, verbose = 5)
        
        if (!"fitIndices" %in% names(sft)) {
            stop("pickSoftThreshold() failure")
        }
        
        softPower = 10
        adjacency = adjacency(IgMexpmat_sub, power = softPower)
        
        TOM = TOMsimilarity(adjacency);
        
        colnames(TOM) <- colnames(adjacency)
        rownames(TOM) <- rownames(adjacency)
        
        
        dissTOM = 1-TOM
        
        
        
        # Call the hierarchical clustering function
        geneTree = hclust(as.dist(dissTOM), method = "average")
        # Plot the resulting clustering tree (dendrogram)
        module_hier1 <- geneTree
        
        
        
        minModuleSize = 30;
        # Module identification using dynamic tree cut:
        dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                    deepSplit = 2, pamRespectsDendro = FALSE,
                                    minClusterSize = minModuleSize);
        table(dynamicMods)
        
        
        
        # Convert numeric lables into colors
        dynamicColors = labels2colors(dynamicMods)
        table(dynamicColors)
        
        .networkTypes = c("unsigned", "signed", "signed hybrid");
        .adjacencyTypes = c(.networkTypes, "distance");
        
        
        .checkAndScaleWeights = function(weights, expr, scaleByMax = TRUE, verbose = 1)
        {
            if (length(weights)==0) return(weights);
            weights = as.matrix(weights);
            if (!isTRUE(all.equal(dim(expr), dim(weights))))
                stop("When 'weights' are given, they must have the same dimensions as 'expr'.")
            if (any(weights<0, na.rm = TRUE))
                stop("Found negative weights. All weights must be non-negative.");
            nf = !is.finite(weights);
            if (any(nf))
            {
                if (verbose > 0)
                    warning("Found non-finite weights. The corresponding data points will be removed.");
                weights[nf] = NA;
            }
            if (scaleByMax)
            {
                maxw = colMaxs(weights, na.rm = TRUE);
                maxw[maxw==0] = 1;
                weights = weights/matrix(maxw, nrow(weights), ncol(weights), byrow = TRUE);
            }
            weights;
        }
        
        
        adjacency = function(datExpr, selectCols=NULL, 
                             type = "unsigned", power = if (type=="distance") 1 else 6,
                             corFnc = "cor", corOptions = list(use = 'p'), weights = NULL, 
                             distFnc = "dist", distOptions = "method = 'euclidean'",
                             weightArgNames = c("weights.x", "weights.y"))
        {
            intType = charmatch(type, .adjacencyTypes)
            if (is.na(intType))
                stop(paste("Unrecognized 'type'. Recognized values are", paste(.adjacencyTypes, collapse = ", ")));
            
            corFnc.fnc = match.fun(corFnc);
            
            .checkAndScaleWeights(weights, datExpr, scaleByMax = FALSE);
            
            if (length(weights) > 0)
            {
                if (is.null(selectCols))
                {
                    if (is.list(corOptions))
                    {
                        weightOpt = list(weights.x = weights);
                        names(weightOpt) = weightArgNames[1];
                    } else weightOpt = spaste(weightArgNames[1], " = weights");
                } else {
                    if (is.list(corOptions))
                    {
                        weightOpt = list(weights.x = weights, weights.y = weights[, selectCols]);
                        names(weightOpt) = weightArgNames[c(1,2)];
                    } else weightOpt = spaste(weightArgNames[1], " = weights, ", weightArgNames[2], " = weights[, selectCols]");
                }
            } else {
                weightOpt = if (is.list(corOptions)) list() else ""
            }
            
            if (intType < 4)
            {
                if (is.null(selectCols))
                {
                    if (is.list(corOptions))
                    {
                        cor_mat = do.call(corFnc.fnc, c(list(x = datExpr), weightOpt, corOptions)) # option "cor" is fast pearson calculation, replacing standard "cor" function from R base package
                    } else {
                        corExpr = parse(text = paste(corFnc, "(datExpr ", prepComma(weightOpt), prepComma(corOptions), ")"));
                        # cor_mat = cor(datExpr, use = "p");
                        cor_mat = eval(corExpr);
                    }
                } else {
                    if (is.list(corOptions))
                    {
                        cor_mat = do.call(corFnc.fnc, c(list(x = datExpr, y = datExpr[, selectCols]), weightOpt, corOptions))
                    } else {
                        corExpr = parse(text = paste(corFnc, "(datExpr, datExpr[, selectCols] ", prepComma(weightOpt), 
                                                     prepComma(corOptions), ")"));
                        #cor_mat = cor(datExpr, datExpr[, selectCols], use="p");
                        cor_mat = eval(corExpr);
                    }
                }
            } else {
                if (!is.null(selectCols)) 
                    stop("The argument 'selectCols' cannot be used for distance adjacency.");
                if (is.list(distOptions))
                {
                    d = do.call(distFnc, c(list(x = t(datExpr)), distOptions));
                } else {
                    corExpr = parse(text = paste(distFnc, "(t(datExpr) ", prepComma(distOptions), ")"));
                    # cor_mat = cor(datExpr, use = "p");
                    d = eval(corExpr);
                }
                if (any(d<0)) 
                    warning("Function WGCNA::adjacency: Distance function returned (some) negative values.");
                cor_mat = 1-as.matrix( (d/max(d, na.rm = TRUE))^2 );
            }
            
            
            if (intType==1)
            {
                cor_mat = abs(cor_mat);
                # print(cor_mat)
                # cor_mat = cor_mat
                print("unsigned adjacency")
                
            } else if (intType==2) # signed
            {
                # cor_mat = (1+cor_mat)/2; # old, makes sure all output values are positive
                cor_mat = cor_mat # new, keep negative values
                print("signed adjacency")
            }
            else if (intType==3)
            { cor_mat[cor_mat < 0] = 0;
            }
            
            cor_mat^power;
            
            
        }
        
        
        plotNetworkHeatmap = function(datExpr,  plotGenes, weights = NULL, useTOM = TRUE, power = 6 , 
                                      networkType = "unsigned", main = "Heatmap of the network") 
        {
            match1=match( plotGenes ,colnames(datExpr) )
            match1=match1[ !is.na(match1)]
            nGenes=length(match1)
            if (  sum( !is.na(match1) )  != length(plotGenes) ) 
            {
                printFlush(paste("Warning: Not all gene names were recognized.", 
                                 "Only the following genes were recognized. "));
                printFlush(paste("   ", colnames(datExpr)[match1], collapse = ", " ))
            }
            if (nGenes< 3 ) 
            { 
                warning(paste("Since you have fewer than 3 genes, the network will not be visualized.\n",
                              "   Hint: please input more genes.")); plot(1,1)
            } else {
                datErest=datExpr[, match1 ]
                # print(datErest)
                # print(nrow((datErest)))
                # print(ncol((datErest)))
                # 
                if (!is.null(weights)) weights = weights[, match1];
                print(networkType)
                ADJ1 = adjacency(datErest, weights = weights, power = power, type = networkType)
                
                if (useTOM) {
                    diss1= 1-TOMsimilarity(ADJ1)   
                } else {
                    diss1 = 1-ADJ1;
                }
                diag(diss1)=NA
                hier1=fastcluster::hclust(as.dist(diss1), method="average" )
                colors1=rep("white", nGenes)
                labeltree = names(data.frame(datErest))
                
                
                labelrow  = names(data.frame(datErest))
                labelrow[hier1$order[length(labeltree):1]]=labelrow[hier1$order]
                
                options(expressions = 10000)
                
                # print(labeltree)
                
                diss1_reordered <- diss1[unlist(hier1["order"]), unlist(hier1["order"])]
                
                if (useTOM) {
                    rownames(diss1) <- labeltree
                    colnames(diss1) <- labeltree
                }
    
                diss1 = 1.0-diss1
               
                
                newList <- list("diss1" = diss1, "hier1" = hier1, "labeltree" = labeltree)
                
                return (newList)
                
                
            } # end of if (nGenes> 2 )
        } # end of function
        
        
        # find most important proteins, based on eigengenes
        datME=moduleEigengenes(IgMexpmat_sub,dynamicColors)$eigengenes
        datKME=signedKME(IgMexpmat_sub, datME, outputColumnName="MM.")
        
        NS1=networkScreening(y=y, datME=datME, datExpr=IgMexpmat_sub,
                             oddPower=3, blockSize=900, minimumSampleSize=4,
                             addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)
        
        
        topList=rank(NS1$p.Weighted,ties.method="first")<=30
        gene.names= colnames(IgMexpmat_sub)[topList]
        
        # The following shows the correlations between the top genes
        newList = plotNetworkHeatmap(IgMexpmat_sub, plotGenes = gene.names,
                                     networkType="signed", useTOM=FALSE,
                                     power=1, main="signed correlations")
        
        
        diss1 <- newList$diss1
        hier1 <- newList$hier1
        labeltree <- newList$labeltree
        
        
        diss1_clustered <- as.data.frame(diss1[unlist(hier1["order"]), unlist(hier1["order"])])
        gene_names <- row.names(diss1_clustered)
        
        # write.csv(diss1_clustered,file="100_diss1.csv",fileEncoding = "UTF-8")
        
        dynamicMods_all <- dynamicMods
        dynamicMods_vis <- dynamicMods[geneTree$order]
        colors_vis <- dynamicColors[geneTree$order]
        
        
        # split proteins by module 
        # list of groups of modules of proteins
        all_protein_names <- colnames(dissTOM)
        all_protein_names <- all_protein_names[geneTree$order]
        
        paired <- mapply(c, all_protein_names, dynamicMods_vis, SIMPLIFY=FALSE)
        l0 <- list()
        l1 <- list()
        l2 <- list()
        mods <- list()
        mods_expmat <- list()
        ns <- list()
        diss_l <- list()
        genenames_l <- list()
        
        
        for (p in paired) {
            val1 <- p[1]
            val2 <- p[2]
            
            mods[[val2]] <- list()
            mods_expmat[[val2]] <- list()
            ns[[val2]] <- list()
        }
        
        for (p in paired) {
            val1 <- p[1]
            val2 <- p[2]
            
            mods[[val2]] <- append(mods[[val2]],val1)
        }
        
        # subset expr mat for each module
        for (key in names(mods)) {
            expmat_sub <- IgMexpmat_sub[, unlist(mods[[key]])]
            mods_expmat[[key]] <- expmat_sub
        }
        
        # redo adjacency, tree cut, and eigengenes for each module
        
        for (key in names(mods_expmat)) {
            exp_mat <- mods_expmat[[key]]
            powers = c(c(1:10), seq(from = 12, to=20, by=2))
            sft = pickSoftThreshold(exp_mat, powerVector = powers, verbose = 5)
            
            if (!"fitIndices" %in% names(sft)) {
                stop("pickSoftThreshold() failure")
            }
            
            softPower = 10
            
            print("TOM adjacency")
            adj = adjacency(exp_mat, power = softPower)
            
            TOM = TOMsimilarity(adj);
            
            colnames(TOM) <- colnames(adj)
            rownames(TOM) <- rownames(adj)
            
            
            dissTOM = 1-TOM
            
            
            # Call the hierarchical clustering function
            geneTree = hclust(as.dist(dissTOM), method = "average")
            
            
            
            minModuleSize = 30;
            # Module identification using dynamic tree cut:
            dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                        deepSplit = 2, pamRespectsDendro = FALSE,
                                        minClusterSize = minModuleSize);
            table(dynamicMods)
            
            dynamicColors = labels2colors(dynamicMods)
            table(dynamicColors)
           
            
            # find most important proteins, based on eigengenes
            datME=moduleEigengenes(exp_mat,dynamicColors)$eigengenes # first PC per sample (module EG per sample)
            datKME=signedKME(exp_mat, datME, outputColumnName="MM.") # pearson cor between exp mat and module EG per sample, for each module
            
            NS1=networkScreening(y=y, datME=datME, datExpr=exp_mat,
                                 oddPower=3, blockSize=900, minimumSampleSize=4,
                                 addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)
            
            
            ns[[key]] <- NS1
    
            
            topList=rank(NS1$p.Weighted,ties.method="first")<=30
            gene.names= colnames(exp_mat)[topList]
            
            
            # The following shows the correlations between the top genes
            newList = plotNetworkHeatmap(exp_mat, plotGenes = gene.names,
                                         networkType="signed", useTOM=FALSE,
                                         power=1, main="signed correlations")
            
            
            diss1 <- newList$diss1
            hier1 <- newList$hier1
            labeltree <- newList$labeltree
            
            
            
            diss1 <- as.data.frame(diss1[unlist(hier1["order"]), unlist(hier1["order"])])
            gene_names <- row.names(diss1)
            
    
            diss_l <- append(diss_l, list(diss1))
            genenames_l <- append(genenames_l, list(gene_names))
            
        }
        
        return (list(diss1_clustered, gene_names, module_hier1$order, module_hier1$merge, module_hier1$height, dynamicMods_all,  dynamicMods_vis, diss_l, genenames_l))
        
    }, error = function(e) {
        stop(e)  
    })
}



wgcna_compute_all <- function(col_names, row_names, array_vst, levels, lfcshrink_type, cur_dataset, exp_mat, data_dir) {

    tryCatch({
        
        cur_dataset_load <- paste(data_dir, cur_dataset, sep="")
        PatientCharacter <- read.csv(file=cur_dataset_load, header=TRUE, sep=",")
        
        # Read in  data and Stromal Score
        Disease=factor(PatientCharacter$Disease)
        y = as.numeric(Disease == levels[2]) + 1
        
        array_vst <- t(array_vst)
        IgMexpmat_sub <- array_vst[unlist(row_names),unlist(col_names)] # col_names comes in from flask as a list
    
        if (!(all(is.finite(IgMexpmat_sub)))) {
            stop("NA/NaN/Inf in your WGCNA selection. Please try with different subset of entities for clustering. ")
        }
        
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        sft = pickSoftThreshold(IgMexpmat_sub, powerVector = powers, verbose = 5)
        
        if (!"fitIndices" %in% names(sft)) {
            stop("pickSoftThreshold() failure")
        }
        
        softPower = 10
        adjacency = adjacency(IgMexpmat_sub, power = softPower)
        
        TOM = TOMsimilarity(adjacency);
        
        colnames(TOM) <- colnames(adjacency)
        rownames(TOM) <- rownames(adjacency)
        
        
        dissTOM = 1-TOM
        
        
        
        # Call the hierarchical clustering function
        geneTree = hclust(as.dist(dissTOM), method = "average")
        
        module_hier1 <- geneTree
        
        
        
        minModuleSize = 30;
        # Module identification using dynamic tree cut:
        dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                    deepSplit = 2, pamRespectsDendro = FALSE,
                                    minClusterSize = minModuleSize);
        table(dynamicMods)
        
        
        
        # Convert numeric lables into colors
        dynamicColors = labels2colors(dynamicMods)
        table(dynamicColors)
       
        .networkTypes = c("unsigned", "signed", "signed hybrid");
        .adjacencyTypes = c(.networkTypes, "distance");
        
        
        .checkAndScaleWeights = function(weights, expr, scaleByMax = TRUE, verbose = 1)
        {
            if (length(weights)==0) return(weights);
            weights = as.matrix(weights);
            if (!isTRUE(all.equal(dim(expr), dim(weights))))
                stop("When 'weights' are given, they must have the same dimensions as 'expr'.")
            if (any(weights<0, na.rm = TRUE))
                stop("Found negative weights. All weights must be non-negative.");
            nf = !is.finite(weights);
            if (any(nf))
            {
                if (verbose > 0)
                    warning("Found non-finite weights. The corresponding data points will be removed.");
                weights[nf] = NA;
            }
            if (scaleByMax)
            {
                maxw = colMaxs(weights, na.rm = TRUE);
                maxw[maxw==0] = 1;
                weights = weights/matrix(maxw, nrow(weights), ncol(weights), byrow = TRUE);
            }
            weights;
        }
        
        
        adjacency = function(datExpr, selectCols=NULL, 
                             type = "unsigned", power = if (type=="distance") 1 else 6,
                             corFnc = "cor", corOptions = list(use = 'p'), weights = NULL, 
                             distFnc = "dist", distOptions = "method = 'euclidean'",
                             weightArgNames = c("weights.x", "weights.y"))
        {
            intType = charmatch(type, .adjacencyTypes)
            if (is.na(intType))
                stop(paste("Unrecognized 'type'. Recognized values are", paste(.adjacencyTypes, collapse = ", ")));
            
            corFnc.fnc = match.fun(corFnc);
            
            .checkAndScaleWeights(weights, datExpr, scaleByMax = FALSE);
            
            if (length(weights) > 0)
            {
                if (is.null(selectCols))
                {
                    if (is.list(corOptions))
                    {
                        weightOpt = list(weights.x = weights);
                        names(weightOpt) = weightArgNames[1];
                    } else weightOpt = spaste(weightArgNames[1], " = weights");
                } else {
                    if (is.list(corOptions))
                    {
                        weightOpt = list(weights.x = weights, weights.y = weights[, selectCols]);
                        names(weightOpt) = weightArgNames[c(1,2)];
                    } else weightOpt = spaste(weightArgNames[1], " = weights, ", weightArgNames[2], " = weights[, selectCols]");
                }
            } else {
                weightOpt = if (is.list(corOptions)) list() else ""
            }
            
            if (intType < 4)
            {
                if (is.null(selectCols))
                {
                    if (is.list(corOptions))
                    {
                        cor_mat = do.call(corFnc.fnc, c(list(x = datExpr), weightOpt, corOptions)) # option "cor" is fast pearson calculation, replacing standard "cor" function from R base package
                    } else {
                        corExpr = parse(text = paste(corFnc, "(datExpr ", prepComma(weightOpt), prepComma(corOptions), ")"));
                        # cor_mat = cor(datExpr, use = "p");
                        cor_mat = eval(corExpr);
                    }
                } else {
                    if (is.list(corOptions))
                    {
                        cor_mat = do.call(corFnc.fnc, c(list(x = datExpr, y = datExpr[, selectCols]), weightOpt, corOptions))
                    } else {
                        corExpr = parse(text = paste(corFnc, "(datExpr, datExpr[, selectCols] ", prepComma(weightOpt), 
                                                     prepComma(corOptions), ")"));
                        cor_mat = eval(corExpr);
                    }
                }
            } else {
                if (!is.null(selectCols)) 
                    stop("The argument 'selectCols' cannot be used for distance adjacency.");
                if (is.list(distOptions))
                {
                    d = do.call(distFnc, c(list(x = t(datExpr)), distOptions));
                } else {
                    corExpr = parse(text = paste(distFnc, "(t(datExpr) ", prepComma(distOptions), ")"));
                    d = eval(corExpr);
                }
                if (any(d<0)) 
                    warning("Function WGCNA::adjacency: Distance function returned (some) negative values.");
                cor_mat = 1-as.matrix( (d/max(d, na.rm = TRUE))^2 );
            }
            
            
            if (intType==1)
            {
                cor_mat = abs(cor_mat);
                print("unsigned adjacency")
                
            } else if (intType==2) # signed
            {
                cor_mat = cor_mat # new, keep negative values
                print("signed adjacency")
            }
            else if (intType==3)
            { cor_mat[cor_mat < 0] = 0;
            }
            
            cor_mat^power;
            
            
        }
        
        
        plotNetworkHeatmap = function(datExpr,  plotGenes, weights = NULL, useTOM = TRUE, power = 6 , 
                                      networkType = "unsigned", main = "Heatmap of the network") 
        {
            match1=match( plotGenes ,colnames(datExpr) )
            match1=match1[ !is.na(match1)]
            nGenes=length(match1)
            if (  sum( !is.na(match1) )  != length(plotGenes) ) 
            {
                printFlush(paste("Warning: Not all gene names were recognized.", 
                                 "Only the following genes were recognized. "));
                printFlush(paste("   ", colnames(datExpr)[match1], collapse = ", " ))
            }
            if (nGenes< 3 ) 
            { 
                warning(paste("Since you have fewer than 3 genes, the network will not be visualized.\n",
                              "   Hint: please input more genes.")); plot(1,1)
            } else {
                datErest=datExpr[, match1 ]
               
                if (!is.null(weights)) weights = weights[, match1];
                print(networkType)
                ADJ1 = adjacency(datErest, weights = weights, power = power, type = networkType)
                
                if (useTOM) {
                    diss1= 1-TOMsimilarity(ADJ1)   
                } else {
                    diss1 = 1-ADJ1;
                }
                diag(diss1)=NA
                hier1=fastcluster::hclust(as.dist(diss1), method="average" )
                colors1=rep("white", nGenes)
                labeltree = names(data.frame(datErest))
                
                
                labelrow  = names(data.frame(datErest))
                labelrow[hier1$order[length(labeltree):1]]=labelrow[hier1$order]
                
                options(expressions = 10000)
                
    
                diss1_reordered <- diss1[unlist(hier1["order"]), unlist(hier1["order"])]
                
                if (useTOM) {
                    rownames(diss1) <- labeltree
                    colnames(diss1) <- labeltree
                }
    
                diss1 = 1.0-diss1
            
                
                newList <- list("diss1" = diss1, "hier1" = hier1, "labeltree" = labeltree)
                
                return (newList)
                
                
            } # end of if (nGenes> 2 )
        } # end of function
        
        
        # find most important proteins, based on eigengenes
        datME=moduleEigengenes(IgMexpmat_sub,dynamicColors)$eigengenes
        datKME=signedKME(IgMexpmat_sub, datME, outputColumnName="MM.")
        
        print(" network screening ")
        NS1=networkScreening(y=y, datME=datME, datExpr=IgMexpmat_sub,
                             oddPower=3, blockSize=900, minimumSampleSize=4,
                             addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)
        
        
        topList=rank(NS1$p.Weighted,ties.method="first")<=30
        gene.names= colnames(IgMexpmat_sub)[topList]
        
        print(" plt network heatmap ")
        # The following shows the correlations between the top genes
        newList = plotNetworkHeatmap(IgMexpmat_sub, plotGenes = gene.names,
                                     networkType="signed", useTOM=FALSE,
                                     power=1, main="signed correlations")
        
        
        diss1 <- newList$diss1
        hier1 <- newList$hier1
        labeltree <- newList$labeltree
        
        
        
        diss1_clustered <- as.data.frame(diss1[unlist(hier1["order"]), unlist(hier1["order"])])
        gene_names <- row.names(diss1_clustered)
        
    
        dynamicMods_all <- dynamicMods
        dynamicMods_vis <- dynamicMods[geneTree$order]
        colors_vis <- dynamicColors[geneTree$order]
        
        
        # split proteins by module 
        # list of groups of modules of proteins
        all_protein_names <- colnames(dissTOM)
        all_protein_names <- all_protein_names[geneTree$order]
        
        paired <- mapply(c, all_protein_names, dynamicMods_vis, SIMPLIFY=FALSE)
        l0 <- list()
        l1 <- list()
        l2 <- list()
        mods <- list()
        mods_expmat <- list()
        ns <- list()
        diss_l <- list()
        genenames_l <- list()
        
        
        for (p in paired) {
            val1 <- p[1]
            val2 <- p[2]
            
            mods[[val2]] <- list()
            mods_expmat[[val2]] <- list()
            ns[[val2]] <- list()
        }
        
        for (p in paired) {
            val1 <- p[1]
            val2 <- p[2]
            
            mods[[val2]] <- append(mods[[val2]],val1)
        }
        
        # subset expr mat for each module
        for (key in names(mods)) {
            expmat_sub <- IgMexpmat_sub[, unlist(mods[[key]])]
            mods_expmat[[key]] <- expmat_sub
        }
        
        # redo adjacency, tree cut, and eigengenes for each module
        
        for (key in names(mods_expmat)) {
            exp_mat <- mods_expmat[[key]]
            powers = c(c(1:10), seq(from = 12, to=20, by=2))
            
            if (!(all(is.finite(exp_mat)))) {
                stop("WGCNA failed: error in generate intermediate value. NA/NaN/Inf detected in intermediate.")
            }
            
            sft <- NULL 
            
            sft = pickSoftThreshold(exp_mat, powerVector = powers, verbose = 5)
            
            
            # sft <<- pickSoftThreshold(exp_mat, powerVector = powers, verbose = 5)
            
            if (!"fitIndices" %in% names(sft)) {
                stop("pickSoftThreshold() failure")
            }
            
            softPower = 10
            
            print("TOM adjacency")
            adj = adjacency(exp_mat, power = softPower)
            
            TOM = TOMsimilarity(adj);
            
            # TODO add row and col names to TOM matrix
            colnames(TOM) <- colnames(adj)
            rownames(TOM) <- rownames(adj)
            
            
            dissTOM = 1-TOM
            
            
            # Call the hierarchical clustering function
            geneTree = hclust(as.dist(dissTOM), method = "average")
            
            
            
            minModuleSize = 30;
            # Module identification using dynamic tree cut:
            dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                        deepSplit = 2, pamRespectsDendro = FALSE,
                                        minClusterSize = minModuleSize);
            table(dynamicMods)
            
            dynamicColors = labels2colors(dynamicMods)
            table(dynamicColors)
            
            # find most important proteins, based on eigengenes
            datME=moduleEigengenes(exp_mat,dynamicColors)$eigengenes # first PC per sample (module EG per sample)
            datKME=signedKME(exp_mat, datME, outputColumnName="MM.") # pearson cor between exp mat and module EG per sample, for each module
            
            NS1=networkScreening(y=y, datME=datME, datExpr=exp_mat,
                                 oddPower=3, blockSize=900, minimumSampleSize=4,
                                 addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)
            
            ns[[key]] <- NS1
    
            
            topList=rank(NS1$p.Weighted,ties.method="first")<=30
            gene.names= colnames(exp_mat)[topList]
            
            
            # The following shows the correlations between the top genes
            newList = plotNetworkHeatmap(exp_mat, plotGenes = gene.names,
                                         networkType="signed", useTOM=FALSE,
                                         power=1, main="signed correlations")
            
            
            diss1 <- newList$diss1
            hier1 <- newList$hier1
            labeltree <- newList$labeltree
            
            
            diss1 <- as.data.frame(diss1[unlist(hier1["order"]), unlist(hier1["order"])])
            gene_names <- row.names(diss1)
            
    
            diss_l <- append(diss_l, list(diss1))
            genenames_l <- append(genenames_l, list(gene_names))
            
        }
        
        return (list(diss1_clustered, gene_names, module_hier1$order, module_hier1$merge, module_hier1$height, dynamicMods_all,  dynamicMods_vis, diss_l, genenames_l))
    
    }, error = function(e) {
        stop(e)
    })

    
}



volcano <- function(data_dir, cur_dataset, exp_mat, levels, lfcshrink_type, folder_name, volcano_fig_count) {
    # cur_dataset <- 'LupusPatientCharacter.csv'
    # data_dir <- "./processed_data/"
    # exp_mat <- "expmat_SNR mean 635.csv"
    # levels <- c("Disease", "Control")
    # lfcshrink_type <- "normal"
    # folder_name<-"LupusPatientCharacter"
    # volcano_fig_count <- 0

    cur_dataset_load <- paste(data_dir, cur_dataset, sep="")
    PatientCharacter <- read.csv(file=cur_dataset_load, header=TRUE, sep=",")
    # load(cur_dataset_load)
    
    expmat_load <- paste(data_dir, exp_mat, sep="")
    IgMexpmat <- read.csv(expmat_load, header=TRUE, row.names=1)
    

    
    # Read in  data and Stromal Score
    Disease=factor(PatientCharacter$Disease)
    #strScore=MNAdsAdvRxns$StromalScore
    
    
    sele=(!is.na(PatientCharacter$Disease))
    Disease=Disease[sele]
    
    
    IgMexpmat <- IgMexpmat %>% mutate(across(where(is.numeric), ~ .x+1))
    count_input=IgMexpmat[,sele]
    
    x=PatientCharacter[sele,]
    output=c()
    for(i in 1:nrow(x)){
        if(x$Disease[i]==levels[2]){output=c(output,levels[2])}
        if(x$Disease[i]==levels[1]){output=c(output,levels[1])}
    }
    
    levels <- c(levels[1], levels[2])
    
    Disease=factor(output, levels = levels)
    
    # Differential analysis
    samples <- data.frame(Disease)
    print(samples)
    rownames(samples) <- PatientCharacter$TypeID[sele]
    print("rownames")
 
    compare_names <- data.frame(
        index     = seq_along(colnames(round(count_input))),
        countData = colnames(round(count_input)),
        colData   = rownames(samples),
        equal     = colnames(round(count_input)) == rownames(samples)
    )
    
    print(head(compare_names, 20))  # or View(compare_names)
    
    dds <- DESeqDataSetFromMatrix(countData = round(count_input),
                                  colData = samples,
                                  design= ~ Disease)
    dds <- DESeq(dds)
    vst_dds <- varianceStabilizingTransformation(dds, blind = FALSE)
    array_vst <- assay(vst_dds)
    
    # Diseaseneg vs. Diseasepos
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
    # cat("Down-regulated:",sum(res2$log2FoldChange<(1.0)&res2$pvalue<0.05,na.rm = TRUE),"\n")
    
    cat("Up-regulated:",sum(res2$log2FoldChange>1.0&res2$padj<0.05,na.rm = TRUE),"\n")
    
    res1name <- paste("./differential_expression_tables/", folder_name, sep = "", collapse = NULL)
    res1name <- paste(res1name, "/res1", sep = "", collapse = NULL)
    res1name <- paste(res1name, folder_name, sep = "_", collapse = NULL)
    res1name <- paste(res1name, "volcano", sep = "_", collapse = NULL)
    res1name <- paste(res1name, lfcshrink_type, sep = "_", collapse = NULL)
    res1name <- paste(res1name, (volcano_fig_count+1), sep = "_", collapse = NULL)
    res1name <- paste(res1name, "csv", sep = ".", collapse = NULL)
    
    res2name <- paste("./differential_expression_tables/", folder_name, sep = "", collapse = NULL)
    res2name <- paste(res2name, "/res2Shrink", sep = "", collapse = NULL)
    res2name <- paste(res2name, folder_name, sep = "_", collapse = NULL)
    res2name <- paste(res2name, "volcano", sep = "_", collapse = NULL)
    res2name <- paste(res2name, lfcshrink_type, sep = "_", collapse = NULL)
    res2name <- paste(res2name, (volcano_fig_count+1), sep = "_", collapse = NULL)
    res2name <- paste(res2name, "csv", sep = ".", collapse = NULL)
    
    # res2name <- paste("./differential_expression_tables/res2Shrink", cur_dataset, sep = "_", collapse = NULL)
    # res2name <- paste(res2name, "csv", sep = ".", collapse = NULL)
    
    write.csv(res1,file=res1name)
    
    
    resShrink=data.frame(res2$log2FoldChange,res2$pvalue, res2$padj)
    colnames(resShrink)=c("Log2Fold_Change","pvalue","adjPvalue")
    rownames(resShrink)=rownames(count_input)
    

    adjpvalue_cutoff_Lupus=0.05
    adjpvalue_cutoff_Control=0.05
    
    # now adjPvalue not padj
    resShrink$label_interest=ifelse((resShrink$pvalue<=adjpvalue_cutoff_Lupus & resShrink$Log2Fold_Change>1.0),rownames(resShrink), ifelse((resShrink$pvalue<=adjpvalue_cutoff_Control & resShrink$Log2Fold_Change<(-1.0)),rownames(resShrink),''))
    
    
    
    resShrink$neglog10pval <- -log10(resShrink$pvalue)
    resShrink$name <- rownames(resShrink)
    
    write.csv(resShrink,file=res2name)
    
    res <- list(resShrink)
}


pathways <- function(data_dir, cur_dataset, exp_mat, levels, pathways_data, lfcshrink_type, folder_name, pathway_fig_count) {
    # pathways_data <- c("B_Cell_Receptor_Signaling.txt", "Cell_Type_Astrocytes.txt")
    # data_dir <- "./processed_data/"
    # cur_dataset <- "PatientCharacter1.csv"
    # lfcshrink_type <- "normal"
    # levels <- c("Control", "Disease")
    # exp_mat <- "ExpMat1.csv"
    
    # cur_dataset <- 'PatientCharacter2.csv'
    # data_dir <- "./processed_data/"
    # exp_mat <- "ExpMat2.csv"
    # levels <- c("Disease", "Control")
    # lfcshrink_type <- "apeglm"
    # folder_name<-"PatientCharacter2"
    # pathway_fig_count <- 0

    
    cur_dataset_load <- paste(data_dir, cur_dataset, sep="")
    PatientCharacter <- read.csv(file=cur_dataset_load, header=TRUE, sep=",")
    
    
    
    expmat_load <- paste(data_dir, exp_mat, sep="")
    IgMexpmat <- read.csv(expmat_load, header=TRUE, row.names=1)
    
   
    # Read in  data and Stromal Score
    Disease=factor(PatientCharacter$Disease)

    sele=(!is.na(PatientCharacter$Disease))
    Disease=Disease[sele]
    
    IgMexpmat <- IgMexpmat %>% mutate(across(where(is.numeric), ~ .x+1))
    count_input=IgMexpmat[,sele]
    
    x=PatientCharacter[sele,]
    output=c()
    for(i in 1:nrow(x)){
        if(x$Disease[i]==levels[2]){output=c(output,levels[2])}
        if(x$Disease[i]==levels[1]){output=c(output,levels[1])}
    }
    
    levels <- c(levels[1], levels[2])
    
    Disease=factor(output, levels = levels)
    
    # Differential analysis
    samples <- data.frame(Disease)
    rownames(samples) <- PatientCharacter$TypeID[sele]
    dds <- DESeqDataSetFromMatrix(countData = round(count_input),
                                  colData = samples,
                                  design= ~ Disease)
    dds <- DESeq(dds)
    vst_dds <- varianceStabilizingTransformation(dds, blind = FALSE)
    array_vst <- assay(vst_dds)
    
    # Diseaseneg vs. Diseasepos
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

    cat("Down-regulated:",sum(res2$log2FoldChange<(-1.0)&res2$padj<0.05,na.rm = TRUE),"\n")

    cat("Up-regulated:",sum(res2$log2FoldChange>1.0&res2$padj<0.05,na.rm = TRUE),"\n")

    res1name <- paste("./differential_expression_tables/", folder_name, sep = "", collapse = NULL)
    res1name <- paste(res1name, "/res1", sep = "", collapse = NULL)
    res1name <- paste(res1name, folder_name, sep = "_", collapse = NULL)
    res1name <- paste(res1name, "pathway", sep = "_", collapse = NULL)
    res1name <- paste(res1name, lfcshrink_type, sep = "_", collapse = NULL)
    res1name <- paste(res1name, (pathway_fig_count+1), sep = "_", collapse = NULL)
    res1name <- paste(res1name, "csv", sep = ".", collapse = NULL)
    
    res2name <- paste("./differential_expression_tables/", folder_name, sep = "", collapse = NULL)
    res2name <- paste(res2name, "/res2Shrink", sep = "", collapse = NULL)
    res2name <- paste(res2name, folder_name, sep = "_", collapse = NULL)
    res2name <- paste(res2name, "pathway", sep = "_", collapse = NULL)
    res2name <- paste(res2name, lfcshrink_type, sep = "_", collapse = NULL)
    res2name <- paste(res2name, (pathway_fig_count+1), sep = "_", collapse = NULL)
    res2name <- paste(res2name, "csv", sep = ".", collapse = NULL)
    
    write.csv(res1,file=res1name)
    
    
    resShrink=data.frame(res2$log2FoldChange,res2$pvalue, res2$padj)
    colnames(resShrink)=c("Log2Fold_Change","pvalue","adjPvalue")
    rownames(resShrink)=rownames(count_input)
    

    adjpvalue_cutoff_Lupus=0.05
    adjpvalue_cutoff_Control=0.05
    
    # now adjPvalue not padj
    resShrink$label_interest=ifelse((resShrink$pvalue<=adjpvalue_cutoff_Lupus & resShrink$Log2Fold_Change>1.0),rownames(resShrink), ifelse((resShrink$pvalue<=adjpvalue_cutoff_Control & resShrink$Log2Fold_Change<(-1.0)),rownames(resShrink),''))
    
    # resShrink$label_interest=ifelse((resShrink$adjPvalue<=adjpvalue_cutoff_Lupus & resShrink$Log2Fold_Change>1.0),rownames(resShrink), ifelse((resShrink$adjPvalue<=adjpvalue_cutoff_Control & resShrink$Log2Fold_Change<(-1.0)),rownames(resShrink),''))
    # 
    
    resShrink$neglog10pval <- -log10(resShrink$pvalue)
    
    # resShrink$neglog10pval <- -log10(resShrink$adjPvalue)
    resShrink$name <- rownames(resShrink)
    
    
    res2=res2
    
    # stat2=res2$stat
    

    stat2 = 0
    if (lfcshrink_type == "normal") {
        stat2=res2$stat
    } else {
        stat2=res2$pvalue
    }
    
    names(stat2)=rownames(count_input)
    
    # fgseaname <- paste("./differential_expression_tables/fgsea", cur_dataset, sep = "_", collapse = NULL)
    # fgseaname <- paste(fgseaname, "tsv", sep = ".", collapse = NULL)
    
    fgseaname <- paste("./differential_expression_tables/", folder_name, sep = "", collapse = NULL)
    fgseaname <- paste(fgseaname, "/fgsea", sep = "", collapse = NULL)
    fgseaname <- paste(fgseaname, folder_name, sep = "_", collapse = NULL)
    fgseaname <- paste(fgseaname, "pathway", sep = "_", collapse = NULL)
    fgseaname <- paste(fgseaname, (pathway_fig_count+1), sep = "_", collapse = NULL)
    fgseaname <- paste(fgseaname, "tsv", sep = ".", collapse = NULL)
    
    
    pathways_list = list()
    
    for (e in pathways_data) {
        pathway_load <- paste(data_dir, e, sep="")
        pathway_name <- tools::file_path_sans_ext(basename(e))
        print(pathway_name)
        pathways_list[[ pathway_name ]] <- readLines(pathway_load)
    }
    
    fgseaRes <- fgsea(pathways = pathways_list, stats =stat2,nperm = 1000000)
    
    fwrite(fgseaRes, file=fgseaname, sep="\t", sep2=c("", " ", ""))
    
    #  Visualize
    fig=data.frame(fgseaRes$pathway,fgseaRes$pval,fgseaRes$NES)
    fig$stars=cut(fig$fgseaRes.pval, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("3", "2", "1", "0"))
    
    
    class(fig$fgseaRes.pathway)
    indexfactor=fig$fgseaRes.pathway[order(fig$fgseaRes.NES,decreasing=F)]
   
    
    write.csv(resShrink,file=res2name)
    
    res <- list(fig, indexfactor)
}


get_gpr_header <- function(gprDirectory) {

    
    # grab any gpr file to get header. all files should have the same column names
    gprFiles <- list.files(path = gprDirectory, pattern = "\\.gpr$", full.names = TRUE)
    fullname <- gprFiles[1]
    
    source <- "genepix.custom"
    
    source2 <- strsplit(source,split=".",fixed=TRUE)[[1]][1]
    
    switch(source2,
           quantarray = {
               firstfield <- scan(fullname,what="",sep="\t",flush=TRUE,quiet=TRUE,blank.lines.skip=FALSE,multi.line=FALSE,allowEscapes=FALSE)
               skip <- grep("Begin Data",firstfield)
               
           }, arrayvision = {
               skip <- 1
           }, 
           # SHERRY we're not using readGenericHeader. it defeats the purpose of extracting column names if we have to give it the column names.
           genepix = {
               h <- readGPRHeader(fullname)
               skip <- h$NHeaderRecords
           }, imagene9 = {
               h <- readImaGeneHeader(fullname)
               skip <- h$NHeaderRecords
           }, smd = {
               skip <- readSMDHeader(fullname)$NHeaderRecords
               nspots <- nrow(obj)
           }
           )
    
    
    res <- scan(fullname,what="list", skip=skip, nlines=1)

}

compute_expmat_snr <- function(snr_type, gprDirectory, data_dir) {
    read.maimages <- function(files=NULL,source="generic",path=NULL,ext=NULL,names=NULL,columns=NULL,other.columns=NULL,annotation=NULL,green.only=FALSE,wt.fun=NULL,verbose=TRUE,sep="\t",quote=NULL,...)
        #	Extracts either an RGList object from a set of two-color image analysis output files
        #	or an EListRaw object from a set of one-color files
        #	Gordon Smyth. 
        #	1 Nov 2002.  Last revised 26 Feb 2019.
    {
        
        #	Determine type of input file
        source <- match.arg(source,c("generic","agilent","agilent.mean","agilent.median","arrayvision","arrayvision.ARM","arrayvision.MTM","bluefuse","genepix","genepix.mean","genepix.median","genepix.custom","imagene","imagene9","quantarray","scanarrayexpress","smd.old","smd","spot","spot.close.open"))
        #	source2 is the source type with qualifications removed
        source2 <- strsplit(source,split=".",fixed=TRUE)[[1]][1]
        
        #	ImaGene is special case
        if(source2=="imagene") return(read.imagene(files=files,path=path,ext=ext,names=names,columns=columns,other.columns=other.columns,wt.fun=wt.fun,verbose=verbose,sep=sep,quote=quote,...))
        
        #	Get list of files and associated targets frame
        if(is.null(files)) {
            if(is.null(ext))
                stop("Must specify input files")
            else {
                extregex <- paste("\\.",ext,"$",sep="")
                files <- dir(path=ifelse(is.null(path),".",path),pattern=extregex)
                files <- sub(extregex,"",files)
            }
        } else
            if(is.data.frame(files)) {
                targets <- files
                files <- files$FileName
                if(is.null(files)) stop("targets frame doesn't contain FileName column")
                if(is.null(names)) names <- targets$Label
            } else {
                targets <- NULL
            }
        slides <- as.vector(as.character(files))
        if(!is.null(ext)) slides <- paste(slides,ext,sep=".")
        nslides <- length(slides)
        
        #	Default sample names
        if(is.null(names)) names <- removeExt(files)
        
        #	Default for quote
        if(is.null(quote)) if(source2=="agilent") quote <- "" else quote <- "\""
        
        if(is.null(columns)) {
            if(source2=="generic") stop("must specify columns for generic input")
            columns <- switch(source,
                              agilent.mean = list(G="gMeanSignal",Gb="gBGMedianSignal",R="rMeanSignal",Rb="rBGMedianSignal"),
                              agilent =,
                              agilent.median = list(G="gMedianSignal",Gb="gBGMedianSignal",R="rMedianSignal",Rb="rBGMedianSignal"),
                              arrayvision=,
                              arrayvision.ARM = list(G="ARM Dens - Levels",Gb="Bkgd",R="ARM Dens - Levels",Rb="Bkgd"),
                              arrayvision.MTM = list(G="MTM Dens - Levels",Gb="Bkgd",R="MTM Dens - Levels",Rb="Bkgd"),
                              bluefuse = list(G="AMPCH1",R="AMPCH2"),
                              genepix =,
                              genepix.mean= list(R="F635 Mean",G="F532 Mean",Rb="B635 Median",Gb="B532 Median"),
                              genepix.median = list(R="F635 Median",G="F532 Median",Rb="B635 Median",Gb="B532 Median"),
                              genepix.custom = list(R="F635 Mean",G="F532 Mean",Rb="B635 Mean",Gb="B532 Mean"), #this is the part in the source code where I made the changes. 
                              quantarray = list(R="ch2 Intensity",G="ch1 Intensity",Rb="ch2 Background",Gb="ch1 Background"),
                              imagene9 = list(R="Signal Mean 2",G="Signal Mean 1",Rb="Background Median 2",Gb="Background Median 1"),
                              scanarrayexpress = list(G="Ch1 Mean",Gb="Ch1 B Median",R="Ch2 Mean",Rb="Ch2 B Median"),
                              smd.old = list(G="CH1I_MEAN",Gb="CH1B_MEDIAN",R="CH2I_MEAN",Rb="CH2B_MEDIAN"),
                              smd = list(G="Ch1 Intensity (Mean)",Gb="Ch1 Background (Median)",R="Ch2 Intensity (Mean)",Rb="Ch2 Background (Median)"),
                              spot = list(R="Rmean",G="Gmean",Rb="morphR",Gb="morphG"),
                              spot.close.open = list(R="Rmean",G="Gmean",Rb="morphR.close.open",Gb="morphG.close.open"),
                              NULL
            )
            if(green.only) {
                columns$R <- columns$Rb <- NULL
                nRG <- 1
                E <- FALSE
            } else {
                nRG <- 2
                E <- FALSE
            }
            cnames <- names(columns)
        } else {
            columns <- as.list(columns)
            #		if(!is.list(columns)) stop("columns must be a list")
            cnames <- names(columns)
            if(is.null(cnames)) {
                if(length(columns)==1) {
                    #				Single channel with no background
                    names(columns) <- "E"
                    E <- TRUE
                    nRG <- 0
                } else {
                    stop("columns needs to be a named list")
                }
            } else {
                names(columns)[cnames=="Gf"] <- "G"
                names(columns)[cnames=="Rf"] <- "R"
                cnames <- names(columns)
                nRG <- sum(c("R","G") %in% cnames)
                E <- ("E" %in% cnames)
                if(E && nRG>0) stop("columns can be R,G for two color data, or E for single channel, but not both")
                if(!E && nRG==0) stop("columns must specify foreground G or R or E")
                if(!all(cnames %in% c("G","R","Gb","Rb","E","Eb"))) warning("non-standard columns specified")
            }
        }
        
        if(is.null(annotation)) annotation <- switch(source2,
                                                     agilent = c("Row","Col","Start","Sequence","SwissProt","GenBank","Primate","GenPept","ProbeUID","ControlType","ProbeName","GeneName","SystematicName","Description"),
                                                     arrayvision = c("Spot labels","ID"),
                                                     bluefuse = c("ROW","COL","SUBGRIDROW","SUBGRIDCOL","BLOCK","NAME","ID"),   
                                                     genepix = c("Block","Row","Column","ID","Name"),
                                                     imagene9 = c("Meta Row","Meta Column","Row","Column","Gene ID"),
                                                     quantarray= c("Array Row","Array Column","Row","Column","Name"),
                                                     scanarrayexpress = c("Array Row","Array Column","Spot Row","Spot Column"), 	
                                                     smd = c("Spot","Clone ID","Gene Symbol","Gene Name","Cluster ID","Accession","Preferred name","Locuslink ID","Name","Sequence Type","X Grid Coordinate (within sector)","Y Grid Coordinate (within sector)","Sector","Failed","Plate Number","Plate Row","Plate Column","Clone Source","Is Verified","Is Contaminated","Luid"),
                                                     NULL
        )
        if(source=="smd.old") annotation <- c("SPOT","NAME","Clone ID","Gene Symbol","Gene Name","Cluster ID","Accession","Preferred name","SUID")
        
        #	End checking input arguments
        
        #	Read first file to get nspots
        fullname <- slides[1]
        if(!is.null(path)) fullname <- file.path(path,fullname)
        required.col <- unique(c(annotation,unlist(columns),other.columns))
        text.to.search <- if(is.null(wt.fun)) "" else deparse(wt.fun)
        switch(source2,
               quantarray = {
                   firstfield <- scan(fullname,what="",sep="\t",flush=TRUE,quiet=TRUE,blank.lines.skip=FALSE,multi.line=FALSE,allowEscapes=FALSE)
                   skip <- grep("Begin Data",firstfield)
                   if(length(skip)==0) stop("Cannot find \"Begin Data\" in image output file")
                   nspots <- grep("End Data",firstfield) - skip -2
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,nrows=nspots,flush=TRUE,...)
               }, arrayvision = {
                   skip <- 1
                   cn <- scan(fullname,what="",sep=sep,quote=quote,skip=1,nlines=1,quiet=TRUE,allowEscapes=FALSE)
                   fg <- grep(" Dens - ",cn)
                   if(length(fg) != 2) stop(paste("Cannot find foreground columns in",fullname))
                   bg <- grep("^Bkgd$",cn)
                   if(length(bg) != 2) stop(paste("Cannot find background columns in",fullname))
                   #		Note that entries for columns for ArrayVision are now numeric
                   columns <- list(R=fg[1],Rb=bg[1],G=fg[2],Gb=bg[2])
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,...)
                   #		obj <- read.table(fullname,skip=skip,header=TRUE,sep=sep,quote=quote,stringsAsFactors=FALSE,check.names=FALSE,fill=TRUE,comment.char="",flush=TRUE,...)
                   fg <- grep(" Dens - ",names(obj))
                   bg <- grep("^Bkgd$",names(obj))
                   columns <- list(R=fg[1],Rb=bg[1],G=fg[2],Gb=bg[2])
                   nspots <- nrow(obj)
               }, bluefuse = {
                   skip <- readGenericHeader(fullname,columns=c(columns$G,columns$R))$NHeaderRecords
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,...)
                   nspots <- nrow(obj)
               }, genepix = {
                   h <- readGPRHeader(fullname)
                   if(verbose && source=="genepix.custom") cat("Custom background:",h$Background,"\n")
                   skip <- h$NHeaderRecords
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,...)
                   nspots <- nrow(obj)
               }, imagene9 = {
                   h <- readImaGeneHeader(fullname)
                   skip <- h$NHeaderRecords
                   FD <- h$"Field Dimensions"
                   if(is.null(FD)) stop("Can't find Field Dimensions in ImaGene header")
                   nspots <- sum(apply(FD,1,prod))
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,nrows=nspots,...)
               }, smd = {
                   skip <- readSMDHeader(fullname)$NHeaderRecords
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,...)
                   nspots <- nrow(obj)
               }, {
                   skip <- readGenericHeader(fullname,columns=columns,sep=sep)$NHeaderRecords
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,...)
                   nspots <- nrow(obj)
               })
        
        #	Initialize RG list object (object.size for matrix of NAs is smaller)
        Y <- matrix(NA,nspots,nslides)
        colnames(Y) <- names
        RG <- columns
        for (a in cnames) RG[[a]] <- Y
        if(!is.null(wt.fun)) RG$weights <- Y
        if(is.data.frame(targets)) {
            rownames(targets) <- names
            RG$targets <- targets
        } else {
            RG$targets <- data.frame(FileName=files,row.names=names,stringsAsFactors=FALSE)
        }
        
        #	Set annotation columns
        if(!is.null(annotation)) {
            j <- match(annotation,colnames(obj),0)
            if(any(j>0)) RG$genes <- data.frame(obj[,j,drop=FALSE],check.names=FALSE)
        }
        
        RG$source <- source
        
        #	Set printer layout, if possible
        if(source2=="agilent") {
            if(!is.null(RG$genes$Row) && !is.null(RG$genes$Col)) {
                nr <- length(unique(RG$genes$Row))
                nc <- length(unique(RG$genes$Col))
                if(nspots==nr*nc) RG$printer <- list(ngrid.r=1,ngrid.c=1,nspot.r=nr,nspot.c=nc)
            }
        }
        if(source2=="genepix") {
            if(!is.null(RG$genes$Block) && !is.null(RG$genes$Row) && !is.null(RG$genes$Column)) {
                RG$printer <- getLayout(RG$genes,guessdups=FALSE)
                nblocks <- RG$printer$ngrid.r*RG$printer$ngrid.c
                if(!is.na(nblocks) && (nblocks>1) && !is.null(obj$X)) {
                    blocksize <- RG$printer$nspot.r*RG$printer$nspot.c
                    i <- (1:(nblocks-1))*blocksize
                    ngrid.r <- sum(obj$X[i] > obj$X[i+1]) + 1
                    if(!is.na(ngrid.r) && nblocks%%ngrid.r==0) {
                        RG$printer$ngrid.r <- ngrid.r
                        RG$printer$ngrid.c <- nblocks/ngrid.r
                    } else {
                        warning("Can't determine number of grid rows")
                        RG$printer$ngrid.r <- RG$printer$ngrid.c <- NA
                    }
                }
            }
        }
        if(source2=="imagene9") {
            printer <- list(ngrid.r=FD[1,"Metarows"],ngrid.c=FD[1,"Metacols"],nspot.r=FD[1,"Rows"],nspot.c=FD[1,"Cols"])
            if(nrow(FD)==1) {
                RG$printer <- printer
            } else {
                printer$ngrid.r <- sum(FD[,"Metarows"])
                if(all(printer$ngrid.c==FD[,"Metacols"]) &&
                   all(printer$nspot.r==FD[,"Rows"]) &&
                   all(printer$nspot.c==FD[,"Cols"]) ) RG$printer <- printer
            }
        }
        
        #	Other columns
        if(!is.null(other.columns)) {
            other.columns <- as.character(other.columns)
            j <- match(other.columns,colnames(obj),0)
            if(any(j>0)) {
                other.columns <- colnames(obj)[j]
                RG$other <- list()
                for (j in other.columns) RG$other[[j]] <- Y 
            } else {
                other.columns <- NULL
            }
        }
        
        #	Read remainder of files
        for (i in 1:nslides) {
            if(i > 1) {
                fullname <- slides[i]
                if(!is.null(path)) fullname <- file.path(path,fullname)
                switch(source2,
                       quantarray = {
                           firstfield <- scan(fullname,what="",sep="\t",flush=TRUE,quiet=TRUE,blank.lines.skip=FALSE,multi.line=FALSE,allowEscapes=FALSE)
                           skip <- grep("Begin Data", firstfield)
                       },  arrayvision = {
                           skip <- 1
                       }, genepix = {
                           skip <- readGPRHeader(fullname)$NHeaderRecords
                       }, smd = {
                           skip <- readSMDHeader(fullname)$NHeaderRecords
                       }, {
                           skip <- readGenericHeader(fullname,columns=columns)$NHeaderRecords
                       })
                if(verbose && source=="genepix.custom") cat("Custom background:",h$Background,"\n")
                obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,stringsAsFactors=FALSE,quote=quote,fill=TRUE,nrows=nspots,flush=TRUE,...)
                if(nrow(obj) < nspots) stop("File ",slides[i]," has fewer rows than files previously read")
            }
            for (a in cnames) RG[[a]][,i] <- obj[,columns[[a]]]
            if(!is.null(wt.fun)) RG$weights[,i] <- wt.fun(obj)
            if(!is.null(other.columns)) for (j in other.columns) {
                RG$other[[j]][,i] <- obj[,j] 
            }
            if(verbose) cat(paste("Read",fullname,"\n"))
        }
        if(nRG==1) {
            n <- names(RG)
            n[n=="G"] <- "E"
            n[n=="Gb"] <- "Eb"
            n[n=="R"] <- "E"
            n[n=="Rb"] <- "Eb"
            names(RG) <- n
        }
        if(E || nRG==1)
            new("EListRaw",RG)
        else
            new("RGList",RG)
    }
    
    gprFiles <- list.files(path = gprDirectory, pattern = "\\.gpr$", full.names = TRUE)
    
    #function to extract the raw values 
    raw_values <- function(element, source) {
        gene_names <- source[["genes"]][["Name"]]
        #Extracting the column number from the RGlist. It will later be used to aggregate values from the same spot. 
        col_num <- source[["genes"]][["Column"]]
        
        # adding an if-else statement for file names that are imported in without the \n name
        if ("\n" %in% substr(gene_names, nchar(gene_names), nchar(gene_names))) {
            rownames <- gsub("\n", "", gene_names)
        } else {
            rownames <- gene_names
        }
        #Extracting the sample id from the gpr files 
        # colname <- sub(".*_(\\w+)\\.gpr", "\\1", gprFiles)
        colname <- sub(".*/(.*)\\.gpr", "\\1", gprFiles)
        df <- as.data.frame(source[[element]], row.names = rownames)
        colnames(df) <- colname
        df$Column <- col_num
        df <- df[, c("Column", colname)]
        
        return(df)
    }
    #function to calculate the raw SNR
    SNR_raw <- function(F, B) {
        rows <- rownames(F)
        cols <- colnames(F)
        num_rows <- nrow(F)
        num_cols <- ncol(F)
        
        col_start <- which(cols == "Column")
        
        if (length(col_start) == 0) {
            stop("The 'Column' column is missing in the input data.")
        }
        
        ratio_cols <- !(cols %in% c("RowName", "Column"))
        
        ratio_df <- data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
        rownames(ratio_df) <- rows
        colnames(ratio_df) <- cols
        
        ratio_df$Column <- F$Column
        
        start_col <- max(1, col_start + 1)
        for (i in start_col:num_cols) {
            if (ratio_cols[i]) {
                ratio_df[, i] <- F[, i] / B[, i]
            }
        }
        
        return(ratio_df)
    }
    
    #function to aggregate the SNR
    averaged_values <- function(df) {
        base_rownames <- gsub("\\.\\d+$", "", rownames(df))
        column <- df$Column  
        #Aggregrating the spots 
        base_and_col_aggregated <- aggregate(. ~ base_rownames + column, df, mean)
        
        #Aggregating the aggregates so that we have one SNR value for one protein 
        aggregated_df <- aggregate(. ~ base_rownames, base_and_col_aggregated, mean)
        
        aggregated_df$Column <- NULL
        aggregated_df$column <- NULL
        
        colnames(aggregated_df)[1] <- "RowNames"
        rownames(aggregated_df) <- aggregated_df$RowNames
        aggregated_df$RowNames <- NULL
        
        
        return(aggregated_df)
        
    }
    
    #Reading and creating median matrices 
    if (grepl( "median", snr_type, fixed = TRUE)) {
        print("median")
        medianlist <- read.maimages(files = gprFiles, source = "genepix.median")
        medianlist[["genes"]][["Name"]] <- gsub("\n", "", medianlist[["genes"]][["Name"]])
        
        #Extracting the raw median values 
        if (grepl( "635", snr_type, fixed = TRUE)) {
            print("635")
            
            R_median <- raw_values('R', medianlist)
            Rb_median <- raw_values('Rb', medianlist)
            Rmedian_SNR_raw <- SNR_raw(R_median, Rb_median)
            median_635 <- averaged_values(Rmedian_SNR_raw)
            
            csvname <- gsub("[^[:alnum:]]", "_", snr_type)
            csvname <- paste("expmat_", snr_type, ".csv", sep = "", collapse = NULL)
            csvname_write <- paste(data_dir, csvname, sep = "", collapse = NULL)
            write.csv(median_635, file = csvname_write, row.names = TRUE, quote=FALSE)
            
            return (csvname)
        }
        
        if (grepl( "532", snr_type, fixed = TRUE)) {
            print("532")
            
            G_median <- raw_values('G', medianlist)
            Gb_median <- raw_values('Gb', medianlist)
            Gmedian_SNR_raw <- SNR_raw(G_median, Gb_median)
            median_532 <- averaged_values(Gmedian_SNR_raw)
            
            csvname <- gsub("[^[:alnum:]]", "_", snr_type)
            csvname <- paste("expmat_", snr_type, ".csv", sep = "", collapse = NULL)
            csvname_write <- paste(data_dir, csvname, sep = "", collapse = NULL)
            write.csv(median_532, file = csvname_write, row.names = TRUE, quote=FALSE)
            
            return (csvname)
        }
        
    }
    
    if (grepl( "mean", snr_type, fixed = TRUE)) {
        print("mean")
        
        #reading and creating mean matrices
        meanlist <- read.maimages(files = gprFiles, source = "genepix.custom")
        meanlist[["genes"]][["Name"]] <- gsub("\n", "", meanlist[["genes"]][["Name"]])
        
        #Extracting the raw median values 
        if (grepl( "635", snr_type, fixed = TRUE)) {
            print("635")
            
            R_mean <- raw_values('R', meanlist)
            Rb_mean <- raw_values('Rb', meanlist)
            Rmean_SNR_raw <- SNR_raw(R_mean, Rb_mean)
            mean_635 <- averaged_values(Rmean_SNR_raw)
            
            csvname <- gsub("[^[:alnum:]]", "_", snr_type)
            csvname <- paste("expmat_", snr_type, ".csv", sep = "", collapse = NULL)
            csvname_write <- paste(data_dir, csvname, sep = "", collapse = NULL)
            write.csv(mean_635, file = csvname_write, row.names = TRUE, quote=FALSE)
            
            return (csvname)
        }
        
        if (grepl( "532", snr_type, fixed = TRUE)) {
            print("532")
            
            G_mean <- raw_values('G', meanlist)
            Gb_mean <- raw_values('Gb', meanlist)
            Gmean_SNR_raw <- SNR_raw(G_mean, Gb_mean)
            mean_532 <- averaged_values(Gmean_SNR_raw)
            
            csvname <- gsub("[^[:alnum:]]", "_", snr_type)
            csvname <- paste("expmat_", snr_type, ".csv", sep = "", collapse = NULL)
            csvname_write <- paste(data_dir, csvname, sep = "", collapse = NULL)
            write.csv(mean_532, file = csvname_write, row.names = TRUE, quote=FALSE)
            
            return (csvname)
        }
    }
}

compute_expmat_column <- function(column_input, gprDirectory, data_dir) {
    
    read.maimages <- function(files=NULL,source="generic",path=NULL,ext=NULL,names=NULL,columns=NULL,other.columns=NULL,annotation=NULL,green.only=FALSE,wt.fun=NULL,verbose=TRUE,sep="\t",quote=NULL,...)
        #	Extracts either an RGList object from a set of two-color image analysis output files
        #	or an EListRaw object from a set of one-color files
        #	Gordon Smyth. 
        #	1 Nov 2002.  Last revised 26 Feb 2019.
    {
        
        #	Determine type of input file
        source <- match.arg(source,c("generic","agilent","agilent.mean","agilent.median","arrayvision","arrayvision.ARM","arrayvision.MTM","bluefuse","genepix","genepix.mean","genepix.median","genepix.custom","imagene","imagene9","quantarray","scanarrayexpress","smd.old","smd","spot","spot.close.open"))
        #	source2 is the source type with qualifications removed
        source2 <- strsplit(source,split=".",fixed=TRUE)[[1]][1]
        
        #	ImaGene is special case
        if(source2=="imagene") return(read.imagene(files=files,path=path,ext=ext,names=names,columns=columns,other.columns=other.columns,wt.fun=wt.fun,verbose=verbose,sep=sep,quote=quote,...))
        
        #	Get list of files and associated targets frame
        if(is.null(files)) {
            if(is.null(ext))
                stop("Must specify input files")
            else {
                extregex <- paste("\\.",ext,"$",sep="")
                files <- dir(path=ifelse(is.null(path),".",path),pattern=extregex)
                files <- sub(extregex,"",files)
            }
        } else
            if(is.data.frame(files)) {
                targets <- files
                files <- files$FileName
                if(is.null(files)) stop("targets frame doesn't contain FileName column")
                if(is.null(names)) names <- targets$Label
            } else {
                targets <- NULL
            }
        slides <- as.vector(as.character(files))
        if(!is.null(ext)) slides <- paste(slides,ext,sep=".")
        nslides <- length(slides)
        
        #	Default sample names
        if(is.null(names)) names <- removeExt(files)
        
        #	Default for quote
        if(is.null(quote)) if(source2=="agilent") quote <- "" else quote <- "\""
        
        if(is.null(columns)) {
            if(source2=="generic") stop("must specify columns for generic input")
            columns <- switch(source,
                              agilent.mean = list(G="gMeanSignal",Gb="gBGMedianSignal",R="rMeanSignal",Rb="rBGMedianSignal"),
                              agilent =,
                              agilent.median = list(G="gMedianSignal",Gb="gBGMedianSignal",R="rMedianSignal",Rb="rBGMedianSignal"),
                              arrayvision=,
                              arrayvision.ARM = list(G="ARM Dens - Levels",Gb="Bkgd",R="ARM Dens - Levels",Rb="Bkgd"),
                              arrayvision.MTM = list(G="MTM Dens - Levels",Gb="Bkgd",R="MTM Dens - Levels",Rb="Bkgd"),
                              bluefuse = list(G="AMPCH1",R="AMPCH2"),
                              genepix =,
                              genepix.mean= list(R="F635 Mean",G="F532 Mean",Rb="B635 Median",Gb="B532 Median"),
                              genepix.median = list(R="F635 Median",G="F532 Median",Rb="B635 Median",Gb="B532 Median"),
                              genepix.custom = list(R="F635 Mean",G="F532 Mean",Rb="B635 Mean",Gb="B532 Mean"), #this is the part in the source code where I made the changes
                              # genepix.select = list(R="F635 Mean",G="F532 Mean",Rb="B635 Mean",Gb="B532 Mean"), #for selected columns 
                              quantarray = list(R="ch2 Intensity",G="ch1 Intensity",Rb="ch2 Background",Gb="ch1 Background"),
                              imagene9 = list(R="Signal Mean 2",G="Signal Mean 1",Rb="Background Median 2",Gb="Background Median 1"),
                              scanarrayexpress = list(G="Ch1 Mean",Gb="Ch1 B Median",R="Ch2 Mean",Rb="Ch2 B Median"),
                              smd.old = list(G="CH1I_MEAN",Gb="CH1B_MEDIAN",R="CH2I_MEAN",Rb="CH2B_MEDIAN"),
                              smd = list(G="Ch1 Intensity (Mean)",Gb="Ch1 Background (Median)",R="Ch2 Intensity (Mean)",Rb="Ch2 Background (Median)"),
                              spot = list(R="Rmean",G="Gmean",Rb="morphR",Gb="morphG"),
                              spot.close.open = list(R="Rmean",G="Gmean",Rb="morphR.close.open",Gb="morphG.close.open"),
                              NULL
            )
            if(green.only) {
                columns$R <- columns$Rb <- NULL
                nRG <- 1
                E <- FALSE
            } else {
                nRG <- 2
                E <- FALSE
            }
            cnames <- names(columns)
        } else {
            columns <- as.list(columns)
            cnames <- names(columns)
            if(is.null(cnames)) {
                if(length(columns)==1) {
                    #				Single channel with no background
                    names(columns) <- "E"
                    E <- TRUE
                    nRG <- 0
                } else {
                    stop("columns needs to be a named list")
                }
            }
            
        }
        
        if(is.null(annotation)) annotation <- switch(source2,
                                                     agilent = c("Row","Col","Start","Sequence","SwissProt","GenBank","Primate","GenPept","ProbeUID","ControlType","ProbeName","GeneName","SystematicName","Description"),
                                                     arrayvision = c("Spot labels","ID"),
                                                     bluefuse = c("ROW","COL","SUBGRIDROW","SUBGRIDCOL","BLOCK","NAME","ID"),   
                                                     genepix = c("Block","Row","Column","ID","Name"),
                                                     imagene9 = c("Meta Row","Meta Column","Row","Column","Gene ID"),
                                                     quantarray= c("Array Row","Array Column","Row","Column","Name"),
                                                     scanarrayexpress = c("Array Row","Array Column","Spot Row","Spot Column"), 	
                                                     smd = c("Spot","Clone ID","Gene Symbol","Gene Name","Cluster ID","Accession","Preferred name","Locuslink ID","Name","Sequence Type","X Grid Coordinate (within sector)","Y Grid Coordinate (within sector)","Sector","Failed","Plate Number","Plate Row","Plate Column","Clone Source","Is Verified","Is Contaminated","Luid"),
                                                     NULL
        )
        if(source=="smd.old") annotation <- c("SPOT","NAME","Clone ID","Gene Symbol","Gene Name","Cluster ID","Accession","Preferred name","SUID")
        
        #	End checking input arguments
        
        #	Read first file to get nspots
        fullname <- slides[1]
        if(!is.null(path)) fullname <- file.path(path,fullname)
        required.col <- unique(c(annotation,unlist(columns),other.columns))
        text.to.search <- if(is.null(wt.fun)) "" else deparse(wt.fun)
        switch(source2,
               quantarray = {
                   firstfield <- scan(fullname,what="",sep="\t",flush=TRUE,quiet=TRUE,blank.lines.skip=FALSE,multi.line=FALSE,allowEscapes=FALSE)
                   skip <- grep("Begin Data",firstfield)
                   if(length(skip)==0) stop("Cannot find \"Begin Data\" in image output file")
                   nspots <- grep("End Data",firstfield) - skip -2
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,nrows=nspots,flush=TRUE,...)
               }, arrayvision = {
                   skip <- 1
                   cn <- scan(fullname,what="",sep=sep,quote=quote,skip=1,nlines=1,quiet=TRUE,allowEscapes=FALSE)
                   fg <- grep(" Dens - ",cn)
                   if(length(fg) != 2) stop(paste("Cannot find foreground columns in",fullname))
                   bg <- grep("^Bkgd$",cn)
                   if(length(bg) != 2) stop(paste("Cannot find background columns in",fullname))
                   #		Note that entries for columns for ArrayVision are now numeric
                   columns <- list(R=fg[1],Rb=bg[1],G=fg[2],Gb=bg[2])
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,...)
                   #		obj <- read.table(fullname,skip=skip,header=TRUE,sep=sep,quote=quote,stringsAsFactors=FALSE,check.names=FALSE,fill=TRUE,comment.char="",flush=TRUE,...)
                   fg <- grep(" Dens - ",names(obj))
                   bg <- grep("^Bkgd$",names(obj))
                   columns <- list(R=fg[1],Rb=bg[1],G=fg[2],Gb=bg[2])
                   nspots <- nrow(obj)
               }, bluefuse = {
                   skip <- readGenericHeader(fullname,columns=c(columns$G,columns$R))$NHeaderRecords
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,...)
                   nspots <- nrow(obj)
               }, genepix = {
                   h <- readGPRHeader(fullname)
                   if(verbose && source=="genepix.custom") cat("Custom background:",h$Background,"\n")
                   skip <- h$NHeaderRecords
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,...)
                   nspots <- nrow(obj)
               }, imagene9 = {
                   h <- readImaGeneHeader(fullname)
                   skip <- h$NHeaderRecords
                   FD <- h$"Field Dimensions"
                   if(is.null(FD)) stop("Can't find Field Dimensions in ImaGene header")
                   nspots <- sum(apply(FD,1,prod))
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,nrows=nspots,...)
               }, smd = {
                   skip <- readSMDHeader(fullname)$NHeaderRecords
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,...)
                   nspots <- nrow(obj)
               }, {
                   skip <- readGenericHeader(fullname,columns=columns,sep=sep)$NHeaderRecords
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,...)
                   nspots <- nrow(obj)
               })
        
        print(columns)
        print("2")
        
        #	Initialize RG list object (object.size for matrix of NAs is smaller)
        Y <- matrix(NA,nspots,nslides)
        colnames(Y) <- names
        RG <- columns
        for (a in cnames) RG[[a]] <- Y
        if(!is.null(wt.fun)) RG$weights <- Y
        if(is.data.frame(targets)) {
            rownames(targets) <- names
            RG$targets <- targets
        } else {
            RG$targets <- data.frame(FileName=files,row.names=names,stringsAsFactors=FALSE)
        }
        
        #	Set annotation columns
        if(!is.null(annotation)) {
            j <- match(annotation,colnames(obj),0)
            if(any(j>0)) RG$genes <- data.frame(obj[,j,drop=FALSE],check.names=FALSE)
        }
        
        RG$source <- source
        
        #	Set printer layout, if possible
        if(source2=="agilent") {
            if(!is.null(RG$genes$Row) && !is.null(RG$genes$Col)) {
                nr <- length(unique(RG$genes$Row))
                nc <- length(unique(RG$genes$Col))
                if(nspots==nr*nc) RG$printer <- list(ngrid.r=1,ngrid.c=1,nspot.r=nr,nspot.c=nc)
            }
        }
        if(source2=="genepix") {
            if(!is.null(RG$genes$Block) && !is.null(RG$genes$Row) && !is.null(RG$genes$Column)) {
                RG$printer <- getLayout(RG$genes,guessdups=FALSE)
                nblocks <- RG$printer$ngrid.r*RG$printer$ngrid.c
                if(!is.na(nblocks) && (nblocks>1) && !is.null(obj$X)) {
                    blocksize <- RG$printer$nspot.r*RG$printer$nspot.c
                    i <- (1:(nblocks-1))*blocksize
                    ngrid.r <- sum(obj$X[i] > obj$X[i+1]) + 1
                    if(!is.na(ngrid.r) && nblocks%%ngrid.r==0) {
                        RG$printer$ngrid.r <- ngrid.r
                        RG$printer$ngrid.c <- nblocks/ngrid.r
                    } else {
                        warning("Can't determine number of grid rows")
                        RG$printer$ngrid.r <- RG$printer$ngrid.c <- NA
                    }
                }
            }
        }
        if(source2=="imagene9") {
            printer <- list(ngrid.r=FD[1,"Metarows"],ngrid.c=FD[1,"Metacols"],nspot.r=FD[1,"Rows"],nspot.c=FD[1,"Cols"])
            if(nrow(FD)==1) {
                RG$printer <- printer
            } else {
                printer$ngrid.r <- sum(FD[,"Metarows"])
                if(all(printer$ngrid.c==FD[,"Metacols"]) &&
                   all(printer$nspot.r==FD[,"Rows"]) &&
                   all(printer$nspot.c==FD[,"Cols"]) ) RG$printer <- printer
            }
        }
        
        #	Other columns
        if(!is.null(other.columns)) {
            other.columns <- as.character(other.columns)
            j <- match(other.columns,colnames(obj),0)
            if(any(j>0)) {
                other.columns <- colnames(obj)[j]
                RG$other <- list()
                for (j in other.columns) RG$other[[j]] <- Y 
            } else {
                other.columns <- NULL
            }
        }
        
        
        #	Read remainder of files
        for (i in 1:nslides) {
            if(i > 1) {
                fullname <- slides[i]
                if(!is.null(path)) fullname <- file.path(path,fullname)
                switch(source2,
                       quantarray = {
                           firstfield <- scan(fullname,what="",sep="\t",flush=TRUE,quiet=TRUE,blank.lines.skip=FALSE,multi.line=FALSE,allowEscapes=FALSE)
                           skip <- grep("Begin Data", firstfield)
                       },  arrayvision = {
                           skip <- 1
                       }, genepix = {
                           skip <- readGPRHeader(fullname)$NHeaderRecords
                       }, smd = {
                           skip <- readSMDHeader(fullname)$NHeaderRecords
                       }, {
                           skip <- readGenericHeader(fullname,columns=columns)$NHeaderRecords
                       })
                if(verbose && source=="genepix.custom") cat("Custom background:",h$Background,"\n")
                obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,stringsAsFactors=FALSE,quote=quote,fill=TRUE,nrows=nspots,flush=TRUE,...)
                if(nrow(obj) < nspots) stop("File ",slides[i]," has fewer rows than files previously read")
            }
            for (a in cnames) RG[[a]][,i] <- obj[,columns[[a]]]
            if(!is.null(wt.fun)) RG$weights[,i] <- wt.fun(obj)
            if(!is.null(other.columns)) for (j in other.columns) {
                RG$other[[j]][,i] <- obj[,j] 
            }
            if(verbose) cat(paste("Read",fullname,"\n"))
        }
        
        new("RGList",RG)
        
    }
    
    gprFiles <- list.files(path = gprDirectory, pattern = "\\.gpr$", full.names = TRUE)
    
    #function to extract the raw values 
    raw_values <- function(element, source) {
        gene_names <- source[["genes"]][["Name"]]
        #Extracting the column number from the RGlist. It will later be used to aggregate values from the same spot. 
        col_num <- source[["genes"]][["Column"]]
        
        # adding an if-else statement for file names that are imported in without the \n name
        if ("\n" %in% substr(gene_names, nchar(gene_names), nchar(gene_names))) {
            rownames <- gsub("\n", "", gene_names)
        } else {
            rownames <- gene_names
        }
        #Extracting the sample id from the gpr files 
        # colname <- sub(".*_(\\w+)\\.gpr", "\\1", gprFiles)
        colname <- sub(".*/(.*)\\.gpr", "\\1", gprFiles)
        print(colname)
        df <- as.data.frame(source[[element]], row.names = rownames)
        colnames(df) <- colname
        df$Column <- col_num
        df <- df[, c("Column", colname)]
        
        return(df)
    }
    #function to calculate the raw SNR
    SNR_raw <- function(F, B) {
        rows <- rownames(F)
        cols <- colnames(F)
        num_rows <- nrow(F)
        num_cols <- ncol(F)
        
        col_start <- which(cols == "Column")
        
        if (length(col_start) == 0) {
            stop("The 'Column' column is missing in the input data.")
        }
        
        ratio_cols <- !(cols %in% c("RowName", "Column"))
        
        ratio_df <- data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
        rownames(ratio_df) <- rows
        colnames(ratio_df) <- cols
        
        ratio_df$Column <- F$Column
        
        start_col <- max(1, col_start + 1)
        for (i in start_col:num_cols) {
            if (ratio_cols[i]) {
                ratio_df[, i] <- F[, i] / B[, i]
            }
        }
        
        return(ratio_df)
    }
    
    #function to aggregate the SNR
    averaged_values <- function(df) {
        base_rownames <- gsub("\\.\\d+$", "", rownames(df))
        column <- df$Column  
        #Aggregrating the spots 
        base_and_col_aggregated <- aggregate(. ~ base_rownames + column, df, mean)
        
        #Aggregating the aggregates so that we have one SNR value for one protein 
        aggregated_df <- aggregate(. ~ base_rownames, base_and_col_aggregated, mean)
        
        aggregated_df$Column <- NULL
        aggregated_df$column <- NULL
        
        colnames(aggregated_df)[1] <- "RowNames"
        rownames(aggregated_df) <- aggregated_df$RowNames
        aggregated_df$RowNames <- NULL
        
        
        return(aggregated_df)
        
    }
    
    col_vals <- read.maimages(files = gprFiles, source = "genepix.median", columns=list(column=column_input))
    col_vals_raw <- raw_values('column', col_vals)

    #Calculating the final median count matrix
    col_vals_avg <- averaged_values(col_vals_raw)

    csvname <- gsub("[^[:alnum:]]", "_", column_input)
    csvname <- paste("expmat_", csvname, ".csv", sep = "", collapse = NULL)
    csvname_write <- paste(data_dir, csvname, sep = "", collapse = NULL)
    write.csv(col_vals_avg, file = csvname_write, row.names = TRUE, quote=FALSE)
    
    return (csvname)
    
}


finalcountmatrices <- function() {
    
    gprDirectory <- "./gpr/"
    
    read.maimages <- function(files=NULL,source="generic",path=NULL,ext=NULL,names=NULL,columns=NULL,other.columns=NULL,annotation=NULL,green.only=FALSE,wt.fun=NULL,verbose=TRUE,sep="\t",quote=NULL,...)
        #	Extracts either an RGList object from a set of two-color image analysis output files
        #	or an EListRaw object from a set of one-color files
        #	Gordon Smyth. 
        #	1 Nov 2002.  Last revised 26 Feb 2019.
    {

        #	Determine type of input file
        source <- match.arg(source,c("generic","agilent","agilent.mean","agilent.median","arrayvision","arrayvision.ARM","arrayvision.MTM","bluefuse","genepix","genepix.mean","genepix.median","genepix.custom","imagene","imagene9","quantarray","scanarrayexpress","smd.old","smd","spot","spot.close.open"))
        #	source2 is the source type with qualifications removed
        source2 <- strsplit(source,split=".",fixed=TRUE)[[1]][1]
        
        #	ImaGene is special case
        if(source2=="imagene") return(read.imagene(files=files,path=path,ext=ext,names=names,columns=columns,other.columns=other.columns,wt.fun=wt.fun,verbose=verbose,sep=sep,quote=quote,...))
        
        #	Get list of files and associated targets frame
        if(is.null(files)) {
            if(is.null(ext))
                stop("Must specify input files")
            else {
                extregex <- paste("\\.",ext,"$",sep="")
                files <- dir(path=ifelse(is.null(path),".",path),pattern=extregex)
                files <- sub(extregex,"",files)
            }
        } else
            if(is.data.frame(files)) {
                targets <- files
                files <- files$FileName
                if(is.null(files)) stop("targets frame doesn't contain FileName column")
                if(is.null(names)) names <- targets$Label
            } else {
                targets <- NULL
            }
        slides <- as.vector(as.character(files))
        if(!is.null(ext)) slides <- paste(slides,ext,sep=".")
        nslides <- length(slides)
        
        #	Default sample names
        if(is.null(names)) names <- removeExt(files)
        
        #	Default for quote
        if(is.null(quote)) if(source2=="agilent") quote <- "" else quote <- "\""
        
        if(is.null(columns)) {
            if(source2=="generic") stop("must specify columns for generic input")
            columns <- switch(source,
                              agilent.mean = list(G="gMeanSignal",Gb="gBGMedianSignal",R="rMeanSignal",Rb="rBGMedianSignal"),
                              agilent =,
                              agilent.median = list(G="gMedianSignal",Gb="gBGMedianSignal",R="rMedianSignal",Rb="rBGMedianSignal"),
                              arrayvision=,
                              arrayvision.ARM = list(G="ARM Dens - Levels",Gb="Bkgd",R="ARM Dens - Levels",Rb="Bkgd"),
                              arrayvision.MTM = list(G="MTM Dens - Levels",Gb="Bkgd",R="MTM Dens - Levels",Rb="Bkgd"),
                              bluefuse = list(G="AMPCH1",R="AMPCH2"),
                              genepix =,
                              genepix.mean= list(R="F635 Mean",G="F532 Mean",Rb="B635 Median",Gb="B532 Median"),
                              genepix.median = list(R="F635 Median",G="F532 Median",Rb="B635 Median",Gb="B532 Median"),
                              genepix.custom = list(R="F635 Mean",G="F532 Mean",Rb="B635 Mean",Gb="B532 Mean"), #this is the part in the source code where I made the changes
                              # genepix.select = list(R="F635 Mean",G="F532 Mean",Rb="B635 Mean",Gb="B532 Mean"), #for selected columns 
                              quantarray = list(R="ch2 Intensity",G="ch1 Intensity",Rb="ch2 Background",Gb="ch1 Background"),
                              imagene9 = list(R="Signal Mean 2",G="Signal Mean 1",Rb="Background Median 2",Gb="Background Median 1"),
                              scanarrayexpress = list(G="Ch1 Mean",Gb="Ch1 B Median",R="Ch2 Mean",Rb="Ch2 B Median"),
                              smd.old = list(G="CH1I_MEAN",Gb="CH1B_MEDIAN",R="CH2I_MEAN",Rb="CH2B_MEDIAN"),
                              smd = list(G="Ch1 Intensity (Mean)",Gb="Ch1 Background (Median)",R="Ch2 Intensity (Mean)",Rb="Ch2 Background (Median)"),
                              spot = list(R="Rmean",G="Gmean",Rb="morphR",Gb="morphG"),
                              spot.close.open = list(R="Rmean",G="Gmean",Rb="morphR.close.open",Gb="morphG.close.open"),
                              NULL
            )
            if(green.only) {
                columns$R <- columns$Rb <- NULL
                nRG <- 1
                E <- FALSE
            } else {
                nRG <- 2
                E <- FALSE
            }
            cnames <- names(columns)
        } else {
            columns <- as.list(columns)
            #		if(!is.list(columns)) stop("columns must be a list")
            cnames <- names(columns)
            if(is.null(cnames)) {
                if(length(columns)==1) {
                    #				Single channel with no background
                    names(columns) <- "E"
                    E <- TRUE
                    nRG <- 0
                } else {
                    stop("columns needs to be a named list")
                }
            }
            
        }
        
        print(columns)
        
        if(is.null(annotation)) annotation <- switch(source2,
                                                     agilent = c("Row","Col","Start","Sequence","SwissProt","GenBank","Primate","GenPept","ProbeUID","ControlType","ProbeName","GeneName","SystematicName","Description"),
                                                     arrayvision = c("Spot labels","ID"),
                                                     bluefuse = c("ROW","COL","SUBGRIDROW","SUBGRIDCOL","BLOCK","NAME","ID"),   
                                                     genepix = c("Block","Row","Column","ID","Name"),
                                                     imagene9 = c("Meta Row","Meta Column","Row","Column","Gene ID"),
                                                     quantarray= c("Array Row","Array Column","Row","Column","Name"),
                                                     scanarrayexpress = c("Array Row","Array Column","Spot Row","Spot Column"), 	
                                                     smd = c("Spot","Clone ID","Gene Symbol","Gene Name","Cluster ID","Accession","Preferred name","Locuslink ID","Name","Sequence Type","X Grid Coordinate (within sector)","Y Grid Coordinate (within sector)","Sector","Failed","Plate Number","Plate Row","Plate Column","Clone Source","Is Verified","Is Contaminated","Luid"),
                                                     NULL
        )
        if(source=="smd.old") annotation <- c("SPOT","NAME","Clone ID","Gene Symbol","Gene Name","Cluster ID","Accession","Preferred name","SUID")
        
        #	End checking input arguments
        
        #	Read first file to get nspots
        fullname <- slides[1]
        if(!is.null(path)) fullname <- file.path(path,fullname)
        required.col <- unique(c(annotation,unlist(columns),other.columns))
        text.to.search <- if(is.null(wt.fun)) "" else deparse(wt.fun)
        switch(source2,
               quantarray = {
                   firstfield <- scan(fullname,what="",sep="\t",flush=TRUE,quiet=TRUE,blank.lines.skip=FALSE,multi.line=FALSE,allowEscapes=FALSE)
                   skip <- grep("Begin Data",firstfield)
                   if(length(skip)==0) stop("Cannot find \"Begin Data\" in image output file")
                   nspots <- grep("End Data",firstfield) - skip -2
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,nrows=nspots,flush=TRUE,...)
               }, arrayvision = {
                   skip <- 1
                   cn <- scan(fullname,what="",sep=sep,quote=quote,skip=1,nlines=1,quiet=TRUE,allowEscapes=FALSE)
                   fg <- grep(" Dens - ",cn)
                   if(length(fg) != 2) stop(paste("Cannot find foreground columns in",fullname))
                   bg <- grep("^Bkgd$",cn)
                   if(length(bg) != 2) stop(paste("Cannot find background columns in",fullname))
                   #		Note that entries for columns for ArrayVision are now numeric
                   columns <- list(R=fg[1],Rb=bg[1],G=fg[2],Gb=bg[2])
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,...)
                   #		obj <- read.table(fullname,skip=skip,header=TRUE,sep=sep,quote=quote,stringsAsFactors=FALSE,check.names=FALSE,fill=TRUE,comment.char="",flush=TRUE,...)
                   fg <- grep(" Dens - ",names(obj))
                   bg <- grep("^Bkgd$",names(obj))
                   columns <- list(R=fg[1],Rb=bg[1],G=fg[2],Gb=bg[2])
                   nspots <- nrow(obj)
               }, bluefuse = {
                   skip <- readGenericHeader(fullname,columns=c(columns$G,columns$R))$NHeaderRecords
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,...)
                   nspots <- nrow(obj)
               }, genepix = {
                   h <- readGPRHeader(fullname)
                   if(verbose && source=="genepix.custom") cat("Custom background:",h$Background,"\n")
                   skip <- h$NHeaderRecords
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,...)
                   nspots <- nrow(obj)
               }, imagene9 = {
                   h <- readImaGeneHeader(fullname)
                   skip <- h$NHeaderRecords
                   FD <- h$"Field Dimensions"
                   if(is.null(FD)) stop("Can't find Field Dimensions in ImaGene header")
                   nspots <- sum(apply(FD,1,prod))
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,nrows=nspots,...)
               }, smd = {
                   skip <- readSMDHeader(fullname)$NHeaderRecords
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,...)
                   nspots <- nrow(obj)
               }, {
                   skip <- readGenericHeader(fullname,columns=columns,sep=sep)$NHeaderRecords
                   obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,quote=quote,stringsAsFactors=FALSE,fill=TRUE,flush=TRUE,...)
                   nspots <- nrow(obj)
               })
        
        print(columns)
        print("2")
        
        #	Initialize RG list object (object.size for matrix of NAs is smaller)
        Y <- matrix(NA,nspots,nslides)
        colnames(Y) <- names
        RG <- columns
        for (a in cnames) RG[[a]] <- Y
        if(!is.null(wt.fun)) RG$weights <- Y
        if(is.data.frame(targets)) {
            rownames(targets) <- names
            RG$targets <- targets
        } else {
            RG$targets <- data.frame(FileName=files,row.names=names,stringsAsFactors=FALSE)
        }
        
        #	Set annotation columns
        if(!is.null(annotation)) {
            j <- match(annotation,colnames(obj),0)
            if(any(j>0)) RG$genes <- data.frame(obj[,j,drop=FALSE],check.names=FALSE)
        }
        
        RG$source <- source
        
        #	Set printer layout, if possible
        if(source2=="agilent") {
            if(!is.null(RG$genes$Row) && !is.null(RG$genes$Col)) {
                nr <- length(unique(RG$genes$Row))
                nc <- length(unique(RG$genes$Col))
                if(nspots==nr*nc) RG$printer <- list(ngrid.r=1,ngrid.c=1,nspot.r=nr,nspot.c=nc)
            }
        }
        if(source2=="genepix") {
            if(!is.null(RG$genes$Block) && !is.null(RG$genes$Row) && !is.null(RG$genes$Column)) {
                RG$printer <- getLayout(RG$genes,guessdups=FALSE)
                nblocks <- RG$printer$ngrid.r*RG$printer$ngrid.c
                if(!is.na(nblocks) && (nblocks>1) && !is.null(obj$X)) {
                    blocksize <- RG$printer$nspot.r*RG$printer$nspot.c
                    i <- (1:(nblocks-1))*blocksize
                    ngrid.r <- sum(obj$X[i] > obj$X[i+1]) + 1
                    if(!is.na(ngrid.r) && nblocks%%ngrid.r==0) {
                        RG$printer$ngrid.r <- ngrid.r
                        RG$printer$ngrid.c <- nblocks/ngrid.r
                    } else {
                        warning("Can't determine number of grid rows")
                        RG$printer$ngrid.r <- RG$printer$ngrid.c <- NA
                    }
                }
            }
        }
        if(source2=="imagene9") {
            printer <- list(ngrid.r=FD[1,"Metarows"],ngrid.c=FD[1,"Metacols"],nspot.r=FD[1,"Rows"],nspot.c=FD[1,"Cols"])
            if(nrow(FD)==1) {
                RG$printer <- printer
            } else {
                printer$ngrid.r <- sum(FD[,"Metarows"])
                if(all(printer$ngrid.c==FD[,"Metacols"]) &&
                   all(printer$nspot.r==FD[,"Rows"]) &&
                   all(printer$nspot.c==FD[,"Cols"]) ) RG$printer <- printer
            }
        }
        
        #	Other columns
        if(!is.null(other.columns)) {
            other.columns <- as.character(other.columns)
            j <- match(other.columns,colnames(obj),0)
            if(any(j>0)) {
                other.columns <- colnames(obj)[j]
                RG$other <- list()
                for (j in other.columns) RG$other[[j]] <- Y 
            } else {
                other.columns <- NULL
            }
        }
        
        print(columns)
        print("3")
        
        #	Read remainder of files
        for (i in 1:nslides) {
            if(i > 1) {
                fullname <- slides[i]
                if(!is.null(path)) fullname <- file.path(path,fullname)
                switch(source2,
                       quantarray = {
                           firstfield <- scan(fullname,what="",sep="\t",flush=TRUE,quiet=TRUE,blank.lines.skip=FALSE,multi.line=FALSE,allowEscapes=FALSE)
                           skip <- grep("Begin Data", firstfield)
                       },  arrayvision = {
                           skip <- 1
                       }, genepix = {
                           skip <- readGPRHeader(fullname)$NHeaderRecords
                       }, smd = {
                           skip <- readSMDHeader(fullname)$NHeaderRecords
                       }, {
                           skip <- readGenericHeader(fullname,columns=columns)$NHeaderRecords
                       })
                if(verbose && source=="genepix.custom") cat("Custom background:",h$Background,"\n")
                obj <- read.columns(fullname,required.col,text.to.search,skip=skip,sep=sep,stringsAsFactors=FALSE,quote=quote,fill=TRUE,nrows=nspots,flush=TRUE,...)
                if(nrow(obj) < nspots) stop("File ",slides[i]," has fewer rows than files previously read")
            }
            for (a in cnames) RG[[a]][,i] <- obj[,columns[[a]]]
            if(!is.null(wt.fun)) RG$weights[,i] <- wt.fun(obj)
            if(!is.null(other.columns)) for (j in other.columns) {
                RG$other[[j]][,i] <- obj[,j] 
            }
            if(verbose) cat(paste("Read",fullname,"\n"))
        }
        
        new("RGList",RG)
    
    }
    
    gprFiles <- list.files(path = gprDirectory, pattern = "\\.gpr$", full.names = TRUE)
    
    #function to extract the raw values 
    raw_values <- function(element, source) {
        gene_names <- source[["genes"]][["Name"]]
        #Extracting the column number from the RGlist. It will later be used to aggregate values from the same spot. 
        col_num <- source[["genes"]][["Column"]]
        
        # adding an if-else statement for file names that are imported in without the \n name
        if ("\n" %in% substr(gene_names, nchar(gene_names), nchar(gene_names))) {
            rownames <- gsub("\n", "", gene_names)
        } else {
            rownames <- gene_names
        }
        #Extracting the sample id from the gpr files 
        # colname <- sub(".*_(\\w+)\\.gpr", "\\1", gprFiles)
        colname <- sub(".*/(.*)\\.gpr", "\\1", gprFiles)
        df <- as.data.frame(source[[element]], row.names = rownames)
        colnames(df) <- colname
        df$Column <- col_num
        df <- df[, c("Column", colname)]
        
        return(df)
    }
    #function to calculate the raw SNR
    SNR_raw <- function(F, B) {
        rows <- rownames(F)
        cols <- colnames(F)
        num_rows <- nrow(F)
        num_cols <- ncol(F)
        
        col_start <- which(cols == "Column")
        
        if (length(col_start) == 0) {
            stop("The 'Column' column is missing in the input data.")
        }
        
        ratio_cols <- !(cols %in% c("RowName", "Column"))
        
        ratio_df <- data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
        rownames(ratio_df) <- rows
        colnames(ratio_df) <- cols
        
        ratio_df$Column <- F$Column
        
        start_col <- max(1, col_start + 1)
        for (i in start_col:num_cols) {
            if (ratio_cols[i]) {
                ratio_df[, i] <- F[, i] / B[, i]
            }
        }
        
        return(ratio_df)
    }
    
    #function to aggregate the SNR
    averaged_values <- function(df) {
        base_rownames <- gsub("\\.\\d+$", "", rownames(df))
        column <- df$Column  
        #Aggregrating the spots 
        base_and_col_aggregated <- aggregate(. ~ base_rownames + column, df, mean)
        
        #Aggregating the aggregates so that we have one SNR value for one protein 
        aggregated_df <- aggregate(. ~ base_rownames, base_and_col_aggregated, mean)
        
        aggregated_df$Column <- NULL
        aggregated_df$column <- NULL
        
        colnames(aggregated_df)[1] <- "RowNames"
        rownames(aggregated_df) <- aggregated_df$RowNames
        aggregated_df$RowNames <- NULL
        
        
        return(aggregated_df)
        
    }
    
    #Reading and creating median matrices 
    medianlist <- read.maimages(files = gprFiles, source = "genepix.median")
    medianlist <- read.maimages(files = gprFiles, columns=list(block="Block", column="Column"))
    
    
    medianlist[["genes"]][["Name"]] <- gsub("\n", "", medianlist[["genes"]][["Name"]])
    
    #Extracting the raw median values 
    R_median <- raw_values('R', medianlist)
    Rb_median <- raw_values('Rb', medianlist)
    G_median <- raw_values('G', medianlist)
    Gb_median <- raw_values('Gb', medianlist)
    
    #Calculating the median SNR
    Rmedian_SNR_raw <- SNR_raw(R_median, Rb_median)
    Gmedian_SNR_raw <- SNR_raw(G_median, Gb_median)
    
    #Calculating the final median count matrix
    median_635 <- averaged_values(Rmedian_SNR_raw)
    median_532 <- averaged_values(Gmedian_SNR_raw)
    
    #reading and creating mean matrices
    meanlist <- read.maimages(files = gprFiles, source = "genepix.custom")
    meanlist[["genes"]][["Name"]] <- gsub("\n", "", meanlist[["genes"]][["Name"]])
    
    #Extracting the raw mean values 
    R_mean <- raw_values('R', meanlist)
    Rb_mean <- raw_values('Rb', meanlist)
    G_mean <- raw_values('G', meanlist)
    Gb_mean <- raw_values('Gb', meanlist)
    
    #Calculating the mean SNR
    Rmean_SNR_raw <- SNR_raw(R_mean, Rb_mean)
    Gmean_SNR_raw <- SNR_raw(G_mean, Gb_mean)
    
    #Calculating the final mean count matrix
    mean_635 <- averaged_values(Rmean_SNR_raw)
    mean_532 <- averaged_values(Gmean_SNR_raw)
    
   
    write.csv(median_532, file = "expmat_median532.csv", row.names = TRUE, quote=FALSE)
    write.csv(median_635, file = "expmat_median635.csv", row.names = TRUE, quote=FALSE)
    write.csv(mean_532, file = "expmat_mean532.csv", row.names = TRUE, quote=FALSE)
    write.csv(mean_635, file = "expmat_mean635.csv", row.names = TRUE, quote=FALSE)

}


# snr_type <- "SNR mean 532"
# gprDirectory <- "./microarray_data/"
# data_dir <- "./processed_data/"
# 
# compute_expmat_snr(snr_type, gprDirectory, data_dir)