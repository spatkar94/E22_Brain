
#steps encapsulated for end-to-end processing of un-normalized readcount data to co-expression tables

process_read_counts <- function(dataset){
    #take mapped mRNA read counts of each gene in each cell and
    #filter out genes that were detected in fewer than 5 cells. These are likely to be dropout genes
    dataset = dataset[rowSums(dataset > 0)>5,]
    
    #nomralization:
    
    #First, normalize for gene length bias that may arise from the process of read mapping.
    #Recent protocols use UMI (unique molecular identifiers) and may be free of this bias.
    #https://www.ncbi.nlm.nih.gov/pubmed/28529717. Nonetheless added as a precautionary step
    
    common = intersect(rownames(mmusculus_genes),rownames(dataset))
    dataset = dataset[common,]
    mmusculus_genes = mmusculus_genes[common,]
    
    coverage = colSums(dataset,na.rm = TRUE)
    dataset = apply(dataset, 2, function(x) (x*1000)/mmusculus_genes$length)
    
    #Second, normalize for differences in cell to cell read coverage
    dataset.norm = t(apply(dataset, 1, function(x) (x*1000000)/coverage))
    
    #finally log transform the RPKM values to get a normal like distribution
    mat = log10(as.matrix(dataset.norm)+1)
    
    #glimpse of the processed data
    print(mat[1:5,1:5])
    return(mat)
}

#find cluster-specific expression (measured as zscores)
cluster_zscores <- function(mat,celltypes){
    clusters = unique(celltypes)
    
    #center the expression data of gene and normalize by standard deviation to get zscores for each gene
    zmat = apply(mat,1,function(x) (x - mean(x))/sd(x))
    
    #compute average zscore per cluster
    zmat = aggregate.data.frame(zmat,by = list(celltypes),FUN = mean)
    
    rownames(zmat) = zmat[,1]
    zmat = zmat[,-1]
    print(t(zmat[1:5,1:5]))
    return(t(zmat))
}

#find genes coexpressed with TCF4
intra_cluster_correlation <- function(mat,celltypes){
    #begin co-expression computation...
    clusters = unique(celltypes)
    
    #create 2 empty matrices (gene x cluster); one for storing correlation values, other for pvalues.
    #store them in a list
    res = list()
    res$cor = matrix(0,nrow = nrow(mat),ncol = length(clusters))
    res$pval = matrix(0,nrow = nrow(mat),ncol = length(clusters))
    colnames(res$cor) = clusters
    rownames(res$cor) = rownames(mat)
    colnames(res$pval) = clusters
    rownames(res$pval) = rownames(mat)
    
    #begin pairwise correlation computation.
    #note to myself:
    #pearson correlation test tests for linear relationship between two variables
    #spearman correlation test tests for monotonic relationship between two variables (not necessarily linear)
    #pearson correlation test is a parametric test (null distribution of rho is assumed normal)
    #spearman correlation is non-parametric test.
    #I have switched to measuring spearman correlation since it is a non-parametric test and does not make
    #any assumptions on the form of the null distribution
    for(x in clusters){
        idx = which(celltypes == x)
        cat("Processing cluster: ",x,"\n")
        start = Sys.time()
        cor.dat = NULL
        
        #set as 0 if less than 5 observations exist
        if(length(idx) < 5){
            cor.dat = matrix(0,nrow = nrow(mat),2)
            next
        }
        
        cor.dat = t(apply(mat,1,function(y) {
            cc = cor.test(as.matrix(mat["Tcf4",idx]),y[idx],method="spearman")
            return(c(cc$estimate,cc$p.value))
        }))
        res$cor[,x] = cor.dat[,1]
        res$pval[,x] = cor.dat[,2]
        cat(format(Sys.time()-start,digits = 3)," seconds taken\n")
    }
    return(res)

}


#main
load('./linarson_adult.rda')
load('./linarson_adult_clusters.rda')
load('./mmusculus_genes.rda')

#process read count data
mat = suppressWarnings(process_read_counts(linarson_adult))

#measure cluster specific relative expression
cl_exp = suppressWarnings(cluster_zscores(mat,linarson_adult_clusters))

#measure cluster specific coexpression
coexp = suppressWarnings(intra_cluster_correlation(mat,linarson_adult_clusters))

#FDR correction of p-values
coexp$pval = apply(coexp$pval,2,function(x) p.adjust(x,method = "fdr"))

#glimpse of results
print(coexp$cor[1:5,1:5])
print(coexp$pval[1:5,1:5])

#Add gene filtering code here...

