ReCalculatePvalue <- function(Example.Go.adjusted.by.exon) {
  pwf=Example.Go.adjusted.by.exon[[2]]
  pwf=unfactor(pwf)
  test.cats=c("GO:CC","GO:BP","GO:MF")
  gene2cat=getgo(rownames(pwf),'hg19','ensGene',fetch.cats=test.cats)
  names(gene2cat)=rownames(pwf)
  cat2gene=reversemapping(gene2cat)
  gene2cat=reversemapping(cat2gene)
  nafrac=(sum(is.na(pwf$pwf))/nrow(pwf))*100
  if(nafrac>50){
    warning(paste("Missing length data for ",round(nafrac),"% of genes.  Accuarcy of GO test will be reduced.",sep=''))
  }
  pwf$pwf[is.na(pwf$pwf)]=pwf$pwf[match(sort(pwf$bias.data[!is.na(pwf$bias.data)])[ceiling(sum(!is.na(pwf$bias.data))/2)],pwf$bias.data)]
  
  unknown_go_terms=nrow(pwf)-length(gene2cat)
  
  if(unknown_go_terms>0 ){
    message(paste("For",unknown_go_terms,"genes, we could not find any categories. These genes will be excluded."))
    message("To force their use, please run with use_genes_without_cat=TRUE (see documentation).")
    message("This was the default behavior for version 1.15.1 and earlier.")
    pwf=pwf[rownames(pwf) %in% names(gene2cat),]
  } 
  #A few variables are always useful so calculate them
  cats=names(cat2gene)
  DE=rownames(pwf)[pwf$DEgenes==1]
  
  #total number of DE genes
  num_de=length(DE)
  
  #total number of genes
  num_genes=nrow(pwf)
  
  pvals=data.frame(category=cats,over_represented_pvalue=NA,under_represented_pvalue=NA,stringsAsFactors=FALSE,numDEInCat=NA,numInCat=NA)
  
  degenesnum=which(pwf$DEgenes==1)
  #Turn all genes into a reference to the pwf object
  cat2genenum=relist(match(unlist(cat2gene),rownames(pwf)),cat2gene)
  #This value is used in every calculation, by storing it we need only calculate it once
  alpha=sum(pwf$pwf)
  
  pvals[,2:3]=t(sapply(cat2genenum,function(u){
    
    print(u)
    
    #The number of DE genes in this category
    num_de_incat=sum(degenesnum%in%u)
    #The total number of genes in this category
    num_incat=length(u)
    #This is just a quick way of calculating weight=avg(PWF within category)/avg(PWF outside of category)
    avg_weight=mean(pwf$pwf[u])
    weight=(avg_weight*(num_genes-num_incat))/(alpha-num_incat*avg_weight)
    if(num_incat==num_genes){ weight=1 } #case for the root GO terms			
    #Now calculate the sum of the tails of the Wallenius distribution (the p-values)
    
    
    c(dWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight)
      +pWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight,lower.tail=FALSE),
      pWNCHypergeo(num_de_incat,num_incat,num_genes-num_incat,num_de,weight))
  }))
  
  return(pvals)
  
}



