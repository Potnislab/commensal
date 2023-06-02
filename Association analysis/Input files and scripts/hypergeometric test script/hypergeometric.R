# (C) Copyright 2017 Isai Salas Gonzalez
#
#    This file is part of gfobap.
#
#    gfobap is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    gfobap is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with gfobap.  If not, see <http://www.gnu.org/licenses/>.

args = commandArgs(trailingOnly=TRUE)
if (length(args)<1) {
  stop("Two argument must be supplied a matrix ", call.=FALSE)
}
mat<-args[1]

######Commensal versus Pathogenic#####
suppressMessages(library(phytools))
print(mat)
print ("Commensal versus Pathogenic")
#Read the Map #Change the path if needed
#Map<-read.table("/path",header=T,sep="\t")
Map<-read.table("metadata_COG_337xantho.txt",header=T,sep="\t")
#Read the matrix file passed
Tab<-read.table(mat,header=T,sep="\t",row.names=1,check.names=F)
Tab<-t(Tab)
#The logic is to subset first using the columns of the matrix
Map<-Map[match(rownames(Tab),Map$taxon_oid),]
#Here we are  removing soil
Map_commensal<-as.character(Map$taxon_oid[which(Map$Classification=="Commensal")])
Map_pathogenic<-as.character(Map$taxon_oid[which(Map$Classification=="Pathogenic")])
myids<-c(Map_commensal,Map_pathogenic)
Map<-Map[match(myids,Map$taxon_oid),]
Map<-droplevels(Map)
Tab<-Tab[match(Map$taxon_oid,rownames(Tab)),]
#Remove zero that sum to zero
Tab<-Tab[,which(colSums(Tab)!=0)]

#Read the phylogenetic tree and root if 
#tree<-read.tree("path/nwk")
tree<-read.tree("new_ortho_tree.nwk")

allgenes_hyp <- function(tree,genes,Map){
  #tree <- Dat$tree
  #Map <- Dat$Map
  #genes <- Dat$genes >=1
  
  N.commensal <- sum(genes[ Map$Classification == "Commensal", ])
  N.pathogenic <- sum(genes[ Map$Classification == "Pathogenic", ])
  #N.genomes <- nrow(Map)
  Res <- NULL
  for(gene in 1:ncol(genes)){
    #gene <- 1
    print (gene)
    gene.id <- colnames(genes)[gene]
    Map$gene <- genes[,gene]
    
    #genes.per.class <- aggregate(gene ~ Classification, data = Map,FUN = sum)
    #genomes.per.class <- aggregate(gene ~ Classification, data = Map,FUN = length)
    #N.genes <- sum(Map$gene)
    
    #ftable(gene ~ Classification , Map)
    
    #dhyper(x = 0, m = 2, n = 2, k = 2)
    #dhyper(x = 1, m = 2, n = 2, k = 2)
    #dhyper(x = 2, m = 2, n = 2, k = 2)
    #phyper(q = 0, m = 2, n = 2, k = 2)
    #phyper(q = 1, m = 2, n = 2, k = 2)
    #phyper(q = 2, m = 2, n = 2, k = 2)
    # Following is wrong, counts genomes for zero and genes for 1+
    # pval <- phyper(q = sum(subset(Map, Classification == "Commensal" & gene > 0)$gene) - 1,
    #                m = sum(Map$gene),
    #                n = nrow(subset(Map,gene == 0)),
    #                k = sum(subset(Map, Classification == "Commensal")$gene))
    #pval <- 1 - pval
    
    # Binary version
    pval <- phyper(q = nrow(subset(Map, Classification == "Commensal" & gene > 0)) - 1,
                   m = sum(Map$gene > 0),
                   n = nrow(subset(Map,gene == 0)),
                   k = nrow(subset(Map, Classification == "Commensal")),lower.tail=FALSE)
    #pval <- 1 - pval
    #pval[ pval == 0 ] <- 1e-16
    score <- -log10(pval)
    
    pval2 <- phyper(q = sum(subset(Map, Classification == "Commensal")$gene) - 1, 
                    m = N.commensal, n = N.pathogenic, k = sum(Map$gene),lower.tail=FALSE)
    #pval2 <- 1 - pval2
    #pval2[ pval2 == 0 ] <- 1e-16 
    
    
    #Test for depletion binary version 
    pval_dep <- phyper(q = nrow(subset(Map, Classification == "Commensal" & gene > 0)),
                   m = sum(Map$gene > 0),
                   n = nrow(subset(Map,gene == 0)),
                   k = nrow(subset(Map, Classification == "Commensal")),
                   lower.tail=TRUE)
    #pval_dep[ pval_dep == 0 ] <- 1e-16
    score_dep<--log10(pval_dep)

    #Test for depletion Raw counts version 
    pval2_dep <- phyper(q = sum(subset(Map, Classification == "Commensal")$gene),
                    m = N.commensal, n = N.pathogenic, k = sum(Map$gene),lower.tail=TRUE)
    #pval2_dep[ pval2_dep == 0 ] <- 1e-16
    score2_dep<--log10(pval2_dep)



    #res <- data.frame(gene.id = gene.id, score = score, p.value = pval,
    #                  full.score = -log10(pval2), full.p.value = pval2)
    res <- data.frame(gene.id = gene.id, score_enriched_binary = score, p.value_enriched_binary = pval,
                      score_depletion_binary=score_dep,p.value_depletion_binary=pval_dep,
                      score_enriched_rawcounts = -log10(pval2), p.value_enriched_rawcounts = pval2,
                      score_depletion_rawcounts=score2_dep,p.value_depletion_rawcounts=pval2_dep)

    Map$gene <- NULL
    Res <- rbind(Res,res)
  }
  
  Res<-data.frame(gene.id=Res$gene.id,score_enriched_binary=Res$score_enriched_binary,
  z.score_enriched_binary= (Res$score_enriched_binary - mean(Res$score_enriched_binary)) / sd(Res$score_enriched_binary),
  p.value_enriched_binary=Res$p.value_enriched_binary,
  score_depletion_binary=Res$score_depletion_binary,
  z.score_depletion_binary=(Res$score_depletion_binary - mean(Res$score_depletion_binary)) / sd(Res$score_depletion_binary),
  p.value_depletion_binary=Res$p.value_depletion_binary,
  score_enriched_rawcounts=Res$score_enriched_rawcounts,
  z.score_enriched_rawcounts=(Res$score_enriched_rawcounts - mean(Res$score_enriched_rawcounts)) / sd(Res$score_enriched_rawcounts),
  p.value_enriched_rawcounts=Res$p.value_enriched_rawcounts,
  score_depletion_rawcounts=Res$score_depletion_rawcounts,
  z.score_depletion_rawcounts=(Res$score_depletion_rawcounts - mean(Res$score_depletion_rawcounts)) / sd(Res$score_depletion_rawcounts),
  p.value_depletion_rawcounts=Res$p.value_depletion_rawcounts)



  #Res <- data.frame(gene.id = Res$gene.id, score = Res$score,
  #                  z.score = (Res$score - mean(Res$score)) / sd(Res$score),
  #                  p.value = Res$p.value,
  #                  full.score = Res$full.score,
  #                  full.z.score = (Res$full.score - mean(Res$full.score)) / sd(Res$full.score),
  #                  full.p.value = Res$full.p.value)
  return(Res)
}

finRes<-allgenes_hyp(tree,Tab,Map)
write.table(finRes,file = "hyp_res_commensal_pathogenic.txt", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE,append=FALSE)

######Weaklypathogenic versus pathogenic######
#Read the Map #Change the path if needed
rm(Map,finRes,allgenes_hyp,Tab,tree)
print("Weaklypathogenic versus pathogenic")
#Map<-read.table("/pine/scr/i/s/isai/3837_march_2017/metadata_3837_genomes_june2016_complete.tsv",header=T,sep="\t")
Map<-read.table("metadata_COG_337xantho.txt",header=T,sep="\t")
#Read the matrix file passed
Tab<-read.table(mat,header=T,sep="\t",row.names=1,check.names=F)
Tab<-t(Tab)
#The logic is to subset first using the columns of the matrix
Map<-Map[match(rownames(Tab),Map$taxon_oid),]
#Subset soil and root at the same time
Map_weaklypathogenic<-as.character(Map$taxon_oid[which(Map$Classification=="Weaklypathogenic")])
Map_pathogenic<-as.character(Map$taxon_oid[which(Map$Classification=="Pathogenic")])
myids<-c(Map_weaklypathogenic,Map_pathogenic)
Map<-Map[match(myids,Map$taxon_oid),]
Map<-droplevels(Map)
Tab<-Tab[match(Map$taxon_oid,rownames(Tab)),]
#Remove zero that sum to zero
Tab<-Tab[,which(colSums(Tab)!=0)]

#Read the phylogenetic tree and root if 
#tree<-read.tree("/pine/scr/i/s/isai/3837_march_2017/3837_genomes_31scg_june2016.newick")
tree<-read.tree("new_ortho_tree.nwk")


allgenes_hyp <- function(tree,genes,Map){
  #tree <- Dat$tree
  #Map <- Dat$Map
  #genes <- Dat$genes >=1

  N.weaklypathogenic <- sum(genes[ Map$Classification == "Weaklypathogenic", ])
  N.pathogenic <- sum(genes[ Map$Classification == "Pathogenic", ])
  #N.genomes <- nrow(Map)
  Res <- NULL
  for(gene in 1:ncol(genes)){
    #gene <- 1
    print (gene)
    gene.id <- colnames(genes)[gene]
    Map$gene <- genes[,gene]
    #genes.per.class <- aggregate(gene ~ Classification, data = Map,FUN = sum)
    #genomes.per.class <- aggregate(gene ~ Classification, data = Map,FUN = length)
    #N.genes <- sum(Map$gene)

    #ftable(gene ~ Classification , Map)

    #dhyper(x = 0, m = 2, n = 2, k = 2)
    #dhyper(x = 1, m = 2, n = 2, k = 2)
    #dhyper(x = 2, m = 2, n = 2, k = 2)
    #phyper(q = 0, m = 2, n = 2, k = 2)
    #phyper(q = 1, m = 2, n = 2, k = 2)
    #phyper(q = 2, m = 2, n = 2, k = 2)
    # Following is wrong, counts genomes for zero and genes for 1+
    # pval <- phyper(q = sum(subset(Map, Classification == "Weaklypathogenic" & gene > 0)$gene) - 1,
    #                m = sum(Map$gene),
    #                n = nrow(subset(Map,gene == 0)),
    #                k = sum(subset(Map, Classification == "Weaklypathogenic")$gene))
    #pval <- 1 - pval

    # Binary version
    #Rest for enrichment
    pval <- phyper(q = nrow(subset(Map, Classification == "Weaklypathogenic" & gene > 0)) - 1,
                   m = sum(Map$gene > 0),
                   n = nrow(subset(Map,gene == 0)),
                   k = nrow(subset(Map, Classification == "Weaklypathogenic")),lower.tail=FALSE)
    #pval <- 1 - pval
    #pval[ pval == 0 ] <- 1e-16
    score <- -log10(pval)
    #Raw Counts version
    #Test for enrichment
    pval2 <- phyper(q = sum(subset(Map, Classification == "Weaklypathogenic")$gene) - 1,
                    m = N.weaklypathogenic, n = N.pathogenic, k = sum(Map$gene),lower.tail=FALSE)
    #pval2 <- 1 - pval2
    #pval2[ pval2 == 0 ] <- 1e-16
    
    #Test for depletion binary version 
    pval_dep <- phyper(q = nrow(subset(Map, Classification == "Weaklypathogenic" & gene > 0)),
                   m = sum(Map$gene > 0),
                   n = nrow(subset(Map,gene == 0)),
                   k = nrow(subset(Map, Classification == "Weaklypathogenic")),
		   lower.tail=TRUE)
    #pval_dep[ pval_dep == 0 ] <- 1e-16
    score_dep<--log10(pval_dep)

    #Test for depletion Raw counts version 
    pval2_dep <- phyper(q = sum(subset(Map, Classification == "Weaklypathogenic")$gene),
                    m = N.weaklypathogenic, n = N.pathogenic, k = sum(Map$gene),lower.tail=TRUE)
    #pval2_dep[ pval2_dep == 0 ] <- 1e-16
    score2_dep<--log10(pval2_dep)
    


    res <- data.frame(gene.id = gene.id, score_enriched_binary = score, p.value_enriched_binary = pval,
		      score_depletion_binary=score_dep,p.value_depletion_binary=pval_dep,
                      score_enriched_rawcounts = -log10(pval2), p.value_enriched_rawcounts = pval2,
		      score_depletion_rawcounts=score2_dep,p.value_depletion_rawcounts=pval2_dep)
    Map$gene <- NULL
    Res <- rbind(Res,res)
  }

  Res<-data.frame(gene.id=Res$gene.id,score_enriched_binary=Res$score_enriched_binary,
  z.score_enriched_binary= (Res$score_enriched_binary - mean(Res$score_enriched_binary)) / sd(Res$score_enriched_binary),
  p.value_enriched_binary=Res$p.value_enriched_binary,
  score_depletion_binary=Res$score_depletion_binary,
  z.score_depletion_binary=(Res$score_depletion_binary - mean(Res$score_depletion_binary)) / sd(Res$score_depletion_binary),
  p.value_depletion_binary=Res$p.value_depletion_binary,
  score_enriched_rawcounts=Res$score_enriched_rawcounts,
  z.score_enriched_rawcounts=(Res$score_enriched_rawcounts - mean(Res$score_enriched_rawcounts)) / sd(Res$score_enriched_rawcounts),
  p.value_enriched_rawcounts=Res$p.value_enriched_rawcounts,
  score_depletion_rawcounts=Res$score_depletion_rawcounts,
  z.score_depletion_rawcounts=(Res$score_depletion_rawcounts - mean(Res$score_depletion_rawcounts)) / sd(Res$score_depletion_rawcounts),
  p.value_depletion_rawcounts=Res$p.value_depletion_rawcounts)

  #Res <- data.frame(gene.id = Res$gene.id, score = Res$score,
  #                  z.score = (Res$score - mean(Res$score)) / sd(Res$score),
  #                  p.value = Res$p.value,
  #                  full.score = Res$full.score,
  #                  full.z.score = (Res$full.score - mean(Res$full.score)) / sd(Res$full.score),
  #                  full.p.value = Res$full.p.value)
  return(Res)
}
finRes<-allgenes_hyp(tree,Tab,Map)
write.table(finRes,file = "hyp_res_weaklypathogenic_pathogenic.txt", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE,append=FALSE)


#######Commensal versus Weaklypathogenic######
#Read the Map #Change the path if needed
rm(Map,finRes,allgenes_hyp,Tab,tree)
print("Commensal versus Weaklypathogenic")
#Map<-read.table("/pine/scr/i/s/isai/3837_march_2017/metadata_3837_genomes_june2016_complete.tsv",header=T,sep="\t")
Map<-read.table("metadata_COG_337xantho.txt",header=T,sep="\t")
#Read the matrix file passed
Tab<-read.table(mat,header=T,sep="\t",row.names=1,check.names=F)
Tab<-t(Tab)
#The logic is to subset first using the columns of the matrix
Map<-Map[match(rownames(Tab),Map$taxon_oid),]
#Subset soil and root at the same time
Map_commensal<-as.character(Map$taxon_oid[which(Map$Classification=="Commensal")])
Map_weaklypathogenic<-as.character(Map$taxon_oid[which(Map$Classification=="Weaklypathogenic")])
myids<-c(Map_commensal, Map_weaklypathogenic)
Map<-Map[match(myids,Map$taxon_oid),]
Map<-droplevels(Map)
Tab<-Tab[match(Map$taxon_oid,rownames(Tab)),]
#Remove zero that sum to zero
Tab<-Tab[,which(colSums(Tab)!=0)]


#Read the phylogenetic tree and root if
#tree<-read.tree("/pine/scr/i/s/isai/3837_march_2017/3837_genomes_31scg_june2016.newick")
tree<-read.tree("new_ortho_tree.nwk")


allgenes_hyp <- function(tree,genes,Map){
  #tree <- Dat$tree
  #Map <- Dat$Map
  #genes <- Dat$genes >=1

  N.commensal <- sum(genes[ Map$Classification == "Commensal", ])
  N.weaklypathogenic <- sum(genes[ Map$Classification == "Weaklypathogenic", ])
  #N.genomes <- nrow(Map)
  Res <- NULL
  for(gene in 1:ncol(genes)){
    #gene <- 1
print (gene)
    gene.id <- colnames(genes)[gene]
    Map$gene <- genes[,gene]
    #genes.per.class <- aggregate(gene ~ Classification, data = Map,FUN = sum)
    #genomes.per.class <- aggregate(gene ~ Classification, data = Map,FUN = length)
    #N.genes <- sum(Map$gene)

    #ftable(gene ~ Classification , Map)

    #dhyper(x = 0, m = 2, n = 2, k = 2)
    #dhyper(x = 1, m = 2, n = 2, k = 2)
    #dhyper(x = 2, m = 2, n = 2, k = 2)
    #phyper(q = 0, m = 2, n = 2, k = 2)
    #phyper(q = 1, m = 2, n = 2, k = 2)
#phyper(q = 2, m = 2, n = 2, k = 2)
    # Following is wrong, counts genomes for zero and genes for 1+
    # pval <- phyper(q = sum(subset(Map, Classification == "Commensal" & gene > 0)$gene) - 1,
    #                m = sum(Map$gene),
    #                n = nrow(subset(Map,gene == 0)),
    #                k = sum(subset(Map, Classification == "Commensal")$gene))
    #pval <- 1 - pval

    # Binary version
    #Rest for enrichment
    pval <- phyper(q = nrow(subset(Map, Classification == "Commensal" & gene > 0)) - 1,
                   m = sum(Map$gene > 0),
                   n = nrow(subset(Map,gene == 0)),
                   k = nrow(subset(Map, Classification == "Commensal")),lower.tail=FALSE)


 #pval_dep[ pval_dep == 0 ] <- 1e-16
    score_dep<--log10(pval_dep)

    #Test for depletion Raw counts version
    pval2_dep <- phyper(q = sum(subset(Map, Classification == "Commensal")$gene),
                    m = N.commensal, n = N.weaklypathogenic, k = sum(Map$gene),lower.tail=TRUE)
    #pval2_dep[ pval2_dep == 0 ] <- 1e-16
    score2_dep<--log10(pval2_dep)



    res <- data.frame(gene.id = gene.id, score_enriched_binary = score, p.value_enriched_binary = pval,
                      score_depletion_binary=score_dep,p.value_depletion_binary=pval_dep,
                      score_enriched_rawcounts = -log10(pval2), p.value_enriched_rawcounts = pval2,
                      score_depletion_rawcounts=score2_dep,p.value_depletion_rawcounts=pval2_dep)
    Map$gene <- NULL
    Res <- rbind(Res,res)
  }
Res<-data.frame(gene.id=Res$gene.id,score_enriched_binary=Res$score_enriched_binary,
  z.score_enriched_binary= (Res$score_enriched_binary - mean(Res$score_enriched_binary)) / sd(Res$score_enriched_binary),
  p.value_enriched_binary=Res$p.value_enriched_binary,
  score_depletion_binary=Res$score_depletion_binary,
  z.score_depletion_binary=(Res$score_depletion_binary - mean(Res$score_depletion_binary)) / sd(Res$score_depletion_binary),
  p.value_depletion_binary=Res$p.value_depletion_binary,
  score_enriched_rawcounts=Res$score_enriched_rawcounts,
  z.score_enriched_rawcounts=(Res$score_enriched_rawcounts - mean(Res$score_enriched_rawcounts)) / sd(Res$score_enriched_rawcounts),
  p.value_enriched_rawcounts=Res$p.value_enriched_rawcounts,
  score_depletion_rawcounts=Res$score_depletion_rawcounts,
  z.score_depletion_rawcounts=(Res$score_depletion_rawcounts - mean(Res$score_depletion_rawcounts)) / sd(Res$score_depletion_rawcounts),
  p.value_depletion_rawcounts=Res$p.value_depletion_rawcounts)

  #Res <- data.frame(gene.id = Res$gene.id, score = Res$score,
  #                  z.score = (Res$score - mean(Res$score)) / sd(Res$score),
  #                  p.value = Res$p.value,
  #                  full.score = Res$full.score,
  #                  full.z.score = (Res$full.score - mean(Res$full.score)) / sd(Res$full.score),
  #                  full.p.value = Res$full.p.value)
  return(Res)
}
finRes<-allgenes_hyp(tree,Tab,Map)
write.table(finRes,file = "hyp_res_commensal_weaklypathogenic.txt", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE,append=FALSE)


