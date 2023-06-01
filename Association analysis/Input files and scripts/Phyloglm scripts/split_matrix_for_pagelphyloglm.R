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
# test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("Two argument must be supplied a matrix and number of chunks to create", call.=FALSE)
}
mat<-args[1]
numsplits<-args[2]
prefix<-args[3]

#######Commensal and Pathogenic######
#Read the Map #Change the path if needed
#Map<-read.table("/pine/scr/i/s/isai/3837_march_2017/metadata_3837_genomes_june2016_complete.tsv",header=T,sep="\t")
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

trait<-c(rep(1,length(Map_pathogenic)),rep(0,length(Map_commensal)))
df_traits<-data.frame(TraitY=trait)
rownames(df_traits)<-Map$taxon_oid

#Create directory to contain the split files
outdir<-paste("split_matrix_pathogenic_commensal_",prefix,"/",sep="")
dir.create(outdir)
#Number of splits
numsplits<-as.numeric(numsplits)
span=numsplits-1
num_paral=0;
start <- seq(from=1,to=ncol(Tab), by = numsplits)
#Loop over the matrix splitting it
for(index in start){
	num_paral=num_paral+1
	outfile=paste(outdir,"splitted_file_",num_paral,".tsv",sep="")
	top<-span + index
	if(top >ncol(Tab)){
		top<-ncol(Tab)
		subTab<-Tab[,index:top]
		subTab<-cbind(df_traits,subTab)
		write.table(file=outfile,x=subTab,quote=F,row.names=T,col.names=T,sep="\t",append=F)
	}else{
		subTab<-Tab[,index:top]
		subTab<-cbind(df_traits,subTab)
		write.table(file=outfile,x=subTab,quote=F,row.names=T,col.names=T,sep="\t",append=F)

	}
}

#######Weaklypathogenic and Pathogenic######
rm(Map,Tab,trait,df_traits,outdir,span,num_paral,start)
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
#Remove cols that sum to zero
Tab<-Tab[,which(colSums(Tab)!=0)]

trait<-c(rep(1,length(Map_pathogenic)),rep(0,length(Map_weaklypathogenic)))
df_traits<-data.frame(TraitY=trait)
rownames(df_traits)<-Map$taxon_oid

#Create directory to contain the split files
outdir<-paste("split_matrix_pathogenic_weaklypathogenic_",prefix,"/",sep="")
dir.create(outdir)
#Number of splits
numsplits=args[2]
numsplits<-as.numeric(numsplits)
span=numsplits-1
num_paral=0;
start <- seq(from=1,to=ncol(Tab), by = numsplits)
#Loop over the matrix splitting it
for(index in start){
        num_paral=num_paral+1
        outfile=paste(outdir,"splitted_file_",num_paral,".tsv",sep="")
        top<-span + index
        if(top >ncol(Tab)){
                top<-ncol(Tab)
                subTab<-Tab[,index:top]
                subTab<-cbind(df_traits,subTab)
                write.table(file=outfile,x=subTab,quote=F,row.names=T,col.names=T,sep="\t",append=F)
        }else{
                subTab<-Tab[,index:top]
                subTab<-cbind(df_traits,subTab)
                write.table(file=outfile,x=subTab,quote=F,row.names=T,col.names=T,sep="\t",append=F)

        }
}

####### commensal vs weaklypathogenic
rm(Map,Tab,trait,df_traits,outdir,span,num_paral,start)
#Read the Map #Change the path if needed
#Map<-read.table("/pine/scr/i/s/isai/3837_march_2017/metadata_3837_genomes_june2016_complete.tsv",header=T,sep="\t")
Map<-read.table("metadata_COG_337xantho.txt",header=T,sep="\t")
#Read the matrix file passed
Tab<-read.table(mat,header=T,sep="\t",row.names=1,check.names=F)
Tab<-t(Tab)
#The logic is to subset first using the columns of the matrix
Map<-Map[match(rownames(Tab),Map$taxon_oid),]
#Here we are  removing soil
Map_commensal<-as.character(Map$taxon_oid[which(Map$Classification=="Commensal")])
Map_weaklypathogenic<-as.character(Map$taxon_oid[which(Map$Classification=="Weaklypathogenic")])
myids<-c(Map_commensal,Map_weaklypathogenic)
Map<-Map[match(myids,Map$taxon_oid),]
Map<-droplevels(Map)
Tab<-Tab[match(Map$taxon_oid,rownames(Tab)),]
#Remove zero that sum to zero
Tab<-Tab[,which(colSums(Tab)!=0)]

trait<-c(rep(1,length(Map_commensal)),rep(0,length(Map_weaklypathogenic)))
df_traits<-data.frame(TraitY=trait)
rownames(df_traits)<-Map$taxon_oid

#Create directory to contain the split files
outdir<-paste("split_matrix_commensal_weaklypathogenic_",prefix,"/",sep="")
dir.create(outdir)
#Number of splits
numsplits<-as.numeric(numsplits)
span=numsplits-1
num_paral=0;
start <- seq(from=1,to=ncol(Tab), by = numsplits)
#Loop over the matrix splitting it
for(index in start){
	num_paral=num_paral+1
	outfile=paste(outdir,"splitted_file_",num_paral,".tsv",sep="")
	top<-span + index
	if(top >ncol(Tab)){
		top<-ncol(Tab)
		subTab<-Tab[,index:top]
		subTab<-cbind(df_traits,subTab)
		write.table(file=outfile,x=subTab,quote=F,row.names=T,col.names=T,sep="\t",append=F)
	}else{
		subTab<-Tab[,index:top]
		subTab<-cbind(df_traits,subTab)
		write.table(file=outfile,x=subTab,quote=F,row.names=T,col.names=T,sep="\t",append=F)

	}
}
