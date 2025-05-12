library(readr)
library(dplyr)

### Please provide 1. RDC bed files: Refer examples from Peggy
###                2. Reference gene 
read_rdc=function(x){ # x is the data path output from rdc calling 
  data <- read_delim(x,col_names = F)
  colnames(data)=c("Chr","Start","End","RDC_ID","jptm","Strand")
  data
}

read_ref=function(x){
  ref=read_delim(x,col_names = F)
  colnames(ref)=c("Chr","Start","End","Gene_ID","None","Strand")
  ref=ref[!grepl("_rev",ref$Gene_ID),]
  ref
}


produce_gtf=function(data,ref,fname,wd){

  rdc_gene_rpkm=data.frame(GeneID=character(),Chr=character(),Start=integer(),End=integer(),Strand=character())
  rdc_rpkm=data.frame(GeneID=character(),Chr=character(),Start=integer(),End=integer(),Strand=character())
  for(i in 1:nrow(data)){
    d=data[i,]
    d$Start=as.numeric(d$Start)
    d$End=as.numeric(d$End)
    corr_gene=ref[which(ref$Chr==d$Chr[1]),] %>% dplyr::filter(.,between(Start,d$Start[1],d$End[1])) %>% 
      dplyr::filter(.,between(End,d$Start[1],d$End[1]))
    
    if(nrow(corr_gene)!=0){
      for(k in 1:nrow(corr_gene)){
        gene=corr_gene$Gene_ID[k]
        id=paste0(d$RDC_ID,"_",gene)
        temp_d=data.frame(GeneID=id,Chr=corr_gene$Chr[k],Start=corr_gene$Start[k],End=corr_gene$End[k],Strand=corr_gene$Strand[k])
        rdc_gene_rpkm=rbind(rdc_gene_rpkm,temp_d)
      }
      
    }else{
      id=paste0(d$RDC_ID,"_NULL")
      temp_d=data.frame(GeneID=id,Chr=d$Chr[1],Start=d$Start[1],End=d$End[1],Strand=".")
      rdc_rpkm=rbind(rdc_rpkm,temp_d)
      
    }
  }
  write_delim(rdc_gene_rpkm,paste0(wd,"result/","Gene_",fname,".gtf"),delim = "\t")
  write_delim(rdc_rpkm,paste0(wd,"result/","RDC_",fname,".gtf"),delim = "\t")
  rdc_gene_rpkm
}


calculate_rpkm=function(bam_path,rdc_input,refPath,wd){
  
  rdc_files=list.files(rdc_input,pattern = "*.bed$")
  samples=list.files(bam_path,pattern = "*.bam$")
  samples=samples[!grepl("pos.bam",samples) & !grepl("neg.bam",samples)]
  
  for(i in 1:length(samples)){
    print(samples[i])
    s=gsub(".bam","",samples[i])
    cmd=paste0("/software/SAMtools/1.20-GCC-14.1.0/bin/samtools view -c ",bam_path,samples[i]," > ","result/",s,"_totalreads.txt")
    system(cmd)
    for(j in 1:length(rdc_files)){
      f=rdc_files[j]
      data=read_rdc(paste0(rdc_input,f))
      ref=read_ref(refPath)

      fname=gsub(".bed","",f)
      produce_gtf(data=data,ref=ref,fname=fname,wd=wd)

      gtf_name=paste0("Gene_",fname,".gtf")
      cmd=paste0("/software/Subread/2.0.6-GCC-14.1.0/bin/featureCounts -O -T 15 -F SAF -t gene -s 1 -a ","result/",gtf_name," -o ","result/Gene_",fname,"_",s,"_featureCounts_output.txt ",bam_path,samples[i])
      system(cmd)
      
      gtf_name=paste0("RDC_",fname,".gtf")
      cmd=paste0("/software/Subread/2.0.6-GCC-14.1.0/bin/featureCounts -O -T 15 -F SAF -t gene -a ","result/",gtf_name," -o ","result/RDC_",fname,"_",s,"_featureCounts_output.txt ",bam_path,samples[i])
      system(cmd)
    }
  } 
  
  files=list.files(paste0(wd,"result/"),pattern = "_featureCounts_output.txt$")
  
  for(i in 1:length(rdc_files)){
    f=rdc_files[i]
    f=gsub(".bed","",f)
    print(paste0("Process: ",f))
    targets=files[grepl(paste0(f,"_"),files)]
    
    for(j in 1:length(targets)){
      r=targets[j]
      t=gsub(paste0(f,"_"),"",r)
      t=gsub("Gene_","",t)
      t=gsub("RDC_","",t)
      t=gsub("_featureCounts_output.txt","",t)
      
      print(paste0("Map to ",t))
      
      total_counts=read_csv(paste0(paste0(wd,"result/"),t,"_totalreads.txt"), col_names = FALSE)$X1
      
      data=read_delim(paste0(wd,"result/",r), delim="\t",comment = "#")
      colnames(data)[7]="read_counts"
      
      
      out=gsub("_featureCounts_output.txt","",r)
      data["RPKM"]=(data$read_counts/data$Length*1000)/(total_counts/1000000)
      print(paste0("Write ",out))
      write.csv(data,paste0(wd,"result/",out,"_rpkm.csv"),row.names = F)
      
    }
    
    
    
  }
  
  print("Done")
}

