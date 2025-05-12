library(optparse)
source("/home/l538g/workingf/GRO_seq/Groseq_result/RDC_RPKM/boyu_hg38/0_utils.r")


option_list = list(
  make_option(c("-o", "--work_dir"), type="character", default=NULL, help="working directory", metavar="character"),
  make_option(c("-r", "--genome"), type="character", default=NULL, help="genome reference", metavar="character"),
  make_option(c("-b", "--bam_files"), type="character", default=NULL, help="the folder of bam files", metavar="character"),
  make_option(c("-i", "--rdc_input"), type="character", default=NULL, help="the folder of rdc bed files", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

setwd(opt$work_dir)
dir.create("result",showWarnings = F)
calculate_rpkm(bam_path=opt$bam_files,rdc_input=opt$rdc_input,refPath=opt$genome,wd=opt$work_dir)
