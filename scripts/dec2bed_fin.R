suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
options(scipen=999)

args <- commandArgs(trailingOnly = T)
if (length(args) < 3) {
  warning("Some required files are missing, please provide: decompositon-raw-file windows2analize-bed-file out-file")
  quit(status = 1)
}

dec_file <- args[1]
bed_file <- args[2]
out_file <- args[3]

dec=read.table(dec_file)

dec <- dec %>% mutate(end=V4+1) %>% select(V1,V2,V3,end)
colnames(dec) <- c("seq","match","rel_start","rel_end")

str2split <- as.character(dec$seq[1])

chr=strsplit(str2split,split = ":")[[1]][1]
pos=strsplit(str2split,split = ":")[[1]][2]

start_piece=as.numeric(strsplit(pos,split = "-")[[1]][1])
end_piece=as.numeric(strsplit(pos,split = "-")[[1]][2]) 

dec <- dec %>% mutate(chr4head=chr) %>%   mutate(start4head=rel_start+start_piece) %>%  
  mutate(end4head=rel_end+start_piece) %>% mutate(name= paste0(chr4head,":",start4head, "-", end4head))

StringDec <- makeGRangesFromDataFrame(dec, start.field = "start4head", end.field = "end4head", seqnames.field = "seq",
                                      ignore.strand = T, keep.extra.columns = T)

bed=read.table(bed_file)
bed <- bed %>% filter(V1==chr) %>% mutate(V1=(unique(dec$seq)))

bed <- makeGRangesFromDataFrame(bed, start.field = "V2", end.field = "V3", seqnames.field = "V1",keep.extra.columns = T)  
 
StringDec_filtered <- subsetByOverlaps(StringDec, bed, type = "within")

StringDec_filtered <-  as.data.frame(StringDec_filtered)

StringDec_filtered <- StringDec_filtered %>% mutate(s= case_when(endsWith(match, "'") ~ "-", .default = "+"))

StringDec_filtered <- StringDec_filtered %>% mutate(score=".") %>% select(seqnames,rel_start,rel_end,name, score, s)

write.table(StringDec_filtered, sep="\t", quote = F, row.names = F, col.names = F,
            file = out_file)
