library(data.table)
library(sequenza)
library(stringr)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 20000)

files.sources = list.files(path = '/mnt/storage2/sequencing/software/WES/scarHRD-master/', pattern='.R$', full.names=T, recursive=T)
none <- capture.output(none <- sapply(files.sources, source))
load('/mnt/storage2/sequencing/software/WES/scarHRD-master/R/sysdata.rda')

args <- commandArgs(trailingOnly = TRUE)
seqz_file = as.character(args[1])

sample_name = basename(seqz_file)
sample_name = sub('\\.bin50.seqz.gz', '', sample_name)
sample_name = sub('\\.bin200.seqz.gz', '', sample_name)

output_dir=paste0(dirname(dirname(seqz_file)),'/Results/', sample_name)
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive=TRUE)}

#usetwd(output_dir)

scar_score(seqz_file,reference = "grch37", seqz=TRUE, outputdir = output_dir)#,estimate_tzg=.85)

