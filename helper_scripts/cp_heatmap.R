library(sequenza)
#library(scarHRD)
library(stringr)


## Get HRD values for all possible sequenza solutions
## Call:
## Rscript /mnt/data/Novaseq/.software/cp_heatmap.R {input.seqz} {output.lpp} | grep -v "HRD" | sed 's/\s\+/ /g' > {output.hrd}

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10000)

# Load custom scarHRD 
files.sources = list.files(path = '/mnt/storage2/sequencing/software/WES/scarHRD-master/', pattern='.R$', full.names=T, recursive=T)
none <- capture.output( none <- sapply(files.sources, source))
load('/mnt/storage2/sequencing/software/WES/scarHRD-master/R/sysdata.rda')

# sequenza input
args = commandArgs(trailingOnly = TRUE)
seqz_file = as.character(args[1])


# Sequenza calls to get lpp scores 

extract<-sequenza.extract(seqz_file, chromosome.list=paste('chr',c(1:24),sep=''),gamma = 60, kmin = 50)#, parallel=10)

cp <- sequenza::sequenza.fit(extract, N.ratio.filter = 10, N.BAF.filter = 1 , segment.filter = 3e6, mufreq.treshold = 0.10, ratio.priority = FALSE,ploidy=seq(1, 5.5, 0.1), mc.cores = 10)

# write lpp values to a file for heatmap plotting
write.table(cp$lpp, file=as.character(args[2]))


output_dir = tempdir()


# Go through each ploidy and cellularity combination
for (ploid in seq(1, 5.5, 0.1)){

    seg.tab <- do.call(rbind, extract$segments[extract$chromosomes])
    seg.len <- (seg.tab$end.pos - seg.tab$start.pos)/1e+06


    for (cell in seq(0.05, 1, 0.05)){
	allele.cn = sequenza:::baf.bayes(Bf = seg.tab$Bf, CNt.max = 20
						 , depth.ratio = seg.tab$depth.ratio
						 , avg.depth.ratio = 1,
					  cellularity = cell
					  , ploidy = ploid
					  , sd.ratio = seg.tab$sd.ratio
					  , weight.ratio = seg.len
					  , sd.Bf = seg.tab$sd.BAF,
					  weight.Bf = 1, ratio.priority = FALSE, CNn = 2)
	allele.cn <- as.data.frame(allele.cn)
	seg <- data.frame(SampleID = "NoID"
		, Chromosome = seg.tab$chromosome
		, Start_position = seg.tab$start.pos
		, End_position = seg.tab$end.pos
		, total_cn = allele.cn$CNt
		, A_cn = allele.cn$B
		,B_cn = allele.cn$A
		, ploidy = ploid)

	seg<-seg[!is.na(seg$A_cn),]
	seg<-seg[!is.na(seg$B_cn),]


        # Write sequenza data bc scar hrd wants to read a file, not be given a data frame
        output.seg = paste0(output_dir,"/sequenza.out")
	write.table(x = seg, file = output.seg, sep = '\t' ) 
        capture.output( r <- scar_score(seg = output.seg, reference='grch37',  seqz=F,outputdir='/tmp'), file = nullfile())

        df = as.data.frame.matrix(r)
        df$cell <- cell
        df$ploidy <- ploid

        print(df,row.names=FALSE) 
    } 
}
