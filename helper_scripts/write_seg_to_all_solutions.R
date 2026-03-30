library(sequenza)
#library(scarHRD)
library(stringr)

options("width"=2000)
	
files.sources = list.files(path = '/mnt/storage2/sequencing/software/WES/scarHRD-master/', pattern='.R$', full.names=T, recursive=T)
none <- capture.output( none <- sapply(files.sources, source))
load('/mnt/storage2/sequencing/software/WES/scarHRD-master/R/sysdata.rda')

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10000)

args = commandArgs(trailingOnly = TRUE)
seqz_file = as.character(args[1])

sample_name = basename(seqz_file)
sample_name = sub('\\.bin50.seqz.gz', '', sample_name)

output_dir=paste0('Results/',sample_name, "_alternative_solutions")
if (!dir.exists(output_dir)) {dir.create(output_dir, recursive=TRUE)}

tryCatch(expr = {
extract<-sequenza.extract(seqz_file, chromosome.list=paste('chr',c(1:24, 'X'),sep=''), gamma = 60, kmin = 50) 
},
error = function(e){
	print("Error in sequenza extract, wrong sequenza version? ")
	print(e)
	quit(status=1)
})

ploidy01 = seq(1, 5.6, 0.1)

extract.fit<-sequenza::sequenza.fit(extract, N.ratio.filter = 10, N.BAF.filter = 1
                                    , segment.filter = 3e6, mufreq.treshold = 0.10
                                    #, ratio.priority = FALSE,ploidy=ploidy01, mc.cores = 10, cellularity = seq(0.8,0.9,0.01))
                                    , ratio.priority = FALSE,ploidy=ploidy01, mc.cores = 10, cellularity = seq(0.1,1,0.01))
temp.dir = tempdir()

sequenza.results(extract, extract.fit , out.dir = temp.dir ,sample.id = 'Temp')

seg.tab <- do.call(rbind, extract$segments[extract$chromosomes])
seg.len <- (seg.tab$end.pos - seg.tab$start.pos)/1e+06
cint <- get.ci(extract.fit)

ploidy <- cint$max.ploidy
altsolutions.file = paste0(temp.dir,'/',"Temp_alternative_solutions.txt")
altsolutions <- read.delim( altsolutions.file , stringsAsFactors=FALSE)

#print(altsolutions)


setwd(output_dir) 
##################################

write_all_seg_files = function(x) {
	allele.cn = sequenza:::baf.bayes(Bf = seg.tab$Bf, CNt.max = 20
						 , depth.ratio = seg.tab$depth.ratio
						 , avg.depth.ratio = 1,
					  cellularity = altsolutions$cellularity[x]
					  , ploidy = altsolutions$ploidy[x],
					  sd.ratio = seg.tab$sd.ratio
					  , weight.ratio = seg.len
					  , sd.Bf = seg.tab$sd.BAF,
					  weight.Bf = 1, ratio.priority = FALSE, CNn = 2)
	allele.cn <- as.data.frame(allele.cn)

	seg <- data.frame(SampleID = as.character(sample_name)
		, Chromosome = seg.tab$chromosome
		, Start_position = seg.tab$start.pos
		, End_position = seg.tab$end.pos
		, total_cn = allele.cn$CNt
		, A_cn = allele.cn$A
		, B_cn = allele.cn$B
		, ploidy = ploidy)

	seg<-seg[!is.na(seg$A_cn),]
	seg<-seg[!is.na(seg$B_cn),]
 

	output.seg = paste(sample_name, altsolutions$cellularity[x],'.txt', sep = '_')
	write.table(x = seg, file = output.seg, sep = '\t' )

	capture.output( r <- scar_score(seg = output.seg, reference='grch37',  seqz=F))#,estimate_tzg=.85)) 
	#scar_score(seg = output.seg, reference='grch37',  seqz=F)

	output.hrd = gsub(x = list.files(pattern='HRDresults')[1], pattern='_HRDresults' , replacement=paste0('__',altsolutions$cellularity[x],'-HRD'))

	if (x == 1){
		output.hrd = gsub(x = list.files(pattern='HRDresults')[1] , pattern='_HRDresults' , replacement=paste0('__main_',altsolutions$cellularity[x],'-HRD')) 
	}

	file.rename(from = list.files(pattern='HRDresults')[1], to = output.hrd)

	df = as.data.frame.matrix(r)

        df$ploidy <- altsolutions[x,]$ploidy
        df$cellularity <- altsolutions[x,]$cellularity
        df$SLPP <- altsolutions[x,]$SLPP

        df$name <- sample_name

	print(df, row.names=FALSE)


}

foo <- sapply(1:nrow(altsolutions), write_all_seg_files)
