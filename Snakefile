
import os
import glob

#configfile: "config_files/config.yaml"

DATASET1=config["dataset1"]
DATASET2=config["dataset2"]

out_dir=config["output_dir"]
barcode=config["barcode"]
bam_normal=config["bam_normal"]

fastq_files = sorted(glob.glob(DATASET1 + "/*.fastq.gz") + glob.glob(DATASET2 + "/*.fastq.gz"))

print("Output dir:"+out_dir)
print(barcode)
print(len(fastq_files))


rule all:
    input:
        out_dir + f"/BAMs/{barcode}.bam", #minimap
        out_dir + f"/BAMs/{barcode}.filtered.sorted.bam", #filtered reads
        out_dir + f"/qc_results/NanoPlot-report.html", #qc 
        out_dir + f"/CuteSV/{barcode}.cuteSV.vcf", #structural variants
        out_dir + f"/Sniffles/{barcode}.sniffles.vcf", #structural variants
        #out_dir + f"/CuteSV/{barcode}.cuteSV.annotated_vep.vcf", #structural variant annotation
        #out_dir + f"/Sniffles/{barcode}.sniffles.annotated_vep.vcf", #structural variant annotation
        out_dir + f"/Seqz_files/{barcode}.bin50.seqz.gz", #sequenza file segmentation
        out_dir + f"/Results/{barcode}/{barcode}-HRD.txt", #hrd
        out_dir + f"/Results/{barcode}_alternative_solutions/{barcode}-HRD_alternative_solutions.txt", #hrd alternative solutions
        out_dir + f"/Results/{barcode}_alternative_solutions/{barcode}_all_hrd_solutions", #hrd heatpmap alternative solutions
        out_dir + f"/Results/{barcode}_alternative_solutions/{barcode}_cp_heatmap.pdf" #heatmap plot alternatives solutions


rule minimap:
    input:
        fastqs=fastq_files
    output:
        bam=out_dir + "/BAMs/{barcode}.bam"
    conda:
        "RNAseq_norm"
    priority:
        15
    threads: 
        32
    shell:
        """
        minimap2 -ax map-ont /home/alejandro/hg19/hg19.fa -t {threads} {input.fastqs} > aligned.sam
        samtools view -@ 16 -Sb aligned.sam > aligned.bam
        samtools sort -@ 16 -o {output.bam} aligned.bam
        samtools index {output.bam}
        rm aligned.sam aligned.bam
        """

rule filter_bam_with_qc:
    input:
        bam=out_dir + f"/BAMs/{barcode}.bam"
    
    output:
        filtered_bam=out_dir + f"/BAMs/{barcode}.filtered.sorted.bam"

    conda:
        "base"

    priority: 
        14

    threads: 
        50

    shell:
        """
        samtools view -@ {threads} -h -q 20 -F 0x900 -b {input.bam} | samtools sort -@ {threads} -o {output.filtered_bam}
        samtools index -@ {threads} {output.filtered_bam}
        """

rule nanoplot:
    input:
        bam=out_dir + f"/BAMs/{barcode}.filtered.sorted.bam"
    output:
        report=out_dir + f"/qc_results/NanoPlot-report.html",
        qc_dir=directory(out_dir + "/qc_results")
    conda:
        "RNAseq_norm"
    priority:
        13
    threads: 
        10
    shell:
        """
        NanoPlot --huge --threads 10 --bam {input.bam} -o {output.qc_dir} 
        """

rule cuteSV:
    input:
        bam=out_dir + f"/BAMs/{barcode}.filtered.sorted.bam"
    output:
        vcf=out_dir + f"/CuteSV/{barcode}.cuteSV.vcf"
    conda:
        "Nanovar_seq"
    threads: 
        10
    priority:
        12
    params:
        folder=out_dir + f"/CuteSV/"
    shell:
        """
        cuteSV {input.bam} /home/alejandro/hg19/hg19.fa {output.vcf} {params.folder} --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3
        """

rule sniffles:
    input:
        bam=out_dir + f"/BAMs/{barcode}.filtered.sorted.bam"
    
    output:
        vcf=out_dir + f"/Sniffles/{barcode}.sniffles.vcf"
    
    conda:
        "Nanovar_seq"

    threads: 
        10
    priority:
        11
    shell:
        """
        sniffles --input {input.bam} --reference /home/alejandro/hg19/hg19.fa -v {output.vcf} --threads {threads}
        """

# rule cuteSV_vep:
#     input:
#         vcf=out_dir + f"/CuteSV/{barcode}.cuteSV.vcf"
#     output:
#         annotated_sv=out_dir + f"/CuteSV/{barcode}.cuteSV.annotated_vep.vcf"
#     threads:
#         12
#     priority:
#         9
#     shell:
#         """
#         vep -i {input.vcf} -o {output.annotated_sv} --offline --fork {threads} --vcf --fasta /staging/hg19.fa --format vcf --hgvs --cache --dir /mnt/storage2/sequencing/software/WES/vep/hg19 --force_overwrite --flag_pick --refseq --use_transcript_ref --no_stats --uniprot --numbers --domains --af_gnomade --pubmed --tsl --variant_class --custom file=/mnt/storage2/sequencing/software/WES/vep/clinvar_20231021.vcf.gz,short_name=ClinVar,format=vcf,type=exact,fields=ID%CLNSIG --custom file=/mnt/storage2/sequencing/software/WES/vep/00-common_all.vcf.gz,short_name=COSMIC,format=vcf,type=exact,fields=ID
#         """

# rule sniffles_vep:
#     input:
#         vcf=out_dir + f"/Sniffles/{barcode}.sniffles.vcf"
#     output:
#         annotated_sv=out_dir + f"/Sniffles/{barcode}.sniffles.annotated_vep.vcf"
#     threads:
#         12
#     priority:
#         8
#     shell:
#         """
#         vep -i {input.vcf} -o {output.annotated_sv} --offline --fork {threads} --vcf --fasta /staging/hg19.fa --format vcf --hgvs --cache --dir /mnt/storage2/sequencing/software/WES/vep/hg19 --force_overwrite --flag_pick --refseq --use_transcript_ref --no_stats --uniprot --numbers --domains --af_gnomade --pubmed --tsl --variant_class --custom file=/mnt/storage2/sequencing/software/WES/vep/clinvar_20231021.vcf.gz,short_name=ClinVar,format=vcf,type=exact,fields=ID%CLNSIG --custom file=/mnt/storage2/sequencing/software/WES/vep/00-common_all.vcf.gz,short_name=COSMIC,format=vcf,type=exact,fields=ID
#        """


rule sequenza_bam2seqz:
    input:
        tumor=out_dir + f"/BAMs/{barcode}.filtered.sorted.bam"
    params:
        normal=bam_normal
    output:
        output=out_dir + f"/Seqz_files/{barcode}.bin50.seqz.gz",
    conda:
        "RNAseq_norm"
    priority:
        7
    threads:
        23
    shell:
        """
        touch BAMs/*.bam.bai
        sequenza-utils bam2seqz -n {params.normal} -t {input.tumor} -gc /mnt/storage2/sequencing/software/WES/genome_gc50.wig.gz -F  /home/alejandro/hg19/hg19.fa -C chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr21 chr22 chrX --output Seqz_files/{barcode}.seqz.gz  --parallel {threads} -T tabix > /dev/null
        /mnt/storage1/malibu/Share/Alejandro/FF_vs_FFPE_pilot_study/Nanopore_analysis/Scripts/bam2seqz.sh {barcode} {output}
        """

rule hrd:
    input:
        seqz=out_dir + f"/Seqz_files/{barcode}.bin50.seqz.gz"
    output:
        out_dir + f"/Results/{barcode}/{barcode}-HRD.txt",
        out_dir + f"/Results/{barcode}/{barcode}_segments.txt",
        out_dir + f"/Results/{barcode}/{barcode}_genome_view.pdf"
    conda:
        "HRD2"
    priority:
        6
    threads:
        10
    shell:
        """
        Rscript /mnt/storage2/sequencing/software/WES/script_HRD.R {input.seqz}
        """

rule hrd_alternative_solutions:
     input:
         seqz=out_dir + f"/Seqz_files/{barcode}.bin50.seqz.gz"
     output:
         out_dir + f"/Results/{barcode}_alternative_solutions/{barcode}-HRD_alternative_solutions.txt"
     conda:
         "RNAseq_norm"
     threads:
         10
     priority:
         5
     shell:
        """
        echo "LOH TAI LST HRDsum ploidy cellularity SLPP name" > {output}
        Rscript /mnt/storage2/sequencing/software/WES/write_seg_to_all_solutions.R {input.seqz} | awk '{{$1=$1;print}}' | grep -v character | grep -v HRD-sum >> {output}
        """

rule heatmap_data:
    input:
        alt_sol=out_dir + f"/Results/{barcode}_alternative_solutions/{barcode}-HRD_alternative_solutions.txt",
        seqz=out_dir + f"/Seqz_files/{barcode}.bin50.seqz.gz"
    output:
        lpp=out_dir + f"/Results/{barcode}_alternative_solutions/{barcode}.lpp",
        hrd=out_dir + f"/Results/{barcode}_alternative_solutions/{barcode}_all_hrd_solutions"
    conda:
        "RNAseq_norm"
    priority:
         4
    shell:
        """
        Rscript /mnt/storage2/sequencing/software/WES/cp_heatmap.R {input.seqz} {output.lpp} | grep -v "HRD-sum" | sed 's/\s\+/ /g' > {output.hrd}
        """

rule heatmap_plot:
    input:
        alt_sol=out_dir + f"/Results/{barcode}_alternative_solutions/{barcode}-HRD_alternative_solutions.txt",
        lpp=out_dir + f"/Results/{barcode}_alternative_solutions/{barcode}.lpp",
        hrd=out_dir + f"/Results/{barcode}_alternative_solutions/{barcode}_all_hrd_solutions"
    conda:
        "RNAseq_norm"
    output:
        pdf=out_dir + f"/Results/{barcode}_alternative_solutions/{barcode}_cp_heatmap.pdf",
    shell:
        """
        python /mnt/storage2/sequencing/software/WES/plot_cp_hrd.py {input.alt_sol} {input.lpp} {input.hrd} {output.pdf}
        """



#Annotate vep hard-filtere somatic variants
#vep -i R-2013-2166_ONT.hard_filterd.vcf.gz -o R-2013-2166_ONT.hard_filterd.annotated.vcf --offline --fork 12 --vcf --fasta /staging/hg19.fa --format vcf --hgvs --cache --dir /mnt/storage2/sequencing/software/WES/vep/hg19 --force_overwrite --flag_pick --refseq --use_transcript_ref --no_stats --uniprot --numbers --domains --af_gnomade --pubmed --tsl --variant_class --custom file=/mnt/storage2/sequencing/software/WES/vep/clinvar_20231021.vcf.gz,short_name=ClinVar,format=vcf,type=exact,fields=ID%CLNSIG --custom file=/mnt/storage2/sequencing/software/WES/vep/00-common_all.vcf.gz,short_name=COSMIC,format=vcf,type=exact,fields=ID
#python /mnt/storage2/sequencing/software/WES/filter_VEP.py  {output.mutect2_annotated_vcf}  {output.mutect2_annotated_csv}

#"zcat {input} | java -jar /staging/software/snpEff/SnpSift.jar intervals /mnt/storage2/sequencing/software/WES/keimbahn_targets.bed | java -jar /staging/software/snpEff/SnpSift.jar filter \" (DP > 30 & FILTER = 'PASS') \" > {output}"
#bcftools view -f  PASS -Oz -o R-2011-3670_ONT.PASSED_variants.vcf.gz deepsomatic_R-2011-3760_ONT.vcf.gz
#bcftools view -i 'FORMAT/VAF>0.05 & FORMAT/DP>30'-f GERMLINE -Oz -o R-2011-3670_ONT.germline_variants.vcf.gz deepsomatic_R-2011-3760_ONT.vcf.gz
#bcftools view -i 'FORMAT/VAF>0.05 & FORMAT/DP>30' -f PASS -Oz -o R-2013-2166_ONT.hard_filterd.vcf.gz deepsomatic_R-2013-2166_ONT.vcf.gz
#nohup Rscript /mnt/storage2/sequencing/software/WES/mutsig.R deepsomatic_test/R-2013-2166_ONT.hard_filterd.vcf.gz MutSig/ > nohup.mutsig 2>&1 & #dragen 2 mutsig conda enviroment

#conda:
#     "mustig" #or RNAseq_norm conda enviroment
#shell:
#        "Rscript /mnt/storage2/sequencing/software/WES/mutsig.R {input} MutSig/"
#Manta SVs calling:
# ./mambaforge/envs/manta/bin/configManta.py --bam={params.normal} --tumorBam={input.bam_tumor} --referenceFasta= --runDir=Manta/{wildcards.sample} --generateEvidenceBam --config /mnt/storage2/sequencing/software/WGS/manta/config.ini
#/mnt/storage2/sequencing/software/WGS/manta/manta-1.6.0.centos6_x86_64/bin/configManta.py --bam={params.normal} --tumorBam={input.bam_tumor} --referenceFasta=/staging/human/reference/hg19/hg19.fa --runDir=Manta/{wildcards.sample} --generateEvidenceBam --config /mnt/storage2/sequencing/software/WGS/manta/config.ini
#Manta/{wildcards.sample}/runWorkflow.py