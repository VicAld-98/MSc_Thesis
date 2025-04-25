#This is the code used in my project. 
#Some of the steps do not have to be done in any particular order though.
#I've done most of it as a general code so that I remember what bits I need to change if I need to come back to it.
#This shows about 10% of the work I have done, the main code.
#I have not included all of the file and data manipulation to get the inputs and extract from the outputs of each section - this made up the vast majority of this project.





######################################## GENOME ASSEMBLY ######################################

#Step 1 - Jellyfish
#conda install -c bioconda -c conda-forge jellyfish

jellyfish count -m 21 -s 1G -t 64 -C -o kmer.jf S1_S4_reads.fastq
jellyfish histo -o 21kmer.histo kmer.jf







#Use this as input for GenomeScope2 (used server but here is code)

#Step 2 - GenomeScope: Tried with different ploidy then compared
#conda install bioconda::genomescope2
genomescope.R -i 21kmer.histo -o GenomeScope -k 21 -p 2
genomescope.R -i 21kmer.histo -o GenomeScope -k 21 -p 3





#Step 3 - Smudgeplot
#conda install bioconda::smudgeplot=0.4.0
#conda install bioconda::fastk

fastK -t 64 -k 21 -T32 S1_S4_reads.fastq
#This will output reads.ktab
smudgeplot.py hetmers -L 10 -t 64 -o smudgeplot reads.ktab
#This will produce smudgeplot.smu
smudgeplot.py all -o smudgeplot_final smudgeplot.smu





#Step 4 - ASSEMBLERS

#Step 4a - HiCanu
#conda install bioconda::canu
canu -p HiCanu -d . genomeSize=3G -pacbio-hifi S1_S4_reads.fastq


#Step 4b - HiFiASM
#conda install bioconda::hifiasm
#Initial run without duplicate purging 
hifiasm -l 0 -o HiFiASM_no_purge -t 64 -k 21 S1_S4_reads.fastq
#Run with intrinsic duplicate purging 
hifiasm -o HiFiASM_purge -t 64 -k 21 S1_S4_reads.fastq


#Step 4c - Flye
flye –pacbio-hifi S1_S4_reads.fastq --genome-size 3G –threads 64 -o .


#Step 4d - Shasta
shasta –-input S1_S4_reads.fastq --assemblyDirectory . --config HiFi-Oct2021 





#REPEAT STEPS 5 - 9 FOR ALL ASSEMBLIES!!!!

#Step 5 - QUAST - repeat for each assembly
#conda install bioconda::quast
quast.py [genome_assembly.fasta] -o quast_assembly

#Step 6 - Bandage (Used the GUI for this inital visualisation)

#Step 7 - BUSCO
#conda install bioconda::busco
busco -i [genome.fasta] -m genome -l insecta_odb10 --cpu 64 -o BUSCO_output

#Step 8 - Merqury
#Download meryl
#conda install bioconda::meryl
#Download Merqury
#conda install bioconda::merqury
#Create meryl database
meryl k=21 count S1_S4_reads.fastq output meryl_database.meryl

#Run Merqury from within MERQURY directory
./merqury.sh meryl_database.meryl [genome.fasta] [out_prefix]

#Step 9 - Mapping quality
#Install minimap2
#conda install bioconda::minimap2
#Install samtools
#conda install bioconda::samtools

#Mapping reads to the reference genome
minimap2 -t 64 -ax map-hifi [assembly.fasta[] S1_S4_reads.fastq | samtools view -b -o [output.bam]
#Sort the BAM file
samtools sort [output.bam] -o [sorted_reads.bam]
#Index the sorted BAM if necessary
samtools index [sorted_reads.bam]
#Look at the statistics
Samtools flagstat [sorted_reads.bam]
Samtools depth [sorted_reads.bam] > depth.txt
#Extract mapping information
awk '{sum+=$3; count++} END {print "Average Coverage:", sum/count}' depth.txt awk '{if($3==0) count++} END {print "Uncovered Bases:", count}' depth.txt





#Step 10 - Duplicate Purging 
#Install purge_dups
#conda install bioconda::purge_dups
#DO NOT USE THE WRAPPER SCRIPT!!!!!!!!
#2 rounds of purge_dups ran - best assembly was chosen for further purge_dups rounds

#Map the reads to the assembly
minimap2 -t 64 -x map-hifi -I 20G [assembly.fasta] S1_S4_reads.fastq | gzip -c - > assembly.paf.gz
#pbcstat 
/pub63/vica/miniconda3/envs/purge_dups/bin/pbcstat assembly.paf.gz
#calcuts 
/pub63/vica/miniconda3/envs/purge_dups/bin/calcuts PB.stat > cutoffs 2>calcuts.log
#Split fasta
/pub63/vica/miniconda3/envs/purge_dups/bin/split_fa [assembly.fasta] > assembly.split
#Self-alignment
minimap2 -t 64 -xasm5 -DP assembly.split assembly.split | gzip -c - > assembly.split.self.paf.gz
#Purge dups
/pub63/vica/miniconda3/envs/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov assembly.split.self.paf.gz > dups.bed 2> purge_dups.log
#Purged fasta
/pub63/vica/miniconda3/envs/purge_dups/bin/get_seqs -e dups.bed [assembly.fasta]





#Step 11 - Quality Control
#Rerun Step 5 - 10 to check assembly quality





#Step 12 - Blobtoolkit
#Install blobtoolkit
#conda install -c conda-forge -c hcc blobtoolkit=4.3.1
#Download taxdump
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzvf taxdump.tar.gz
#Create blobdirectory
blobtools create --fasta [purged_genome.fasta] --taxid 1133654 --taxdump [path to taxdump] /pub63/vica/miniconda3/envs/blobtoolkit/ .
#Add coverage files - MAKE SURE THE CONTIG NAMES ARE THE SAME AS IN THE FASTA FILE!!
samtools coverage [Assembly]_sorted_mapped.bam > [Assembly].coverage.txt
blobtools add --text [Assembler].coverage.txt --text-header --text-cols '#rname=identifier,meandepth=coverage' .
blobtools add --key plot.y=coverage .
#Add BUSCO from previous steps
blobtools add --busco [path to BUSCO output busco.tsv file] .
#Add BLAST hits
blastn -task megablast -query [assembly.fasta] -db /pub65/markw/CGR/progs/ncbi_dbs/nt/nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -culling_limit 5 -num_threads 64 -evalue 1e-25 -out [megablast.out]
blobtools add --hits [output].ncbi.blastn.out --hits-cols 1=qseqid,2=staxids,3=bitscore,5=sseqid,10=qstart,11=qend,14=evalue --taxrule bestsumorder --taxdump /pub63/vica/miniconda3/envs/blobtoolkit .
#View plots
blobtools view --plot .
blobtools view --view snail --plot . 
#If the coverage does not work then go into the json file and in the coverage
#section change scaleLog to scaleLinear.
#To look at the contig phylum results. You can see which contigs had which hits by opening the bestsumoder_phylum.json file.
#I looked at the coverage of the unusual contigs using the coverage.json file.




#Step 13 - BLAST the unusual contigs
#This was the Ascomycota and Mollusca contigs. 
#Download BLAST
#conda install -c bioconda blast
#BLAST the unusual contigs
blastn -query [contig.fasta] -db /pub65/markw/CGR/progs/ncbi_dbs/nt/nt -out [output.out] -outfmt 6 -max_target_seqs 10
#I also inspected the coverage of the unusual contigs by looking at the coverage.json file.
#Took the BLAST results from blobtools. See Teams convo with Al for the full explanation of Ascomycota contig.




#Step 14 - Looked at general genome coverage
#Overall coverage
#Mapping reads to the reference genome
minimap2 -t 64 -ax map-hifi FINAL_GENOME_HIFIASM3.fasta S1_S4_reads.fastq| samtools view -b -o [output.bam]
#Sort the BAM file
samtools sort [output.bam] -o [sorted_reads.bam]
#Index the sorted BAM if necessary
samtools index [sorted_reads.bam]
#Look at the statistics
Samtools flagstat [sorted_reads.bam]
Samtools depth [sorted_reads.bam] > [depth.txt]
#Extract mapping information
awk '{sum+=$3; count++} END {print "Average Coverage:", sum/count}' [coverage.txt] awk '{if($3==0) count++} END {print "Uncovered Bases:", count}' [coverage.txt]






################################### CHOSEN ASSEMBLY WAS HIFIASM AFTER 3 ROUNDS OF PURGE_DUPS ##################################

#Step 14 - RepeatModeller
#Install RepeatModeler
wget https://github.com/Dfam-consortium/RepeatModeler/archive/refs/tags/2.0.6.tar.gz
tar -zxvf RepeatModeler-2.0.6.tar.gz
#Build database
buildDatabase -name HiFiASM3 [final_genome.fasta]
repeatModeler -database [HiFiASM3] -threads 64 -LTRStruct





#Step 15 - RepeatMasker
#Install RepeatMasker
wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.7-p1.tar.gz
tar -zxvf RepeatMasker-open-4.1.7-p1.tar.gz
tar -zxvf RepBaseRepeatMaskerEdition-20181026.tar.gz
cd ./path/to/RepeatMasker
perl ./configure

#Install TRF
git clone https://github.com/Benson-Genomics-Lab/TRF.git
cd TRF
tar -xzvf trf-4.10.0.tar.gz
cd trf-4.10.0
mkdir build
cd build
../configure
make

#Install RMBlast
wget https://www.repeatmasker.org/rmblast/rmblast-2.14.1+-x64-linux.tar.gz
tar -zxvf rmblast-2.14.1-x64-linux.tar.gz
./configure

#Run
RepeatMasker -lib HiFiASM3-families.fa -pa 64 --xsmall final_assembly.fasta





#Step 16 - cdhit
#Used this to see what groups of repeats are present.
#conda install bioconda::cdhit
conda install bioconda::cdhit
cd-hit -i [input.fasta] -o [output.fasta] -c 0.95 -n 11
#-c is seq. identity min. -n = word size





#Step 17 - Plotting Repeats ################RUN IN R - check the R file too ###################
data <-read.csv("HiFiASM3_repeats2.csv")

#Count the number of Unknown class repeats
count_unknown <- sum(data$repeat_class == "Unknown")
print(count_unknown)

#Find the most common repeat family
most_common_string <- names(sort(table(data$matching_repeat), decreasing = TRUE))[1]
print(most_common_string)

#This was the most common repeat family - find the count
count_family1 <- sum(data$matching_repeat == "rnd-5_family-571")
print(count_repeat)

#Find the second most common repeat
freq_table <- table(data$matching_repeat)
sorted_freq <- sort(freq_table, decreasing = TRUE)
second_most_common <- names(sorted_freq)[2]
print(second_most_common)

#second most common
second_most_common <- names(sorted_freq)[2]
count_family2 <- sum(data$matching_repeat == "rnd-1_family-330")


#third most common
third_most_common <- names(sorted_freq)[3]
count_family3 <- sum(data$matching_repeat == "rnd-2_family-120")


########MOST COMMON REPEAT THAT IS NOT UNKNOWN########
filtered_data <- data[data$repeat_class == "Unknown", ]
dim(filtered_data)
#Find the most common value in 'matching_repeats' column of the filtered data
most_common_known_value <- names(sort(table(filtered_data$matching_repeat), decreasing = TRUE))[1]

#Print the most common value
print(most_common_known_value)
#Then count it
count_known <- sum(data$matching_repeat == "ltr-1_family-30")
print(count_known)

library(stringr)
#Count the number of repeats that are in each LTR family
count_ltr_1 <- sum(str_detect(filtered_data$matching_repeat, "^rnd-1_family-"))
count_ltr_2 <- sum(str_detect(filtered_data$matching_repeat, "^rnd-2_family-"))
count_ltr_3 <- sum(str_detect(filtered_data$matching_repeat, "^rnd-3_family-"))
count_ltr_4 <- sum(str_detect(filtered_data$matching_repeat, "^rnd-4_family-"))
count_ltr_5 <- sum(str_detect(filtered_data$matching_repeat, "^rnd-5_family-"))

#Create a summary table of the above data
summary_table <- data.frame(
  Family = c("rnd-1_family-", "rnd-2_family-", "rnd-3_family-", "rnd-4_family-", "rnd-5_family-"),
  Count = c(count_ltr_1, count_ltr_2, count_ltr_3, count_ltr_4, count_ltr_5)
)

#Print the summary table
print(summary_table)

#Data frame 
data <- data.frame(
  category = c("LINEs", "LTR elements", "Transposons", "Rolling Circles", "Simple Repeats", "Low Complexity Repeats"),
  percentage = c(2.06, 16.46, 8.56, 1, 1, 1),
  subgroup = NA # No subgroups here
)

#Data for the Unclassified subgroup %s
unclassified_data <- data.frame(
  category = rep("Unclassified", 5),
  subgroup = c("rnd1", "rnd2", "rnd3", "rnd4", "rnd5"),
  percentage = c(22 * 0.36, 22 * 0.075, 22 * 0.10, 22 * 0.189, 22 * 0.271)
)

#Ensure both data frames have the same columns (with NA for missing ones)
data$subgroup <- NA  
#Add a 'subgroup' column for the main data, to match unclassified_data

#combin the Unclassified data with the main data
combined_data <- rbind(data, unclassified_data)

#Calculate the total percentage for each category, and then list the categories in descending order (I think? I annotated retrospectively
#which I regret)
category_totals <- combined_data %>%
  group_by(category) %>%
  summarize(total_percentage = sum(percentage)) %>%
  arrange(desc(total_percentage))

combined_data$category <- factor(combined_data$category, levels = category_totals$category)

#Create the bar plot so each category has a bar, but the LTR bar has subgroups.
ggplot(combined_data, aes(x = category, y = percentage, fill = ifelse(is.na(subgroup), category, subgroup))) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(
    title = "Repeat composition of an assembled Phyllium philippinicum genome",
    x = "Category",
    y = "Percentage (%)",
    fill = "Subgroup"
  ) +
  scale_fill_manual(values = c("Unclassified" = "grey95", "LINEs" = "grey95", "LTR elements" = "green", 
                               "Transposons" = "red", "Rolling Circles" = "purple", "Simple Repeats" = "orange",
                               "Low Complexity Repeats" = "yellow",
                               "rnd1" = "darkblue", "rnd2" = "darkgreen", "rnd3" = "darkred", 
                               "rnd4" = "darkorange", "rnd5" = "lightgoldenrod")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))







###################################### Step 18 - ANNOTATION #######################################

#Attempt 1 - BRAKER3
#Issues with this due to the novelty of the genome - RNA-seq would help here
#Download BRAKER3
singularity build braker3.sif docker://teambraker/braker3:latest

#Download newest GeneMark from https://exon.gatech.edu/GeneMark/license_download.cgi
#NB!!!!!!!! This has to be version GeneMark-ES/ET/EP+ ver 4.72_lic and the Linux 64 option. 
#Also download the key. Unzip what I need to. 
#Put the GeneMark (rename GeneMarkES) in the same file as the genome file and the proteins database fasta file. 
#Make sure the key file is also in here. Try not to cry whilst doing this for the 100th time.
#(I downloaded the earlier version by accident and it created an output file which was not readable further down the pipeline) 

#Run from HOME
singularity exec braker3.sif braker.pl --genome=[genome.fa] --prot_seq=[OrthoDB_proteins.fa] --threads=48

#This did not work - come back to this - see portfolio for more details if you dare







#Attempt 2.1 - GALBA with OrthoDB
#Comparison against the closest Phasmid with a genome.
#Even then, this genome is of a stick insect which isn't really THAT close so this is not the best annotation.

#Download GALBA
singularity build galba.sif docker://katharinahoff/galba-notebook:latest

#Download the OrthoDB protein set Arthropoda
wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Arthropoda.fa.gz

#Run
singularity exec galba.sif galba.pl --species=[myspecies] --genome=[assembly.fasta] --prot_seq=[proteins.fa] #In this case OrthoDB


#1st attempt was alright. I used the full Arthropoda database.

######Attempt 2.2 - Using the genome of D. australis. 
#Get the protein database from NCBI. unzip etc to get the protein file.
datasets download genome accession GCA_029891345.1 --include protein
singularity exec galba.sif galba.pl --species=[myspecies] --genome=[assembly.fasta] --prot_seq=[proteins.fa] #In this case D. australis

#I then re-read the GALBA Doc and it says that if the protein ref is distantly related or not trusted 
#then disable the DIAMOND filter
singularity exec galba.sif galba.pl --species=[myspecies] --genome=[assembly.fasta] --prot_seq=[proteins.fa] --disable_diamond_filter=true #In this case D. australis

#This attempt gave me a tonne of proteins, the vast majority false hits.
#To sort this, I then put all of these annotations through EggNOG to find the ones with true functions.

######### See Portfolio for the full testing code - redundant testing of ProtHint and GeneMark to try identify where the error was occuring



#Step 19 - EggNOG Mapper
#I used the server. 
http://eggnog-mapper.embl.de/

#This gave a csv file of annotated genes.





#Step 20 - ColabFold
#I then looked in the GALBA annotation file to extract the protein sequence of these genes and mapped it with ColabFold
#See the portfolio for the issues I had with fragmentation.
https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb#scrollTo=kOblAo-xetgx


#Step x - Comparison of GC Content, Coverage, Repeats, Errors before and after purge_dups on HiFiASM.

############# Error rate #############
#For this I need the Merqury_HiFiASM.HiFiASM.qv file for BEFORE AND AFTER PURGE_DUPS

data <- read.delim("HiFiASM_merqury.HiFiASM.qv", header = FALSE)
head(data)
library(ggplot2)
#This is to add the contig number so I can plot it as the x-axis
data$X <- 1:nrow(data)

#Plot the contig index (from the previous line) against the Error Rate for the corresponding contig (Column 5)
ggplot(data, aes(x = X, y = V5)) +
  geom_bar(stat = "identity", fill = "gray", width =1) +  #Bar chart
  theme_minimal() +
  labs(title = "Error rate across contigs of HiFiASM assembly of Phyllium philippinicum genome",
       x = "Contig Index",
       y = "Error Rate") +
  scale_x_continuous(breaks = seq(0, nrow(data), by = 50)) +  #Label every 1000th contig
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


############ GC CONTENT #################
#Calculate the GC content % per contig - this was written with ChatGPT
seqkit fx2tab HiFiASM.fasta > hifi_contigs.txt 
gcdata <- read.csv("gc_content_per_contig.csv", header = FALSE)
awk '{gc=0; total=length($2); for(i=1;i<=total;i++){if(substr($2,i,1)=="G"||substr($2,i,1)=="C") gc++} gc_content=gc/total*100; print $1, gc_content}' hifi_contigs.txt > gc_content_per_contig.txt

#Add the contig index number
gcdata$X <- 1:nrow(gcdata)

#Plot the same as the previous code.
ggplot(gcdata, aes(x = X, y = V2)) +
  geom_bar(stat = "identity", fill = "gray", width = 1) +  #Bar chart
  theme_minimal() +
  labs(title = "GC content across contigs of HiFiASM assembly of Phyllium philippinicum genome",
       x = "Contig Index",
       y = "GC (%)") +
  scale_x_continuous(breaks = seq(0, nrow(data), by = 50)) +  #Label every 1000th contig
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylim(0,100)


########### Coverage #############
#Use the samtools depth data
awk '{sum[$1]+=$3; count[$1]++} END {for (contig in sum) print contig, sum[contig]/count[contig]}' [coverage.txt] > [contig_coverage.txt]
#REMEMBER FOR EVERYTHING TO SORT THE CONTIGS BY NUMBER SO THEY ARE IN THE ORDER ASSEMBLED - THIS DOES NOT MEAN IT IS FINAL ORDER!!

cvdata <- read.csv("coverage_per_contig_SORTED.csv", header = FALSE)
cvdata$X <- 1:nrow(cvdata)

ggplot(cvdata, aes(x = X, y = V2)) +
  geom_bar(stat = "identity", fill = "gray", width = 1) +  
  theme_minimal() +
  labs(title = "Coverage across contigs of HiFiASM assembly of Phyllium philippinicum genome",
       x = "Contig Index",
       y = "Coverage") +
  scale_x_continuous(breaks = seq(0, nrow(data), by = 50)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
which.max(cvdata$V2)

############# REPEATS ################## 
#Using this code I extracted the start, end and length of repeat region for each contig from the output table from RepeatMasker:

awk '{print $5, $6, $7}' repeatmasker.out > repeat_length_hifi.txt

#Then I merged the repeats - this was written with ChatGPT and even then, it took a few days of faffing to get it to work.

sort -k1,1 -k2,2n repeat_positions_HIFIASM3.txt | awk '
BEGIN {contig=""; start=0; end=0}
{
  if ($1 == contig) {
    if ($2 <= end) {
      #Overlapping or adjacent regions, extend the end
      if ($3 > end) {
        end = $3
      }
    } else {
      #No overlap, print the previous region and reset start/end
      print contig, start, end
      start = $2
      end = $3
    }
  } else {
    #New contig, print the previous one if it exists
    if (contig != "") {
      print contig, start, end
    }
    contig = $1
    start = $2
    end = $3
  }
}
END {print contig, start, end}
' > merged_HIFIASM3.txt

#Then calculated the length of each repeat:
awk '{print $1, $2, $3, $3-$2+1}' repeat_length_hifi.txt > repeat_length_hifi.txt

#I then had to sum the length of the repeats for each contig:
awk '{length = $3 - $2; sum[$1] += length} END {for (contig in sum) print contig, sum[contig]}' input_file > output_file

#Then had to find the length of each contig, I am sure they are in a file somewhere but I used:
awk '/^>/ {if (seqlen) print seqlen; print $0; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' HiFiASM.fasta > contig_lengths.txt

#I joined the contig length and repeat length file:
join repeats_sum.txt contig_lengths.txt > joined_output.txt

#Finally I calculated the percentage of repeats for each contig (Note that there are only 585 contigs with repeats out of the total 620):
awk '{print $1, ($2 / $3) * 100}' joined_output.txt > repeats_percentage.txt

#THEN I PLOTTED IT IN R LIKE PREVIOUS
rpdata <- read.csv("repeat_percentage.csv", header = FALSE)
rpdata$X <- 1:nrow(rpdata)

ggplot(rpdata, aes(x = X, y = V2)) +
  geom_bar(stat = "identity", fill = "gray", width = 1) +
  theme_minimal() +
  labs(title = "Repeat content of contigs across HiFiASM assembly of Phyllium philippinicum genome",
       x = "Contig Index",
       y = "Repeat Content (%)") +
  scale_x_continuous(breaks = seq(0, nrow(rpdata), by = 50)) +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 



###################################### MITOGENOMICS ####################################

# Step 21 - MitoHiFi
#Install MitoHiFi from the docker. It isn’t the one on GitHub.
singularity pull docker://ghcr.io/marcelauliano/mitohifi:master
#Download the nearest mitogenome
singularity exec mitohifi_master.sif findMitoReference.py --species "Phyllium philippinicum" --outfolder ./mitohifi --min_length 13000
#Run from HOME!
singularity exec mitohifi_master.sif mitohifi.py -c [contig_file.fasta] -f [Related_mitogenome.fasta] -g [related_mito.gbk] -t [threads]
#This gave an error that the parse txt files were empty. 
#Tried again but using the reads this time instead of the assembly.
#Turns out the contigs are too long so had to use it on the reads instead of contigs
#Run from HOME! THIS ONE!
singularity exec mitohifi_master.sif mitohifi.py -r [reads.fasta] -f [Related_mitogenome.fasta] -g [related_mito.gb] -t [threads]




#Step 22 - Mitos2 
#Used the server
https://usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fmitos2%2Fmitos2%2F2.1.7%20galaxy0
#Used the Invertebrate Mito code (5). Also got the structure, bed, gff, fasta file etc.




#Step 23 - Clustal Omega for alignment of full 13 CDS from mitogenomes from Phasmatodea.
#Full seqs. obtained from NCBI and then concatated manually in NotePad++
https://www.ebi.ac.uk/jdispatcher/msa/clustalo





#Step 24 - IQTree
http://iqtree.cibiv.univie.ac.at/
#See the results in the portfolio/paper for info on the model selection etc etc. 


#Step 25 - Visualisation of the circular mitogenome.
#Used a server for this too - GeSeq
https://chlorobox.mpimp-golm.mpg.de/geseq.html
#Uploaded the CDS and tRNA sequences for annotation.


#Comparison of mitogenomes
#To calculate the length of PCGs, tRNAs, rRNAs and CR, I used Excel. I also used this for means and averages.
#I then created a boxplot of all of this information.

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)


df <- read.xlsx("mito_stats.xlsx")


df_long <- df %>%
  pivot_longer(cols = c(rrnL, rrnS, tRNA, Control_Region, Total_Length, PCG), 
               names_to = "Region", values_to = "Length")

#Replace underscores with spaces in Region labels
df_long$Region <- str_replace_all(df_long$Region, "_", " ")

#Create boxplot!
ggplot(df_long, aes(x = Region, y = Length / 1000)) +
  geom_boxplot(fill = "grey", color = "black") +
  geom_point(aes(color = Color), size = 3) +
  #Add labels only for Control Region and Total Length for P. philippinicum because these are the only features where P. phil is an outlier
  geom_text(
    data = df_long %>%
      filter(Species == "P. philippinicum", Region %in% c("Control Region", "Total Length")),
    aes(label = Species),
    vjust = -0.5,
    size = 3,
    color = "darkgreen"
  ) +
  scale_color_identity() +
  scale_y_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 2)
  ) +
  labs(
    y = "Length (Kb)",
    x = "Feature",
    title = "Boxplot of Length of the Mitochondrial Features in Euphasmatodea"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )



#Created a PCA plot to see what features contribute to the differences (but tbh it is pretty obvious from the data)

library(openxlsx)
df <- read.xlsx("mito_stats.xlsx")
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggrepel)
library(factoextra)

df_scaled <- scale(df[,1:6], scale = T, center = T)

#Do PCA
pca_result <- prcomp(df_scaled, center = TRUE, scale = TRUE)

#Get PCA results into a dataframe
pca_df <- as.data.frame(pca_result$x)
pca_df

#Add the "Group" column back to the df :)
pca_df$Group <- df$Group

#Plot PCA results with ggplot
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "darkgreen")) +  # Customize colors
  stat_ellipse(aes(group = Group), level = 0.95, linetype = "dashed", linewidth = 1) +  # Add ellipses
  labs(title = "PCA of Length of Mitochondrial Features:Total, PCG, tRNA, rrnS, rrnL, Control Region", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(legend.title = element_blank())+
  geom_text_repel(aes(label = ifelse(df$Species == "P. philippinicum", "P. philippinicum", "")),
                  size = 3, color = "black")

#Check what PC loadings are
summary(pca_result)


pca_result$rotation[, 1]
#Define PC1 loadings - probably couldve done this a better way
loadings1 <- data.frame(
  Variable = c("PCG", "rrnL", "rrnS", "tRNA", "Control_Region", "Total_Length"),
  PC1 = c(0.20552986, -0.07613023, -0.21105527, 0.25623724, -0.65267376, -0.64480755)
)


#Create loading plot for PC1
ggplot(loadings1, aes(x = Variable, y = PC1, label = Variable)) +
  geom_bar(stat = "identity", fill = "darkgray") +
  geom_text(aes(label = round(PC1, 3)), vjust = -0.5) +
  theme_minimal() +
  labs(title = "PC1 Loadings", x = "Variables", y = "Loading Value")


#PC2 Loadings
pca_result$rotation[,2]
loadings2 <- data.frame(
  Variable = c("PCG", "rrnL", "rrnS", "tRNA", "Control_Region", "Total_Length"),
  PC2 = c(-0.48770640, -0.66637362, -0.33158341, -0.444467221, -0.07629180, -0.06772958))


#Create loading plot for PC2
ggplot(loadings2, aes(x = Variable, y = PC2, label = Variable)) +
  geom_bar(stat = "identity", fill = "darkgray") +
  geom_text(aes(label = round(PC2, 3)), vjust = -0.5) +
  theme_minimal() +
  labs(title = "PC2 Loadings", x = "Variables", y = "Loading Value")

#####Cumulative Proportions#####

cumulative_prop <- c(0.3553, 0.5836, 0.7851, 0.9133, 0.99775, 1.00000)
pcs <- paste("PC", 1:6)

#Create my df
prop_df <- data.frame(PC = pcs, Cumulative_Prop = cumulative_prop)

#Plot the chart
ggplot(prop_df, aes(x = PC, y = Cumulative_Prop)) +
  geom_bar(stat = "identity", fill = "darkgray", color = "black") +
  labs(title = "Cumulative Proportion of Variance by Principal Component", 
       x = "Principal Components", 
       y = "Cumulative Proportion") +
  theme_minimal()


#2 of the mitos2 structure plots looked unusual - stem bulges. I re-modelled them using RNAFold server, but decided to stick
#with mitos2 - see portfolio for explanation.

#Used https://jamiemcgowan.ie/bioinf/rscu.html for RSCU calculations
#Then visualised this using R.

data <- read.csv("Codon.csv")
library(ggplot2)
library(ggpattern)

#Codon frequency plot
#Order the data by 'Count' in descending ord
data$Symbol <- factor(data$Symbol, levels = data$Symbol[order(-data$Frequency)])
library(ggplot2)
library(ggpattern)

ggplot(data, aes(x = reorder(Symbol, -Frequency), y = Frequency, pattern = Symbol)) +
  geom_bar_pattern(stat = "identity", 
                   aes(fill = ifelse(Symbol %in% c("STOP", "Ile", "Met"), "white", "grey")), 
                   color = "black",  
                   pattern_density = 0.1, 
                   pattern_spacing = 0.02,
                   pattern_fill = "black") +  
  scale_fill_identity() +  
  scale_pattern_manual(values = c("STOP" = "stripe", "Ile" = "circle", "Met" = "circle", 
                                  setNames(rep("none", length(unique(data$Symbol)) - 3), 
                                           setdiff(unique(data$Symbol), c("STOP", "Ile", "Met")))), 
                       na.translate = FALSE) +
  theme_minimal() +
  labs(x = "Codon Family", y = "Frequency", title = "Codon Family Usage Within the PCGs of P. philippinicum mitogenome") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


###RSCU Plot

data2 <- read.csv("RSCU.csv")

ggplot(data2, aes(x = reorder(Codon, -RSCU, na.rm = TRUE), y = RSCU)) +
  geom_bar_pattern(stat = "identity", 
                   aes(fill = ifelse(Codon %in% c("ATT", "ATC", "ATA", "ATG", "TAA"), "white", "grey"),
                       pattern = ifelse(Codon %in% c("TAA", "TAG"), "stripe", 
                                        ifelse(Codon %in% c("ATT", "ATC", "ATA", "ATG"), "circle", "none"))),  
                   color = "black",  
                   pattern_density = 0.1, 
                   pattern_spacing = 0.02,
                   pattern_fill = "black") +  
  scale_fill_identity(guide = "none") +  #Removes legend
  scale_pattern_manual(
    values = c("stripe" = "stripe", "circle" = "circle", "none" = "none"),
    guide = "none"  #Removes  legend
  ) +
  theme_minimal() +
  labs(x = "Codon", y = "RSCU", title = "Codon RSCU within the PCGs of P. philippinicum mitogenome") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.1, vjust = 0.5))

#CDspT plot
#Used the NCBI mito data for the PCGs for each phylliidae species, then used the above codon usage calculator to get the count.
#CDspT is then (codon count / total mito length)*1000. Did this for each codon for each species. Then plotted it. 
#CDSPT data
codon_df <- read.csv("cdspt.csv", header = T)

library(ggplot2)
library(tidyr)

#Long format
codon_long <- pivot_longer(codon_df, 
                           cols = -Codon, 
                           names_to = "Sample", 
                           values_to = "Value")

codon_long$Species <- gsub("_", ". ", codon_long$Sample)

#Plot for the CDspT
ggplot(codon_long, aes(x = Codon, y = Value, color = Sample, group = Sample)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = "Codon Frequency per 1000 Codons Across Phylliidae Mitogenomes",
       x = "Codon",
       y = "Frequency per Thousand (CDspT)",
       color = "Species") +
  scale_color_discrete(
    labels = c(
      "C_oyae" = "C. oyae",
      "C_tibetense" = "C. tibetense",
      "C_westwoodii" = "C. westwoodii",
      "P_bioculatum" = "P. bioculatum",
      "P_giganteum" = "P. giganteum",
      "P_philippinicum" = "P. philippinicum"
    )
  ) +
  theme(
    axis.text.x = element_text(angle = 90, size = 8),
    axis.line = element_line(color = "black", linewidth = 0.8),
    panel.grid = element_blank()
  )






################################## COMPARATIVE GENOMICS ###########################

#Step  - Mummer
#Due to there being no close genomes for alignment this was not very successful.
#However it did work for mtDNA comparison 
./nucmer --mum -c 100 -l 30 -g 500 -b 1000 -t 64 -p HOOKERI [final_genome.fasta] [Hookeri.fasta] 

#Filter high identity and length alignmenets
delta-filter -i 90 -l 1000 -m [alignment.delta] > [filtered.delta]

#Generate the plot
./mummerplot  [filtered.delta] --png




#Step 22 - DGenie
#Used the server for genome comparison here
https://dgenies.toulouse.inra.fr/run
#Created plots, they weren't great either due to the lack of close genomes.
#The '''''best''''''' one was the comparison to C. hookeri.





########################### ATTEMPT 1 AT PANGENOME #######################################

#Step 23 - ProgressiveMauve ############### DID NOT DO IN THE END!!!!!!!! ###########################
#Install
conda install -c bioconda progressivemauve
#Run comparing my genome to the other genomes
progressiveMauve --output=[output.xmfa] genome1.fasta genome2.fasta genome3.fasta genome4.fasta genome5.fasta
#This ran for WEEKS with no end in sight. I cancelled it. 

#ARTEMIS and ProgressiveCactus were also attempted without success. 





#Step 24 - odgi ####################    Did not use this either  ########################
conda install -c bioconda odgi
#Map the new reads to the old reads using minimap2
minimap2 -x asm20 old_assembly.fasta new_assembly.fasta > mappings.paf
#Convert the old gaf file to ODGI format:
odgi build -g old_assembly.gfa -o old_graph.og
#Create an image of the graph
odgi viz -i old_graph.og -o old_graph.png
#Highlight the mappings on the old graph.
odgi depth -i old_graph.og -p mappings.paf -o overlay.og
#Create the png image
odgi viz -i overlay.og -o overlay.png
#Also looked awful.





#Step 25 - Looking to see if there are beginning of chromosomes
#So to do this, I opened up the final assembly fasta.
#I manually looked for repeating regions at the beginning of contigs. 
#I found the repeat "TTAGG" and the rev comp - "CCTAA" - the common insect telomere.
#This was found at the end of 17 contigs:
#Contig 7
#Contig 42 
#Contig 53 
#Contig 55
#Contig 64 
#Contig 69 
#Contig 76 
#Contig 82
#Contig 85 
#Contig 92 
#Contig 95 
#Contig 104
#Contig 124
#Contig 130
#Contig 154
#Contig 178
#Contig 199

#It may be that the other ends of the chromsomes were not accurately sequenced
#due to the repeating sections, however, there were some sub-telomere regions that were present. 
#Other Euphasmids have been shown to have approx. 17 - 18 chromosomes


############################# MISC. ####################################

#FastQ to Fasta
sed -n '1~4s/^@/>/p;2~4p' input.fastq > output.fasta


#Making a BLASTDB
makeblastdb -in input.fasta -dbtype nucl -out reads_db -max_file_sz 900MB


#Kill screen
screen -X -S [screen_name] quit

##############THIS IS A BRANCH TEST ######################