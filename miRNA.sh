#***************************************************************************************
#Basic sequence statistics. Print total number of reads, total number unique reads, percentage of 
#unique reads, most abundant sequence, its frequency, and percentage of total in file.fq:
#***************************************************************************************
cat RNA_Barcode_None_001.R_2011_10_24_03_13_31_user_IIC-100_212T.fastq  | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'

#for removing reads smaller than 17 nt and if base quality value is less than 20#
fastq_quality_trimmer -t 20 -l 17 -Q 33 -i  RNA_Barcode_None_001.R_2011_10_24_03_13_31_user_IIC-100_212T.fastq -o  without_barcode_trimmed.fastq

#for making uniform read length of 35bp#
fastx_trimmer -l 35 -Q 33 -i  without_barcode_trimmed.fastq -o  without_barcode_trimmed_35bp.fastq

#After pre processing total number of reads#
cat without_barcode_trimmed_35bp.fastq  | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'

# mapping with shrimp#
/data/results/pijush_data/DGSample56mirna/Alignment_Test/SHRiMP_2_2_3/bin/gmapper-ls -N 10 -o 1 -E -M mirna --qv-offset 33 without_barcode_trimmed_35bp.fastq /data/results/all_genomes/Mirbase19/mirbase_20/mirbase_20_with_mature_seq.fa >  mapped_reads.sam

#Mapping  Work from share apps
/share/apps/SHRiMP_2_2_3/bin/gmapper-ls -N 10 -o 1 -E -M mirna --qv-offset 33 without_barcode_trimmed_35bp.fastq /data/results/all_genomes/Mirbase19/mirbase_20/mirbase_20_with_mature_seq.fa >  mapped_reads.sam

#Sam to text file#
awk -F"\t" '{print $1, "\t", $3}' mapped_reads.sam >  mapped_reads_id.txt
	
#Separate each element#
less mapped_reads_id.txt | grep "hsa"  >  hsa_reads.txt
less mapped_reads_id.txt | grep "NR"  >  rrna_reads.txt
less mapped_reads_id.txt | grep "chr"  >  trna_reads.txt
less mapped_reads_id.txt | grep "PGM"  >  adaptor_reads.txt

#Find out the unique element which will be use in perl program#
less hsa_reads.txt | cut -f 2 | sort | uniq >  uniq_hsa_reads.txt
less rrna_reads.txt | cut -f 2 | sort | uniq >  uniq_rRNA_reads.txt
less trna_reads.txt | cut -f 2 | sort | uniq >  uniq_tRNA_reads.txt
less adaptor_reads.txt | cut -f 2 | sort | uniq >  uniq_adaptor_reads.txt

#Run perl asript#
perl mirna_count1.pl > mirna_read_count.txt
perl mirna_count3.pl > rrna_read_count.txt
perl mirna_count2.pl > trna_read_count.txt
perl mirna_count4.pl > adaptor_read_count.txt

#Some another read count can be found there (3107) which are not miRNA, tRNA, rRNA or adaptor reads.
#When calculate the summary it is required to discard those number of reads.
