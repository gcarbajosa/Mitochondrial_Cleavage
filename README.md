# Mitochondrial RNA Cleavage

The scripts in this repository together with other freely available bionformatics tools listed in the workflow below will let you calculate mitochondrial RNA cleavage ratios from any set of RNA-Seq fastq files. It is possible to adapt the pipeline to run with single end reads but the following workflow applies to paired end read data only.

These are the steps needed to go from fastq files to cleavage ratios:

1 - Trim sequence adapters from the fastq files with Trimgalore. We recommend using stringency 3. We recommend keeping quality as 0 because we are interested in the detection the position where the RNA fragments end and not the base calling. The command would be something like the following:

    ./trim_galore --quality 0 --stringency 3 --output_dir 3 --paired fq_file1 fq_file2
 
2- Trim PolyA reads using prinseq. We recommend trimming 5 bases. This step uses as input the fastq files obtained on the previous step.This would be an example command: 

    perl prinseq-lite.pl -fastq fq_file1_from_trimgalore -fastq2 fq_file2_from_trimgalore -out_good fq_out_filename -min_len 20 -trim_tail_left 5 -trim_tail_right 5

3- Map the trimmed reads from the previous steps with STAR using the EndToEnd option to avoid soft-clipping and therefore detect more accurately the positions of the original RNA fragment ends. Please note tha STAR requires a large amount of RAM memory (up to 32G) so this step would have to be run in a suitably powered computer. The --outSAMunmapped Within option can be specificed to keep unmapped reads for QC purposes. The --limitBAMsortRAM 100000000000 option is also recommended to avoid computer crashes due to RAM overload. We recommend reading the STAR manual if in doubt about these or other options (see https://github.com/alexdobin/STAR) An example command would be:

    ./STAR --runThreadN 18 --alignEndsType EndToEnd --outSAMtype BAM SortedByCoordinate --genomeDir Path_2_STAR_genome_dir --readFilesIn $readFilesInPath/fq_file1_from_trimgalore_and_prinseq fq_file2_from_trimgalore_and_prinseq --outFileNamePrefix outFileNamePrefix 

4 - Keep only mitochondrial reads (specify MT or chrM depending on the genome reference used) reads in bam file using samtools (http://samtools.sourceforge.net/):
     samtools index outFileNamePrefix.bam 
     samtools view -hb  outFileNamePrefix.bam chrM > outFileNamePrefix.MT.bam

5 - Generate cleavage ratios running get_linear_cleavage_ratio.pl. The script requires 2 arguments:
    - The name of the bam file name obtained using the previous steps
    - A prefix to be used by the script to generate the output

    perl get_linear_cleavage_ratio.pl outFileNamePrefix.MT.bam my_out_cleavage >my_out_cleavage.log 2>my_out_cleavage.err


Once the ratios have been generated it is possible to determine the presence of significantly high ratios across samples choosing clevage ratio threshold. It is also possible to generate a matrix with the chosen thershold to use as input for quantitive trait loci analysis (for example with PLINK). These are the steps required: 


6 - Filter peaks running get_himalaya_ratios.pl. The script uses a cleavage ratio threshold (we recommend 0.1) and keeps the highest ratio above the threshold. When there is a cluster of peaks (Himalaya) within a bp range (we recommend 5bp) the script selects the peak with the highest ratio (everest) to represent the cluster. It can be run using a command like: 

    perl get_himalaya_ratios.checkCov.pl --ratiosFile=rat_filename --everestOutfile=eve_out_filename --himalayaOutfile=him_out_filename >him_out_filename.log 2>him_out_filename.err

7 - Using the filtered peaks from the previous step, you can filter across your samples for the presence of cleavage ratios across a proportion of your samples running generate_all_himalayas_peaks_presence.pl

     - Generate a list of files to process. For example:
        
        ls *eve.txt > my_list_of_everest_files.txt

     - Calculate in which proportion of the files each peak is present (e.g. if a peak is present in 2 out 4 files it will have a proportion of 0.5)

        perl generate_all_himalayas_peaks_presence.pl  my_list_of_everest_files.txt

8 - Finally you can generate a matrix of cleavage ratios per sample to use as input for, for example, quantitative loci analysis with PLINK (https://zzz.bwh.harvard.edu/plink/). Youy need as input:
    - Filename of the peaks presence positions table (required)
    - Filename with a list of ratio files you want to use as input (required)
    - Threshold for the proportion you want to use as threshold (optional, default is 0.5)
    - Filename for the ouput file (optional, default is peaks_presence_matrix.05.txt)
    - Read coverage (optional, default is 20)

    perl generate_all_himalayas_peaks_matrix.pl peaks_presence_full_table.txt ratio_filenames_list.txt 0.5 peaks_presence_matrix.05.txt 20


