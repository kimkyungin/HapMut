# HapMut

HapMut is a somatic point mutation caller for cancer sequencing data. HapMut first determines haplotype of sequence reads probabilistically and use the haplotype information to call somatic mutations. Matched normal and tumor sequencing samples are required. HapMut can be run on a commandline as, for example
  
  $ lua hapmut.lua ref.fa 21 normal.bam tumor.bam

Here, normal.bam and tumor.bam represent sorted BAM files for matched cancer samples. ref.fa is the reference genome used to align the normal and the tumor sample. "21" represents the chromosome 21. Output of the HapMut consists of 

  pos ref_allele alt_allele prob1 ... prob16
  
where pos is the zero-based genomic position of the found somatic mutation. ref_allele and alt_allele are the reference allele and alternate allele for the genomic position respectively. prob1 to prob16 are 16 posterior probabilities that correspond to 16 possible somatic mutation status for two haplotypes of the genomic site; each haplotype has four different DNA configuration, A, C, G, and T.   

Some preliminary steps are required to run HapMut.HapMut uses two existing libraries: lua scripting language v5.2 and samtools v.0.1.19. Those two softwares should be preinstalled and those two library information (location of header file and library file) should be informed to Makefile. To install, one can run "make" on the folder having the source code. It was tested on Mac OS X 10.9.5 and Ubuntu 14.04. After running "make", the following command should be performed before running HapMut to get reference genome information.

  $ lua get-chr-size.lua ref.fa
