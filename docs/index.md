# Germline Variant Detection

## Method

1. BWA-MEM (Li et al 2013) is used to align sequence reads to reference genome. The default genomes is GRCh38, the lastest human build.  GRCm38 (mm10) is also available.
2. Post-BAM processing is done with BWAKit, Samtools 1.4, Sambamba and GATK 3.7  (DePristo et al 2011, McKenna et al 2010).
3. Variants are detected using GATK 3.7 (DePristo et al 2011, McKenna et al 2010), Platypus (Rimmer et al 2014), Samtools version 1.4 and FreeBayes version 0.9.7 (Garrison and Marth 2012).
4. A union VCF file is created with the results from each individual caller for filtering by scientist.
5. Effect of SNPs and INDELs on genes was predicted using snpEff (Cingolani et al 2012).
6. For Human Samples aligned to GRCh38: discovered variants are annotated using SnpSift (Cingolani et al 2012) using the following database:
   a. ExAC to determine general population AF
   b. dbSNP
   c. COSMIC (Forbes et al. 2009)
   d. CLINVAR (Landrum et al 2014)
   e. GWAS Catalog (Welter et al 2014)
   f. dbNSFP (Liu et al 2011) 


## Key Output Files
This workflow package provides the following key output files:

1. VCF file — SNPs/Indels for each sample labeled SampleID.annot.vcf.gz
2. Coverage Histogram for each sample labeled SampleID.coverage_histogram.png
3. Cumulative Distribution Plot for all samples labeled coverage_cdf.png
4. QC for all samples labeled sequence.stats.txt
5. Structural Variants (unfiltered) labeled SampleID.sssv.sv.vcf.gz.annot.txt

## Recommended Filtering for Germline Testing
1. ExAC POPMAX AF (0.01-0.05) -  depends on rarity of the phenotype of the proband
2. Depth >10
3. LOF or Misssense (Coding Changes)
4 . Alt Read Ct > 3
5. Mutation Allele Frequency (MAF) > 0.15
6. If novel: Called by 2+ callers

## Workflow Parameters

1. fastq
   a. Choose one or more FastQ read files to process.
2. genome
   a. Choose a genomic reference (genome).
3. pairs
   a. Choose if pair-ended or single-end sequences
4. capturebed
   a. Choose the bedfile that corresponds to the capture kit used for the sequencing.
   b. This is used to calculate average sequence coverage and filter SNPs
5. design
   This file matches the fastq files to data about the sample

## File Format For Design File

[https://cloud.biohpc.swmed.edu/index.php/s/fpWTgU5UikzUldb](See Workflow Template Design Files)

The following columns are necessary, must be named as in template and can be in any order:

    SampleID
        This ID is a unique identifier for the sample and will be used to name all of the sample files and act as the same names in the VCF
    FqR1
    	Name of the fastq file R1
    FqR2
    	Name of the fastq file R2


## Test Data
[Available with BioHPC Account](https://lamella.biohpc.swmed.edu/index.php/s/rJYLrm96VGg7DsR)

## Credits
[BICF](http://www.utsouthwestern.edu/labs/bioinformatics/)

## References
1.	Andy Rimmer, Hang Phan, Iain Mathieson, Zamin Iqbal, Stephen R. F. Twigg, WGS500 Consortium, Andrew O. M. Wilkie, Gil McVean, Gerton Lunter. Integrating mapping-, assembly- and haplotype-based approaches for calling variants in clinical sequencing applications. Nature Genetics (2014) doi:10.1038/ng.3036
2.	Bernstein, B. E., Birney, E., Dunham, I., Green, E. D., Gunter, C., & Snyder, M. (2012). An integrated encyclopedia of DNA elements in the human genome. Nature, 489, 57–74. doi:10.1038/nature11247
3.	Cantarel, B. L., Weaver, D., McNeill, N., Zhang, J., Mackey, A. J., & Reese, J. (2014). BAYSIC: a Bayesian method for combining sets of genome variants with improved specificity and sensitivity. BMC Bioinformatics, 15, 104. doi:10.1186/1471-2105-15-104
4.	Cingolani, P., Platts, A., Wang, L. L., Coon, M., Nguyen, T., Wang, L., ? Ruden, D. M. (2012). A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. Fly. doi:10.4161/fly.19695
5.	Cingolani, P., Patel, V. M., Coon, M., Nguyen, T., Land, S. J., Ruden, D. M., & Lu, X. (2012). Using Drosophila melanogaster as a Model for Genotoxic Chemical Mutational Studies with a New Program, SnpSift. Frontiers in Genetics. doi:10.3389/fgene.2012.00035
6.	Challis, D., Yu, J., Evani, U. S., Jackson, A. R., Paithankar, S., Coarfa, C., ? Yu, F. (2012). An integrative variant analysis suite for whole exome next-generation sequencing data. BMC Bioinformatics. doi:10.1186/1471-2105-13-8
7.	DePristo, M. A., Banks, E., Poplin, R., Garimella, K. V, Maguire, J. R., Hartl, C., ? Daly, M. J. (2011). A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nature Genetics, 43, 491–498. doi:10.1038/ng.806
8.	Exome Variant Server, NHLBI GO Exome Sequencing Project (ESP), Seattle, WA (URL: http://evs.gs.washington.edu/EVS/) [01 (01, 2016) accessed].
9.	Forbes, S. A., Tang, G., Bindal, N., Bamford, S., Dawson, E., Cole, C., ? Futreal, P. A. (2009). COSMIC (the Catalogue of Somatic Mutations In Cancer): A resource to investigate acquired mutations in human cancer. Nucleic Acids Research, 38. doi:10.1093/nar/gkp995
10.	Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. arXiv preprint arXiv:1207.3907 [q-bio.GN] 2012
11.	Hansen NF, Gartner JJ, Mei L, Samuels Y, Mullikin JC. Shimmer: detection of genetic alterations in tumors using next-generation sequence data. Bioinformatics. 2013 Jun 15;29(12):1498-503. doi: 10.1093/bioinformatics/btt183. Epub 2013 Apr 24. PubMed PMID: 23620360; PubMed Central PMCID: PMC3673219.
12.	Kim S, Jeong K, Bhutani K, Lee J, Patel A, Scott E, Nam H, Lee H, Gleeson JG, Bafna V. Virmid: accurate detection of somatic mutations with sample impurity inference. Genome Biol. 2013 Aug 29;14(8):R90. doi: 10.1186/gb-2013-14-8-r90. PubMed PMID: 23987214; PubMed Central PMCID: PMC4054681.
13.	Koboldt DC, Zhang Q, Larson DE, Shen D, McLellan MD, Lin L, Miller CA, Mardis ER, Ding L, Wilson RK. VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing. Genome Res. 2012 Mar;22(3):568-76. doi: 10.1101/gr.129684.111. Epub 2012 Feb 2. PubMed PMID: 22300766; PubMed Central PMCID: PMC3290792.
14.	Landrum, M. J., Lee, J. M., Riley, G. R., Jang, W., Rubinstein, W. S., Church, D. M., & Maglott, D. R. (2014). ClinVar: Public archive of relationships among sequence variation and human phenotype. Nucleic Acids Research, 42. doi:10.1093/nar/gkt1113
15.	Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv Preprint arXiv, 00, 3. doi:arXiv:1303.3997 [q-bio.GN]
16.	Liu X, Jian X, Boerwinkle E. dbNSFP: a lightweight database of human nonsynonymous SNPs and their functional predictions. Hum Mutat. 2011 Aug;32(8):894-9. doi: 10.1002/humu.21517. PubMed PMID: 21520341
17.	McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., ? DePristo, M. A. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20, 1297–1303. doi:10.1101/gr.107524.110
18.	The 1000 Genome Consortium. An integrated map of genetic variation from 1,092 human genomes. (2012). Nature, 491(7422), 56–65. Retrieved from http://dx.doi.org/10.1038/nature11632.
19.	Saunders CT, Wong WS, Swamy S, Becq J, Murray LJ, Cheetham RK. Strelka: accurate somatic small-variant calling from sequenced tumor-normal sample pairs. Bioinformatics. 2012 Jul 15;28(14):1811-7. doi: 10.1093/bioinformatics/bts271. Epub 2012 May 10. PubMed PMID: 22581179.
20.	Welter D, MacArthur J, Morales J, Burdett T, Hall P, Junkins H, Klemm A, Flicek P, Manolio T, Hindorff L, and Parkinson H. The NHGRI GWAS Catalog, a curated resource of SNP-trait associations. Nucleic Acids Research, 2014, Vol. 42 (Database issue): D1001-D1006.

