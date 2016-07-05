# Astrocyte RNASeq Workflow Package

This workflow carries out a simple RNASeq Analysis, including simple pairwise differential expression analysis and alternative splicing primary analysis. One or more FASTQ files containing reads from a RNASeq experiment can be selected as input. For each file this workflow:

    1) Reads are trimmed by quality (threshold of 25) and adapter sequence at the ends and reads shorter than 35bps are removed.

    2) Trimmed Fastq files are aligned to the selected reference genome using HiSAT2 (Kim et al. 2015)

    3) Marks Duplicates using PICARD
 
    4) Features (genes, transcrips and exone) are counted using featureCounts (Liao et al 2014) and StringTie (Pertea et al. 2015) using the Gencode feature table(Harrow et al. 2012)

    5) Basic pairwise differential expression analysis is performed using EdgeR (Robinson et al. 2010) and DESeq

    6) Abundances of transcripts are calculated using ballgown (Frazee et al. 2014)

##Workflow Parameters

    fastq - Choose one or more RNASeq read files to process.

    index - Choose a genomic reference (genome).

    design - This file matches the fastq files to data abot the sample

 The following columns are necessary, must be named as in template and can be in any order:

    SampleID
        This ID should match the name in the fastq file ie S0001.R1.fastq.gz the sample ID is S0001
    SampleName
        This ID can be the identifier of the researcher or clinician
    SubjectID
        Used in order to link samples from the same patient
    SampleGroup
	This is the group that will be used for pairwise differential expression analysis
    FullPathToFqR1
	Name of the fastq file R1
    FullPathToFqR2
	Name of the fastq file R2

There are some optional columns that might help with the analysis:
      Tissue
      Gender
      CultureDate
      SequenceRun
      Organism
      CellPopulation
      Treatment
      GeneticFeature (WT or KO)
      Race
      Ethnicity
      Age
      

### Test Data


### Credits
This example worklow is derived from original scripts kindly contributed by the Bioinformatic Core Facility (BICF), Department of Bioinformatics and Clinical Sciences.

### References
    Kim D, Langmead B, Salzberg SL (2015). HISAT: a fast spliced aligner with low memory requirements. Nat Methods. Apr;12(4):357-60. doi: 10.1038/nmeth.3317. Epub 2015 Mar 9. PubMed PMID: 25751142.
    Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009). The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9.PMID: 19505943.
    Liao Y, Smyth GK and Shi W (2014). featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30.
    Harrow, J., Frankish, A., Gonzalez, J. M., Tapanari, E., Diekhans, M., Kokocinski, F., Hubbard, T. J. (2012). GENCODE: The reference human genome annotation for The ENCODE Project. Genome Research. doi:10.1101/gr.135350.111.
    Robinson, M. D., McCarthy, D. J., Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics (Oxford, England), 26, 139-140. doi:10.1093/bioinformatics/btp616.
    Yaari G, Bolen CR, Thakar J, Kleinstein SH (2013). Quantitative set analysis for gene expression: a method to quantify gene set differential expression including gene-gene correlations. Nucleic Acids Res. 2013 Oct;41(18):e170. doi: 10.1093/nar/gkt660. Epub Aug 5. PubMed PMID: 23921631; PubMed Central PMCID: PMC3794608.
    Pertea M, Pertea GM, Antonescu CM, Chang TC, Mendell JT, Salzberg SL (2015). StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat Biotechnol. Mar;33(3):290-5. doi: 10.1038/nbt.3122. Epub 2015 Feb 18. PubMed PMID: 25690850.
    Frazee AC, Collado-Torres L, Jaffe AE and Leek JT. ballgown: Flexible, isoform-level differential expression analysis. R package version 2.0.0.

