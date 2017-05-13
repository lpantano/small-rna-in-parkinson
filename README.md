[![Project Status: Inactive - The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows.](http://www.repostatus.org/badges/latest/inactive.svg)](http://www.repostatus.org/#inactive)


Scripts and raw data and reanalysis of data related to:

[Specific small-RNA signatures in the amygdala at premotor and motor stages of Parkinson's disease revealed by deep sequencing analysis](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btv632). Pantano L, Friedländer MR, Escaramís G, Lizano E, Pallarès-Albanell J, Ferrer I, Estivill X, Martí E. Bioinformatics. 2015 Nov 2. pii: btv632. [Epub ahead of print] PMID: 26530722

Data is available at: [GSE97285](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97285)


Data was re-analyze with [bcbio-nextgen](https://bcbio-nextgen.readthedocs.io/en/latest/) following these methods:


Data was analyzed with bcbio-nextgen (https://github.com/chapmanb/bcbio-nextgen)
using piDNA to detect the adapter, cutadapt to remove it, STAR/bowtie to align against
the genome and seqcluster to detect small RNA transcripts. miRNAs were detected using
miraligner tool with miRBase as the reference miRNA database. tRNA profiles were
detected using tdrmapper tool. mirdeep2 was used for discovery of novel miRNAs. FastQC
was used for QC metrics and multiqc for reporting.

See citation at: https://bcbio-nextgen.readthedocs.io/en/latest/contents/citations.html
