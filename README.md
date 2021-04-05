
# Biogeography of Bacterial Communities and Specialized Metabolism in Human Aerodigestive Tract Microbiomes

Reed M. Stubbendieck<sup>1</sup> & Cameron R. Currie<sup>1,2,3</sup>

<sup>1</sup>Department of Bacteriology, University of Wisconsin-Madison, Madison, WI 53706
<sup>2</sup>Department of Energy Great Lakes Bioenergy Research Center, University of Wisconsin-Madison, Madison, WI 53706
<sup>3</sup>Laboratory of Genetics, University of Wisconsin-Madison, Madison, WI 53706

Contact: stubbendieck@wisc.edu

## Introduction

This repository contains the code necessary to replicate the results of our study on the biosynthetic gene clusters (BGCs) in aerodigestive tract (ADT) microbiomes.

## Raw Datasets

The raw genome sequences used in this study were accessed from the expanded Human Oral Microbiome Database ([eHOMD](http://www.ehomd.org/)) version 9.03 ([FTP link](http://www.homd.org/ftp/HOMD_prokka_genomes/fna/ALL_genomes.fna)).

The raw metagenome sequencing reads used in this study were accessed from the NIH integrative Human Microbiome Project ([iHMP](https://www.hmpdacc.org/ihmp/)).

## Derived Datasets

The derived datasets can be downloaded from FigShare [here](https://doi.org/10.6084/m9.figshare.14217326).

## Prerequisites

* [ANCOM version 2.1](https://github.com/FrederickHuangLin/ANCOM) - located in ./scripts/r/
* [antiSMASH version 4.2.0](https://anaconda.org/bioconda/antismash/4.2.0/download/linux-64/antismash-4.2.0-py27hb89731f_1.tar.bz2)
* [BiG-SCAPE](https://git.wageningenur.nl/medema-group/BiG-SCAPE)
* [Cytoscape version 3.8.0](https://cytoscape.org/)
* [genbank_to_fasta.py](https://rocaplab.ocean.washington.edu/tools/genbank_to_fasta/) - in path
* [kallisto version 0.46.0](https://github.com/pachterlab/kallisto)
* [hmmscan](http://hmmer.org/) - in path
* [MetaPhlAn version 3.0](https://huttenhower.sph.harvard.edu/metaphlan)
* [Pfam HMM Databases](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam30.0/Pfam-A.full.gz)

The conda_envs/ directory contains .yml files for replicating the conda environments used in our analyses.

The following R packages are required for analysis and figure generation:
* [circlize](https://jokergoo.github.io/circlize_book/book/)
* [cowplot](https://cran.r-project.org/web/packages/cowplot/index.html)
* [ComplexHeatmap](https://jokergoo.github.io/ComplexHeatmap-reference/book/)
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
* [EnhancedVolcano](https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html)
* [FSA](https://cran.r-project.org/web/packages/FSA/index.html)
* [ggpmisc](https://cran.r-project.org/web/packages/ggpmisc/index.html)
* [ggpubr](https://cran.r-project.org/web/packages/ggpubr/index.html)
* [ggrepel](https://github.com/slowkow/ggrepel)
* [nlme](https://cran.r-project.org/web/packages/nlme/index.html)
* [rcompanion](https://cran.r-project.org/web/packages/rcompanion/index.html)
* [readr](https://cran.r-project.org/web/packages/readr/index.html)
* [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
* [rstatix](https://cran.r-project.org/web/packages/rstatix/index.html)
* [tidyverse](https://www.tidyverse.org/)
* [vegan](https://cran.r-project.org/web/packages/vegan/index.html)

## Workflow

### Running antiSMASH

1. Download genomes from the eHOMD:

	    wget http://www.homd.org/ftp/HOMD_prokka_genomes/fna/ALL_genomes.fna -P ./rawData/

2. Rename the contigs before running antiSMASH

	    python3 ./scripts/python/fasta_rename.py/ ./rawData/ALL_genomes.fna ./derivedData/all_genomes_renamed_contigs.fna

3. Run antiSMASH
	
	Change thread_num to run

		mkdir ./derivedData/eHOMD_as4_out/

	    antismash -c thread_num --taxon bacteria --clusterblast --knownclusterblast --smcogs --outputfolder ./derivedData/eHOMD_as4_out/ ./derivedData/all_genomes_renamed_contigs.fna

4. Remove final.gbk to avoid future problems when scripts look for *.gbk

		rm ./derivedData/eHOMD_as4_out/*.final.gbk

### Running kallisto

1. Convert the antiSMASH GenBank-formatted output to .faa format and identify non-biosynthethic ORFs in pfam_exclude.txt for removal
	
	Change /path/to/Pfam-A.hmm to run
	
	    python3 ./scripts/python/as_gbk_convert.py ./derivedData/eHOMD_as4_out/ ./rawData/pfam_exclude.txt /path/to/Pfam-A.hmm

2. Concatenate the biosynthetic ORFs into a single file

	    for f in ./derivedData/eHOMD_as4_out/*_biosynthetic.fna; do cat $f >> ./derivedData/eHOMD_as4_out/biosynthetic_orfs_bgcs.fna; done

3. Build the kallisto index from the biosynthetic ORFs

	    kallisto index -i ./derivedData/kallisto/eHOMD_kallisto_bgc_orfs_index ./derivedData/eHOMD_as4_out/biosynthetic_orfs_bgcs.fna
    
4. Download raw sequencing reads from the iHMP
	There are multiple methods available for downloading raw sequencing reads from the iHMP. This includes the using manifests and the [portal_client](https://github.com/IGS/portal_client) or using wget to download the reads. My approach was to split the manifests into smaller chunks, extract the paths into a text file named urls.txt, and use the following code to download the reads:

		mkdir ./derivedData/reads/
		while read p; do wget -P ./derivedData/reads/ "$p"; done < urls.txt

	Most reads are compressed as .tar.bz2 files, but there are some reads that are .tar compressed files. 

	Note: SRS043422.tar.bz2 has an error in the sequencing file that required manual correction. 

5. Run fastp and kallisto to process raw reads and pseudoalign to biosynthetic ORFs
	
	Change thread_num to run
	
	Note: these scripts will remove the raw read files when done

	For .tar.bz2 files:

		python3 ./scripts/python/metagenome_read_align_bzip.py ./derivedData/reads/ ./rawData/hmp_hmp_adt_metagenome_manifest_bzip.tsv ./derivedData/kallisto/eHOMD_kallisto_bgc_orfs_index ./derivedData/kallisto/quant_output/ thread_num

	For .tar files:

		python3 ./scripts/python/metagenome_read_align_tar.py ./derivedData/reads/ ./rawData/hmp_adt_metagenome_manifest_tar.tsv ./derivedData/kallisto/eHOMD_kallisto_bgc_orfs_index ./derivedData/kallisto/quant_output/ thread_num

6. Process the completed kallisto quant output into a single count per cluster and convert to tidy format for R. Also removes the "SEQF1003_c2." header from each cluster prefix.

	    python3 ./scripts/python/kallisto_quant_processing.py ./derivedData/kallisto/quant_output/ ./derivedData/

### Running MetaPhlAn

1. Use MetaPhlAn to characterize the composition of ADT microbiomes
	
	Change thread_num and specify (g)enus or (s)pecies with single letter flag to run
	
	Note: these scripts will remove the raw read files when done

	For .tar.bz2 files:

		python3 ./scripts/python/metaphlan_classify_bzip.py ./derivedData/reads/ ./rawData/hmp_adt_metagenome_manifest_bzip.tsv ./derivedData/metaphlan/ [g/s] thread_num

	For .tar files:

		python3 ./scripts/python/metaphlan_classify_tar.py ./derivedData/reads/ ./rawData/hmp_adt_metagenome_manifest_tar.tsv ./derivedData/metaphlan/ [(g)enus/(s)pecies] thread_num

2. Process the completed MetaPhlAn output and convert into tidy format for R.
	Specify genus or species using the full word to run
	
		python3 ./scripts/python/metaphlan_table_processing.py ./derivedData/metaphlan/ ./derivedData/ [genus/species]

### Processing and extracting information from antiSMASH GenBank output

1. Process the file names for the antiSMASH GenBank-formatted files, then extract cluster number and SEQFID to associate genomes and clusters together.

		mkdir ./derivedData/eHOMD_as4_out/cluster_without_prefix_gbks/
		
		for f in ./derivedData/eHOMD_as4_out/*.gbk; do cp $f ./derivedData/eHOMD_as4_out/cluster_without_prefix_gbks/${f:40}; done

		python3 ./scripts/python/gbk_definition_extract.py ./derivedData/eHOMD_as4_out/cluster_without_prefix_gbks/ ./rawData/SEQFID_info.txt ./derivedData/cluster_SEQFIDs.tsv

2. Substitute the locus, definition, accession, and version given by antiSMASH with the cluster numbers

		python3 ./scripts/python/gbk_rename.py ./derivedData/eHOMD_as4_out/cluster_without_prefix_gbks/ ./derivedData/eHOMD_as4_out/cluster_gbks/

3. Extract BGC product type, contig edge, and domain counts from antiSMASH GenBank-formatted files

		python3 ./scripts/python/as_gbk_extract.py ./derivedData/eHOMD_as4_out/cluster_gbks/ ./derivedData/eHOMD_cluster_domain_counts.tsv

### Running BiG-SCAPE on eHOMD BGCs

1. Run BiG-Scape on all of the BGCs detected in the eHOMD

		python3 /path/to/bigscape.py --pfam_dir /path/to/pfam/databases/ -c thread_num --mibig --include_singletons --mix --hybrids-off --verbose -i ./derivedData/eHOMD_as4_out/cluster_gbks/ -o ./derivedData/bigscape/

### Identify RiPP precursor peptides from BGCs

1. Use NLPPrecursor to identify ORFs that encode RiPP precursors
	Change paths to run the following code. The models are available from the GitHub repository for NLPPrecursor ([here](https://github.com/magarveylab/NLPPrecursor/tree/master/training_data))
	
		python3 ./scripts/python/ripp_extract.py /path/to/ripp_bgc_sequence.faa /path/to/output/dir/ /path/to/models/

### Running BiG-SCAPE on oral <i>Streptococcus</i> genomes

1. Download oral <i>Streptococcus</i> genomes from GenBank using preferred method (see rawData/streptococcus_genomes_table.tsv for the list used in this study) into derivedData/streptococcus_genomes/

2. Run antiSMASH on oral <i>Streptococcus</i> genomes

	Change thread_num to run
	
		antismash -c thread_num --taxon bacteria --clusterblast --knownclusterblast --smcogs --outputfolder ./derivedData/streptococcus_genomes_antismash_output/ ./derivedData/streptococcus_genomes/

3. Run BiG-SCAPE on oral <i>Streptococcus</i> genomes

	Change thread_num to run
	
		python3 /path/to/bigscape.py --pfam_dir /path/to/pfam/databases/ -c thread_num --mibig --include_singletons --mix --hybrids-off --verbose -i ./derivedData/streptococcus_genomes_antismash_output/ -o ./derivedData/streptococcus_bigscape/

### Extract read counts for ANCOM

1. Extract read counts from kallisto run_info.json files for use in ANCOM

		python3 ./scripts/python/extract_read_counts.py ./derivedData/kallisto/quant_output/ ./derivedData/