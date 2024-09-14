# PEChyron
A repository for PE Chyron code base

# To analyze peCHYRON loci by nanopore sequencing:

Start with demultiplexed fastq data.

We use the Maple pipeline (github.com/gordonrix/maple) for demultiplexing, then convert the outputted .bam files to .fastq using the bamtofastq utility in bedtools.

First, run extract_insertion_sequences_and_report_indel_sizes_and_signatures.py.

Then, if analyzing overall efficiency of editing is necessary, run calculate_number_edits_from_nanopore_data.py.

# To analyze peCHYRON loci by Illumina sequencing:

First, we demultiplex using the separate_by_barcodes_modified.py script from the liusynevolab/CHYRON-NGS repository. Then, if desired we combine forward and reverse reads using PEAR, as described in the CHYRON-NGS repository. However, we normally use forward reads only, as the repetitive nature of the peCHYRON locus can lead to errors in the combination process.

Next, we analyze

# Please find in other repositories at github.com/liusynevolab:

Plasmid maps.

Flow cytometry data.

Analysis pipeline underlying Figure 5c.

Analysis pipeline underlying Figure 5b and Extended Data Figures 7-9.
