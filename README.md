# PEChyron
A repository for PE Chyron code base

# To analyze peCHYRON loci by nanopore sequencing:

Start with demultiplexed fastq data.

We use the Maple pipeline (github.com/gordonrix/maple) for demultiplexing, then convert the outputted .bam files to .fastq using the bamtofastq utility in bedtools.

First, run extract_insertion_sequences_and_report_indel_sizes_and_signatures.py.

Then, if analyzing overall efficiency of editing is necessary, run calculate_number_edits_from_nanopore_data.py.

# Please find in other repositories at github.com/liusynevolab:

Plasmid maps.

Flow cytometry data.

Analysis pipeline underlying Figure 5c.

Analysis pipeline underlying Figure 5b and Extended Data Figures 7-9.
