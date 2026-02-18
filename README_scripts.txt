Respiratory mNGS Pipeline Scripts
==================================

This directory contains all Python and bash scripts used by the respiratory mNGS pipeline.
These scripts are self-contained copies from the main kirsty_scripts directories.

Scripts:
--------
1. add_taxon.sh - Adds taxonomic annotation to BLAST results
2. best_ref_mergedblastn.py - Determines best reference sequences for viral contigs
3. blastn_filter.py - Filters BLASTN results
4. diamond_m8_filter.py - Filters DIAMOND m8 output files
5. pileup2consensus.py - Generates consensus sequences from pileup summaries
6. pileup_summary.py - Summarizes samtools mpileup output
7. scaffold_nucmer_coords.py - Scaffolds contigs using nucmer alignments
8. summarise_contigs_per_taxon.py - Summarizes viral contigs per taxonomic group
9. trim_contig.py - Trims contigs based on read depth
10. trim_viral_contigs.py - Trims viral contigs (initial processing)

Configuration:
--------------
All pipeline workflows now use a single config parameter to locate these scripts:

    scripts_dir: /home3/2893911k/kirsty_scripts/smk/respiratory_mNGS/scripts

This replaces the previous individual script path parameters.

Maintenance:
------------
If scripts need to be updated, update the files in this directory.
The pipeline is now self-contained and independent of the main kirsty_scripts directories.

Last updated: 2025-10-28
