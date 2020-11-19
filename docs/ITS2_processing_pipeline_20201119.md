# ITS2_ProcessingPipeline
## Germaine Yong
### November 18, 2020

This workflow processes Miseq raw sequence data into a non-normalized ITS2 operational taxonomic unit (OTU) table for further downstream analyses. Taxanomy is assigned to the output ITS2 OTUs using the UNITE database. Note that samples were processed and analyzed on lynchserver2; feel free to adapt to DADA2's ITS2 processing pipeline instead which is newer: https://benjjneb.github.io/dada2/ITS_workflow.html

#### Download data from Miseq and move to directory with raw fastq reads

```
sf166265@lynchserver2:/data$ cp NextSeq_data/MiSeq_Data/160930_M01869_0166_000000000-AW1MA/Data/Intensities/BaseCalls/*.fastq Users/sf166265/MiSeq_20160930
cd Users/sf166265/MiSeq_20160930
```

#### Trim the R1 and R2 reads of Illumina adapter sequence. Since this is a 290bp x 290bp run some of the sequences will read completely through the amplicon and into the adapter. R1 and R2 have different target sequences for trimming as noted below.

```
~/.local/bin/cutadapt -a GCATATCAATAAGC -q 25 -e .2 -n 3 -O 1 Undetermined_S0_L001_R1_001.fastq > trimmed_R1_v2.fastq
~/.local/bin/cutadapt -a GCATATCAATAAGC -q 25 -e .2 -n 3 -O 1 Undetermined_S0_L001_R2_001.fastq > trimmed_R2_v2.fastq
```

#### Assemble adapter-trimmed R1 and R2 sequences with a minimum overlap of 25bp.
```
flash -m 25 -M 290 -r 290 trimmed_R1.fastq trimmed_R2.fastq -o pool1
```

#### Create a .txt file containing only the names of the sequences that assembled during FLASh
```
sed -n '1~4'p pool1.extendedFrags.fastq > FLAShReads_pool1.txt
```

#### Filter the index file to include only the barcodes for the sequences that assembled. Rename this file (FLAShReads_pool1.txt) if you have more than one pool of amplicons.
```
python /data/Users/dfadrosh/scripts/filter_fasta_keep_order.py Index.fastq FLAShReads_pool1.txt
```

#### Demultiplex samples without quality filtering and output in .fastq format; I previously used QIIME1 on lynchserver2 to do so, you can switch to QIIME2 if you prefer
```
split_libraries_fastq.py -i pool1.extendedFrags.fastq -b Index_filtered_ordered.fastq -o splibs -m OBESITY_ITS2_Run1_MappingFile_wBMI.txt --rev_comp_mapping_barcodes --store_demultiplexed_fastq -r 0 -q 0 -n 100 -p 0.1 --barcode_type 12
```

#### If you have more than one pool, remove adapters and demultiplex each sequence run separately then concatenate the non-quality filtered fastq files here.

#### Quality filter based on expected errors and remove any sequence expected to have 2 or more errors.
```
usearch -fastq_filter splibs/seqs.fastq -fastq_maxee 2 –-eeout -fastaout quality_seqs.fasta
```

#### Rename the demultiplexed quality filtered sequences to tidy up the header line.
```
perl /data/Users/dfadrosh/scripts/drive5/bmp_uparse_pipeline.pl –i quality_seqs.fasta -o reads_uparse.fa
```

#### Full length dereplicate sequences and append sequence counts to header line
```
usearch --derep_fulllength reads_uparse.fa -sizeout -output derep.fa
```

#### Sort the dereplicated sequences based on sequence count with the most counts first and remove any representatives that have a count of less than 2 (i.e. singletons) 
```
/data/Software/usearch8.0.1623_i86linux32 -sortbysize derep.fa –fastaout sorted.fa -minsize 2
```

#### Cluster sequences into OTUs using the uparse pipeline. A sequence similarity based chimera checkin (not vs. a database) is done during this step as an initial chimera check pass.
```
/data/Software/usearch8.0.1623_i86linux32 -cluster_otus sorted.fa –sizein -sizeout -otus otus.fa -uparseout otu_clustering.txt
```

#### Check for chimeras against the Unite database. I had used the UNITE 2015 database, as the newer database had special characters that were messing with taxonomy assignments.
```
/data/Software/usearch8.0.1623_i86linux32 -uchime_ref otus.fa -db /data/Users/dfadrosh/scripts/UNITE/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta -strand plus -nonchimeras otus_nonchimeras.fa -chimeras otus_chimeras.fa
```

#### Rename OTUs OTU_1, OTU_2, etc.
```
python /data/Users/dfadrosh/scripts/drive5/fasta_number.py otus_nonchimeras.fa OTU_ > otus_nonchimeras_numbered.fa
```

#### Map demultiplexed quality filtered reads to the OTU database.
```
usearch -usearch_global reads_uparse.fa -db otus_nonchimeras_numbered.fa -strand plus -id 0.97 -uc map.uc
```

#### Construct tab delimited OTU table and convert txt OTU table to biom format
```
python /data/Users/dfadrosh/scripts/drive5/uc2otutab.py map.uc > otu_table.txt
biom convert -i otu_table.txt -o otu_table_uparse.biom --to-hdf5 --table-type="OTU table"
```

#### Assign taxonomy based on May 2015 UNITE database. There was a more current version of the UNITE database (2016; check for futher updates since?) but that needed to be modified due to special characters. It did not work for me following modification/removal of special characters.
```
assign_taxonomy.py -m rdp -c 0.8 -i otus_nonchimeras_numbered.fa -o otus_nonchimeras_numbered_tax_01_08_2015_UNITE -r /data/Users/dfadrosh/scripts/UNITE/01_08_2015_release/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta -t
/data/Users/dfadrosh/scripts/UNITE/01_08_2015_release/sh_taxonomy_qiime_ver7_dynamic_01.08.2015.txt
```

#### Add taxonomy information to the .biom OTU table. This table contains all ITS2 data including non-fungal data. This may be useful for some applications. Downstream processing cleans up the dataset to include only OTUs know to be of fungal origin.
```
biom add-metadata -i otu_table_uparse.biom -o otu_table_tax.biom --observation-metadata-fp otus_nonchimeras_numbered_tax_01_08_2015_UNITE/otus_nonchimeras_numbered_tax_assignments.txt --observation-header OTU,rdp_taxonomy,rdp_confidence --float-fields confidence
```

#### Use ITSx extractor to identify ITS2 regions belonging to fungi. Note: -t F option only returns fungal reads.
```
perl /data/Users/dfadrosh/scripts/ITSx_1.0.11/ITSx -i otus_nonchimeras_numbered.fa --reset --cpu 60 -t F
```

#### Use text manipulator “sed” to create a list of OTU ids that were determined to be of fungal origin.
```
sed "s/[|].*//" ITSx_out.ITS2.fasta > ITSx_out.ITS2_OTU.fasta
sed 's/>//' ITSx_out.ITS2_OTU.fasta > ITS2.otu.list.txt
```

#### Filter .biom OTU table to include only OTUs of fungal origin and discard very low abundance OTUs (less than 1/1000th of a percent of the total read count for all samples).
```{r}
filter_otus_from_otu_table.py -i otu_table_uparse.biom -o otu_table_ITS2.biom -e ITS2.otu.list.txt --negate_ids_to_exclude --min_count_fraction 0.00001
```

#### Summarize OTU table to get information on sequence depth per sample, which will be used for rarefaction later.
```{r}
biom summarize-table -i otu_table_ITS2.biom -o otu_table_tax_ITS2_SUMMARY.txt
```

#### Reassign taxonomy so that OTUs are numerical in order, append taxonomy information to ITS2-only OTU table
```{r}
assign_taxonomy.py -m rdp -c 0.8 -i ITSx_out.ITS2.fasta -o otus_nonchimeras_numbered_tax_01_08_2015_UNITE_ITS2_ONLY -r /data/Users/dfadrosh/scripts/UNITE/01_08_2015_release/sh_refs_qiime_ver7_dynamic_01.08.2015.fasta -t
/data/Users/dfadrosh/scripts/UNITE/01_08_2015_release/sh_taxonomy_qiime_ver7_dynamic_01.08.2015.txt

biom add-metadata -i otu_table_ITS2.biom -o otu_table_ITS2_tax_2015.biom --observation-metadata-fp otus_nonchimeras_numbered_tax_01_08_2015_UNITE_ITS2_ONLY/ITSx_out.ITS2_tax_assignments.txt --observation-header OTU,ITS2_rdp_taxonomy,ITS2_rdp_confidence --float-fields confidence
```

#### Convert .biom OTU table with taxonomy info to tab delimited text format for downstream analysis
```{r}
biom convert -i otu_table_ITS2_tax_2015.biom -o otu_table_ITS2_tax_2015.txt --to-tsv --header-key taxonomy
```

Proceed to normalizing for differences in read depth across samples. This can be done through either Multiply Rarefying or Variance Stabilized Transformations. Refer to 16S_processing_pipeline section "Normalizing for Differences in Read Depth" for more information: https://lynchlab-ucsf.github.io/docs/16s_processing_pipeline_06Nov20.html 
