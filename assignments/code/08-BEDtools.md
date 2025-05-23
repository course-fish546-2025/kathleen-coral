08-BEDtools
================
Kathleen Durkin
2025-05-23

- <a href="#1-convert-bam-to-bed" id="toc-1-convert-bam-to-bed">1 Convert
  bam to bed</a>
- <a href="#2-get-coverage-of-sequence-reads-on-gene-regions"
  id="toc-2-get-coverage-of-sequence-reads-on-gene-regions">2 Get coverage
  of sequence reads on gene regions</a>
- <a href="#3-intersect" id="toc-3-intersect">3 Intersect</a>
- <a href="#4-closest" id="toc-4-closest">4 Closest</a>

See assignment details
[here](https://sr320.github.io/course-fish546-2025/assignments/08-bedtools.html)

Re-download Week 6 data if necessary

``` bash
cd ../data 
curl -O https://gannet.fish.washington.edu/seashell/bu-mox/scrubbed/120321-cvBS/19F_R1_val_1_bismark_bt2_pe.deduplicated.sorted.bam 
curl -O https://gannet.fish.washington.edu/seashell/bu-mox/scrubbed/120321-cvBS/19F_R1_val_1_bismark_bt2_pe.deduplicated.sorted.bam.bai
```

We’ll also need a bed file with gene information.

``` bash
cd ../data 
curl -O https://eagle.fish.washington.edu/Cvirg_tracks/C_virginica-3.0_Gnomon_genes.bed
```

# 1 Convert bam to bed

``` bash
/home/shared/bedtools2/bin/bedtools bamtobed \
-i ../data/19F_R1_val_1_bismark_bt2_pe.deduplicated.sorted.bam \
> ../output/08-BEDtools/08-19F.bed
```

Warning This is big file, so it should *not* be staged or pushed to
Github. Instead, add it to your .gitignore

# 2 Get coverage of sequence reads on gene regions

``` bash
/home/shared/bedtools2/bin/bedtools coverage \
-a ../data/C_virginica-3.0_Gnomon_genes.bed \
-b ../output/08-BEDtools/08-19F.bed \
> ../output/08-BEDtools/08-gene-19F-coverage.out
```

``` bash
head -2 ../output/08-BEDtools/08-gene-19F-coverage.out
```

    ## NC_035780.1  13578   14594   gene-LOC111116054   0   +   65  1008    1016    0.9921260
    ## NC_035780.1  28961   33324   gene-LOC111126949   0   +   2679    4363    4363    1.0000000

# 3 Intersect

Lets grab a bed file of Transposable Elements and lncRNAs

``` bash
cd ../data
curl -O http://owl.fish.washington.edu/halfshell/genomic-databank/cgigas_uk_roslin_v1_gene.gff 
curl -O http://owl.fish.washington.edu/halfshell/genomic-databank/cgigas_uk_roslin_v1_rm.te.bed
curl -O http://owl.fish.washington.edu/halfshell/genomic-databank/cgigas_uk_roslin_v1_lncRNA.gff
```

``` bash
/home/shared/bedtools2/bin/bedtools intersect \
-a ../data/cgigas_uk_roslin_v1_gene.gff \
-b ../data/cgigas_uk_roslin_v1_rm.te.bed \
> ../output/08-BEDtools/08-gene-TE-intersect.out
```

``` bash
head -2 ../output/08-BEDtools/08-gene-TE-intersect.out
```

    ## NC_047559.1  Gnomon  gene    15715   15759   .   +   .   ID=gene-LOC109621113;Dbxref=GeneID:109621113;Name=LOC109621113;gbkey=Gene;gene=LOC109621113;gene_biotype=protein_coding
    ## NC_047559.1  Gnomon  gene    19138   19160   .   -   .   ID=gene-LOC117687066;Dbxref=GeneID:117687066;Name=LOC117687066;gbkey=Gene;gene=LOC117687066;gene_biotype=protein_coding

# 4 Closest

``` bash
/home/shared/bedtools2/bin/bedtools closest \
-a ../data/cgigas_uk_roslin_v1_lncRNA.gff \
-b ../data/cgigas_uk_roslin_v1_gene.gff \
> ../output/08-BEDtools/08-lnc-gene-closet.out
```

``` bash
head -2 ../output/08-BEDtools/08-lnc-gene-closet.out
```

    ## NC_047559.1  Gnomon  lnc_RNA 9839    11386   .   +   .   ID=rna-XR_004604272.1;Parent=gene-LOC117693020;Dbxref=GeneID:117693020,Genbank:XR_004604272.1;Name=XR_004604272.1;gbkey=ncRNA;gene=LOC117693020;model_evidence=Supporting evidence includes similarity to: 1 EST%2C and 100%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 15 samples with support for all annotated introns;product=uncharacterized LOC117693020;transcript_id=XR_004604272.1  NC_047559.1 Gnomon  gene    9839    11386   .   +   .   ID=gene-LOC117693020;Dbxref=GeneID:117693020;Name=LOC117693020;gbkey=Gene;gene=LOC117693020;gene_biotype=lncRNA
    ## NC_047559.1  Gnomon  lnc_RNA 167270  168430  .   -   .   ID=rna-XR_004601744.1;Parent=gene-LOC117689460;Dbxref=GeneID:117689460,Genbank:XR_004601744.1;Name=XR_004601744.1;gbkey=ncRNA;gene=LOC117689460;model_evidence=Supporting evidence includes similarity to: 3 long SRA reads%2C and 98%25 coverage of the annotated genomic feature by RNAseq alignments%2C including 52 samples with support for all annotated introns;product=uncharacterized LOC117689460;transcript_id=XR_004601744.1    NC_047559.1 Gnomon  gene    151758  185673  .   +   .   ID=gene-LOC117687070;Dbxref=GeneID:117687070;Name=LOC117687070;gbkey=Gene;gene=LOC117687070;gene_biotype=protein_coding
