# Region and BED Filtering

<!-- I18N:START -->

**English** | [中文](./region-bed-filtering.zh.md)

<!-- I18N:END -->

## -r / --region

Syntax: `-r CHR:START-END`

- 1-based inclusive coordinates
- Keeps records overlapping the region on the same seqid

Sample data (`demo.gff3`):

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr2	src	gene	200	600	.	-	.	ID=gene02
```

```bash
./gffsub demo.gff3 -r chr1:200-600
```

Output:

```
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
```

## -b / --bed

Syntax: `-b FILE.bed`

- BED is 0-based half-open (start included, end excluded)
- Tab-separated, at least 3 columns (chrom, start, end)
- Optional columns 4-6 (name, score, strand) parsed but not used for filtering
- Keeps records overlapping any BED interval

Example BED file (`regions.bed`):

```
chr1	99	800
chr2	199	700
```

```bash
./gffsub demo.gff3 -b regions.bed
```

## --seqid

Syntax: `--seqid LIST`

- Comma-separated seqid list
- `^` prefix excludes: `--seqid ^chr1` keeps everything except chr1
- AND with region/BED filters

```bash
./gffsub demo.gff3 --seqid chr1
./gffsub demo.gff3 --seqid chr1,chr2
./gffsub demo.gff3 --seqid ^chr1
```

## Coordinate systems

| Format | Base | Convention |
|--------|------|------------|
| GFF3/GTF | 1-based | start and end inclusive |
| BED | 0-based | start included, end excluded (half-open) |
| -r region query | 1-based | start and end inclusive |
