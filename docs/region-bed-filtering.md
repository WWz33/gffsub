# Region and BED Filtering

## Region filter: -r / --region

Syntax: `-r CHR:START-END`

- 1-based inclusive coordinates
- Keeps records overlapping the region on the same seqid
- Overlap is any base-level intersection between the record span and the query region
- Combine with `-f` (feature type) to restrict further

Sample data (save as demo.gff3):
```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr2	src	gene	200	600	.	-	.	ID=gene02
```

Commands:
1. `./gffsub demo.gff3 -r chr1:200-600` — all chr1 records overlapping 200-600
2. `./gffsub demo.gff3 -r chr1:200-600 -f gene` — only gene records in that region
3. `./gffsub demo.gff3 -r chr2:1-1000` — all chr2 records

## BED file filter: -b / --bed

Syntax: `-b FILE.bed`

- BED is 0-based half-open (start inclusive, end exclusive)
- File format: tab-separated, minimum 3 columns (chrom, start, end)
- Optional columns: name (4th), score (5th), strand (6th), parsed but not used for filtering
- Keeps records overlapping any interval in the BED file
- Multiple intervals in one file are OR'd together

Example BED file (save as regions.bed):
```
chr1	99	800
chr2	199	700
```

Commands:
1. `./gffsub demo.gff3 -b regions.bed` — records overlapping any BED interval
2. `./gffsub demo.gff3 -b regions.bed -f gene` — only gene records

## Seqid filter: --seqid

Syntax: `--seqid LIST`

- Comma-separated seqid list
- Keeps records only on listed seqids
- Prefix `^` to exclude: `--seqid ^chr1` keeps everything except chr1
- AND logic with region/BED filters

Commands:
1. `./gffsub demo.gff3 --seqid chr1` — only chr1 records
2. `./gffsub demo.gff3 --seqid chr1,chr2` — chr1 and chr2 records
3. `./gffsub demo.gff3 --seqid ^chr1` — everything except chr1

## Coordinate systems

| Format | Base | Coordinate convention |
|--------|------|-----------------------|
| GFF3/GTF | 1-based | start and end inclusive |
| BED | 0-based | start inclusive, end exclusive (half-open) |
| Region query (-r) | 1-based | start and end inclusive |

gffsub converts BED coordinates internally. Provide coordinates in the native format of each input type.
