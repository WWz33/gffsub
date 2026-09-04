# 输出格式

<!-- I18N:START -->

[English](./output-formats.md) | **中文**

<!-- I18N:END -->

## -t / --format

语法：`-t FMT` 或 `--format FMT`（别名 `--output-format FMT`）。

| 值      | 输出 |
|---------|------|
| `gff3`  | GFF3，保留原始第 9 列（默认）|
| `gtf`   | GTF2.2，gene_id + transcript_id |
| `gtf2`  | `gtf` 的别名 |
| `gtf3`  | GTF2.2.1，mRNA 改名为 transcript，应用 type 白名单 |
| `bed`   | 6 列 BED：chrom、start、end、name、score、strand |

## GFF3

- 头部：`##gff-version 3`
- 第 9 列保留输入原样。
- GTF 输入时：第 9 列重写为 `tag=value`，带 `ID=`/`Parent=`，并做 URL 转义。

## GTF

- 头部：`##gtf-version 2`
- 每行都带 `gene_id`。
- `transcript_id` 仅出现在非 gene 行。
- `gtf`/`gtf2` 中 `mRNA` 保持为 `mRNA`。
- 值转义：反斜杠、引号、tab、换行、CR。
- `gene_id` 无法解析时输出空值（`gene_id "";`）。
- 其余 GFF3 属性保留在必备对之后，URL 解码并加引号（GTF2.2 允许额外属性；AGAT 同样保留）。

## GTF3

- 头部：`##gtf-version 2.2.1`
- type 白名单：`gene`、`transcript`、`exon`、`CDS`、`start_codon`、`stop_codon`、`five_prime_utr`、`three_prime_utr`、`Selenocysteine`、`mRNA`、`5UTR`、`3UTR`、`inter`、`inter_CNS`、`intron_CNS`。
- `mRNA` 改名为 `transcript`。
- UTR 统一为 `five_prime_utr`/`three_prime_utr`。

## BED

- 6 列：chrom、start（0-based）、end、name、score、strand。
- `start` = 记录起始 - 1。
- `score` 取第 6 列，缺失或 `.` 时为 `0`。
- `name` 取 `ID`，缺失时取 type。

## 示例

输入样例（demo.gff3）：
```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
```

GFF3（默认）：
```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
...
```

GTF：
```
##gtf-version 2
chr1	src	gene	100	1000	.	+	.	gene_id "gene01";
chr1	src	mRNA	100	1000	.	+	.	gene_id "gene01"; transcript_id "tx01";
chr1	src	exon	100	250	.	+	.	gene_id "gene01"; transcript_id "tx01";
chr1	src	CDS	100	250	.	+	0	gene_id "gene01"; transcript_id "tx01";
```

BED：
```
chr1	99	1000	gene01	0	+
chr1	99	1000	tx01	0	+
chr1	99	250	ex01	0	+
chr1	99	250	cds01	0	+
```

```bash
./gffsub demo.gff3 --format gtf
./gffsub demo.gff3 --format gtf3
./gffsub demo.gff3 --format bed
./gffsub input.gtf --format gff3
./gffsub demo.gff3 --format gtf -o output.gtf
```
