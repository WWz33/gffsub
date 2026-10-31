# 输出格式

<!-- I18N:START -->

[English](./output-formats.md) | **中文**

<!-- I18N:END -->

## -f / --format

语法：`-f FMT` 或 `--format FMT`（别名 `--output-format FMT`）。

| 值 | 输出 |
|----|------|
| `gff3` | GFF3，保留原始第 9 列（默认）|
| `gtf`、`gtf3` | GTF3 词表：`gene`、`transcript`、`exon`、`CDS`、`start_codon`、`stop_codon`、`five_prime_utr`、`three_prime_utr`、`Selenocysteine` |
| `gtf25` | GTF2.5 词表：同 GTF3，但两类 UTR 合并为 `UTR` |
| `gtf22` | GTF2.2 词表：`CDS`、`exon`、`start_codon`、`stop_codon`、`5UTR`、`3UTR`、`inter`、`inter_CNS`、`intron_CNS`；无 gene 和 transcript 行 |
| `bed` | 6 列 BED |

GTF 版本沿革与使用方：

| 版本 | 类型 | 首次发布 |
|------|------|----------|
| GTF2.2 | CDS、exon、start_codon、stop_codon、5UTR、3UTR、inter、inter_CNS、intron_CNS | 2007 年正式规范 |
| GTF2.5 | gene、transcript、exon、CDS、Selenocysteine、start_codon、stop_codon、UTR | Ensembl release 75（非正式）|
| GTF3 | gene、transcript、exon、CDS、Selenocysteine、start_codon、stop_codon、five_prime_utr、three_prime_utr | Ensembl release 82（非正式）|

Ensembl 与 GENCODE 现在发布 GTF3，StringTie、featureCounts、ballgown 都读它。GTF2.2 是唯一有正式规范的版本。

## GTF 输出

三种 GTF 共用以下规则：

- 每行都有 `gene_id`，非 gene 行还有 `transcript_id`，两者位于第 9 列最前。
- `gene_id` 与 `transcript_id` 优先取同名的第 9 列属性，没有则取 `ID` / `Parent`。Ensembl GFF3 把基因标为 `gene:gene-X`，另带 `gene_id=gene-X`，GTF 输出用 `gene-X`。
- 其余属性跟在必备两项之后，从 GFF3 做 URL 解码并加引号。
- 值转义：`"`、`\`、tab、换行、CR。
- 终止密码子紧贴同一转录本 CDS 的 3' 端时，末端 CDS 结束在密码子之前。
- 词表之外的类型丢弃。
- GTF3 与 GTF2.5 把 `mRNA` 改名为 `transcript`。

抬头按版本分别写 `##gtf-version 3`、`##gtf-version 2.5`、`##gtf-version 2.2`。

## GFF3

- 头部：`##gff-version 3`
- 第 9 列保留输入原样。
- GTF 输入时：第 9 列重写为 `tag=value`，带 `ID=`/`Parent=`，并做 URL 转义。

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
chr1	src	stop_codon	248	250	.	+	0	ID=stop01;Parent=tx01
```

GTF3（`-f gtf`）：

```
##gtf-version 3
chr1	src	gene	100	1000	.	+	.	gene_id "gene01"; Name "BRCA1";
chr1	src	transcript	100	1000	.	+	.	gene_id "gene01"; transcript_id "tx01";
chr1	src	exon	100	250	.	+	.	gene_id "gene01"; transcript_id "tx01";
chr1	src	CDS	100	247	.	+	0	gene_id "gene01"; transcript_id "tx01";
chr1	src	stop_codon	248	250	.	+	0	gene_id "gene01"; transcript_id "tx01";
```

GTF2.2（`-f gtf22`）丢掉 gene 和 transcript 行，CDS 照样裁剪：

```
##gtf-version 2.2
chr1	src	exon	100	250	.	+	.	gene_id "gene01"; transcript_id "tx01";
chr1	src	CDS	100	247	.	+	0	gene_id "gene01"; transcript_id "tx01";
chr1	src	stop_codon	248	250	.	+	0	gene_id "gene01"; transcript_id "tx01";
```

BED：

```
chr1	99	1000	gene01	0	+
chr1	99	1000	tx01	0	+
chr1	99	250	ex01	0	+
chr1	99	250	cds01	0	+
chr1	247	250	stop01	0	+
```

```bash
./gffsub demo.gff3 --format gtf
./gffsub demo.gff3 --format gtf25
./gffsub demo.gff3 --format gtf22
./gffsub demo.gff3 --format bed
./gffsub input.gtf --format gff3
./gffsub demo.gff3 --format gtf -o output.gtf
```

## --out-attrs

只输出列出的第 9 列 tag（逗号分隔）。`ID` 与 `Parent` 始终保留，保证 feature 树一致。

```bash
./gffsub demo.gff3 --out-attrs Name
./gffsub demo.gff3 --out-attrs Name,gene_id --format gtf
```

第 9 列为 `.` 的记录不受影响。
