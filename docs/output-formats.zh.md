# 输出格式

## -t / --format 标志

语法：`-t FMT` 或 `--format FMT`（别名：`--output-format FMT`）

| 格式 | 值 | 说明 |
|--------|-------|-------------|
| GFF3 | `gff3` | GFF3，保留原始第 9 列属性（默认） |
| GTF | `gtf` | GTF2.2，含 gene_id 和 transcript_id |
| GTF2 | `gtf2` | gtf 的别名 |
| GTF3 | `gtf3` | GTF2.2.1，带类型重命名（mRNA 改为 transcript） |
| BED | `bed` | 6 列 BED（chrom、start、end、name、score、strand） |

## GFF3 输出

- 头部：`##gff-version 3`
- 第 9 列原样保留
- 按输入顺序输出记录

## GTF 输出

- 头部：`##gtf-version 2`
- 每条 feature 行含 gene_id（GTF2.2 要求）
- transcript_id 仅出现在非 gene feature 上（GTF2.2：所有行含 gene_id，除 gene 外所有行含 transcript_id）
- gene_id 通过 Parent 链解析：子 feature（exon、CDS）经内部 mRNA-to-gene 映射查其 mRNA 对应的 gene
- 无法解析 gene_id 时输出空值（`gene_id "";`），不丢弃该 feature
- gtf/gtf2 中 mRNA 类型保持 "mRNA"，gtf3 中重命名为 "transcript"
- 属性值做转义（反斜杠和引号）

## GTF3 输出

- 头部：`##gtf-version 2.2.1`
- 只通过 GTF3 兼容的 feature 类型：gene、transcript、exon、CDS、start_codon、stop_codon、five_prime_utr、three_prime_utr、Selenocysteine、mRNA
- mRNA 重命名为 transcript

## BED 输出

- 6 列：chrom、start（0-based）、end、name、score、strand
- start = 记录起始 - 1（GFF3 为 1-based，BED 为 0-based）
- score 取第 6 列，缺失或为 `.` 时取 `0`
- name 取 ID 属性，缺失时取 feature 类型

## 格式转换示例

示例输入（demo.gff3）：
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

## 命令

```bash
# GFF3 转 GTF
./gffsub demo.gff3 -t gtf

# GFF3 转 GTF3（mRNA 重命名为 transcript）
./gffsub demo.gff3 -t gtf3

# GFF3 转 BED
./gffsub demo.gff3 -t bed

# GTF 转 GFF3
./gffsub input.gtf -t gff3

# 输出到文件
./gffsub demo.gff3 -t gtf -o output.gtf
```
