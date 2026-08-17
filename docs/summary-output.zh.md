# summary 输出

<!-- I18N:START -->

[English](./summary-output.md) | **中文**

<!-- I18N:END -->

## -s, --summary

输出过滤后记录的聚合 TSV 统计，而非 feature 记录。可用于主子集路径和 `query` 子命令。

## TSV 列

`seqid`、`type`、`count`、`sum_len`、`min_len`、`avg_len`、`max_len`、`Q1`、`Q2`、`Q3`、`coverage`。

每个不同 (seqid, type) 一行。`length` = `end - start + 1`。`avg_len` 截断到一位小数。`Q1`/`Q2`/`Q3` 为 Tukey hinges（下/上半区中位数）。`coverage` 为去重覆盖碱基数：组内区间重叠或首尾相接时合并（bedtools merge 语义），重叠记录共享的碱基（如两条 isoform 共用的 exon）只计一次。按 seqid 再按 type 排序。

## all 行

过滤后记录涉及多于一个不同 `seqid` 时，每个 type 追加一行 `all`：`count`/`sum_len` 为总和，`min_len`/`max_len` 取全局极值，`avg_len` 为整体均值，`Q1`/`Q2`/`Q3` 基于全部长度计算，`coverage` 为各 seqid union 之和（不同 seqid 的区间不互相合并）。若所有记录同一 `seqid`，不输出 `all` 行。

## 示例

输入样例（demo.gff3）：
```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
chr2	src	gene	500	600	.	-	.	ID=gene02;Name=XYZ1
```

```bash
./gffsub demo.gff3 -s
```

输出：
```tsv
seqid	type	count	sum_len	min_len	avg_len	max_len	Q1	Q2	Q3	coverage
chr1	CDS	2	402	151	201	251	151	201	251	402
chr1	exon	2	402	151	201	251	151	201	251	402
chr1	gene	1	901	901	901	901	901	901	901	901
chr1	mRNA	1	901	901	901	901	901	901	901	901
chr2	gene	1	101	101	101	101	101	101	101	101
all	CDS	2	402	151	201	251	151	201	251	402
all	exon	2	402	151	201	251	151	201	251	402
all	gene	2	1002	101	501	901	101	501	901	1002
all	mRNA	1	901	901	901	901	901	901	901	901
```
