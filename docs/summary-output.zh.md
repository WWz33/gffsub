# summary 输出

<!-- I18N:START -->

[English](./summary-output.md) | **中文**

<!-- I18N:END -->

## -s, --summary

输出过滤后记录的 TSV 统计摘要，而非 feature 记录。可用于主子集路径和 `query` 子命令。计数（`child_count`、`transcript_count`、`exon_count`、`cds_length`）基于完整注释索引计算，不受过滤影响，父子关系保持完整。

## TSV 列

`seqid`、`start`、`end`、`strand`、`type`、`length`、`child_count`、`transcript_count`、`exon_count`、`cds_length`。

- `length` = `end - start + 1`
- `child_count` = 直接子记录数（`Parent` 指向本记录 ID 的记录数）
- `transcript_count` = mRNA/transcript 后代数
- `exon_count` = exon 后代数
- `cds_length` = CDS 后代长度之和

值中的 tab、换行、CR 转义为 `\t`、`\n`、`\r`。

## all 行

过滤后记录涉及多于一个不同 `seqid` 时，末尾输出 `all` 行，数值列求和。其 `start`、`end`、`strand`、`type` 为 `NA`，以区分逐记录的 GFF 行。若所有记录同一 `seqid`，不输出 `all` 行。

## 示例

输入样例（demo.gff3）：
```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1;biotype=protein_coding
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
chr2	src	gene	500	600	.	-	.	ID=gene02;Name=XYZ1
```

```bash
./gffsub demo.gff3 -f gene -s
```

输出：
```tsv
seqid	start	end	strand	type	length	child_count	transcript_count	exon_count	cds_length
chr1	100	1000	+	gene	901	1	1	2	402
chr2	500	600	-	gene	101	0	0	0	0
all	NA	NA	NA	NA	1002	1	1	2	402
```
