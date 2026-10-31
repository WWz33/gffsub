# 基因模型扩展

<!-- I18N:START -->

[English](./gene-model-expansion.md) | **中文**

<!-- I18N:END -->

gffsub 沿 `Parent=` 链向下扩展选择器命中。GTF 通过 `gene_id` / `transcript_id` 合成链。

示例数据（demo.gff3）：

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
chr1	src	mRNA	100	800	.	+	.	ID=tx02;Parent=gene01
chr1	src	exon	100	800	.	+	.	ID=ex03;Parent=tx02
chr2	src	gene	200	600	.	-	.	ID=gene02;Name=XYZ1
chr2	src	mRNA	200	600	.	-	.	ID=tx03;Parent=gene02
chr2	src	exon	200	400	.	-	.	ID=ex04;Parent=tx03
```

## -C / --children

包含选择器命中及其所有后代。

```bash
./gffsub demo.gff3 -i gene01 -C
```

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
chr1	src	mRNA	100	800	.	+	.	ID=tx02;Parent=gene01
chr1	src	exon	100	800	.	+	.	ID=ex03;Parent=tx02
```

选择 transcript 时到该 transcript 为止，不再向下：

```bash
./gffsub demo.gff3 -i tx01 -C
```

```
##gff-version 3
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
```

## -p / --parents

包含选择器命中及其所有祖先，上溯到 gene。

```bash
./gffsub demo.gff3 -i ex01 -p
```

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
```

## -m / --model

包含完整 gene 模型: gene、选择器命中、所有 sibling transcript 及其所有子代。

```bash
./gffsub demo.gff3 -i tx01 -m
```

输出与上面的 `-i gene01 -C` 相同，因为 tx01 属于 gene01，而 gene01 没有其他子代。

## 对比

| Flag | 方向 | 包含 |
|------|------|------|
| `-C` / `--children` | 向下 | 命中 + 所有后代 |
| `-p` | 向上 | 命中 + 所有祖先到 gene |
| `-m` | 双向 | 完整 gene 模型 |

## --drop-orphans

删除所有 Parent 引用都已不在保留集内的记录。在其他过滤之后执行，迭代到不动点，因此删除 transcript 会连带删除其 exon。无 Parent 的记录不会删除。

```bash
./gffsub demo.gff3 -t exon --drop-orphans
```

```
##gff-version 3
```

`-t exon` 删掉了 transcript 行，所有 exon 因此失去 Parent，不动点把整批删光。
