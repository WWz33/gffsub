# 基因模型扩展

<!-- I18N:START -->

[English](./gene-model-expansion.md) | **中文**

<!-- I18N:END -->

gffsub 沿 `Parent=` 链向下扩展选择器命中。GTF 通过 `gene_id` / `transcript_id` 合成链。

示例数据 (demo.gff3):

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

## -C / --children

包含选择器命中及其所有后代。

```bash
# gene01 + 所有 mRNA、exon、CDS
./gffsub demo.gff3 -i gene01 -C

# tx01 + 其 exon 和 CDS
./gffsub demo.gff3 -i tx01 -C
```

`-i gene01 -C` 返回: gene01、tx01、ex01、cds01、ex02、cds02、tx02、ex03。

`-i tx01 -C` 返回: tx01、ex01、cds01、ex02、cds02。

## -p / -p

包含选择器命中及其所有祖先，上溯到 gene。

```bash
./gffsub demo.gff3 -i ex01 -p
```

返回: ex01、tx01、gene01。

## -m / -m

包含完整 gene 模型: gene、选择器命中、所有 sibling transcript 及其所有子代。

```bash
./gffsub demo.gff3 -i tx01 -m
```

返回: gene01、tx01、ex01、cds01、ex02、cds02、tx02、ex03。

## 对比

| Flag | 方向 | 包含 |
|------|------|------|
| `-C` / `-C` | 向下 | 命中 + 所有后代 |
| `-p` | 向上 | 命中 + 所有祖先到 gene |
| `-m` | 双向 | 完整 gene 模型 |
