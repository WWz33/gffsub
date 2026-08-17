# Gene Model 展开

gffsub 根据 ID 和 Parent 属性构建父/子索引。三个 flag 沿层级展开选择器命中的记录。

## 层级

```
gene
  └── mRNA / transcript
        ├── exon
        └── CDS
```

索引沿 `Parent=` 链追溯。GTF 输入使用由 gene_id 和 transcript_id 合成的父/子链。

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
```

## --children / -C：包含后代

递归包含选择器命中的全部子 feature。

```bash
# gene01 加其全部 mRNA、exon、CDS
./gffsub demo.gff3 --id gene01 -C

# tx01 加其 exon 和 CDS
./gffsub demo.gff3 --id tx01 -C
```
--id tx01 -C 返回：tx01、ex01、cds01、ex02、cds02。
--id gene01 -C 返回：gene01、tx01、ex01、cds01、ex02、cds02、tx02、ex03。

## --parents：包含祖先

沿 Parent 链向上走到根（gene）。

```bash
# exon 及其父 mRNA 和祖父 gene
./gffsub demo.gff3 --id ex01 --parents
```
返回：ex01、tx01、gene01。

## --model：完整 gene model

包含完整 gene model：gene、命中的选择器记录、全部兄弟（同 gene 的其他 transcript）及它们的所有子 feature。

```bash
# 包含 tx01 的完整 model
./gffsub demo.gff3 --id tx01 --model
```
返回：gene01、tx01、ex01、cds01、ex02、cds02、tx02、ex03。

## 对比

| Flag | 方向 | 包含 |
|------|-----------|----------|
| `-C` / `--children` | 向下 | 命中记录 + 全部后代 |
| `--parents` | 向上 | 命中记录 + 全部祖先到 gene |
| `--model` | 双向 | 完整 gene model（gene + 全部 transcript + 全部子 feature） |

## 与多个选择器组合

与 --name 或多个 --id 组合时：
- 每个选择器命中独立展开
- 结果按行索引去重

```bash
# 两个 gene 的完整 model
./gffsub demo.gff3 --id gene01 --id gene02 --model
```

## 与 region 筛选组合

Region、seqid 和 feature 类型筛选先执行。Gene model 展开在筛选之后进行，因此可能拉入被 region 筛选排除的记录。
