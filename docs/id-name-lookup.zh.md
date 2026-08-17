# Feature ID 与 Gene Name 查找

<!-- I18N:START -->

[English](./id-name-lookup.md) | **中文**

<!-- I18N:END -->

示例数据（`demo.gff3`）：

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1;Alias=BRCA-1;Dbxref=GeneID:672
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	mRNA	100	800	.	+	.	ID=tx02;Parent=gene01
chr1	src	exon	100	800	.	+	.	ID=ex02;Parent=tx02
```

## -i / --id

语法：`--id ID`（可重复）

- feature ID 精确匹配
- 多个 `--id` 之间为 OR 关系

```bash
./gffsub demo.gff3 -i gene01
./gffsub demo.gff3 -i tx01 -i tx02
./gffsub demo.gff3 -i gene01 -C
```

## --ids

语法：`--ids FILE`

- 每行一个 ID
- 空行与以 `#` 开头的行跳过

示例文件（`ids.txt`）：

```
gene01
tx01
```

```bash
./gffsub demo.gff3 --ids ids.txt -C
```

## -n / --name

语法：`-n NAME`

- 按多个 key 查找 gene：`ID`、`gene_id`、`Name`、`locus_tag`、`Alias`、`Dbxref`
- 仅匹配 gene 类型 feature

```bash
./gffsub demo.gff3 -n BRCA1
./gffsub demo.gff3 -n BRCA-1
./gffsub demo.gff3 -n GeneID:672
./gffsub demo.gff3 -n gene01 -C
```

## -N / -N REGION

语法：`--nearest CHR:START-END`

- 在同一 seqid 上查找距给定 1-based region 最近的 gene

```bash
./gffsub demo.gff3 -N chr1:1500-2000
```

另见：[Gene Model 展开](gene-model-expansion.zh.md)、[Summary 输出](summary-output.zh.md)。
