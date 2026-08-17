# Feature ID 与基因名查找

## 按 ID：--id

语法：`--id ID`（可重复）

- 精确匹配 feature 的 ID 属性（GFF3 第 9 列 `ID=`，GTF 中由 gene_id/transcript_id 合成）
- 多个 --id 为 OR 关系：匹配任一 ID 的记录均保留
- 不加 -C/--children/--parents/--model 时仅返回精确匹配的 ID

示例数据（demo.gff3）：
```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1;Alias=BRCA-1;Dbxref=GeneID:672
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01;Name=BRCA1.1
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	mRNA	100	800	.	+	.	ID=tx02;Parent=gene01
chr1	src	exon	100	800	.	+	.	ID=ex02;Parent=tx02
```

命令：
1. `./gffsub demo.gff3 --id gene01` — 仅 gene 记录
2. `./gffsub demo.gff3 --id tx01 --id tx02` — 两条 mRNA 记录
3. `./gffsub demo.gff3 --id gene01 -C` — gene 加全部后代（见 [Gene Model 展开](gene-model-expansion.md)）

## 按 ID 文件：--ids

语法：`--ids FILE`

- 每个非空行一个 ID
- 等价于多个 --id
- 空行和以 # 开头的行跳过

示例文件（ids.txt）：
```
gene01
tx01
```
命令：`./gffsub demo.gff3 --ids ids.txt -C`

## 按基因名：--name

语法：`--name NAME`

- 按 gene 记录的多个命名键查找，不只是 ID 属性
- 查找键：ID、gene_id、Name、locus_tag、Alias、Dbxref（任一匹配即可；冲突时取文件顺序中第一个 gene）
- 仅匹配 gene 类型 feature（type == "gene"）
- 返回 gene 记录；加 -C 或 --model 展开

命令：
1. `./gffsub demo.gff3 --name BRCA1` — Name=BRCA1 的 gene
2. `./gffsub demo.gff3 --name BRCA-1` — Alias=BRCA-1 的 gene
3. `./gffsub demo.gff3 --name GeneID:672` — Dbxref=GeneID:672 的 gene
4. `./gffsub demo.gff3 --name gene01 -C` — 按 ID 查 gene，带子 feature

## query 子命令

`query` 子命令是 ID/name 查找的首选接口，支持汇总输出。

```
./gffsub query demo.gff3 --id gene01 -C
./gffsub query demo.gff3 --name BRCA1 --summary tsv
./gffsub query demo.gff3 --ids ids.txt --summary json
```

另见：[Gene Model 展开](gene-model-expansion.md)、[汇总输出](summary-output.md)。
