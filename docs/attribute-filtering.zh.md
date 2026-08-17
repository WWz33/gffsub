# 属性筛选

<!-- I18N:START -->

[English](./attribute-filtering.md) | **中文**

<!-- I18N:END -->

## --where：属性精确匹配

语法：`--where KEY=VALUE`（可重复）

- 精确匹配 GFF3 第 9 列属性
- GTF 输入下可访问 gene_id 和 transcript_id（由解析器合成）
- 多个 --where 之间是 OR 逻辑（任一匹配即保留该记录）
- 属性值在比较前做 URL 解码（GFF3 规范）：文件中的 `Note=a%2Cb` 用解码后的形式 `--where Note=a,b` 匹配

示例数据（demo.gff3）：
```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1;biotype=protein_coding;Dbxref=GeneID:672
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01;biotype=protein_coding
chr2	src	gene	200	600	.	-	.	ID=gene02;Name=TP53;biotype=tumor_suppressor
```

命令：
1. `./gffsub demo.gff3 --where biotype=protein_coding -f gene` — biotype 为 protein_coding 的 gene
2. `./gffsub demo.gff3 --where biotype=tumor_suppressor` — 该 biotype 的所有记录
3. `./gffsub demo.gff3 --where biotype=protein_coding --where Name=BRCA1` — 任一条件匹配即保留

## --grep：子串搜索

语法：`--grep FIELD:PATTERN`（可重复）

- 在指定 FIELD 中搜索子串 PATTERN
- 字段：列字段（seqid、source、type、start、end、score、strand、phase）和第 9 列属性键（ID、Name、Parent、gene_id、transcript_id、biotype 等）
- 不在内置列表中的键需加 `attr.` 前缀
- 多个 --grep 之间是 OR 逻辑（任一匹配即保留）
- `-v` / `--invert-match` 取反：保留不匹配的记录

命令：
1. `./gffsub demo.gff3 --grep type:gene` — gene（type 列含子串 "gene"）
2. `./gffsub demo.gff3 --grep Name:BRCA` — Name 属性含 "BRCA" 的记录
3. `./gffsub demo.gff3 --grep type:mRNA -v` — 排除 mRNA 记录外的所有记录

## --grep-regex：正则搜索

语法：`--grep-regex FIELD:REGEX`（可重复）

- ECMAScript 正则语法（std::regex）
- 字段列表同 --grep
- 多个标志之间是 OR 逻辑

命令：
1. `./gffsub demo.gff3 --grep-regex ID:ex0.` — ID 匹配该正则的记录
2. `./gffsub demo.gff3 --grep-regex Name:^BRCA` — 以 BRCA 开头的 Name

## --grep-file：批量模式

语法：`--grep-field FIELD --grep-file FILE [--grep-file-regex]`

- 从 FILE 读取模式（每行一个）
- 对 --grep-field 指定的单个 FIELD 匹配
- 不加 --grep-file-regex：子串匹配
- 加 --grep-file-regex：正则匹配

示例文件（patterns.txt）：
```
BRCA
TP53
```
命令：`./gffsub demo.gff3 --grep-field Name --grep-file patterns.txt`

## --ignore-case

语法：`--ignore-case`

- 使 --grep、--grep-regex、--grep-file、-I、-E 不区分大小写

命令：`./gffsub demo.gff3 --grep type:mrna --ignore-case`

## 字段参考

| 字段 | 来源 |
|-------|--------|
| seqid | 第 1 列 |
| source | 第 2 列 |
| type | 第 3 列 |
| start | 第 4 列 |
| end | 第 5 列 |
| score | 第 6 列 |
| strand | 第 7 列 |
| phase | 第 8 列 |
| ID | 第 9 列属性 |
| Name | 第 9 列属性 |
| Parent | 第 9 列属性 |
| Alias | 第 9 列属性 |
| Dbxref | 第 9 列属性 |
| Note | 第 9 列属性 |
| gene_id | 第 9 列属性 |
| transcript_id | 第 9 列属性 |
| biotype | 第 9 列属性 |
| locus_tag | 第 9 列属性 |
| attr.KEY | 第 9 列属性（任意键） |

参见 [表达式筛选](expression-filtering.zh.md) 了解数值比较、逻辑运算符和复合谓词。
