# summary 输出

<!-- I18N:START -->

[English](./summary-output.md) | **中文**

<!-- I18N:END -->

## --summary FMT

取值：`tsv`、`json`。输出表格行而非 feature 记录。仅可用于查询式 selector（`--id`、`--ids`、`--name`、`--where`、`--nearest`、query 子命令 `--region`）。不可与 `--grep`、`-I`、`-E`、`--seqid`、`--source`、`--score`、`--strand`、`--phase`、`--bed`、`--longest`、`--threads`、`--format`/`--output-format`、`--output` 同时使用。

## --out-attrs KEYS

逗号分隔的 attribute key，作为额外列追加（别名 `--output-attrs`、`--attrs`）。先从记录字段解析 `gene_id`/`transcript_id`，再查 GFF3 attribute，再查 GTF 引号值。

## TSV 列

`query_id`、`matched_id`、`matched_by`、`seqid`、`start`、`end`、`strand`、`type`、`parent_id`、`child_count`、`transcript_count`、`exon_count`、`cds_length`、`status`，外加每个 `--out-attrs` key 一列。值中的 tab、换行、CR 转义为 `\t`、`\n`、`\r`。

## JSON 格式

对象数组，字段同 TSV，另加 `attrs` 对象承载 `--out-attrs` 的 key。

## status 字段

`found` 或 `not_found`。未命中时，`matched_by` 保留查询 key，数值计数为 0，其余字段为空。

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
```

```bash
./gffsub demo.gff3 --id gene01 --summary tsv
./gffsub demo.gff3 --id gene01 --summary json --out-attrs Name,biotype
```

TSV 输出：
```tsv
query_id	matched_id	matched_by	seqid	start	end	strand	type	parent_id	child_count	transcript_count	exon_count	cds_length	status
gene01	gene01	ID	chr1	100	1000	+	gene		1	1	2	402	found
```

JSON 输出：
```json
[
  {"query_id":"gene01","matched_id":"gene01","matched_by":"ID","seqid":"chr1","start":100,"end":1000,"strand":"+","type":"gene","parent_id":"","child_count":1,"transcript_count":1,"exon_count":2,"cds_length":402,"status":"found","attrs":{"Name":"BRCA1","biotype":"protein_coding"}}
]
```
