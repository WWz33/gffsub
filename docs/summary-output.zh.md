# 摘要输出

## 概述

gffsub 不输出注释记录，而是输出 selector 命中结果的表格摘要。便于检查 gene 结构细节。

## --summary 选项

语法：`--summary FMT`

| 格式 | 值 |
|--------|-------|
| TSV | `tsv` |
| JSON | `json` |

默认：GFF3 注释输出（不输出摘要）。

## --out-attrs：选定的属性

语法：`--out-attrs KEYS`（逗号分隔）

别名：`--output-attrs`、`--attrs`

在摘要输出中把选定的第 9 列属性作为额外列追加。KEYS 是逗号分隔的属性名列表。

## TSV 列

| 列 | 含义 |
|--------|-------------|
| query_id | 查询的 ID 或 name |
| matched_id | 命中的 feature ID |
| matched_by | 命中所用的 lookup key（ID、Name、gene_id 等） |
| seqid | 第 1 列 |
| start | 第 4 列 |
| end | 第 5 列 |
| strand | 第 7 列 |
| type | 第 3 列（feature type） |
| parent_id | parent 的 feature ID |
| child_count | 直接子节点数 |
| transcript_count | 命中 feature 下 transcript 类型的后代数 |
| exon_count | 命中 feature 下的 exon 后代数 |
| cds_length | 命中 feature 子树的 CDS 总长度 |
| status | found 或 not_found |
| (extra) | 每个 --out-attrs key 对应一列 |

## JSON 输出

对象数组。每个对象的字段与 TSV 列相同，另外包含一个 `attrs` 对象，存放请求的属性 key。

字段值中的控制字符按 RFC 8259 转义（backspace、formfeed、newline、carriage return、tab 使用短转义；其他 0x20 以下字符使用 \uXXXX）。

## 示例数据

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1;biotype=protein_coding
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
```

## 命令

```bash
# 单个 gene 的 TSV 摘要
./gffsub query demo.gff3 --id gene01 --summary tsv

# 带选定属性的 JSON 摘要
./gffsub query demo.gff3 --id gene01 --summary json --out-attrs ID,Name,biotype

# 带属性的 TSV
./gffsub query demo.gff3 --id gene01 --summary tsv --out-attrs Name,biotype
```

## TSV 输出示例

```tsv
query_id	matched_id	matched_by	seqid	start	end	strand	type	parent_id	child_count	transcript_count	exon_count	cds_length	status
gene01	gene01	ID	chr1	100	1000	+	gene		1	1	2	402	found
```

## JSON 输出示例

```json
[
  {"query_id":"gene01","matched_id":"gene01","matched_by":"ID","seqid":"chr1","start":100,"end":1000,"strand":"+","type":"gene","parent_id":"","child_count":1,"transcript_count":1,"exon_count":2,"cds_length":402,"status":"found","attrs":{"Name":"BRCA1","biotype":"protein_coding"}}
]
```

## not_found 状态

查询 ID 没有命中任何 feature 时，摘要会输出一行 status 为 `not_found` 的记录：`matched_by` 保留 lookup key，数值计数为 0，其余字段为空。批量查询会报告缺失的 ID。

## 限制

--summary 和 --out-attrs 只能与 query 风格的 selector 配合使用（--id、--ids、--name、--where、--nearest、query 子命令 --region）。不能与 --grep、-I、-E、--seqid、--source、--score、--strand、--phase、--bed、--longest、--threads、--format/--output-format、--output 组合使用。
