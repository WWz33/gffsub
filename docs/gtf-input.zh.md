# GTF 输入处理

<!-- I18N:START -->

[English](./gtf-input.md) | **中文**

<!-- I18N:END -->

## 格式识别

gffsub 看第一条数据行：第 9 列第一个片段决定格式。`tag=value` 是 GFF3，`tag "value"` 是 GTF，3 到 12 列且 start/end 为整数是 BED，默认 GFF3。

只看第一个片段，所以值里的标点不会翻转判定：`ID=g1;Name="A"` 仍是 GFF3，`note "a=b"` 仍是 GTF。

## GTF 与 GFF3 attribute

- GFF3：`key=value;` 对。
- GTF：`key "value";` 对。
- `gene_id` 和 `transcript_id` 提取到记录字段。

GTF 第 9 列细则：

- 引号外的 `#` 起注释，到行尾结束；引号内的 `#` 是数据。
- `\"`、`\\`、`\t`、`\n`、`\r` 会被解码。
- 未加引号的值也接受。
- GFF3 风格的裸 `ID=x;Parent=y` 片段也接受；`key="value"` 的 tag 名取 `key`。

## Parent/child 合成

GTF 没有 `ID=`/`Parent=`。gffsub 合成如下：

| Feature 类型        | 合成 ID        | 合成 Parent       |
|---------------------|----------------|--------------------|
| gene                | gene_id        | 无                 |
| transcript / mRNA   | transcript_id  | gene_id            |
| exon、CDS 等         | 无             | transcript_id      |

由此 GTF 可用 `-C`、`-p`、`-m`。

## attribute 访问

对 `--where`、`--grep`、`-I`：`gene_id` 和 `transcript_id` 从记录字段解析。其他 GTF attribute 用 `attr.` 前缀（如 `attr.gene_name`）。

```bash
./gffsub input.gtf -w gene_id=ENSG000001
./gffsub input.gtf --grep gene_id:ENSG
./gffsub input.gtf -I 'gene_id == "ENSG000001"'
```

## GTF 上的 --longest

自动检测 isoform 类型：取最多的 transcript 类类型（GTF 为 `transcript`）。无需 `-t transcript`。

## GTF 转 GFF3

第 9 列重写为 `tag=value`。`ID=`/`Parent=` 由 `gene_id`/`transcript_id` 合成。其他 attribute 做 URL 转义。

```bash
./gffsub input.gtf --format gff3
```

## GTF 转 GTF

每行带 `gene_id`。`transcript_id` 仅出现在非 gene 行。值做转义。

输入样例（demo.gtf）：
```
chr1	src	gene	100	1000	.	+	.	gene_id "g1"; gene_name "BRCA1";
chr1	src	transcript	100	1000	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	src	exon	100	250	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	src	exon	500	750	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	src	CDS	100	250	.	+	0	gene_id "g1"; transcript_id "t1";
```

```bash
./gffsub demo.gtf -w gene_id=g1
./gffsub demo.gtf -i t1 -C
./gffsub demo.gtf --longest
./gffsub demo.gtf -i g1 -C --format gff3
```
