# GTF 输入处理

<!-- I18N:START -->

[English](./gtf-input.md) | **中文**

<!-- I18N:END -->

## 格式识别

gffsub 嗅探文件内容：第 9 列含 `=` 和 `;` 判为 GFF3；第 9 列含 `"` 判为 GTF；3 到 12 列且 start/end 为整数判为 BED；默认 GFF3。

## GTF 与 GFF3 attribute

- GFF3：`key=value;` 对。
- GTF：`key "value";` 对。
- `gene_id` 和 `transcript_id` 提取到记录字段。

## Parent/child 合成

GTF 没有 `ID=`/`Parent=`。gffsub 合成如下：

| Feature 类型        | 合成 ID        | 合成 Parent       |
|---------------------|----------------|--------------------|
| gene                | gene_id        | 无                 |
| transcript / mRNA   | transcript_id  | gene_id            |
| exon、CDS 等         | 无             | transcript_id      |

由此 GTF 可用 `--children`、`--parents`、`--model`。

## attribute 访问

对 `--where`、`--grep`、`-I`：`gene_id` 和 `transcript_id` 从记录字段解析。其他 GTF attribute 用 `attr.` 前缀（如 `attr.gene_name`）。

```bash
./gffsub input.gtf --where gene_id=ENSG000001
./gffsub input.gtf --grep gene_id:ENSG
./gffsub input.gtf -I 'gene_id == "ENSG000001"'
```

## GTF 上的 --longest

自动检测 isoform 类型：不存在 `mRNA` 记录时使用 `transcript`。无需 `-f transcript`。

## GTF 转 GFF3

第 9 列重写为 `tag=value`。`ID=`/`Parent=` 由 `gene_id`/`transcript_id` 合成。其他 attribute 做 URL 转义。

```bash
./gffsub input.gtf -t gff3
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
./gffsub demo.gtf --where gene_id=g1
./gffsub demo.gtf --id t1 -C
./gffsub demo.gtf --longest
./gffsub demo.gtf --id g1 -C -t gff3
```
