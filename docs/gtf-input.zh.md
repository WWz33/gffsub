# GTF 输入处理

## 格式检测

gffsub 根据文件扩展名推断输入格式：
- `.gtf` 或 `.GTF` -> GTF 解析
- `.bed` 或 `.BED` -> BED 解析
- 其余 -> GFF3 解析

## GTF 属性模型

GTF 的属性格式与 GFF3 不同：
- GFF3：`key=value;`（key 与 value 以 `=` 分隔，键值对以 `;` 分隔）
- GTF：`key "value";`（key 与带引号的 value 以空格分隔，键值对以 `;` 分隔）

gffsub 的 `parse_attributes` 函数按 GFF3 风格解析属性。对于 GTF，解析器使用一个感知 GTF 的提取器，把 `gene_id` 和 `transcript_id` 抽到专用的 record 字段中，该提取器能处理带引号的值和转义引号。

## Parent/child 合成

GTF 没有显式的 `ID=` 和 `Parent=` 属性。gffsub 从 `gene_id` 和 `transcript_id` 合成这些属性：

| Feature 类型 | 合成 ID | 合成 Parent |
|--------------|---------------|-------------------|
| gene | gene_id | 无（无 parent） |
| transcript / mRNA | transcript_id | gene_id |
| exon、CDS 等 | 无 | transcript_id |

合成后会构建 parent/child 索引，因此 `--children`、`--parents` 和 `--model` 可用于 GTF 输入。exon 和 CDS 记录不会合成 ID，只通过 parent 字段链接到其 transcript。这样可避免多个 exon 共享同一 transcript_id 时出现 ID 冲突。

## GTF 专属解析

- 带引号的值：`gene_id "ENSG000001";` 提取引号之间的字符串。
- 转义引号：`gene_id "test\"name";` 处理值中反斜杠转义的引号。
- 空值：`gene_id "";` 视为缺失（该 record 没有 gene_id）。

## GTF 的属性访问

对于 `--where`、`--grep` 和 `-I` 表达式：
- `gene_id` 和 `transcript_id` 从合成的 record 字段解析，不是从第 9 列解析。
- 其他 GTF 属性（如 `gene_name`、`gene_biotype`）可通过 `attr.` 前缀访问，但走 GTF 解析路径，直接提取带引号的值。

在 GTF 中按 gene_id 筛选：
```bash
./gffsub input.gtf --where gene_id=ENSG000001
./gffsub input.gtf --grep gene_id:ENSG
./gffsub input.gtf -I 'gene_id == "ENSG000001"'
```

## GTF 上的 --longest

isoform 类型自动检测。GTF 通常用 `transcript` 而非 `mRNA`，gffsub 会检查数据中存在哪种类型：
- 若存在 `mRNA` 记录，用 `mRNA`。
- 否则若存在 `transcript` 记录，用 `transcript`。

`--longest` 在 GTF 上无需 `-f transcript` 即可工作。

## GTF 输出

GTF 转 GTF 时，gffsub 读取 gene_id 和 transcript_id，重建层级，输出符合 GTF2.2 的属性：
- 每行都有 gene_id（GTF2.2 要求）。
- 非 gene 行才有 transcript_id。
- gene 行不写 transcript_id（符合 GTF2.2 和 GENCODE 约定）。
- 属性值做转义处理（反斜杠和引号）。

详见 [输出格式](output-formats.md)。

## 示例

示例 GTF（demo.gtf）：
```
chr1	src	gene	100	1000	.	+	.	gene_id "g1"; gene_name "BRCA1";
chr1	src	transcript	100	1000	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	src	exon	100	250	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	src	exon	500	750	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	src	CDS	100	250	.	+	0	gene_id "g1"; transcript_id "t1";
```

```bash
# 按 gene_id 选择
./gffsub demo.gtf --where gene_id=g1

# transcript 的子节点
./gffsub demo.gtf --id t1 -C

# 最长 isoform（自动检测 transcript 类型）
./gffsub demo.gtf --longest

# GTF 转 GFF3：第 9 列改写为 tag=value，并合成
# ID=/Parent=（gene 行 ID=<gene_id>，transcript 行 ID=<transcript_id>
# +Parent=<gene_id>，子节点 Parent=<transcript_id>）；其余属性
# 转为 key=value 形式并做 URL 转义
./gffsub demo.gtf --id g1 -C -t gff3
```
