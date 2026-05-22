# gffsub

<!-- README-I18N:START -->

[English](./README.md) | **汉语**

<!-- README-I18N:END -->

`gffsub` 是一个 C++17 编写的小型命令行工具，用于提取、查询、汇总和检查 GFF3/GTF 风格的基因组注释文件。它面向日常注释处理场景：按区域取子集、用多种 ID 找基因、提取基因上下游窗口、保留最长转录本，或在下游流程前快速做 QC。

## 从你的任务开始

| 我想要... | 使用 |
|----------|------|
| 提取某个基因组区间内的 feature | `gffsub annotation.gff3 -r chr1:1-100000 -f gene` |
| 用 BED 区间作为输入 | `gffsub annotation.gff3 -b regions.bed -f exon` |
| 按精确 ID 提取一个 feature | `gffsub annotation.gff3 --id GeneA` |
| 批量提取精确 ID | `gffsub annotation.gff3 --id-list genes.txt` |
| 按 name 找一个基因或 feature | `gffsub annotation.gff3 --name GeneA` |
| 按属性值提取 feature | `gffsub annotation.gff3 --attr biotype=protein_coding` |
| 批量提取 ID，并带上子 feature | `gffsub annotation.gff3 --id-list genes.txt --include-children` |
| 输出适合 pipeline 读取的 summary | `gffsub annotation.gff3 --id GeneA --summary-format tsv` |
| 提取指定属性值 | `gffsub annotation.gff3 --id GeneA --attrs ID,Name,Parent` |
| 提取基因上下游背景区域 | `gffsub annotation.gff3 --id GeneA --upstream 2000 --downstream 500 --strand-aware` |
| 每个基因只保留最长转录本 | `gffsub annotation.gff3 --longest` |
| 检查注释 graph 问题 | `gffsub annotation.gff3 --qc` |

## 简化命令速查

这些是常见工作流的最短写法：

```bash
# 按区间取子集
gffsub annotation.gff3 -r chr1:1-100000

# 按区间取子集，并限制 feature 类型
gffsub annotation.gff3 -r chr1:1-100000 -f gene

# 用 BED 区间取子集
gffsub annotation.gff3 -b regions.bed

# 转为 BED 输出
gffsub annotation.gff3 -r chr1:1-100000 -t bed

# 写入文件
gffsub annotation.gff3 -r chr1:1-100000 -o subset.gff3

# 每个基因保留最长转录本
gffsub annotation.gff3 -L

# 提取一个精确 ID
gffsub annotation.gff3 --id GeneA

# 查询一个名称或基因查找键
gffsub annotation.gff3 --name ABC1

# 提取精确属性匹配
gffsub annotation.gff3 --attr biotype=protein_coding

# 批量提取 ID
gffsub annotation.gff3 --id-list genes.txt

# 提取 ID 及其后代 feature
gffsub annotation.gff3 --id GeneA --include-children

# 提取一个区间
gffsub annotation.gff3 --region chr1:1-100000

# 查询并输出 summary
gffsub annotation.gff3 --id GeneA --summary-format tsv

# 查询并提取指定属性
gffsub annotation.gff3 --id GeneA --attrs ID,Name,Parent

# 提取上下游窗口
gffsub annotation.gff3 --id GeneA --upstream 2000 --downstream 500

# 输出 QC 表
gffsub annotation.gff3 --qc
```

## 安装

### 环境要求

- C++17 编译器（`g++` 或 `clang++`）
- CMake 可选；仓库也提供 `Makefile`

### 使用 Make 编译

```bash
cd gffsub
make
```

### 使用 CMake 编译并测试

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

## 场景：按区域提取 Feature

当你已经知道目标区间，或手里有 BED 文件时，使用经典模式。

```bash
# 提取 1-based 闭合 GFF 区间内的基因
./gffsub annotation.gff3 -r chr1:1-100000 -f gene

# 提取与 BED 区间重叠的 exon
./gffsub annotation.gff3 -b regions.bed -f exon

# 将重叠 feature 转为 BED 输出
./gffsub annotation.gff3 -r chr1:1-100000 -t bed
```

坐标规则是显式的：

| 输入或输出 | 坐标系统 |
|-----------|----------|
| GFF3/GTF records | 1-based 闭合区间 |
| `CHR:START-END` regions | 1-based 闭合区间 |
| BED input/output | 0-based 半开区间 |

## 场景：查找基因和 Gene Model

当问题从标识符、名称、属性或区间开始时，使用 `query`。

```bash
# 精确 feature ID
./gffsub annotation.gff3 --id Glyma.01G000100

# 按 ID、Name、gene_id、locus_tag、Alias 或完整 Dbxref 值查基因
./gffsub annotation.gff3 --name ABC1

# 属性选择器
./gffsub annotation.gff3 --attr biotype=protein_coding

# 包含 transcript、exon、CDS、UTR 等后代记录
./gffsub annotation.gff3 --id Glyma.01G000100 --include-children
```

在批处理流程中，可以输出 summary，而不是原始 GFF3：

```bash
./gffsub annotation.gff3 --id-list genes.txt --summary-format tsv
./gffsub annotation.gff3 --id GeneA --summary-format json
```

summary 字段包括 query ID、matched ID、匹配字段、坐标、链方向、feature 类型、parent ID、child 数量、transcript 数量、exon 数量、CDS 长度和状态。

## 场景：提取上下游窗口

当你需要某个基因或 feature 周围的局部注释背景时，使用 `window`。

```bash
# Genomic expansion：按参考序列左右扩展
./gffsub annotation.gff3 --id GeneA --upstream 2000 --downstream 500

# Biological upstream/downstream：按 strand 解释上下游
./gffsub annotation.gff3 --id GeneA --upstream 2000 --downstream 500 --strand-aware
```

不加 `--strand-aware` 时，upstream 表示更小的基因组坐标，downstream 表示更大的基因组坐标。加上 `--strand-aware` 后，上下游按 feature 的链方向解释。

## 场景：每个基因只保留一个转录本

当下游工具要求每个基因只有一个代表转录本时，使用 `--longest`。

```bash
./gffsub annotation.gff3 --longest
./gffsub annotation.gff3 --longest -@ 6
```

最长 isoform 逻辑遵循本项目现有的 AGAT 风格规则：如果基因存在 CDS isoform，则比较 CDS 长度；否则比较 exon 长度。

## 场景：检查注释质量

在把注释送入 graph-aware 工作流前，可以先运行 `--qc`。

```bash
./gffsub annotation.gff3 --qc
```

当前检查会报告重复 ID、非法坐标范围、Parent 缺失，以及 child feature 超出 parent 坐标范围。输出为 TSV，便于用常规命令行工具继续过滤。

## CLI 参数

`gffsub` 有四个命令入口：经典子集模式、`query`、`window` 和 `qc`。

### 经典子集模式

```bash
gffsub <input.gff3> [options]
```

| 参数 | 值 | 含义 |
|------|----|------|
| `<input.gff3>` | 文件 | 输入 GFF3/GTF 风格注释文件。 |
| `--id` | ID | 保留精确 feature `ID`。该参数可以重复使用。默认 GFF3 输出等价于 `gffsub query <input.gff3> --id ID`。 |
| `--id-list` | 文件 | 每个非空行读取一个精确 feature ID。默认 GFF3 输出等价于 `gffsub query <input.gff3> --id-list FILE`。 |
| `--name` | key | 保留一个按 `ID`、`Name`、`gene_id`、`locus_tag`、`Alias` 或完整 `Dbxref` 值找到的基因。默认 GFF3 输出等价于 `gffsub query <input.gff3> --name NAME`。 |
| `--attr` | `KEY=VALUE` | 保留精确 GFF3 属性值匹配的 feature。该参数可以重复使用。默认 GFF3 输出等价于 `gffsub query <input.gff3> --attr KEY=VALUE`。 |
| `--include-children` | 标志 | 包含由 `--id`、`--id-list`、`--name` 或 `--attr` 匹配记录的后代。 |
| `--attrs` | `KEY1,KEY2,...` | 将指定属性值作为额外 TSV/JSON 字段输出。该输出复用 `query` 的 summary 语义。 |
| `--summary-format` | `tsv`, `json` | 输出 summary 行，而不是 GFF3 记录。该输出复用 `query` 的 summary 语义。 |
| `--upstream` | 整数 | 与 `--id` 配合，提取与目标上游扩展窗口重叠的记录。输出等价于 `gffsub window <input.gff3> --id ID --upstream N`。 |
| `--downstream` | 整数 | 与 `--id` 配合，提取与目标下游扩展窗口重叠的记录。输出等价于 `window` 命令。 |
| `--strand-aware` | 标志 | 窗口提取时，按 feature 链方向解释 biological upstream/downstream。 |
| `--qc` | 标志 | 运行注释 QC。输出等价于 `gffsub qc <input.gff3>`。 |
| `-r`, `--region` | `CHR:START-END` | 保留与 1-based 闭合区间重叠的 feature。 |
| `-b`, `--bed` | 文件 | 保留与 BED 区间重叠的 feature；BED 按 0-based 半开区间读取。 |
| `-f`, `--feature`, `--type` | 类型 | 只保留第三列等于该类型的记录，例如 `gene`、`mRNA`、`transcript`、`exon` 或 `CDS`。 |
| `-L`, `--longest` | 标志 | 每个基因只保留最长转录本；存在 CDS isoform 时按 CDS 长度比较，否则按 exon 长度比较。 |
| `-@`, `--threads` | 整数 | 设置 `--longest` 使用的线程数；超过 256 会被限制为 256。需要固定资源使用时建议显式设置。 |
| `-t`, `--output-format` | `gff3`, `gtf`, `gtf2`, `gtf3`, `bed` | 选择输出格式。`gtf` 会按 `gtf2` 处理。默认输出 `gff3`。 |
| `-o`, `--output` | 文件 | 写入文件，而不是 stdout。 |
| `-h`, `--help` | 标志 | 显示经典子集模式帮助。 |

参数可以组合。例如 `-r chr1:1-100000 -f gene -t bed` 会先按区间过滤，再按 feature 类型过滤，最后输出 BED 坐标。

### Query 模式

```bash
gffsub query <input.gff3> [options]
```

| 参数 | 值 | 含义 |
|------|----|------|
| `<input.gff3>` | 文件 | 输入注释文件。 |
| `--id` | ID | 按精确 feature `ID` 查询。该参数可以重复使用。 |
| `--name` | key | 按 `ID`、`Name`、`gene_id`、`locus_tag`、`Alias` 或完整 `Dbxref` 值查找基因。 |
| `--id-list` | 文件 | 每个非空行读取一个精确 feature ID。 |
| `--region` | `CHR:START-END` | 查询与 1-based 闭合区间重叠的 feature。 |
| `--type` | 类型 | 将查询输出限制为一个 feature 类型。 |
| `--attr` | `KEY=VALUE` | 按精确属性值查询。该参数可以重复使用。 |
| `--attrs` | `KEY1,KEY2,...` | 将指定属性值追加为额外 TSV/JSON 字段。 |
| `--include-children` | 标志 | 包含匹配记录的后代，例如 transcript、exon、CDS 和 UTR feature。 |
| `--summary-format` | `tsv`, `json` | 输出 summary 行，而不是 GFF3 记录。 |
| `-h`, `--help` | 标志 | 显示 query 模式帮助。 |

使用 `--summary-format` 时，summary 字段为 `query_id`、`matched_id`、`matched_by`、`seqid`、`start`、`end`、`strand`、`type`、`parent_id`、`child_count`、`transcript_count`、`exon_count`、`cds_length` 和 `status`。如果使用 `--attrs`，这些属性键会追加为额外 TSV 列，或在 JSON 中输出为 `attrs` 对象。不加 `--summary-format` 时，`--attrs` 默认输出 TSV。

### Window 模式

```bash
gffsub window <input.gff3> --id ID [options]
```

| 参数 | 值 | 含义 |
|------|----|------|
| `<input.gff3>` | 文件 | 输入注释文件。 |
| `--id` | ID 或基因查找键 | 必填目标。命令会先查精确 `ID`，再做基因查找。 |
| `--upstream` | 整数 | 加到目标上游的碱基数；必须非负。默认 `0`。 |
| `--downstream` | 整数 | 加到目标下游的碱基数；必须非负。默认 `0`。 |
| `--strand-aware` | 标志 | 按 feature 链方向解释 biological upstream/downstream。不加时，upstream 表示更小的基因组坐标，downstream 表示更大的基因组坐标。 |
| `-h`, `--help` | 标志 | 显示 window 模式帮助。 |

输出为与扩展窗口重叠的 GFF3 记录。

### QC 模式

```bash
gffsub qc <input.gff3>
```

`qc` 只接受输入文件。它输出 TSV 表，字段为 `severity`、`code`、`line_idx`、`id` 和 `message`。当前检查代码包括 `duplicate_id`、`invalid_range`、`missing_parent` 和 `child_outside_parent`。

## 输出格式

| 格式 | Header | 坐标 |
|------|--------|------|
| gff3 | `##gff-version 3` | 1-based 闭合区间 |
| gtf2 | `##gtf-version 2` | 1-based 闭合区间 |
| gtf3 | `##gtf-version 2.2.1` | 1-based 闭合区间 |
| bed | 无 header | 0-based 半开区间 |

## C++ API

公开 C++ API 以 `AnnotationIndex` 为中心：

```cpp
auto index = gffsub::AnnotationIndex::from_gff3("annotation.gff3");
auto gene = index.find_gene("GeneA");
auto model = index.gene_model("GeneA");
auto hits = index.overlap("chr1", 1000, 2000);
```

可用查询包括 `find_by_id`、`find_gene`、`parents_of`、`children_of`、`descendants_of`、`gene_model`、`overlap`、`nearest_gene` 和 `with_attribute`。

## 分发

`gffsub` 会构建为单个命令行二进制文件。将二进制复制到 glibc-based Linux x86_64 机器后，即可配合注释文件运行。

## 许可证

MIT License
