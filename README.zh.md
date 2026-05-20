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
| 按 ID/name/attribute 找一个基因或 feature | `gffsub query annotation.gff3 --id GeneA` |
| 批量查询 ID，并带上子 feature | `gffsub query annotation.gff3 --id-list genes.txt --include-children` |
| 输出适合 pipeline 读取的 summary | `gffsub query annotation.gff3 --id GeneA --summary-format tsv` |
| 提取基因上下游背景区域 | `gffsub window annotation.gff3 --id GeneA --upstream 2000 --downstream 500 --strand-aware` |
| 每个基因只保留最长转录本 | `gffsub annotation.gff3 --longest` |
| 检查注释 graph 问题 | `gffsub qc annotation.gff3` |

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
./gffsub query annotation.gff3 --id Glyma.01G000100

# 按 ID、Name、gene_id、locus_tag、Alias 或完整 Dbxref 值查基因
./gffsub query annotation.gff3 --name ABC1

# 属性选择器
./gffsub query annotation.gff3 --attr biotype=protein_coding

# 包含 transcript、exon、CDS、UTR 等后代记录
./gffsub query annotation.gff3 --id Glyma.01G000100 --include-children
```

在批处理流程中，可以输出 summary，而不是原始 GFF3：

```bash
./gffsub query annotation.gff3 --id-list genes.txt --summary-format tsv
./gffsub query annotation.gff3 --id GeneA --summary-format json
```

summary 字段包括 query ID、matched ID、匹配字段、坐标、链方向、feature 类型、parent ID、child 数量、transcript 数量、exon 数量、CDS 长度和状态。

## 场景：提取上下游窗口

当你需要某个基因或 feature 周围的局部注释背景时，使用 `window`。

```bash
# Genomic expansion：按参考序列左右扩展
./gffsub window annotation.gff3 --id GeneA --upstream 2000 --downstream 500

# Biological upstream/downstream：按 strand 解释上下游
./gffsub window annotation.gff3 --id GeneA --upstream 2000 --downstream 500 --strand-aware
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

在把注释送入 graph-aware 工作流前，可以先运行 `qc`。

```bash
./gffsub qc annotation.gff3
```

当前检查会报告重复 ID、非法坐标范围、Parent 缺失，以及 child feature 超出 parent 坐标范围。输出为 TSV，便于用常规命令行工具继续过滤。

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
