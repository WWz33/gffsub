# gffsub

<!-- README-I18N:START -->

[English](./README.md) | **汉语**

<!-- README-I18N:END -->

`gffsub` 是一个 C++17 编写的小型命令行工具，用于提取、查询、汇总和检查 GFF3/GTF 风格的基因组注释文件。它面向日常注释处理场景：按区域取子集、用多种 ID 找基因、提取基因上下游窗口、保留最长转录本，或在下游流程前快速做 QC。

## 从你的任务开始

| 我想要... | 使用 |
|----------|------|
| 提取某个基因组区间内的 feature | `gffsub annotation.gff3 -r chr1:1-100000 -f gene` |
| 提取某个染色体或 contig | `gffsub annotation.gff3 --seqid chr1` |
| 提取某个注释来源的记录 | `gffsub annotation.gff3 --source Gnomon` |
| 提取某个 score 的记录 | `gffsub annotation.gff3 --score 42.5` |
| 提取某条链上的记录 | `gffsub annotation.gff3 --strand -` |
| 提取某个 phase 的记录 | `gffsub annotation.gff3 --phase 0 -f CDS` |
| 用 BED 区间作为输入 | `gffsub annotation.gff3 -b regions.bed -f exon` |
| 按精确 ID 提取一个 feature | `gffsub annotation.gff3 --id GeneA` |
| 批量提取精确 ID | `gffsub annotation.gff3 --ids genes.txt` |
| 按 name 找一个基因或 feature | `gffsub annotation.gff3 --name GeneA` |
| 按属性值提取 feature | `gffsub annotation.gff3 --where biotype=protein_coding` |
| 查找离某个区间最近的基因 | `gffsub annotation.gff3 --nearest chr1:1000-2000` |
| 批量提取 ID，并带上子 feature | `gffsub annotation.gff3 --ids genes.txt -C` |
| 提取一个 feature 及其父级 | `gffsub annotation.gff3 --id ExonA --parents` |
| 从任意 feature 提取完整 gene model | `gffsub annotation.gff3 --id ExonA --model` |
| 输出适合 pipeline 读取的 summary | `gffsub annotation.gff3 --id GeneA --summary tsv` |
| 提取指定属性值 | `gffsub annotation.gff3 --id GeneA --out-attrs ID,Name,Parent` |
| 提取基因上下游背景区域 | `gffsub annotation.gff3 --id GeneA --up 2000 --down 500 --strand-aware` |
| 每个基因只保留最长转录本 | `gffsub annotation.gff3 --longest` |
| 检查注释 graph 问题 | `gffsub annotation.gff3 --qc` |

## 简化命令速查

这些是常见工作流的最短写法：

```bash
# 按区间取子集
gffsub annotation.gff3 -r chr1:1-100000

# 按第 1 列 seqid 取子集
gffsub annotation.gff3 --seqid chr1

# 按第 2 列 source 取子集
gffsub annotation.gff3 --source Gnomon

# 按第 6 列 score 取子集
gffsub annotation.gff3 --score 42.5

# 按第 7 列 strand 取子集
gffsub annotation.gff3 --strand -

# 按第 8 列 phase 取子集
gffsub annotation.gff3 --phase 0

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
gffsub annotation.gff3 --where biotype=protein_coding

# 查找离 1-based 闭合区间最近的基因
gffsub annotation.gff3 --nearest chr1:1000-2000

# 批量提取 ID
gffsub annotation.gff3 --ids genes.txt

# 提取 ID 及其后代 feature
gffsub annotation.gff3 --id GeneA -C

# 提取 ID 及其祖先 feature
gffsub annotation.gff3 --id ExonA --parents

# 提取包含该 ID 的完整 gene model
gffsub annotation.gff3 --id ExonA --model

# 提取一个区间
gffsub annotation.gff3 --region chr1:1-100000

# 查询并输出 summary
gffsub annotation.gff3 --id GeneA --summary tsv

# 查询并提取指定属性
gffsub annotation.gff3 --id GeneA --out-attrs ID,Name,Parent

# 提取上下游窗口
gffsub annotation.gff3 --id GeneA --up 2000 --down 500

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

# 提取某个第 1 列 seqid 下的全部记录
./gffsub annotation.gff3 --seqid chr1

# 提取某个第 2 列 source 下的全部记录
./gffsub annotation.gff3 --source Gnomon

# 提取某个第 6 列 score 下的全部记录
./gffsub annotation.gff3 --score 42.5

# 提取某个第 7 列 strand 下的全部记录
./gffsub annotation.gff3 --strand -

# 提取某个第 8 列 phase 下的全部记录
./gffsub annotation.gff3 --phase 0

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

当问题从标识符、名称、属性或区间开始时，优先使用顶层 selector 写法。

```bash
# 精确 feature ID
./gffsub annotation.gff3 --id Glyma.01G000100

# 按 ID、Name、gene_id、locus_tag、Alias 或完整 Dbxref 值查基因
./gffsub annotation.gff3 --name ABC1

# 属性选择器
./gffsub annotation.gff3 --where biotype=protein_coding

# 查找离 1-based 闭合区间最近的基因
./gffsub annotation.gff3 --nearest chr1:1000-2000

# 包含 transcript、exon、CDS、UTR 等后代记录
./gffsub annotation.gff3 --id Glyma.01G000100 -C

# 包含通过 Parent 连接的祖先记录
./gffsub annotation.gff3 --id ExonA --parents

# 返回包含该 feature 的完整 gene model
./gffsub annotation.gff3 --id ExonA --model
```

在批处理流程中，可以输出 summary，而不是原始 GFF3：

```bash
./gffsub annotation.gff3 --ids genes.txt --summary tsv
./gffsub annotation.gff3 --id GeneA --summary json
```

summary 字段包括 query ID、matched ID、匹配字段、坐标、链方向、feature 类型、parent ID、child 数量、transcript 数量、exon 数量、CDS 长度和状态。

非 region 查询常用这些第 9 列属性键：

| 任务 | 命令 | 使用的键 |
|------|------|----------|
| 精确 feature 查询 | `--id gene0001` | `ID` |
| 批量精确 feature 查询 | `--ids genes.txt` | 每行一个 `ID` 值 |
| 基因名称查询 | `--name ABC1` | gene 记录上的 `ID`、`gene_id`、`Name`、`locus_tag`、`Alias` 或完整 `Dbxref` 值 |
| 任意精确属性过滤 | `--where Parent=gene0001` | 任意第 9 列 `KEY=VALUE`，包括 `ID`、`Name`、`Alias`、`Parent`、`Dbxref`、`Accession` 或 `Parent_Accession` |
| 最近基因查询 | `--nearest chr1:1000-2000` | 同一 seqid 上距离 1-based 闭合区间最近的 gene；与区间重叠的 gene 距离为 0；并列时按输入文件顺序 |
| 包含匹配记录后代 | `-C`, `--children` | 通过 `Parent` 连接的 child 记录；从 `--id`、`--ids`、`--name`、`--where` 或 `--nearest` 匹配记录出发；`--include-children` 是较长兼容别名 |
| 包含匹配记录祖先 | `--parents` | 沿 `Parent` 链向上找到的 parent 记录；从 `--id`、`--ids`、`--name`、`--where` 或 `--nearest` 匹配记录出发；`--include-parents` 是较长兼容别名 |
| 提取完整 gene model | `--model`, `--gene-model` | 从 `--id`、`--ids`、`--name`、`--where` 或 `--nearest` 匹配的记录出发，返回所属 gene 及其后代记录 |
| 打印指定属性 | `--out-attrs ID,Name,Parent` | 记录匹配后，将指定第 9 列键输出为 summary 字段 |

## 场景：提取上下游窗口

当你需要某个基因或 feature 周围的局部注释背景时，优先使用顶层 window 参数。

```bash
# Genomic expansion：按参考序列左右扩展
./gffsub annotation.gff3 --id GeneA --up 2000 --down 500

# Biological upstream/downstream：按 strand 解释上下游
./gffsub annotation.gff3 --id GeneA --up 2000 --down 500 --strand-aware
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

当前检查会报告重复 ID、非法坐标范围、非法 CDS phase、Parent 缺失，以及 child feature 超出 parent 坐标范围。输出为 TSV，便于用常规命令行工具继续过滤。

## CLI 参数

`gffsub` 以顶层入口为主：常用 GFF3 工作流从 `gffsub <input.gff3> [options]` 开始。`query`、`window` 和 `qc` 子命令作为兼容和进阶入口保留，并共享同一套输出语义。

### 顶层模式

```bash
gffsub <input.gff3> [options]
```

| 参数 | 值 | 含义 |
|------|----|------|
| `<input.gff3>` | 文件 | 输入 GFF3/GTF 风格注释文件。 |
| `--id` | ID | 保留精确 feature `ID`。该参数可以重复使用。这是精确 ID 提取的顶层 selector。 |
| `--ids`, `--id-list` | 文件 | 每个非空行读取一个精确 feature ID。这是批量精确 ID 提取的顶层 selector。`--id-list` 是较长兼容别名。 |
| `--name` | key | 保留一个按 `ID`、`Name`、`gene_id`、`locus_tag`、`Alias` 或完整 `Dbxref` 值找到的基因。这是常见基因命名键的顶层 selector。 |
| `--where`, `--attr` | `KEY=VALUE` | 保留精确 GFF3 属性值匹配的 feature。该参数可以重复使用。这是精确第 9 列 `KEY=VALUE` 匹配的顶层 selector。`--attr` 是兼容别名。 |
| `-C`, `--children`, `--include-children` | 标志 | 包含由 `--id`、`--ids`、`--name`、`--where` 或 `--nearest` 匹配记录的后代。`--include-children` 是较长兼容别名。 |
| `--parents`, `--include-parents` | 标志 | 包含由 `--id`、`--ids`、`--name`、`--where` 或 `--nearest` 匹配记录的祖先。`--include-parents` 是较长兼容别名。 |
| `--model`, `--gene-model` | 标志 | 包含由 `--id`、`--ids`、`--name`、`--where` 或 `--nearest` 匹配记录所属的完整 gene model；输出 gene 以及 transcript/exon/CDS/UTR 等后代。`--gene-model` 是较长兼容别名。 |
| `--nearest`, `--nearest-gene` | `CHR:START-END` | 保留同一 seqid 上离 1-based 闭合区间最近的 gene；并列时按输入文件顺序。`--nearest-gene` 是较长兼容别名。 |
| `--out-attrs`, `--output-attrs` | `KEY1,KEY2,...` | 将指定属性值作为额外 TSV/JSON 字段输出；只与 query-style selector 组合。 |
| `--attrs` | `KEY1,KEY2,...` | `--out-attrs` 的兼容别名，已不推荐使用。 |
| `--summary`, `--summary-format` | `tsv`, `json` | 输出 summary 行，而不是 GFF3 记录；只与 query-style selector 组合。`--summary-format` 是较长兼容别名。 |
| `--up`, `--upstream` | 整数 | 与 `--id` 配合，提取与目标上游扩展窗口重叠的记录。 |
| `--down`, `--downstream` | 整数 | 与 `--id` 配合，提取与目标下游扩展窗口重叠的记录。 |
| `--strand-aware` | 标志 | 窗口提取时，按 feature 链方向解释 biological upstream/downstream。 |
| `--qc` | 标志 | 运行注释 QC。 |
| `--seqid` | seqid | 只保留第 1 列 GFF3 seqid 精确等于该值的记录。 |
| `--source` | source | 只保留第 2 列 GFF3 source 精确等于该值的记录。 |
| `--score` | 数字, `.` | 只保留第 6 列 GFF3 score 匹配该数值，或第 6 列为 `.` 的记录。 |
| `--strand` | `+`, `-`, `.`, `?` | 只保留第 7 列 GFF3 strand 精确等于该值的记录。这是过滤器；`--strand-aware` 只改变 window 上下游解释。 |
| `--phase` | `0`, `1`, `2`, `.` | 只保留第 8 列 GFF3 phase 精确等于该值的记录。对 CDS 记录，GFF3 phase 通常是 `0`、`1` 或 `2`；`.` 会匹配第 8 列为 `.` 的任何记录。 |
| `-r`, `--region` | `CHR:START-END` | 保留与 1-based 闭合区间重叠的 feature。 |
| `-b`, `--bed` | 文件 | 保留与 BED 区间重叠的 feature；BED 按 0-based 半开区间读取。 |
| `-f`, `--feature`, `--type` | 类型 | 只保留第三列等于该类型的记录，例如 `gene`、`mRNA`、`transcript`、`exon` 或 `CDS`。 |
| `-L`, `--longest` | 标志 | 每个基因只保留最长转录本；存在 CDS isoform 时按 CDS 长度比较，否则按 exon 长度比较。 |
| `-@`, `--threads` | 整数 | 设置 `--longest` 使用的线程数；超过 256 会被限制为 256。需要固定资源使用时建议显式设置。 |
| `-t`, `--format`, `--output-format` | `gff3`, `gtf`, `gtf2`, `gtf3`, `bed` | 选择输出格式。`gtf` 会按 `gtf2` 处理。默认输出 `gff3`。`--output-format` 是较长兼容别名。 |
| `-o`, `--output` | 文件 | 写入文件，而不是 stdout。 |
| `-h`, `--help` | 标志 | 显示顶层模式帮助。 |

参数可以组合。例如 `-r chr1:1-100000 -f gene -t bed` 会先按区间过滤，再按 feature 类型过滤，最后输出 BED 坐标。

### Query 兼容模式

```bash
gffsub query <input.gff3> [options]
```

大多数默认 GFF3 和 summary 工作流都可以不写 `query` 子命令。该兼容模式保留给已有脚本和显式 query 风格命令行。

| 参数 | 值 | 含义 |
|------|----|------|
| `<input.gff3>` | 文件 | 输入注释文件。 |
| `--id` | ID | 按精确 feature `ID` 查询。该参数可以重复使用。 |
| `--name` | key | 按 `ID`、`Name`、`gene_id`、`locus_tag`、`Alias` 或完整 `Dbxref` 值查找基因。 |
| `--ids`, `--id-list` | 文件 | 每个非空行读取一个精确 feature ID。`--id-list` 是较长兼容别名。 |
| `--region` | `CHR:START-END` | 查询与 1-based 闭合区间重叠的 feature。 |
| `--nearest`, `--nearest-gene` | `CHR:START-END` | 查询同一 seqid 上离 1-based 闭合区间最近的 gene；并列时按输入文件顺序。`--nearest-gene` 是较长兼容别名。 |
| `--type` | 类型 | 将查询输出限制为一个 feature 类型。 |
| `--where`, `--attr` | `KEY=VALUE` | 按精确属性值查询。该参数可以重复使用。`--attr` 是兼容别名。 |
| `--out-attrs`, `--output-attrs` | `KEY1,KEY2,...` | 将指定属性值追加为额外 TSV/JSON 字段。 |
| `--attrs` | `KEY1,KEY2,...` | `--out-attrs` 的兼容别名，已不推荐使用。 |
| `-C`, `--children`, `--include-children` | 标志 | 包含匹配记录的后代，例如 transcript、exon、CDS 和 UTR feature。`--include-children` 是较长兼容别名。 |
| `--parents`, `--include-parents` | 标志 | 沿 GFF3 `Parent` 链向上包含匹配记录的祖先。`--include-parents` 是较长兼容别名。 |
| `--model`, `--gene-model` | 标志 | 包含匹配记录所属的完整 gene model。`--gene-model` 是较长兼容别名。 |
| `--summary`, `--summary-format` | `tsv`, `json` | 输出 summary 行，而不是 GFF3 记录。`--summary-format` 是较长兼容别名。 |
| `-h`, `--help` | 标志 | 显示 query 模式帮助。 |

使用 `--summary` 时，summary 字段为 `query_id`、`matched_id`、`matched_by`、`seqid`、`start`、`end`、`strand`、`type`、`parent_id`、`child_count`、`transcript_count`、`exon_count`、`cds_length` 和 `status`。如果使用 `--out-attrs`，这些属性键会追加为额外 TSV 列，或在 JSON 中输出为 `attrs` 对象。不加 `--summary` 时，`--out-attrs` 默认输出 TSV。

### Window 兼容模式

```bash
gffsub window <input.gff3> --id ID [options]
```

同一工作流也可以写成顶层形式：`gffsub <input.gff3> --id ID --up N --down N`。

| 参数 | 值 | 含义 |
|------|----|------|
| `<input.gff3>` | 文件 | 输入注释文件。 |
| `--id` | ID 或基因查找键 | 必填目标。命令会先查精确 `ID`，再做基因查找。 |
| `--up`, `--upstream` | 整数 | 加到目标上游的碱基数；必须非负。默认 `0`。 |
| `--down`, `--downstream` | 整数 | 加到目标下游的碱基数；必须非负。默认 `0`。 |
| `--strand-aware` | 标志 | 按 feature 链方向解释 biological upstream/downstream。不加时，upstream 表示更小的基因组坐标，downstream 表示更大的基因组坐标。 |
| `-h`, `--help` | 标志 | 显示 window 模式帮助。 |

输出为与扩展窗口重叠的 GFF3 记录。

### QC 兼容模式

```bash
gffsub qc <input.gff3>
```

同一工作流也可以写成顶层形式：`gffsub <input.gff3> --qc`。QC 输出 TSV 表，字段为 `severity`、`code`、`line_idx`、`id` 和 `message`。当前检查代码包括 `duplicate_id`、`invalid_range`、`invalid_strand`、`invalid_cds_phase`、`missing_parent` 和 `child_outside_parent`。

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
