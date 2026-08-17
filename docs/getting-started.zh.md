# 快速上手

<!-- I18N:START -->

[English](./getting-started.md) | **中文**

<!-- I18N:END -->

gffsub 按 region、feature ID、属性或 gene model 从 GFF3 和 GTF 注释文件中提取子集。

## 安装

```bash
git clone https://github.com/WWz33/gffsub.git
cd gffsub && make -j
```

生成 `./gffsub`。需要 C++17 编译器（g++ 9+、clang 10+），无外部依赖。

## 快速上手

示例数据（保存为 `demo.gff3`）：

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
```

### 按 region 提取

保留与基因组区间重叠的记录。

```bash
./gffsub demo.gff3 -r chr1:200-600
```

### 按 ID 选取并包含子 feature

保留一个 feature 及其全部后代。

```bash
./gffsub demo.gff3 --id tx01 -C
```

### 按 feature 类型筛选

```bash
./gffsub demo.gff3 -f exon
```

### 转换为 GTF

```bash
./gffsub demo.gff3 --id gene01 -C -t gtf
```

### 最长 isoform

每个 gene 保留最长的 mRNA。

```bash
./gffsub demo.gff3 --longest
```

## 命令结构

```
gffsub <input.gff3> [options]
gffsub query <input.gff3> [options]
gffsub window <input.gff3> [options]
```

- 默认模式：提取并筛选记录，输出注释。
- `query`：按 ID、name 或属性查找 feature，输出汇总。
- `window`：某个 feature 周围窗口内重叠的记录。

## 帮助

```bash
./gffsub -h          # 完整 flag 列表
./gffsub query -h    # query 专用选项
./gffsub window -h   # window 专用选项
```
