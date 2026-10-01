# 入门

<!-- I18N:START -->

[English](./getting-started.md) | **中文**

<!-- I18N:END -->

gffsub 按 region、feature ID、attribute 或 gene model 从 GFF3、GTF 注释文件中提取子集。

## 编译

```bash
git clone https://github.com/WWz33/gffsub.git
cd gffsub && make -j
```

生成 `./gffsub`。需 C++17 编译器（g++ 9+、clang 10+）与 zlib（`-lz`；macOS 自带，Debian/Ubuntu 安装 `zlib1g-dev`）。

## 输入

- GFF3、GTF 或 BED，按内容识别（不依赖扩展名）
- 纯文本或 gzip 压缩文件均可直接读取（`ann.gtf.gz`）
- `-` 读取 stdin，纯文本或 gzip 均可

## 示例数据

存为 `demo.gff3`：

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
```

## 示例

按 region 提取：

```bash
./gffsub demo.gff3 -r chr1:200-600
```

按 ID 选择并带子 feature：

```bash
./gffsub demo.gff3 -i tx01 -C
```

按 feature 类型过滤：

```bash
./gffsub demo.gff3 -t exon
```

转 GTF：

```bash
./gffsub demo.gff3 --format gtf
```

每个 gene 取最长 isoform：

```bash
./gffsub demo.gff3 --longest
```

## 长命令

一行一个阶段。selector 挑行，`--longest`/`--drop-orphans` 整理层级，`--format` 在输出时转换格式。

蛋白编码品系，每 gene 一个 isoform，完整 model，输出 GTF：

```bash
./gffsub ann.gff3 -w biotype=protein_coding -L -m -f gtf -o coding.gtf
```

`-w` 保留 gene 行；`-m` 把选中的 gene 重新展开成完整 model，transcript 得以在
`--longest` 中参与竞争。

区域去掉重复区，层级清理：

```bash
./gffsub ann.gff3 -S chr1 --exclude-region chr1:500000-800000 --drop-orphans
```

按列表取 gene，只留 model，精简第 9 列：

```bash
./gffsub ann.gff3 --ids genes.txt -m --out-attrs Name,biotype
```

离某个位点最近的 gene 及其 model：

```bash
./gffsub ann.gff3 -N chr1:1000000-1000500 -m
```

排序、gzip 输入、BED 输出：

```bash
./gffsub ann.gff3.gz -k seqid,natural-seqid,start -t gene -f bed
```

各阶段也可通过 stdin 串接：

```bash
./gffsub ann.gff3 -w biotype=protein_coding -C | ./gffsub - -L -t mRNA
```

## 命令结构

```
gffsub <input> [options]
gffsub query <input> [options]
gffsub window <input> [options]
```

## 帮助

```bash
./gffsub -h
./gffsub query -h
./gffsub window -h
```
