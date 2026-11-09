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

- GFF3 或 GTF，按内容识别（不依赖扩展名）
- BED 文件只用于 BED 输出（`-f bed`）和统计输出（`-s`）；BED 不转换为 GFF3 或 GTF
- 纯文本或 gzip 压缩文件均可直接读取（`ann.gtf.gz`）
- `-` 读取 stdin，纯文本或 gzip 均可

## 解析告警

格式有问题的输入不会中断运行。每个文件结束后在 stderr 报一次：

- 因列数不足被跳过的行
- 因坐标非法被跳过的行
- phase 不是 0/1/2 的 CDS 行
- 构成环的 Parent 链

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

```bash
# 与 chr1:200-600 重叠的记录
./gffsub demo.gff3 -r chr1:200-600

# tx01 及其 exon、CDS
./gffsub demo.gff3 -i tx01 -C

# 仅 exon 行
./gffsub demo.gff3 -t exon

# 输出 GTF3
./gffsub demo.gff3 --format gtf

# 每个 gene 只留一个转录本
./gffsub demo.gff3 --longest
```

多个过滤条件写在同一次调用里。对完整注释文件：

```bash
# 蛋白编码基因，最长 isoform，完整 gene 模型，输出 GTF
./gffsub ann.gff3 -w biotype=protein_coding -L -m -f gtf -o coding.gtf

# chr1 去掉 chr1:500000-800000，输出不含悬空 Parent 引用
./gffsub ann.gff3 -S chr1 --exclude-region chr1:500000-800000 --drop-orphans

# ID 来自文件，取完整 gene 模型，第 9 列只留 Name 和 biotype
./gffsub ann.gff3 --ids genes.txt -m --out-attrs Name,biotype

# 距 500 bp 区间最近的 gene，取其完整模型
./gffsub ann.gff3 -N chr1:1000000-1000500 -m

# gzip 输入，按 seqid 再按起点排序，gene 行输出 BED
./gffsub ann.gff3.gz -k seqid,natural-seqid,start -t gene -f bed

# 第二遍从 stdin 读第一遍的结果
./gffsub ann.gff3 -w biotype=protein_coding -C | ./gffsub - -L -t mRNA
```

## 选择条件组合

默认命令中，选择器（`-i`、`--ids`、`-n`、`-w`、`-N`）给出候选记录，过滤器（`-r`、`-b`、`-t` 等）再收窄，两边必须同时满足。

```bash
# ex01 位于 chr1:100-250，与 chr1:600-650 不重叠
./gffsub demo.gff3 -i ex01 -r chr1:600-650
```

```
##gff-version 3
```

`query` 子命令把每个选择器命中的记录相加。同样的参数返回按 ID 命中的 ex01，以及与该区间重叠的记录：

```bash
./gffsub query demo.gff3 -i ex01 -r chr1:600-650
```

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
```

## 输出

- `-o FILE` 写入文件；不写 `-o` 或写 `-o -` 时写 stdout。
- gzip 输入透明读取（`.gz` 路径或 gzip 压缩的 stdin）。gzip 输出未实现：`-o out.gz` 会报错，请写纯文本后另行压缩。
- `query` 与 `window` 同样接受 `-o`、`-f`、`--out-attrs`。

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
