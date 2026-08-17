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

生成 `./gffsub`。需 C++17 编译器（g++ 9+、clang 10+），无外部依赖。

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
