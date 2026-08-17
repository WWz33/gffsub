# 最长 Isoform 选择

<!-- I18N:START -->

[English](./longest-isoform.md) | **中文**

<!-- I18N:END -->

## -L / --longest

每个 gene 保留一个 transcript。长度度量:

- gene 含 CDS 时用 CDS 长度 (CDS 片段之和)，否则用 exon 长度 (exon 片段之和)。
- 不连续 CDS (同 ID 多行): 片段累加。
- 同一 transcript 下多个 CDS 变体 (不同 ID): 取最长变体，不累加。
- 长度相同取第一个遇到的。
- 自动检测 isoform 类型: 取文件中最多的 transcript 类类型（`mRNA`、`transcript`、`ncRNA`、`tRNA` 等）。数量并列时优先 `mRNA`，其次 `transcript`。
- 单 isoform 的 gene 保持不变。

示例数据 (demo.gff3):

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
chr1	src	mRNA	100	800	.	+	.	ID=tx02;Parent=gene01
chr1	src	CDS	100	400	.	+	0	ID=cds03;Parent=tx02
chr2	src	gene	200	600	.	-	.	ID=gene02
chr2	src	mRNA	200	600	.	-	.	ID=tx03;Parent=gene02
chr2	src	exon	200	400	.	-	.	ID=ex04;Parent=tx03
chr2	src	exon	450	600	.	-	.	ID=ex05;Parent=tx03
```

- tx01 CDS 长度: 151 + 251 = 402 bp。
- tx02 CDS 长度: 301 bp。
- tx03 exon 长度: 201 + 151 = 352 bp (无 CDS)。

```bash
# gene01 保留 tx01 (402 bp)，gene02 保留 tx03 (唯一 isoform)
./gffsub demo.gff3 --longest

# 只输出最长 transcript 行 (去掉 gene/exon)
./gffsub demo.gff3 --longest -t mRNA
```

## --longest-type TYPE

`--longest` 选择 isoform 时用的类型。与 `-t` 相互独立: `--longest-type` 决定哪些记录参与最长竞争，`-t` 过滤输出行。需要 `--longest`。

```bash
# ncRNA 文件: 每 gene 选最长 ncRNA，保留层级
./gffsub lnc.gff3 --longest --longest-type ncRNA

# 同上，但只输出 ncRNA 行
./gffsub lnc.gff3 --longest --longest-type ncRNA -t ncRNA
```

## -@ / --threads

按染色体并行。默认 1，最大 256。

```bash
./gffsub demo.gff3 --longest --threads 4
```
