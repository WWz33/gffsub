# window 子命令

<!-- I18N:START -->

[English](./window-subcommand.md) | **中文**

<!-- I18N:END -->

## 用法

```
gffsub window <input> -i ID [-u N] [-D N] [-a]
```

快捷形式（默认模式）：

```
gffsub <input> -i ID -u N -D N
```

## -i ID

目标 feature 的 ID 或 gene name。优先匹配 ID；未命中时按 gene name 查找（Name、gene_id、locus_tag、Alias、Dbxref）。必填。

## -u / -u N / --upstream N

feature 起始位置上游碱基数。默认 0。须为非负整数；尾部非法字符（如 `50abc`）会被拒绝。

## -D / -D N / --downstream N

feature 结束位置下游碱基数。默认 0。须为非负整数。

## -a / -a

- 正链：上游向 start 左侧延伸，下游向 end 右侧延伸。
- 负链：上游向 end 右侧延伸，下游向 start 左侧延伸。
- 不带此 flag：上游始终向左，下游始终向右。
- window 起始小于 1 时截断为 1。

## 不连续 feature

多行 CDS 共享同一 ID 时，window 覆盖所有片段的整体跨度。

## 示例

输入样例（demo.gff3）：
```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr2	src	gene	200	600	.	-	.	ID=gene02
```

```bash
# gene01 上游 200 bp，下游 500 bp
# window chr1:1-1500（start 100-200=-100，截断为 1；end 1000+500=1500）
./gffsub window demo.gff3 -i gene01 -u 200 -D 500

# 负链 gene02 的 strand-aware window
# window chr2:1-800（下游 500 向 start 左侧：200-500=-300，截断为 1；上游 200 向 end 右侧：600+200=800）
./gffsub window demo.gff3 -i gene02 -u 200 -D 500 -a

# 按 gene name 查找
./gffsub window demo.gff3 -i BRCA1 -u 1000 -D 1000

# 默认模式下的快捷形式
./gffsub demo.gff3 -i gene01 -u 200 -D 500
```

快捷形式要求恰好一个 `--id`，加上 `--up`、`--down`、`--strand-aware` 任意组合。不接受其他过滤 flag；传入会被拒绝并报错，提示允许的 flag 集合。
