# window 子命令

## 语法

```
gffsub window <input.gff3> --id ID [--up N] [--down N] [--strand-aware]
```

返回与一个 feature 周基因组窗口重叠的所有记录。

## 参数

| 标志 | 默认值 | 说明 |
|------|---------|-------------|
| `--id ID` | 必填 | 按 ID 或 gene 名指定目标 feature |
| `--up N` | 0 | feature 起始上游碱基数 |
| `--down N` | 0 | feature 终止下游碱基数 |
| `--strand-aware` | 关 | 窗口跟随 feature 链方向 |

`--up` 有长别名 `--upstream`；`--down` 有长别名 `--downstream`。

## 窗口构造

不加 --strand-aware（默认）：
- 窗口起始 = feature 起始 - 上游
- 窗口终止 = feature 终止 + 下游
- 与方向无关：窗口向左向右延伸

加 --strand-aware：
- 正链：上游向起始左侧延伸，下游向终止右侧延伸
- 负链：两侧延伸互换。上游向终止右侧延伸，下游向起始左侧延伸

计算起始低于 1 时在重叠查找前钳为 1。

## --id 查找

--id 接受 feature ID 和 gene 名键，使用与 --name 相同的查找方式。ID 匹配优先。无 ID 匹配时再用 gene 名查找。索引的 gene 名属性为 Name、gene_id、locus_tag、Alias、Dbxref。

## 示例数据

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr2	src	gene	200	600	.	-	.	ID=gene02
```

## 命令

```bash
# gene01 上游 200 bp、下游 500 bp
./gffsub window demo.gff3 --id gene01 --up 200 --down 500
# 窗口：chr1:1 到 1500（起始 100 - 200 = -100，钳为 1；终止 1000 + 500 = 1500）

# 负链 gene 的链方向窗口
./gffsub window demo.gff3 --id gene02 --up 200 --down 500 --strand-aware
# 窗口：chr2:1 到 800（下游 500 向起始左侧延伸：200 - 500 = -300，钳为 1；上游 200 向终止右侧延伸：600 + 200 = 800）

# 按 gene 名查找
./gffsub window demo.gff3 --id BRCA1 --up 1000 --down 1000
```

## 默认模式下的窗口快捷方式

默认（非子命令）模式下，--up 和 --down 配合恰好一个 --id 时内部触发窗口查询：

```bash
./gffsub demo.gff3 --id gene01 --up 200 --down 500
```

要求恰好一个 --id，外加 --up、--down 或 --strand-aware 中任意一个。不接受其他筛选标志。传入其他筛选标志时报错中止，并列出允许集合：--id、--up/--upstream、--down/--downstream、--strand-aware。

## 校验

- --up 和 --down 必须为非负整数
- 整数后的尾部杂字符（如 "50abc"）被拒绝
- --id 必填
