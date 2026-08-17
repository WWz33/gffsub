# 表达式筛选

<!-- I18N:START -->

[English](./expression-filtering.md) | **中文**

<!-- I18N:END -->

gffsub 通过 `-I` / `--include-expr`（保留匹配记录）和 `-E` / `--exclude-expr`（丢弃匹配记录）按字段值、数值比较和属性匹配筛选记录。

## 语法

```
-I 'FIELD OP VALUE'
```

多个 `-I` 和 `-E` 之间是 AND 逻辑。单个表达式内支持 `&&`（与）、`||`（或）、`!`（非）和括号 `()`。

### 操作符

| 操作符 | 含义 | 示例 |
|--------|------|------|
| `==` | 字符串相等 | `-I 'type == "exon"'` |
| `!=` | 字符串不等 | `-I 'type != "gene"'` |
| `~` | 正则匹配 | `-I 'ID ~ "exon.*"'` |
| `!~` | 正则不匹配 | `-I 'Name !~ "test.*"'` |
| `<` | 数值小于 | `-I 'length < 100'` |
| `<=` | 数值小于等于 | `-I 'length <= 200'` |
| `>` | 数值大于 | `-I 'length > 500'` |
| `>=` | 数值大于等于 | `-I 'length >= 1000'` |

值含空格或特殊字符时用双引号包裹。

### 逻辑运算符

在单个表达式中组合条件：

```bash
# 长度 >= 100 且正链
-I 'length >= 100 && strand == "+"'

# gene 或 transcript
-I 'type == "gene" || type == "mRNA"'

# 不是 gene
-I '!(type == "gene")'
```

## 字段

### 内置字段

| 字段 | 类型 | 说明 |
|------|------|------|
| `seqid` | 字符串 | 第 1 列（序列名） |
| `source` | 字符串 | 第 2 列（来源） |
| `type` | 字符串 | 第 3 列（feature 类型） |
| `start` | 数值 | 第 4 列（1-based 起始） |
| `end` | 数值 | 第 5 列（1-based 终止） |
| `length` | 数值 | `end - start + 1` |
| `score` | 数值 | 第 6 列 |
| `strand` | 字符串 | 第 7 列（`+ - . ?`） |
| `phase` | 字符串 | 第 8 列（`0 1 2 .`） |
| `attrs` | 字符串 | 原始第 9 列 |

### 属性字段

GFF3 第 9 列属性按键名访问。不在内置列表中的键用 `attr.` 前缀。

| 字段 | GFF3 键 | 说明 |
|------|---------|------|
| `ID` | `ID` | feature 标识 |
| `Name` | `Name` | 显示名 |
| `Parent` | `Parent` | 父 feature ID |
| `Alias` | `Alias` | 别名 |
| `Dbxref` | `Dbxref` | 数据库交叉引用 |
| `Note` | `Note` | 自由文本注释 |
| `gene_id` | `gene_id` | 基因 ID（GTF / GFF3） |
| `transcript_id` | `transcript_id` | 转录本 ID（GTF / GFF3） |
| `biotype` | `biotype` | 生物类型 |
| `locus_tag` | `locus_tag` | locus 标签 |
| `attr.KEY` | `KEY` | 任意第 9 列属性 |

GTF 输入：`gene_id`、`transcript_id`、`ID`、`Parent` 从解析器合成的字段获取。其他 GTF 属性用 `attr.` 前缀加 GTF 键名。

## 数值比较

`start`、`end`、`length`、`score` 字段在 `<`、`<=`、`>`、`>=`、`==`、`!=` 中按数值处理。非数值内容比较失败（不崩溃，记录不匹配）。

## 缺失属性

当记录缺少引用的属性时：

- `==`、`~`、`<`、`<=`、`>`、`>=` 不匹配（`-I` 排除该记录）。
- `!=` 和 `!~` 匹配（`-I` 保留该记录）。

即 `-I 'attr.biotype != "lncRNA"'` 保留没有 `biotype` 属性的记录，将缺失视为"不等于"。

## 按长度筛选示例

按基因组跨度筛选 feature：

```bash
# 长度 >= 100 bp 的 exon
gffsub ann.gff3 -f exon -I 'length >= 100'

# 长度 < 200 bp 的 CDS
gffsub ann.gff3 -f CDS -I 'length < 200'

# 长度 > 1 kb 的 mRNA
gffsub ann.gff3 -f mRNA -I 'length > 1000'

# 长度在 50 到 200 bp 之间的 exon
gffsub ann.gff3 -f exon -I 'length >= 50 && length <= 200'

# 丢弃长度 < 10 bp 的 feature
gffsub ann.gff3 -E 'length < 10'
```

## 与 feature 类型组合

`-f` / `--feature` 在表达式求值前按第 3 列类型筛选：

```bash
# 正链上的长 exon（>= 200 bp）
gffsub ann.gff3 -f exon -I 'length >= 200 && strand == "+"'
```
