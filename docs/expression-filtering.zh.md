# 表达式过滤

<!-- I18N:START -->

[English](./expression-filtering.md) | **中文**

<!-- I18N:END -->

## 语法

```
-I 'FIELD OP VALUE'   # 保留匹配的记录
-E 'FIELD OP VALUE'   # 丢弃匹配的记录
```

多个 `-I` 和 `-E` flag 之间为 AND。单个表达式内: `&&` (与)、`||` (或)、`!` (非)、`()`。

含空格或特殊字符的值用双引号括起。

## 运算符

| 运算符 | 含义 |
|--------|------|
| `==` | 字符串相等 |
| `!=` | 字符串不等 |
| `~` | 正则匹配 |
| `!~` | 正则不匹配 |
| `<` | 数值小于 |
| `<=` | 数值小于等于 |
| `>` | 数值大于 |
| `>=` | 数值大于等于 |

## 逻辑示例

```bash
-I 'length >= 100 && strand == "+"'
-I 'type == "gene" || type == "mRNA"'
-I '!(type == "gene")'
```

## 内置字段

| 字段 | 来源 |
|------|------|
| `seqid` | 第 1 列 |
| `source` | 第 2 列 |
| `type` | 第 3 列 |
| `start` | 第 4 列 |
| `end` | 第 5 列 |
| `length` | `end - start + 1` |
| `score` | 第 6 列 |
| `strand` | 第 7 列 |
| `phase` | 第 8 列 |
| `attrs` | 第 9 列 (原始) |

## 属性字段

| 字段 | GFF3 key |
|------|----------|
| `ID` | `ID` |
| `Name` | `Name` |
| `Parent` | `Parent` |
| `Alias` | `Alias` |
| `Dbxref` | `Dbxref` |
| `Note` | `Note` |
| `biotype` | `biotype` |
| `gene_id` | `gene_id` |
| `transcript_id` | `transcript_id` |
| `locus_tag` | `locus_tag` |
| `attr.KEY` | 第 9 列任意属性 |

GTF 输入: `gene_id` / `transcript_id` 从记录字段解析。其他属性用 `attr.` 前缀。

## 数值比较

`start`、`end`、`length`、`score` 作为数值参与 `<`、`<=`、`>`、`>=`、`==`、`!=`。

## 缺失属性

- `==`、`~`、`<`、`<=`、`>`、`>=` 不匹配 (记录被 `-I` 排除)。
- `!=` 和 `!~` 匹配 (记录被 `-I` 保留)。

## 示例

```bash
gffsub ann.gff3 -f exon -I 'length >= 100'
gffsub ann.gff3 -E 'length < 10'
gffsub ann.gff3 -f exon -I 'length >= 200 && strand == "+"'
```
