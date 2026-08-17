# 属性过滤

<!-- I18N:START -->

[English](./attribute-filtering.md) | **中文**

<!-- I18N:END -->

## --where KEY=VALUE

对第 9 列属性做精确匹配。可重复；多个 flag 之间为 OR。

```bash
./gffsub demo.gff3 --where biotype=protein_coding
```

## --grep FIELD:PATTERN

对字段或属性做子串匹配。可重复；多个 flag 之间为 OR。

列字段: `seqid`、`source`、`type`、`start`、`end`、`length`、`score`、`strand`、`phase`、`attrs`。

属性字段: `ID`、`Name`、`Parent`、`Alias`、`Dbxref`、`Note`、`biotype`、`gene_id`、`transcript_id`、`locus_tag`。其他属性用 `attr.KEY` 形式。

```bash
./gffsub demo.gff3 --grep Name:ABC
```

## --grep-regex FIELD:REGEX

正则匹配 (ECMAScript)。可重复；多个 flag 之间为 OR。

```bash
./gffsub demo.gff3 --grep-regex ID:exon.*
```

## --grep-file FILE

从 FILE 读取 pattern，每行一个。需配合 `--grep-field FIELD` 指定目标字段。

- 不带 `--grep-file-regex`: 子串匹配。
- 带 `--grep-file-regex`: 正则匹配。

pattern 文件 (patterns.txt):

```
BRCA
TP53
```

```bash
./gffsub demo.gff3 --grep-field Name --grep-file patterns.txt
```

## -v / --invert-match

反转 grep 匹配: 保留不匹配的记录。

```bash
./gffsub demo.gff3 --grep type:mRNA -v
```

## --ignore-case

对 `--grep`、`--grep-regex`、`--grep-file`、`-I`、`-E` 启用大小写不敏感匹配。

```bash
./gffsub demo.gff3 --grep type:mrna --ignore-case
```

参见: [表达式过滤](./expression-filtering.zh.md) 了解数值比较和逻辑运算符。
