#!/usr/bin/env bash
set -euo pipefail

# === 输入与输出 ===
GFF="NM_genes.renamed.gff3"                     # 你的注释文件路径
IPR_IMMUNE_CSV="interproscan_immune_results.csv" # 第一步筛好的免疫相关 InterPro 结果
OUT_PREFIX="immune"                              # 输出前缀，可改

# === 基本检查 ===
[[ -s "$GFF" ]] || { echo "ERROR: GFF not found: $GFF" >&2; exit 1; }
[[ -s "$IPR_IMMUNE_CSV" ]] || { echo "ERROR: CSV not found: $IPR_IMMUNE_CSV" >&2; exit 1; }

# === 临时文件 ===
WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

LIST_MRNA="$WORKDIR/listmRNA"
ATTR_MRNA="$WORKDIR/last_column_only_mRNA"
FINAL_MRNA="$WORKDIR/final_mRNA"
IMMUNE_TXT_IDS="$WORKDIR/immune_transcript_ids"
IMMUNE_MRNA_GENE_IDS="${OUT_PREFIX}_mRNA_gene_ids"
IMMUNE_GFF="${OUT_PREFIX}_genes.gff"

echo "[1/6] Extract mRNA features from GFF..."
# 更严格：按第3列等于 mRNA 过滤，避免误匹配
awk -F'\t' '$3=="mRNA"' "$GFF" > "$LIST_MRNA"

echo "[2/6] Keep only attributes (column 9)..."
# GFF3 规范中的 attributes 在第9列；用 $9 比 $NF 更稳
awk -F'\t' '{print $9}' "$LIST_MRNA" > "$ATTR_MRNA"

echo "[3/6] Keep first two attributes (usually ID=...;Parent=...)..."
awk -F';' '{print $1 ";" $2}' "$ATTR_MRNA" > "$FINAL_MRNA"

echo "[4/6] Get transcript IDs from InterPro immune CSV (skip header)..."
# 取第一列为转录本/蛋白ID；跳过表头，并去重
awk -F',' 'NR>1 && $1!="" {print $1}' "$IPR_IMMUNE_CSV" | sort -u > "$IMMUNE_TXT_IDS"

echo "[5/6] Filter mRNA: keep only those whose transcript ID is in InterPro list..."
# 拆分分隔符为 ; 或 = ： ID=transcript ; Parent=gene  -> $2 是 transcript 值
awk -F'[;=]' 'NR==FNR{ids[$1]; next} ($2 in ids)' \
  "$IMMUNE_TXT_IDS" "$FINAL_MRNA" > "$IMMUNE_MRNA_GENE_IDS"

echo "[6/6] Recover full mRNA lines from GFF..."
# 再提取转录本 ID（$2），到原始 mRNA 行里做精确字面匹配
grep -F -f <(awk -F'[=;]' '{print $2}' "$IMMUNE_MRNA_GENE_IDS") \
  "$LIST_MRNA" > "$IMMUNE_GFF"

# 将关键产物拷贝到当前目录
cp "$IMMUNE_MRNA_GENE_IDS" .
cp "$IMMUNE_GFF" .

echo "Done."
echo "  - 转录本-基因对应: $(pwd)/$(basename "$IMMUNE_MRNA_GENE_IDS")"
echo "  - 免疫相关 mRNA 的完整 GFF 行: $(pwd)/$(basename "$IMMUNE_GFF")"
