import pandas as pd
import os
import sys

print("开始处理 GFF 文件转换为 BED 格式...")

# 检查输入文件是否存在
print("检查输入文件...")
if not os.path.exists('immune_genes.gff'):
    print("错误: immune_genes.gff 文件不存在")
    sys.exit(1)

# 读取 GFF 文件
print("读取 GFF 文件...")
try:
    # 使用标准的 GFF3 列名
    column_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    df = pd.read_csv('immune_genes.gff', sep='\t', names=column_names, comment='#')
    print(f"成功读取 {len(df)} 行数据")
except Exception as e:
    print(f"错误: 读取文件失败 - {e}")
    sys.exit(1)

# 提取基因ID
print("提取基因ID...")
df['gene_id'] = df['attribute'].str.extract(r'ID=(transcript:[^;]+)')
df['parent_gene'] = df['attribute'].str.extract(r'Parent=(gene:[^;]+)')

# 使用parent_gene作为主要ID，如果没有则使用gene_id
df['final_id'] = df['parent_gene'].fillna(df['gene_id'])

# 清理数据
print("清理和格式化数据...")
df = df.dropna(subset=['final_id'])  # 移除没有有效ID的行

# 准备 BED 格式数据
print("转换为 BED 格式...")
bed_df = pd.DataFrame()
bed_df['chrom'] = df['seqname']
bed_df['chromStart'] = df['start'] - 1  # GFF是1-based，BED是0-based
bed_df['chromEnd'] = df['end']
bed_df['name'] = df['final_id']
bed_df['score'] = df['score'].replace('.', '0')
bed_df['strand'] = df['strand']

# 计算转录本长度
bed_df['length'] = bed_df['chromEnd'] - bed_df['chromStart']

print("筛选最长的转录本...")
# 对于每个基因，保留最长的转录本
idx = bed_df.groupby(['name'])['length'].transform(max) == bed_df['length']
genes_filtered = bed_df[idx].drop_duplicates(subset=['name'])

print(f"筛选后保留 {len(genes_filtered)} 个基因")

# 重新排列列顺序（标准BED格式）
bed_output = genes_filtered[['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']].copy()

# 保存结果
print("保存 BED 文件...")
try:
    bed_output.to_csv("interproscan_genes.bed", sep='\t', index=False, header=False)
    print(f"成功保存 {len(bed_output)} 行数据到 'interproscan_genes.bed'")
    
    # 显示前几行作为预览
    print("\n前5行预览:")
    print(bed_output.head().to_string(index=False))
    
    # 输出统计信息
    print(f"\n处理统计:")
    print(f"   染色体数量: {bed_output['chrom'].nunique()}")
    print(f"   基因数量: {len(bed_output)}")
    print(f"   正链基因: {len(bed_output[bed_output['strand'] == '+'])}")
    print(f"   负链基因: {len(bed_output[bed_output['strand'] == '-'])}")
    
except Exception as e:
    print(f"错误: 保存文件失败 - {e}")
    sys.exit(1)

print("处理完成!")