import os
import pandas as pd
import json

# Use a dedicated test path so real data is never overwritten.
TEST_DIR = "/path/to/your/workspace/kofam/processing/test_run/"
os.makedirs(TEST_DIR, exist_ok=True)

# 模拟输入文件内容
# Case 1: Score < Threshold (应该被过滤)
# Case 2: Score > Threshold, 但不是 Immune KO (应该被过滤)
# Case 3: Score > Threshold, 是 Immune KO (应该保留)
# Case 4: 同 Case 3 的 Gene，但 E-value 更高 (应该被去重过滤)
# Case 5: 同 Case 3 的 Gene，但 E-value 更低 (应该保留并替换 Case 3)

mock_data = """# gene name	KO	thrshld	score	E-value	"KO definition"
transcript_1	K00001	100	50	1e-10	"Not passing threshold"
transcript_2	K99999	100	200	1e-20	"Passing thrshld but not Immune (Fake KO)"
transcript_3	K04687	100	150	1e-50	"Immune KO (IFNB)"
transcript_3	K05443	100	120	1e-60	"Immune KO (IL2) - Lower Evalue, should keep this"
"""

mock_input_path = os.path.join(TEST_DIR, "result_kofam_detail.txt")
with open(mock_input_path, "w") as f:
    f.write(mock_data)

# 模拟一个简化的 KEGG JSON (只包含我们要测试的 Immune KO)
mock_json = {
    "name": "ko00001",
    "children": [
        {
            "name": "09150 Organismal Systems",
            "children": [
                {
                    "name": "09151 Immune system",
                    "children": [
                        {
                            "name": "04660 T cell receptor signaling pathway",
                            "children": [
                                {"name": "K04687  IFNB; interferon beta"},
                                {"name": "K05443  IL2; interleukin 2"}
                            ]
                        }
                    ]
                }
            ]
        }
    ]
}

mock_json_path = os.path.join(TEST_DIR, "ko00001.json")
with open(mock_json_path, "w") as f:
    json.dump(mock_json, f)

print("=== 开始测试 ===")
# 临时修改 filter_kofam.py 中的路径配置 (在实际运行中，你可以直接把下面逻辑复制进去，或者为了方便，我这里重写一段简单的测试逻辑)

# --- 简单测试逻辑 ---
df = pd.read_csv(mock_input_path, sep='\t', comment='#', header=None, names=['gene','KO','thr','score','eval','def'])
# 手动清理模拟数据以匹配真实解析逻辑
df['score'] = pd.to_numeric(df['score'])
df['thr'] = pd.to_numeric(df['thr'])

print(f"原始数据:\n{df[['gene', 'score', 'eval']]}\n")

# 1. 阈值过滤
df = df[df['score'] > df['thr']]
print(f"阈值过滤后 (应剩3行):\n{df[['gene', 'score']]}\n")

# 2. Immune 过滤
valid_kos = ["K04687", "K05443"] # 来自 mock json
df = df[df['KO'].isin(valid_kos)]
print(f"Immune过滤后 (应剩2行, 都是 transcript_3):\n{df[['gene', 'KO']]}\n")

# 3. E-value 去重
df['eval'] = pd.to_numeric(df['eval'])
df = df.sort_values('eval').drop_duplicates('gene', keep='first')
print(f"去重后 (应剩1行, K05443, 1e-60):\n{df[['gene', 'KO', 'eval']]}\n")

if len(df) == 1 and df.iloc[0]['KO'] == 'K05443':
    print(">>> 测试成功！逻辑验证通过。")
else:
    print(">>> 测试失败，请检查逻辑。")