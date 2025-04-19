import openai
from openai import OpenAI

# 设置 API 密钥
aizex_api_key = "your_api_key_here"  # 替换为您的 API 密钥

# 创建 OpenAI 客户端
client = OpenAI(
    api_key=aizex_api_key,
    base_url="https://aizex.top/v1"  # 或其他 Aizex 端点
)

# 测试论文摘要
test_abstract = """
[在这里粘贴一个论文摘要用于测试]
"""

# 使用 API 获取评分
response = client.chat.completions.create(
    model="gpt-4o-mini",  # 或其他可用模型
    messages=[
        {"role": "system", "content": "你是一位生物信息学、肿瘤基因组学和演化生物学领域的专家研究员。你擅长筛选有趣、新颖且重要的研究成果。"},
        {"role": "user", "content": f"基于以下文本'{test_abstract}'，请为该论文提供两个评分：\n"
                                   "1. 研究分数 (0-100)：基于研究创新性、方法严谨性、数据可靠性和在该领域的重要性\n"
                                   "2. 社会影响分数 (0-100)：基于公众关注度、临床应用潜力、政策相关性和社会影响\n"
                                   "请以以下格式提供评分：\n"
                                   "研究分数: <分数>\n"
                                   "社会影响分数: <分数>"}
    ],
    max_tokens=100,
    temperature=0.5
)

# 打印结果
print(response.choices[0].message.content)