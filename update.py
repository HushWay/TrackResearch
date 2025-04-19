import feedparser
from datetime import datetime, timedelta, timezone
import json
import requests
import os
import openai
from openai import OpenAI

# 生物信息学、肿瘤基因组与演化领域的PubMed RSS订阅
# 以下URL需要替换为您定制的PubMed搜索RSS
# ((cancer) AND (evolution)) 
# AND ((("last month"[EDAT]) NOT "last week"[EDAT])) 
# AND (
#     Nature[journal] OR 
#     Science[journal] OR 
#     Cell[journal] OR 
#     "Nature Genetics"[journal] OR 
#     "Cancer Cell"[journal] OR 
#     "Genome Biology"[journal] OR 
#     "Molecular Biology and Evolution"[journal] OR
#     "NEW ENGLAND JOURNAL OF MEDICINE"[journal] OR 
#     LANCET[journal] OR 
#     JAMA[journal] OR 
#     "JOURNAL OF CLINICAL ONCOLOGY"[journal] OR 
#     "Cancer Discovery"[journal]
# )


rss_url = 'https://pubmed.ncbi.nlm.nih.gov/rss/search/1-7FJ7kkHKrLcV1C5C04Tw1hGZlKKoM-Rd4u97RQk_RR4dggc7/?limit=100&utm_campaign=pubmed-2&fc=20250418230358'

access_token = os.getenv('GITHUB_TOKEN')
aizex_api_key = os.getenv('AIZEX_API_KEY')

client = OpenAI(
    api_key=aizex_api_key,
    base_url="https://a1.aizex.me/v1"  # 尝试 https://a1.aizex.me/v1（国内优化），默认：https://aizex.top/v1
)

def extract_scores(text):
# Aizex 提供多种模型选择，可以根据您的需求和预算选择不同的模型：
#   经济型选择：
#       gpt-4o-mini：适合基础文献评分，性价比高
#       deepseek-v3：国产大模型，能力不错且价格较低
#   高质量选择：
#       claude-3-7-sonnet：更全面的文献分析能力
#       gpt-4-1106-preview：深度理解科学文献的能力强
    response = client.chat.completions.create(
        model="gpt-4o-mini",  # 根据 Aizex 支持的模型选择
        messages=[
            {"role": "system", "content": "你是一位生物信息学、肿瘤基因组学和演化生物学领域的专家研究员。你擅长筛选有趣、新颖且重要的研究成果。"},
            {"role": "user", "content": f"基于以下文本'{text}'，请为该论文提供两个评分：\n"
                                        "1. 研究分数 (0-100)：基于研究创新性、方法严谨性、数据可靠性和在该领域的重要性\n"
                                        "2. 社会影响分数 (0-100)：基于公众关注度、临床应用潜力、政策相关性和社会影响\n"
                                        "请以以下格式提供评分：\n"
                                        "研究分数: <分数>\n"
                                        "社会影响分数: <分数>"}
        ],
        max_tokens=100,
        temperature=0.5
    )

    generated_text = response.choices[0].message.content.strip()  

    # 提取研究分数
    research_score_start = generated_text.find("研究分数:")
    if research_score_start == -1:
        research_score_start = generated_text.find("Research Score:")
    
    if research_score_start != -1:
        research_score = generated_text[research_score_start+len("研究分数:"):].split("\n")[0].strip()
        if not research_score:
            research_score = generated_text[research_score_start+len("Research Score:"):].split("\n")[0].strip()
    else:
        research_score = "N/A"

    # 提取社会影响分数
    social_impact_score_start = generated_text.find("社会影响分数:")
    if social_impact_score_start == -1:
        social_impact_score_start = generated_text.find("Social Impact Score:")
    
    if social_impact_score_start != -1:
        social_impact_score = generated_text[social_impact_score_start+len("社会影响分数:"):].strip()
        if not social_impact_score:
            social_impact_score = generated_text[social_impact_score_start+len("Social Impact Score:"):].strip()
    else:
        social_impact_score = "N/A"

    return research_score, social_impact_score

def get_pubmed_abstracts(rss_url):
    abstracts_with_urls = []

    # 解析PubMed RSS feed
    feed = feedparser.parse(rss_url)

    # 计算一周前的日期
    one_week_ago = datetime.now(timezone.utc) - timedelta(weeks=1)

    # 遍历PubMed RSS feed中的条目并提取摘要和URL
    for entry in feed.entries:
        try:
            # 获取条目的发布日期
            published_date = datetime.strptime(entry.published, '%a, %d %b %Y %H:%M:%S %z')

            # 如果发布日期在一周内，提取摘要和URL
            if published_date >= one_week_ago:
                # 获取条目的标题、摘要和DOI
                title = entry.title
                abstract = entry.get('summary', '') if hasattr(entry, 'summary') else entry.get('description', '')
                if not abstract and hasattr(entry, 'content'):
                    abstract = entry.content[0].value
                
                doi = entry.get('dc_identifier', 'No DOI available')
                if not isinstance(doi, str):
                    doi = str(doi)
                
                # 获取PubMed链接
                link = entry.link
                
                abstracts_with_urls.append({"title": title, "abstract": abstract, "doi": doi, "link": link})
        except Exception as e:
            print(f"处理条目时出错: {e}")
            continue

    return abstracts_with_urls

# 从PubMed RSS feed获取摘要
pubmed_abstracts = get_pubmed_abstracts(rss_url)

# 创建一个空列表来存储每个摘要及其评分
new_articles_data = []

for abstract_data in pubmed_abstracts:
    title = abstract_data["title"]
    abstract = abstract_data.get("abstract", "")
    research_score, social_impact_score = extract_scores(abstract)
    doi = abstract_data["doi"]
    link = abstract_data.get("link", "")

    new_articles_data.append({
        "title": title,
        "research_score": research_score,
        "social_impact_score": social_impact_score,
        "doi": doi,
        "link": link
    })
    
# 创建issue标题和内容
issue_title = f"生物信息学/肿瘤基因组/演化领域周报 - {datetime.now().strftime('%Y-%m-%d')}"
issue_body = "## 本周生物信息学、肿瘤基因组学和演化生物学领域文献评分结果：\n\n"

# 按研究分数排序
sorted_articles = sorted(new_articles_data, key=lambda x: int(x["research_score"].split()[0]) if x["research_score"] not in ["N/A", ""] and x["research_score"].split()[0].isdigit() else 0, reverse=True)

for article_data in sorted_articles:
    title = article_data["title"]
    research_score = article_data["research_score"]
    social_impact_score = article_data["social_impact_score"]
    doi = article_data.get("doi", "无DOI信息")
    link = article_data.get("link", "")

    issue_body += f"### {title}\n"
    issue_body += f"- **研究分数**: {research_score}\n"
    issue_body += f"- **社会影响分数**: {social_impact_score}\n"
    issue_body += f"- **DOI**: {doi}\n"
    if link:
        issue_body += f"- **链接**: {link}\n"
    issue_body += "\n---\n\n"

def create_github_issue(title, body, access_token):
    # 修改为您自己的仓库信息
    url = f"https://api.github.com/repos/HushWay/TrackResearch/issues"
    headers = {
        "Authorization": f"token {access_token}",
        "Accept": "application/vnd.github.v3+json"
    }
    payload = {
        "title": title,
        "body": body
    }

    response = requests.post(url, headers=headers, data=json.dumps(payload))

    if response.status_code == 201:
        print("Issue创建成功！")
    else:
        print("创建Issue失败。状态码:", response.status_code)
        print("响应:", response.text)

# 创建issue
create_github_issue(issue_title, issue_body, access_token)