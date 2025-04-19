import feedparser
from datetime import datetime, timedelta, timezone
import json
import requests
import os
import openai
from openai import OpenAI

# 肿瘤基因组与演化领域的PubMed RSS订阅
# (cancer AND evolution AND clone)
rss_url = 'https://pubmed.ncbi.nlm.nih.gov/rss/search/1fYfMzmpBkC9YVrZoBY9tdpF7fpHB74272RN58ypVcgf2sokfz/?limit=100&utm_campaign=pubmed-2&fc=20250419020659'

access_token = os.getenv('GITHUB_TOKEN')
aizex_api_key = os.getenv('AIZEX_API_KEY')

client = OpenAI(
    api_key=aizex_api_key,
    base_url="https://a1.aizex.me/v1"  # 尝试 https://a1.aizex.me/v1（国内优化），默认：https://aizex.top/v1
)

def extract_relevance_score(text, title):
    """使用英文提示词评估文章与肿瘤基因组克隆演化领域的相关性"""
    response = client.chat.completions.create(
        model="gpt-4o-mini",  # 根据 Aizex 支持的模型选择
        messages=[
            {"role": "system", "content": "You are an expert researcher in cancer genomics, tumor evolution, and clonal dynamics."},
            {"role": "user", "content": f"Based on the following title and abstract, rate the relevance of this paper to cancer genomics and clonal evolution research on a scale of 0-100:\n\nTitle: {title}\n\nAbstract: {text}\n\nProvide only a single numerical score (0-100) with no additional text."}
        ],
        max_tokens=10,
        temperature=0.3
    )

    generated_text = response.choices[0].message.content.strip()
    
    # 提取数字分数
    try:
        relevance_score = int(''.join(filter(str.isdigit, generated_text)))
        if relevance_score > 100:  # 防止异常值
            relevance_score = 100
    except:
        relevance_score = 0
        
    return relevance_score

def get_pubmed_abstracts(rss_url, days=7):
    """从PubMed RSS获取最近几天的文章摘要和元数据"""
    abstracts_with_urls = []
    
    # 计算日期阈值 (今天减去指定天数)
    cutoff_date = datetime.now(timezone.utc) - timedelta(days=days)
    
    # 解析PubMed RSS feed
    feed = feedparser.parse(rss_url)
    
    # 遍历PubMed RSS feed中的条目并提取信息
    for entry in feed.entries:
        try:
            # 获取条目的发布日期
            pub_date = None
            if hasattr(entry, 'published_parsed'):
                pub_date = datetime(*entry.published_parsed[:6], tzinfo=timezone.utc)
            elif hasattr(entry, 'updated_parsed'):
                pub_date = datetime(*entry.updated_parsed[:6], tzinfo=timezone.utc)
            
            # 如果没有日期或日期早于阈值，则跳过
            if not pub_date or pub_date < cutoff_date:
                continue
                
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
            
            # 添加发布日期到返回的数据中
            abstracts_with_urls.append({
                "title": title, 
                "abstract": abstract, 
                "doi": doi, 
                "link": link,
                "pub_date": pub_date.strftime('%Y-%m-%d')
            })
        except Exception as e:
            print(f"处理条目时出错: {e}")
            continue

    return abstracts_with_urls

def create_github_issue(title, body, access_token):
    """创建GitHub Issue"""
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
        return True
    else:
        print("创建Issue失败。状态码:", response.status_code)
        print("响应:", response.text)
        return False

def main():
    """主函数：获取本周文章，评分，筛选，并创建GitHub Issue"""
    # 从PubMed RSS feed获取本周摘要
    print("正在获取本周PubMed文章...")
    pubmed_abstracts = get_pubmed_abstracts(rss_url, days=7)  # 指定获取最近7天的文章
    print(f"获取到 {len(pubmed_abstracts)} 篇本周文章")
    
    # 如果没有获取到本周文章，则提前结束
    if not pubmed_abstracts:
        print("本周没有新的相关文章")
        return {
            "total_articles": 0,
            "selected_articles": 0,
            "top_articles": [],
            "github_issue_created": False
        }
    
    # 创建一个空列表来存储每个摘要及其评分
    articles_data = []
    
    print("正在评估文章相关性...")
    for i, abstract_data in enumerate(pubmed_abstracts):
        title = abstract_data["title"]
        abstract = abstract_data.get("abstract", "")
        print(f"正在处理第 {i+1}/{len(pubmed_abstracts)} 篇文章: {title[:50]}...")
        
        relevance_score = extract_relevance_score(abstract, title)
        doi = abstract_data["doi"]
        link = abstract_data.get("link", "")
        pub_date = abstract_data.get("pub_date", "")
    
        articles_data.append({
            "title": title,
            "relevance_score": relevance_score,
            "doi": doi,
            "link": link,
            "abstract": abstract,
            "pub_date": pub_date
        })
    
    # 按相关性分数排序并只保留前10篇
    sorted_articles = sorted(articles_data, key=lambda x: x["relevance_score"], reverse=True)
    top_articles = sorted_articles[:10]
    
    # 创建issue标题和内容
    current_date = datetime.now().strftime('%Y-%m-%d')
    one_week_ago = (datetime.now() - timedelta(days=7)).strftime('%Y-%m-%d')
    issue_title = f"肿瘤基因组克隆演化研究周报 ({one_week_ago} 至 {current_date})"
    issue_body = f"## {one_week_ago} 至 {current_date} 期间最相关的肿瘤基因组克隆演化研究：\n\n"
    
    for article in top_articles:
        title = article["title"]
        score = article["relevance_score"]
        doi = article.get("doi", "无DOI信息")
        link = article.get("link", "")
        abstract_excerpt = article.get("abstract", "")
        pub_date = article.get("pub_date", "")
    
        issue_body += f"### {title}\n"
        if pub_date:
            issue_body += f"- **发布日期**: {pub_date}\n"
        issue_body += f"- **相关性分数**: {score}/100\n"
        issue_body += f"- **DOI**: {doi}\n"
        if link:
            issue_body += f"- **链接**: {link}\n"
        if abstract_excerpt:
            issue_body += f"- **摘要**: {abstract_excerpt}\n"
        issue_body += "\n---\n\n"
    
    # 创建issue
    success = create_github_issue(issue_title, issue_body, access_token)
    
    # 输出结果到控制台
    print("\n===== 筛选结果 =====")
    print(f"共分析了 {len(pubmed_abstracts)} 篇本周文章，筛选出相关度最高的 {len(top_articles)} 篇")
    for i, article in enumerate(top_articles):
        print(f"{i+1}. {article['title']} (发布日期: {article['pub_date']}, 分数: {article['relevance_score']})")
    
    # 返回结果
    return {
        "total_articles": len(pubmed_abstracts),
        "selected_articles": len(top_articles),
        "top_articles": top_articles,
        "github_issue_created": success
    }

# 执行主函数并获取结果
if __name__ == "__main__":
    result = main()
    print(f"\n脚本执行完毕，结果: {json.dumps(result, ensure_ascii=False, indent=2)}")