import os
import json
import time
import requests
import feedparser
from datetime import datetime, timedelta, timezone
from xml.etree import ElementTree as ET

from Bio import Entrez
from openai import OpenAI

# ========== 配置 ==========
Entrez.email = "shehuizhuyitese@gmail.com"              # 必填
Entrez.api_key = os.getenv("NCBI_API_KEY")           # 可选，加快并发
AIZEX_API_KEY = os.getenv("AIZEX_API_KEY")
GITHUB_TOKEN = os.getenv("GITHUB_TOKEN")

client = OpenAI(
    api_key=AIZEX_API_KEY,
    base_url="https://a1.aizex.me/v1"
)

# ========== 检索配置 ==========
# TIAB 关键词和 MeSH 词列表
TIAB_TERMS = [
    "cancer[tiab]", "tumor[tiab]",
    "evolution[tiab]", "clonal heterogeneity[tiab]",
    "clonal expansion[tiab]", "intra-tumor heterogeneity[tiab]"
]
MESH_TERMS = [
    "Neoplasm Proteomics[MeSH Terms]",
    "Genetic Heterogeneity[MeSH Terms]"
]

def build_pubmed_query():
    """构造 PubMed E-utilities 查询语句"""
    tiab = " OR ".join(TIAB_TERMS)
    mesh = " OR ".join(MESH_TERMS)
    return f"({tiab}) OR ({mesh})"

# ========== 文献抓取 ==========
def fetch_pubmed(months_back=1, retmax=100):
    """深度检索：过去 months_back 个月内所有匹配条目"""
    end = datetime.utcnow()
    start = end - timedelta(days=30*months_back)
    query = build_pubmed_query()
    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=retmax,
        datetype="pdat",
        mindate=start.strftime("%Y/%m/%d"),
        maxdate=end.strftime("%Y/%m/%d")
    )
    record = Entrez.read(handle)
    pmids = record["IdList"]
    if not pmids:
        return []
    # 批量 efetch
    fetch = Entrez.efetch(
        db="pubmed",
        id=",".join(pmids),
        rettype="xml"
    )
    root = ET.fromstring(fetch.read())
    results = []
    for article in root.findall(".//PubmedArticle"):
        art = {}
        art["title"] = article.findtext(".//ArticleTitle")
        art["abstract"] = "".join([t.text or "" for t in article.findall(".//AbstractText")])
        art["doi"] = article.findtext(".//ArticleId[@IdType='doi']") or ""
        art["link"] = f"https://pubmed.ncbi.nlm.nih.gov/{article.findtext('.//PMID')}/"
        art["pub_date"] = article.findtext(".//PubDate/Year") or ""
        results.append(art)
    return results

def fetch_trending_pubmed(days=7, retmax=100):
    """寻新检索：过去 days 天内最新条目（RSS 或 E-utilities）"""
    end = datetime.utcnow()
    start = end - timedelta(days=days)
    query = build_pubmed_query()
    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=retmax,
        datetype="pdat",
        mindate=start.strftime("%Y/%m/%d"),
        maxdate=end.strftime("%Y/%m/%d")
    )
    record = Entrez.read(handle)
    pmids = record["IdList"]
    return fetch_pubmed_from_pmids(pmids)

def fetch_pubmed_from_pmids(pmids):
    """辅助：给定 pmid 列表，批量获取详细信息"""
    if not pmids: return []
    fetch = Entrez.efetch(
        db="pubmed",
        id=",".join(pmids),
        rettype="xml"
    )
    root = ET.fromstring(fetch.read())
    results = []
    for article in root.findall(".//PubmedArticle"):
        art = {}
        art["title"] = article.findtext(".//ArticleTitle")
        art["abstract"] = "".join([t.text or "" for t in article.findall(".//AbstractText")])
        art["doi"] = article.findtext(".//ArticleId[@IdType='doi']") or ""
        art["link"] = f"https://pubmed.ncbi.nlm.nih.gov/{article.findtext('.//PMID')}/"
        art["pub_date"] = article.findtext(".//PubDate/Year") or ""
        results.append(art)
    return results

def fetch_europe_pmc(days=7, pageSize=100):
    """Europe PMC 抓取 bioRxiv/medRxiv 等预印本和已发表文章"""
    end = datetime.utcnow().strftime("%Y-%m-%d")
    start = (datetime.utcnow() - timedelta(days=days)).strftime("%Y-%m-%d")
    query = build_pubmed_query()
    url = (
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
        f"?query={requests.utils.quote(query)}"
        f"+PUB_TYPE:preprint+OR+PUB_TYPE:journal"
        f"&dateFilter=PUB_CITED_TS:[{start}%20TO%20{end}]"
        "&format=json"
        f"&pageSize={pageSize}"
    )
    resp = requests.get(url)
    data = resp.json().get("resultList", {}).get("result", [])
    results = []
    for item in data:
        art = {
            "title": item.get("title"),
            "abstract": item.get("abstractText", ""),
            "doi": item.get("doi", ""),
            "link": item.get("fullTextUrlList", {}).get("fullTextUrl", [{}])[0].get("url", ""),
            "pub_date": item.get("firstPublicationDate", "")
        }
        results.append(art)
    return results

def fetch_arxiv(max_results=100):
    """ArXiv 抓取 q-bio 相关预印本"""
    # 先获取关键词列表，去掉[tiab]标记
    terms = [t.split("[")[0] for t in TIAB_TERMS]
    
    # 使用 requests.utils.quote 对整个查询进行编码
    raw_query = "cat:q-bio*+AND+(" + "+OR+".join([requests.utils.quote(term) for term in terms]) + ")"
    
    url = (
        "http://export.arxiv.org/api/query"
        f"?search_query={raw_query}"
        f"&start=0&max_results={max_results}"
        "&sortBy=submittedDate&sortOrder=descending"
    )
    feed = feedparser.parse(url)
    results = []
    for entry in feed.entries:
        art = {
            "title": entry.title,
            "abstract": entry.summary,
            "doi": next((l.href for l in entry.links if l.href.startswith("http://dx.doi.org/")), ""),
            "link": entry.link,
            "pub_date": entry.published[:10]
        }
        results.append(art)
    return results

# ========== 评分 & 筛选 ==========
def extract_relevance_score(text, title):
    response = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[
            {"role": "system", "content": "You are an expert in cancer genomics and clonal evolution."},
            {"role": "user", "content": (
                "Rate 0-100 relevance to cancer clonal evolution.\n\n"
                f"Title: {title}\nAbstract: {text}\n\n"
                "Respond with only the number."
            )}
        ],
        max_tokens=5,
        temperature=0
    )
    txt = response.choices[0].message.content.strip()
    try:
        return min(100, int(''.join(filter(str.isdigit, txt))))
    except:
        return 0

# ========== GitHub Issue ==========
def create_github_issue(title, body):
    url = "https://api.github.com/repos/HushWay/TrackResearch/issues"
    headers = {
        "Authorization": f"token {GITHUB_TOKEN}",
        "Accept": "application/vnd.github.v3+json"
    }
    data = {"title": title, "body": body}
    r = requests.post(url, headers=headers, json=data)
    return r.status_code == 201

# ========== 主流程 ==========
def main():
    # 1) 深度检索过去 12 个月经典文献
    deep = fetch_pubmed(months_back=12, retmax=200)
    # 2) 寻新检索过去 7 天最新进展（PubMed + Europe PMC + arXiv）
    new_pm = fetch_trending_pubmed(days=7, retmax=100)
    new_ep = fetch_europe_pmc(days=7, pageSize=100)
    new_ax = fetch_arxiv(max_results=100)

    all_articles = deep + new_pm + new_ep + new_ax
    print(f"总抓取文章数: {len(all_articles)}")

    # 3) 打分并排序
    for art in all_articles:
        art["score"] = extract_relevance_score(art["abstract"], art["title"])
    top10 = sorted(all_articles, key=lambda x: x["score"], reverse=True)[:10]

    # 4) 准备 GitHub Issue
    today = datetime.now().strftime("%Y-%m-%d")
    week_ago = (datetime.now() - timedelta(days=7)).strftime("%Y-%m-%d")
    title = f"肿瘤克隆演化综述 ({week_ago} 至 {today})"
    body = f"## {week_ago} – {today} 最相关 Top10：\n\n"
    for art in top10:
        body += (
            f"### {art['title']}\n"
            f"- 发布: {art.get('pub_date','')}\n"
            f"- 分数: {art['score']}/100\n"
            f"- DOI: {art.get('doi','无')}\n"
            f"- 链接: {art.get('link','')}\n\n"
        )
    success = create_github_issue(title, body)
    print("Issue 已创建" if success else "Issue 创建失败")

if __name__ == "__main__":
    main()
