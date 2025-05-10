import os
import json
import time
import requests
import feedparser
from datetime import datetime, timedelta, timezone
from xml.etree import ElementTree as ET

from Bio import Entrez
from openai import OpenAI

# ========== 用户配置部分 ==========
# 以下部分为可配置项，用户可以根据自己的研究方向进行修改

# 研究方向名称（用于GitHub Issue标题）
RESEARCH_AREA = "肿瘤基因组克隆结构和系统发生树分析"

# API密钥配置 - 通过环境变量传入
ENTREZ_EMAIL = "shehuizhuyitese@gmail.com"  # 必填
ENTREZ_API_KEY = os.getenv("NCBI_API_KEY")  # 可选，但建议配置
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")  # 如果使用OpenAI
AIZEX_API_KEY = os.getenv("AIZEX_API_KEY")    # 如果使用AIZEX
GITHUB_TOKEN = os.getenv("GITHUB_TOKEN")      # 必填，用于创建Issue
GITHUB_REPO_OWNER = "HushWay"      # GitHub仓库所有者
GITHUB_REPO_NAME = "TrackResearch"             # GitHub仓库名称

# API服务配置
API_BASE_URL = "https://a1.aizex.me/v1"  # 可以替换为OpenAI或其他兼容服务
AI_MODEL = "gpt-4o-mini"  # 使用的模型名称

# PubMed/PMC检索关键词配置 - 根据研究方向修改
TIAB_TERMS = [
    "clonal structure[tiab]", "clone structure[tiab]", 
    "phylogenetic tree[tiab]", "clonal phylogeny[tiab]",
    "clonal evolution[tiab]", "clonal architecture[tiab]",
    "clonal reconstruction[tiab]", "subclonal reconstruction[tiab]",
    "tumor phylogeny[tiab]", "cancer phylogeny[tiab]",
    "subclonal architecture[tiab]", "phylogenetic analysis[tiab]",
    "lineage tracing[tiab]", "cancer lineage[tiab]",
    "tumor lineage[tiab]", "clonality analysis[tiab]"
]

MESH_TERMS = [
    "Phylogeny[MeSH Terms]", 
    "Genetic Heterogeneity[MeSH Terms]",
    "Neoplasms/genetics[MeSH Terms]",
    "Genomics/methods[MeSH Terms]",
    "Computational Biology/methods[MeSH Terms]"
]

# 高影响力杂志列表 - 根据研究领域修改
HIGH_IMPACT_JOURNALS = [
    "Nature", "Science", "Cell", "Nature Genetics", 
    "Nature Methods", "Nature Biotechnology", "Cancer Cell", 
    "Nature Communications", "Science Advances", 
    "Genome Biology", "Genome Research", "Cancer Discovery",
    "Nature Reviews Cancer", "Cell Reports", "PNAS", 
    "Molecular Cell", "Bioinformatics", "Nucleic Acids Research"
]

# 高影响力杂志加分
JOURNAL_IMPACT_BONUS = 10  # 高影响力杂志文章的额外分数

# 预过滤相关度关键词 - 根据研究方向修改
HIGH_RELEVANCE_KEYWORDS = [
    'clonal structure', 'phylogenetic tree', 'clonal phylogeny',
    'tumor phylogeny', 'cancer phylogeny', 'phylogenetic analysis',
    'clone structure', 'clonal architecture', 'subclonal architecture',
    'clonal reconstruction', 'subclonal reconstruction',
    'clone evolution', 'tumor lineage', 'cancer lineage',
    'lineage tracing', 'clonality analysis', 'ancestral clone',
    'clonal dynamics', 'mutation tree', 'evolutionary tree',
    'phylogenetic inference', 'clonal inference', 'treeomics',
    'clonevol', 'clone tracking', 'clonal complexity',
    'subclone mapping', 'clonal tracking'
]

# 预过滤匹配阈值
KEYWORD_MATCH_THRESHOLD = 2  # 匹配多少个关键词才保留

# LLM评分提示词模板 - 根据研究方向修改
SCORING_SYSTEM_PROMPT = "你是{research_area}专家。你需要识别与{research_area}相关的研究。"

SCORING_USER_PROMPT = """对下面的文章进行0-100分打分，评估其与{research_area}的相关性。

评分标准：
- 90-100：直接研究{research_area}的核心文章
- 70-89：与{research_area}有密切关系的方法学文章
- 40-69：提到{research_area}但主要研究其他方面的文章
- 0-39：基本不相关的文章

标题: {title}
摘要: {abstract}

只返回分数数字，不要解释。"""

# 批量评分提示词模板
BATCH_SCORING_USER_PROMPT = """对以下多篇文章进行打分（0-100），评估它们与{research_area}的相关性：

评分标准：
- 90-100分：直接研究{research_area}的核心文章
- 70-89分：与{research_area}有密切关系的方法学文章
- 40-69分：提到{research_area}但主要研究其他方面的文章
- 0-39分：基本不相关的文章

{articles}

请用JSON格式返回结果，键为文章ID，值为分数，例如：
{{"article_0": 85, "article_1": 45, ...}}"""

# 检索配置
DEEP_SEARCH_MONTHS = 12     # 深度检索的时间范围（月）
DEEP_SEARCH_RETMAX = 100    # 深度检索的最大文章数
RECENT_SEARCH_DAYS = 7      # 最近检索的时间范围（天）
RECENT_SEARCH_RETMAX = 50   # 最近检索的最大文章数
ARXIV_MAX_RESULTS = 50      # arXiv检索的最大文章数

# 缓存文件配置
CACHE_FILE = "article_cache.json"  # 缓存文件路径

# 日志配置
LOG_FILE = "research_tracker.log"  # 日志文件路径
LOG_LEVEL = "INFO"                 # 日志级别 (DEBUG, INFO, WARNING, ERROR)

# LLM API配置
BATCH_SIZE = 5              # 批量评分的批次大小
MAX_RETRIES = 3             # API调用最大重试次数
BACKOFF_FACTOR = 2          # 重试间隔递增因子

# ========== 初始化部分 ==========

# 初始化API客户端
client = OpenAI(
    api_key=AIZEX_API_KEY or OPENAI_API_KEY,
    base_url=API_BASE_URL
)

# 设置Entrez邮箱和API密钥
Entrez.email = ENTREZ_EMAIL
if ENTREZ_API_KEY:
    Entrez.api_key = ENTREZ_API_KEY

# 设置日志
import logging
logging.basicConfig(
    level=getattr(logging, LOG_LEVEL),
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_FILE, mode='a'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# ========== 功能实现部分 ==========

def build_pubmed_query():
    """构建PubMed查询语句"""
    tiab = " OR ".join(TIAB_TERMS)
    mesh = " OR ".join(MESH_TERMS)
    return f"({tiab}) OR ({mesh})"

# === 缓存管理 ===
def normalize_string(s):
    """简单的字符串规范化，便于比较"""
    if not s:
        return ""
    return "".join(s.lower().split())

def load_article_cache(cache_file=CACHE_FILE):
    """加载文章缓存"""
    if os.path.exists(cache_file):
        try:
            with open(cache_file, "r", encoding="utf-8") as f:
                cache = json.load(f)
            logger.info(f"已加载 {len(cache.get('articles', []))} 篇缓存文章")
            return cache
        except Exception as e:
            logger.error(f"加载缓存失败: {e}")
    return {"last_update": "", "articles": []}

def save_article_cache(cache, cache_file=CACHE_FILE):
    """保存文章缓存"""
    cache["last_update"] = datetime.now().isoformat()
    try:
        with open(cache_file, "w", encoding="utf-8") as f:
            json.dump(cache, f, ensure_ascii=False, indent=2)
        logger.info(f"已保存 {len(cache.get('articles', []))} 篇文章到缓存")
    except Exception as e:
        logger.error(f"保存缓存失败: {e}")

def merge_with_cache(new_articles, cache_file=CACHE_FILE):
    """将新抓取的文章与缓存合并，避免重复请求"""
    cache = load_article_cache(cache_file)
    cached_articles = cache.get("articles", [])
    
    # 构建标题和DOI索引，用于快速查找
    cached_titles = {normalize_string(art.get("title", "")): art for art in cached_articles if art.get("title")}
    cached_dois = {art.get("doi", ""): art for art in cached_articles if art.get("doi")}
    
    merged_articles = []
    new_count = 0
    
    for article in new_articles:
        title = normalize_string(article.get("title", ""))
        doi = article.get("doi", "")
        
        # 检查是否已经在缓存中
        cached_art = None
        if doi and doi in cached_dois:
            cached_art = cached_dois[doi]
        elif title and title in cached_titles:
            cached_art = cached_titles[title]
            
        if cached_art and "score" in cached_art:
            # 使用缓存的分数
            article["score"] = cached_art["score"]
            merged_articles.append(article)
        else:
            # 新文章，需要打分
            merged_articles.append(article)
            new_count += 1
    
    logger.info(f"合并后共 {len(merged_articles)} 篇文章，其中 {new_count} 篇新文章需要打分")
    return merged_articles, new_count

def update_cache_with_scores(articles, cache_file=CACHE_FILE):
    """更新缓存，添加新打分的文章"""
    cache = load_article_cache(cache_file)
    cached_articles = cache.get("articles", [])
    
    # 构建标题和DOI索引
    cached_titles = {normalize_string(art.get("title", "")): i for i, art in enumerate(cached_articles) if art.get("title")}
    cached_dois = {art.get("doi", ""): i for i, art in enumerate(cached_articles) if art.get("doi")}
    
    # 更新现有文章或添加新文章
    for article in articles:
        if "score" not in article:
            continue  # 跳过未打分的文章
            
        title = normalize_string(article.get("title", ""))
        doi = article.get("doi", "")
        
        if doi and doi in cached_dois:
            # 更新已有文章
            cached_articles[cached_dois[doi]]["score"] = article["score"]
        elif title and title in cached_titles:
            # 更新已有文章
            cached_articles[cached_titles[title]]["score"] = article["score"]
        else:
            # 添加新文章
            cached_articles.append(article)
    
    cache["articles"] = cached_articles
    save_article_cache(cache, cache_file)
    return cache

# === 文献抓取 ===
def fetch_pubmed(months_back=DEEP_SEARCH_MONTHS, retmax=DEEP_SEARCH_RETMAX):
    """深度检索PubMed，获取过去几个月的相关文献"""
    logger.info(f"开始深度检索PubMed, 时间范围: {months_back}个月, 最大数量: {retmax}")
    try:
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
            logger.warning("深度检索未找到结果")
            return []

        results = fetch_pubmed_from_pmids(pmids)
        logger.info(f"深度检索完成，获取了 {len(results)} 篇文章")
        return results
    except Exception as e:
        logger.error(f"深度检索出错: {str(e)}")
        return []

def fetch_pubmed_from_pmids(pmids):
    """根据PMID列表获取文章详情"""
    if not pmids:
        return []
    try:
        fetch = Entrez.efetch(
            db="pubmed",
            id=",".join(pmids),
            rettype="xml"
        )
        root = ET.fromstring(fetch.read())
        results = []
        for article in root.findall(".//PubmedArticle"):
            art = {}
            art["title"] = article.findtext(".//ArticleTitle") or ""
            art["abstract"] = "".join([t.text or "" for t in article.findall(".//AbstractText")])
            art["doi"] = article.findtext(".//ArticleId[@IdType='doi']") or ""
            pmid = article.findtext(".//PMID")
            art["link"] = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ""
            art["journal"] = (
                article.findtext(".//Journal/Title")
                or article.findtext(".//MedlineJournalInfo/JournalTitle")
                or ""
            )
            # 提取完整日期
            pubdate = article.find(".//PubDate")
            if pubdate is not None:
                year  = pubdate.findtext("Year") or ""
                month = pubdate.findtext("Month") or ""
                day   = pubdate.findtext("Day") or ""
                dt = None
                try:
                    dt = datetime.strptime(f"{year} {month} {day}", "%Y %b %d")
                except:
                    try:
                        dt = datetime.strptime(f"{year}-{month}-{day}", "%Y-%m-%d")
                    except:
                        dt = None
                art["pub_date"] = dt.strftime("%Y-%m-%d") if dt else year
            else:
                art["pub_date"] = ""
            results.append(art)
        return results
    except Exception as e:
        logger.error(f"从PMID获取文章详情失败: {str(e)}")
        return []

def fetch_trending_pubmed(days=RECENT_SEARCH_DAYS, retmax=RECENT_SEARCH_RETMAX):
    """获取最近几天的PubMed文章"""
    logger.info(f"开始寻新检索PubMed, 时间范围: {days}天, 最大数量: {retmax}")
    try:
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
        results = fetch_pubmed_from_pmids(pmids)
        logger.info(f"寻新检索完成，获取了 {len(results)} 篇文章")
        return results
    except Exception as e:
        logger.error(f"寻新检索出错: {str(e)}")
        return []

def fetch_europe_pmc(days=RECENT_SEARCH_DAYS, pageSize=RECENT_SEARCH_RETMAX):
    """获取Europe PMC最近文章"""
    logger.info(f"开始检索Europe PMC, 时间范围: {days}天, 最大数量: {pageSize}")
    try:
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
        resp = requests.get(url, timeout=30)
        data = resp.json().get("resultList", {}).get("result", [])
        results = []
        for item in data:
            art = {
                "title": item.get("title", ""),
                "abstract": item.get("abstractText", ""),
                "doi": item.get("doi", ""),
                "link": item.get("fullTextUrlList", {}).get("fullTextUrl", [{}])[0].get("url", ""),
                # Europe PMC 本身已提供 YYYY-MM-DD 格式
                "pub_date": item.get("firstPublicationDate", ""),
                "journal": item.get("journalTitle", "")
            }
            results.append(art)
        logger.info(f"Europe PMC检索完成，获取了 {len(results)} 篇文章")
        return results
    except Exception as e:
        logger.error(f"Europe PMC检索出错: {str(e)}")
        return []

def fetch_arxiv(max_results=ARXIV_MAX_RESULTS):
    """获取arXiv相关预印本"""
    logger.info(f"开始检索ArXiv, 最大数量: {max_results}")
    try:
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
            journal_ref = entry.get('arxiv_journal_ref', '')
            art = {
                "title": entry.title,
                "abstract": entry.summary,
                "doi": next((l.href for l in entry.links if l.href.startswith("http://dx.doi.org/")), ""),
                "link": entry.link,
                "pub_date": entry.published[:10],
                "journal": journal_ref or "arXiv"
            }
            results.append(art)
        logger.info(f"ArXiv检索完成，获取了 {len(results)} 篇文章")
        return results
    except Exception as e:
        logger.error(f"ArXiv检索出错: {str(e)}")
        return []

# === 预过滤 ===
def simple_keyword_filter(articles, threshold=KEYWORD_MATCH_THRESHOLD):
    """使用关键词匹配进行预过滤，减少需要API评分的文章数量"""
    filtered_articles = []
    for art in articles:
        text = (art.get('title', '') + ' ' + art.get('abstract', '')).lower()
        # 计算匹配的关键词数量
        matches = sum(1 for keyword in HIGH_RELEVANCE_KEYWORDS if keyword.lower() in text)
        if matches >= threshold:
            filtered_articles.append(art)
    
    logger.info(f"预过滤: 从 {len(articles)} 篇文章中筛选出 {len(filtered_articles)} 篇")
    return filtered_articles

# === 评分 & 筛选 ===
def extract_relevance_score(article, max_retries=MAX_RETRIES, backoff_factor=BACKOFF_FACTOR):
    """单篇文章评分，添加重试机制"""
    text = article.get("abstract", "")
    title = article.get("title", "")
    
    for attempt in range(max_retries):
        try:
            system_prompt = SCORING_SYSTEM_PROMPT.format(research_area=RESEARCH_AREA)
            user_prompt = SCORING_USER_PROMPT.format(
                research_area=RESEARCH_AREA,
                title=title,
                abstract=text
            )
            
            response = client.chat.completions.create(
                model=AI_MODEL,
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_prompt}
                ],
                max_tokens=5,
                temperature=0
            )
            txt = response.choices[0].message.content.strip()
            try:
                score = min(100, int(''.join(filter(str.isdigit, txt))))
                article["score"] = score
                return score
            except:
                article["score"] = 0
                return 0
        except Exception as e:
            logger.warning(f"评分尝试 {attempt+1}/{max_retries} 失败: {str(e)}")
            if attempt < max_retries - 1:
                sleep_time = backoff_factor ** attempt
                logger.info(f"等待 {sleep_time} 秒后重试...")
                time.sleep(sleep_time)
            else:
                logger.error("达到最大重试次数，返回默认分数")
                article["score"] = 30
                return 30  # 返回默认中等分数
    
    article["score"] = 30
    return 30  # 以防万一，确保返回一个默认值

def batch_extract_relevance_scores(articles, batch_size=BATCH_SIZE, max_retries=MAX_RETRIES, backoff_factor=BACKOFF_FACTOR):
    """批量处理文章评分，每次API调用处理多篇文章"""
    all_scores = {}
    
    # 将文章分成批次
    for i in range(0, len(articles), batch_size):
        batch = articles[i:i+batch_size]
        batch_ids = [f"article_{j}" for j in range(i, i+len(batch))]
        
        # 构建批量请求的文章内容
        articles_content = ""
        for idx, art in zip(batch_ids, batch):
            articles_content += f"--- {idx} ---\n"
            articles_content += f"标题: {art.get('title', '')}\n"
            # 摘要可能过长，取前1000字符
            abstract = art.get('abstract', '')
            articles_content += f"摘要: {abstract[:1000]}{'...' if len(abstract) > 1000 else ''}\n\n"
        
        # 构建批量请求的完整提示
        system_prompt = SCORING_SYSTEM_PROMPT.format(research_area=RESEARCH_AREA)
        user_prompt = BATCH_SCORING_USER_PROMPT.format(
            research_area=RESEARCH_AREA,
            articles=articles_content
        )
        
        # 添加重试逻辑
        for attempt in range(max_retries):
            try:
                logger.info(f"处理批次 {i//batch_size + 1}/{(len(articles)-1)//batch_size + 1}")
                response = client.chat.completions.create(
                    model=AI_MODEL,
                    messages=[
                        {"role": "system", "content": system_prompt},
                        {"role": "user", "content": user_prompt}
                    ],
                    max_tokens=500,  # 增加token以容纳多篇文章的评分结果
                    temperature=0
                )
                
                # 解析结果
                response_text = response.choices[0].message.content.strip()
                try:
                    # 提取JSON部分 (可能会有其他文本)
                    import re
                    json_str = re.search(r'{.*}', response_text, re.DOTALL)
                    if json_str:
                        import json
                        scores_dict = json.loads(json_str.group(0))
                        
                        # 更新总分数字典
                        for j, art_id in enumerate(batch_ids):
                            if art_id in scores_dict:
                                article_index = i + j
                                if article_index < len(articles):
                                    all_scores[article_index] = min(100, int(scores_dict[art_id]))
                            else:
                                logger.warning(f"警告: 未找到文章ID {art_id}的分数")
                        
                        # 成功处理，跳出重试循环
                        break
                    else:
                        raise ValueError("响应中未找到JSON格式数据")
                except Exception as e:
                    logger.error(f"解析批次结果失败: {e}")
                    logger.debug(f"原始响应: {response_text}")
                    if attempt < max_retries - 1:
                        continue
                    else:  # 最后一次尝试，为批次中文章分配默认分数
                        for j in range(len(batch)):
                            article_index = i + j
                            if article_index < len(articles):
                                all_scores[article_index] = 30  # 默认中等分数
            
            except Exception as e:
                logger.error(f"批次 {i//batch_size + 1} 尝试 {attempt+1}/{max_retries} 失败: {str(e)}")
                if attempt < max_retries - 1:
                    sleep_time = backoff_factor ** attempt
                    logger.info(f"等待 {sleep_time} 秒后重试...")
                    time.sleep(sleep_time)
                else:
                    logger.error("达到最大重试次数，为批次中所有文章分配默认分数")
                    for j in range(len(batch)):
                        article_index = i + j
                        if article_index < len(articles):
                            all_scores[article_index] = 30
        
        # 在批次之间等待一小段时间，避免API速率限制
        time.sleep(1)
    
    # 将分数应用到文章上
    for i, art in enumerate(articles):
        art["score"] = all_scores.get(i, 30)  # 如果没有分数，使用默认值
    
    return articles

def journal_impact_bonus(article):
    """为高影响力杂志的文章增加额外分数"""
    journal = article.get("journal", "").lower()
    
    # 精确匹配高影响力杂志
    for high_impact in HIGH_IMPACT_JOURNALS:
        if high_impact.lower() in journal:
            return JOURNAL_IMPACT_BONUS
    
    return 0

# === GitHub Issue ===
def create_github_issue(title, body):
    """创建GitHub Issue来展示结果"""
    try:
        url = f"https://api.github.com/repos/{GITHUB_REPO_OWNER}/{GITHUB_REPO_NAME}/issues"
        headers = {
            "Authorization": f"token {GITHUB_TOKEN}",
            "Accept": "application/vnd.github.v3+json"
        }
        data = {"title": title, "body": body}
        r = requests.post(url, headers=headers, json=data, timeout=30)
        if r.status_code == 201:
            logger.info(f"成功创建Issue: {title}")
            return True
        else:
            logger.error(f"创建Issue失败，状态码: {r.status_code}, 响应: {r.text}")
            return False
    except Exception as e:
        logger.error(f"创建Issue过程中发生异常: {str(e)}")
        return False

# ========== 主流程 ==========
def main():
    try:
        # 1) 深度检索过去几个月经典文献
        deep = fetch_pubmed(months_back=DEEP_SEARCH_MONTHS, retmax=DEEP_SEARCH_RETMAX)
        
        # 2) 寻新检索过去几天最新进展（PubMed + Europe PMC + arXiv）
        new_pm = fetch_trending_pubmed(days=RECENT_SEARCH_DAYS, retmax=RECENT_SEARCH_RETMAX)
        new_ep = fetch_europe_pmc(days=RECENT_SEARCH_DAYS, pageSize=RECENT_SEARCH_RETMAX)
        new_ax = fetch_arxiv(max_results=ARXIV_MAX_RESULTS)

        all_articles = deep + new_pm + new_ep + new_ax
        logger.info(f"总抓取文章数: {len(all_articles)}")

        # 3) 去重
        unique_articles = []
        seen_titles = set()
        for art in all_articles:
            norm = normalize_string(art.get("title", ""))
            if norm and norm not in seen_titles:
                seen_titles.add(norm)
                unique_articles.append(art)
        logger.info(f"去重后文章数: {len(unique_articles)}")

        # 4) 与缓存合并
        merged_articles, new_count = merge_with_cache(unique_articles)

        # 5) 对新文章评分
        if new_count > 0:
            new_articles = [a for a in merged_articles if "score" not in a]
            # 预过滤减少需要评分的文章数量
            filtered = simple_keyword_filter(new_articles, threshold=KEYWORD_MATCH_THRESHOLD)
            # 如果过滤后文章太少，降低阈值
            if len(filtered) < 10:
                filtered = simple_keyword_filter(new_articles, threshold=1)
            # 如果还是太少，使用所有新文章
            if len(filtered) < 10:
                filtered = new_articles
                
            logger.info(f"需要打分的新文章: {len(filtered)}")
            
            if filtered:
                # 批量评分
                batch_extract_relevance_scores(filtered, batch_size=BATCH_SIZE, 
                                              max_retries=MAX_RETRIES, 
                                              backoff_factor=BACKOFF_FACTOR)
                                              
                # 将分数更新到合并文章列表
                for art in filtered:
                    if "score" in art:
                        for m in merged_articles:
                            if normalize_string(m.get("title","")) == normalize_string(art.get("title","")):
                                m["score"] = art["score"]
                                
                # 更新缓存
                update_cache_with_scores(filtered)

        # 6) 为没有分数的文章添加默认分数
        for art in merged_articles:
            if "score" not in art:
                art["score"] = 0

        # 7) 计算最终分数（含高影响力杂志加分）
        for art in merged_articles:
            impact_bonus = journal_impact_bonus(art)
            art["final_score"] = art.get("score", 0) + impact_bonus

        # 8) 排序获取Top10
        top10 = sorted(merged_articles, key=lambda x: x.get("final_score", 0), reverse=True)[:10]

        # 9) 创建GitHub Issue
        today    = datetime.now().strftime("%Y-%m-%d")
        week_ago = (datetime.now() - timedelta(days=RECENT_SEARCH_DAYS)).strftime("%Y-%m-%d")
        title = f"{RESEARCH_AREA}研究周报 ({week_ago} 至 {today})"
        body  = f"## {week_ago} – {today} 最相关 Top10：\n\n"
        
        for art in top10:
            is_high_impact = journal_impact_bonus(art) > 0
            body += (
                f"### {art['title']}\n"
                f"- 杂志: {art.get('journal','未知')}{' 🌟' if is_high_impact else ''}\n"
                f"- 发表日期: {art.get('pub_date','未知')}\n"
                f"- 相关性分数: {art.get('score',0)}/100\n"
                f"- DOI: {art.get('doi','无')}\n"
                f"- 链接: {art.get('link','')}\n\n"
            )
        
        success = create_github_issue(title, body)
        logger.info("Issue 已创建" if success else "Issue 创建失败")

    except Exception as e:
        logger.error(f"主流程执行出错: {str(e)}", exc_info=True)
        error_title = f"{RESEARCH_AREA}研究周报执行出错 - {datetime.now().strftime('%Y-%m-%d')}"
        error_body  = f"## 执行过程中出现错误\n\n```\n{str(e)}\n```\n\n请检查日志获取详细信息。"
        create_github_issue(error_title, error_body)

if __name__ == "__main__":
    main()