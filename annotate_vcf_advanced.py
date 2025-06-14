import pymongo
import requests
import vcf # PyVCF 庫
import time
import json
import logging
import gzip
from datetime import datetime # Added for report generation

# 設定日誌
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- 設定區塊 ---
# 請替換為您的 MongoDB Atlas 連接字串
# 連接字串範例: "mongodb+srv://<username>:<password>@cluster0.abcde.mongodb.net/?retryWrites=true&w=majority"
MONGO_URI = "mongodb+srv://<username>:<password>@cluster0.abcde.mongodb.net/?retryWrites=true&w=majority"
DB_NAME = "genomic_data"
COLLECTION_NAME = "na12878_lung_cancer_variants" # 更改 Collection 名稱以區分

# VCF 檔案路徑
VCF_FILE_PATH = "HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" # 請確保此檔案存在並可讀取

# MyVariant.info API 端點
# MyVariant.info 更多資訊: https://myvariant.info/
MYVARIANT_INFO_API_URL = "https://myvariant.info/v1/variant/"

# API 請求間的延遲 (秒) - 避免觸發 API 速率限制
API_DELAY_SECONDS = 0.1 # MyVariant.info 通常較快，但仍建議延遲

# VCF 讀取分塊大小
VCF_CHUNK_SIZE = 100 # 每批次處理 100 個變異

# 肺腺癌相關設定
LUNG_ADENOCARCINOMA_GENE_PANEL = [
    "EGFR", "KRAS", "TP53", "ALK", "ROS1", "BRAF", "MET", "RET", "ERBB2", "NF1", "STK11", "KEAP1"
]
# 肺腺癌相關的 HPO 術語 (概念性使用，這裡不直接用於過濾，而是作為臨床佐證的參考)
LUNG_ADENOCARCINOMA_HPO_TERMS = [
    "HP:0002665",  # 肺腺癌 (Lung adenocarcinoma)
    "HP:0002829",  # 腫瘤生長 (Neoplasm growth)
    "HP:0002828",  # 惡性腫瘤 (Malignant neoplasm)
    "HP:0100526",  # 肺腫瘤 (Lung neoplasm)
    "HP:0002104",  # 呼吸困難 (Dyspnea)
    "HP:0002204",  # 咳嗽 (Cough)
    "HP:0002861",  # 胸痛 (Chest pain)
    "HP:0002090",  # 血痰 (Hemoptysis)
    "HP:0002883",  # 肺功能異常 (Abnormal pulmonary function)
    "HP:0004305",  # 肺浸潤 (Pulmonary infiltrates)
    "HP:0006536",  # 肺纖維化 (Pulmonary fibrosis)
    "HP:0002206",  # 呼吸道感染 (Respiratory tract infection)
    "HP:0003002",  # 肺不張 (Atelectasis)
    "HP:0002105",  # 呼吸衰竭 (Respiratory failure)
    "HP:0002207",  # 慢性阻塞性肺疾病 (Chronic obstructive pulmonary disease)
]

# ACMG 準則頻率閾值 (簡化版，實際應用中需要更詳細的閾值和人群細分)
# 對於顯性遺傳，如果變異頻率在一般人群中高於此閾值，則可能是良性。
# 對於隱性遺傳，閾值可以更高。這裡作為通用閾值使用。
GNOMAD_PATHOGENICITY_THRESHOLD_AF = 0.0001 # 0.01%

# --- MongoDB 相關函數 ---
def get_mongo_collection():
    """
    連接 MongoDB Atlas 並返回指定的 Collection。
    """
    try:
        client = pymongo.MongoClient(MONGO_URI)
        db = client[DB_NAME]
        collection = db[COLLECTION_NAME]
        logging.info("成功連接到 MongoDB Atlas.")
        return collection
    except pymongo.errors.ConnectionFailure as e:
        logging.error(f"無法連接到 MongoDB Atlas: {e}")
        return None

# --- VCF 解析函數 (帶有分塊功能) ---
def parse_vcf_in_chunks(file_path, chunk_size):
    """
    解析 VCF 文件並以指定大小的列表 (chunk) 產生變異記錄。
    支援壓縮的 VCF 文件 (.vcf.gz) by letting PyVCF handle it.
    """
    try:
        # PyVCF's Reader can automatically handle .gz files when given a filename.
        vcf_reader = vcf.Reader(filename=file_path)
        
        chunk = []
        for record in vcf_reader:
            chunk.append(record)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
        if chunk:  # 產生最後一個不完整的 chunk
            yield chunk
        logging.info(f"成功處理 VCF 檔案: {file_path}。") # Adjusted log after full processing
    except FileNotFoundError:
        logging.error(f"VCF 檔案未找到: {file_path}")
        yield []  # 返回空列表
    except gzip.BadGzipFile as e:
        logging.error(f"VCF 檔案 '{file_path}' 不是有效的 GZIP 檔案或已損壞: {e}")
        logging.error("這可能發生在當一個以 .gz 結尾的檔案實際上是純文字檔時 (例如，如果測試檔案生成邏輯未正確創建 GZIP 檔案)。")
        logging.error("請檢查 VCF 檔案是否為正確的 GZIP 格式。")
        yield []
    except Exception as e:
        # Log the type of exception as well for better debugging
        logging.error(f"解析 VCF 檔案 '{file_path}' 時發生錯誤 ({type(e).__name__}): {e}")
        yield []  # 返回空列表

# --- 註釋函數 (使用 MyVariant.info) ---
def annotate_myvariant_info(chrom, pos, ref, alt):
    """
    使用 MyVariant.info API 查詢變異的綜合註釋資訊。
    使用 HGVS ID 格式: "chr7:g.55241707G>T"
    """
    # Construct the HGVS ID format for the variant
    variant_id = f"chr{chrom}:g.{pos}{ref}>{alt}"

    url = f"{MYVARIANT_INFO_API_URL}{variant_id}?fields=clinvar,gnomad,dbnsfp,dbsnp,ensembl.gene"
    try:
        response = requests.get(url)
        response.raise_for_status()  # Check for HTTP errors
        data = response.json()
        if "_id" in data:  # If the variant is successfully found
            return data
        else:
            return {"status": "not_found", "query": variant_id, "message": data.get("error", "未知錯誤")}
    except requests.exceptions.RequestException as e:
        logging.warning(f"MyVariant.info API 請求失敗 (變異 {variant_id}): {e}")
        return {"status": "api_error", "query": variant_id, "message": str(e)}
    except Exception as e:
        logging.warning(f"MyVariant.info 註釋解析失敗 (變異 {variant_id}): {e}")
        return {"status": "parsing_error", "query": variant_id, "message": str(e)}
    finally:
        time.sleep(API_DELAY_SECONDS)  # Respect API rate limits

# --- 致病性評估函數 (簡化 ACMG 準則) ---
def assess_pathogenicity(variant_doc, gene_panel, hpo_terms):
    """
    根據簡化版的 ACMG 準則、ClinVar 評級、gnomAD 頻率以及基因面板來評估變異致病性。
    同時提供相關臨床佐證。
    """
    pathogenicity = "VUS (意義不明變異)"
    evidence = []
    
    myvariant_anno = variant_doc.get("annotation_myvariant_info", {})
    
    # 1. 檢查 ClinVar 臨床意義 (PVS1_Strong / PS1_Strong)
    clinvar_data = myvariant_anno.get("clinvar", {})
    if clinvar_data:
        # ClinVar 可能有多個結果，取第一個或聚合結果
        # 修正: ClinVar rcv 字段可能是一個 dict 或 list of dicts
        if isinstance(clinvar_data.get("rcv"), list) and clinvar_data["rcv"]:
            clin_sig = clinvar_data["rcv"][0].get("clinical_significance", "").lower()
        elif isinstance(clinvar_data.get("rcv"), dict):
            clin_sig = clinvar_data["rcv"].get("clinical_significance", "").lower()
        else:
            clin_sig = "" # 沒有 rcv 資訊

        if "pathogenic" in clin_sig:
            pathogenicity = "Pathogenic (致病性)"
            evidence.append(f"ClinVar 報告為 '{clin_sig}'。")
        elif "likely pathogenic" in clin_sig:
            pathogenicity = "Likely Pathogenic (可能致病性)"
            evidence.append(f"ClinVar 報告為 '{clin_sig}'。")
        elif "benign" in clin_sig:
            pathogenicity = "Benign (良性)"
            evidence.append(f"ClinVar 報告為 '{clin_sig}'。")
        elif "likely benign" in clin_sig:
            pathogenicity = "Likely Benign (可能良性)"
            evidence.append(f"ClinVar 報告為 '{clin_sig}'。")
        elif "uncertain significance" in clin_sig:
            evidence.append(f"ClinVar 報告為 '{clin_sig}'。") # 仍是 VUS
            
    # 2. 檢查 gnomAD 等位基因頻率 (PM2_Moderate / BA1_StandAlone)
    gnomad_exome_af = myvariant_anno.get("gnomad_exome", {}).get("af")
    gnomad_genome_af = myvariant_anno.get("gnomad_genome", {}).get("af")

    if gnomad_exome_af is not None and gnomad_exome_af > GNOMAD_PATHOGENICITY_THRESHOLD_AF:
        evidence.append(f"gnomAD Exome 頻率 ({gnomad_exome_af:.4f}) 高於致病閾值 ({GNOMAD_PATHOGENICITY_THRESHOLD_AF:.4f})。")
        if gnomad_exome_af > 0.05: # 簡化 BA1 (>5% in controls)
            pathogenicity = "Benign (良性)"
            evidence.append(f"gnomAD Exome 頻率 ({gnomad_exome_af:.4f}) 非常高，高度提示良性。")
    elif gnomad_genome_af is not None and gnomad_genome_af > GNOMAD_PATHOGENICITY_THRESHOLD_AF:
        evidence.append(f"gnomAD Genome 頻率 ({gnomad_genome_af:.4f}) 高於致病閾值 ({GNOMAD_PATHOGENICITY_THRESHOLD_AF:.4f})。")
        if gnomad_genome_af > 0.05:
            pathogenicity = "Benign (良性)"
            evidence.append(f"gnomAD Genome 頻率 ({gnomad_genome_af:.4f}) 非常高，高度提示良性。")
    else: # 頻率非常低 (PM2_Moderate)
        if pathogenicity == "VUS (意義不明變異)": # 只有當 ClinVar 沒有明確指示時才考慮
            evidence.append(f"gnomAD Exome 頻率 ({gnomad_exome_af if gnomad_exome_af is not None else 'N/A'}) 和 Genome 頻率 ({gnomad_genome_af if gnomad_genome_af is not None else 'N/A'}) 低，支持致病性。")

    # 3. 檢查變異類型和基因面板 (PVS1_Strong / PM1_Moderate)
    # 這裡假設變異是單一基因的，MyVariant.info 的 ensembl.gene 應提供。
    gene_symbol = None
    if "ensembl" in myvariant_anno and myvariant_anno["ensembl"] and "gene" in myvariant_anno["ensembl"]:
        gene_info = myvariant_anno["ensembl"]["gene"]
        # 修正: ensembl.gene 可能是一個 dict (單個基因) 或 list (多個基因)
        if isinstance(gene_info, list) and gene_info:
            gene_symbol = gene_info[0].get("symbol") # 取第一個基因
        elif isinstance(gene_info, dict):
             gene_symbol = gene_info.get("symbol")

    if gene_symbol and gene_symbol.upper() in [g.upper() for g in gene_panel]:
        evidence.append(f"變異位於肺腺癌相關基因面板中的 '{gene_symbol}' 基因。")
        variant_doc["relevant_to_lung_adenocarcinoma"] = True
        
        # 檢查是否為截斷性變異 (nonsense, frameshift, splice site) - 模擬 PVS1
        dbnsfp_data = myvariant_anno.get("dbnsfp", {})
        if dbnsfp_data and isinstance(dbnsfp_data, list) and dbnsfp_data:
            # dbnsfp 可能為列表，因為一個變異可能有多個轉錄本影響
            consequences = []
            for item in dbnsfp_data:
                if "interpro_domain" in item and item["interpro_domain"] and isinstance(item["interpro_domain"], list):
                     for domain in item["interpro_domain"]:
                        if "description" in domain:
                            evidence.append(f"變異影響功能性蛋白質結構域: {domain['description']}。") # PM1_Moderate
                # 從 dbnsfp 獲取變異影響類型，通常在 `genecode` 或 `variant_effect` 等字段
                # MyVariant.info 的 dbnsfp 字段結構可能因版本而異，這裡嘗試獲取 'consequence'
                if "genecode" in item and item["genecode"]: 
                    if isinstance(item["genecode"], list):
                        for genecode_entry in item["genecode"]:
                            if "consequence" in genecode_entry:
                                consequences.append(genecode_entry["consequence"])
                    elif isinstance(item["genecode"], dict): # 單個 genecode entry
                        if "consequence" in item["genecode"]:
                            consequences.append(item["genecode"]["consequence"])
            
            # 簡化判斷截斷性變異
            is_truncating = False
            for cons in consequences:
                if any(c in cons.lower() for c in ["stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant"]):
                    is_truncating = True
                    break
            
            if is_truncating:
                evidence.append("變異為截斷性變異 (無義突變/移碼突變/剪接位點變異)，強烈支持致病性。")
                if pathogenicity not in ["Pathogenic (致病性)", "Likely Pathogenic (可能致病性)"]:
                    pathogenicity = "Likely Pathogenic (可能致病性)" # 提升致病性評級
        
    else:
        variant_doc["relevant_to_lung_adenocarcinoma"] = False
        if gene_symbol:
            evidence.append(f"變異位於 '{gene_symbol}' 基因，但此基因不在肺腺癌相關基因面板中。")
        else:
            evidence.append("無法識別變異所在的基因符號。")

    # 4. 計算證據分數或提升評級 (這裡簡化為基於關鍵證據提升)
    # 如果 ClinVar 已經是 Pathogenic/Likely Pathogenic，則優先採用。
    # 否則，綜合其他證據。

    # 如果是 VUS，且有低頻率和截斷性/功能影響，可以考慮提升為 Likely Pathogenic
    if pathogenicity == "VUS (意義不明變異)" and \
       any("低頻率" in s for s in evidence) and \
       any("截斷性變異" in s for s in evidence) and \
       variant_doc.get("relevant_to_lung_adenocarcinoma", False): # 確保基因與疾病相關
        pathogenicity = "Likely Pathogenic (可能致病性)"
        evidence.append("綜合低頻率和截斷性變異等證據，評級提升至可能致病性。")
        
    # 5. 添加 HPO 術語的連結 (概念性佐證)
    # 這裡不直接用 HPO 術語篩選變異，而是強調疾病上下文
    evidence.append(f"此分析以肺腺癌 ({', '.join(hpo_terms)}) 為主要疾病上下文進行。")

    # 6. 整合預測工具分數 (PP3_Supporting)
    if "dbnsfp" in myvariant_anno and myvariant_anno["dbnsfp"] and isinstance(myvariant_anno["dbnsfp"], list):
        for item in myvariant_anno["dbnsfp"]:
            sift_score = item.get("sift", {}).get("score")
            polyphen_score = item.get("polyphen", {}).get("score")
            
            if sift_score is not None and sift_score < 0.05: # SIFT 越低越有害 (通常是 < 0.05)
                evidence.append(f"SIFT 預測有害 (分數: {sift_score:.2f})。")
            if polyphen_score is not None and polyphen_score > 0.9: # PolyPhen-2 越高越有害 (通常是 > 0.9)
                evidence.append(f"PolyPhen-2 預測可能有害 (分數: {polyphen_score:.2f})。")
            # 可以添加更多預測工具的判斷

    return pathogenicity, evidence

# --- 模擬 MedGemma 報告生成函數 ---
def generate_report_with_medgemma(variants, vcf_file_name, hpo_terms_list):
    """
    模擬使用類似 MedGemma 的服務生成基因報告。
    接收已篩選的變異列表。
    """
    report_parts = []
    report_parts.append("<h3>Patient Gene Report (Simulated by MedGemma via Python Backend)</h3>")
    report_parts.append(f"<p><strong>Date:</strong> {datetime.now().strftime('%Y-%m-%d')}</p>")
    report_parts.append(f"<p><strong>Patient VCF:</strong> {vcf_file_name}</p>")
    report_parts.append(f"<p><strong>Associated HPO Terms:</strong> {', '.join(hpo_terms_list)}</p>")
    report_parts.append("<hr>")
    report_parts.append("<h4>Summary of Findings:</h4>")
    report_parts.append(
        f"<p>Analysis of the provided VCF file in the context of HPO terms: {', '.join(hpo_terms_list)}, "
        "has identified the following clinically significant variants. These findings should be "
        "correlated with clinical presentation and other diagnostic information by a qualified "
        "healthcare professional.</p>"
    )
    report_parts.append("<h4>Variant Details:</h4>")

    if not variants:
        report_parts.append("<p>No significant pathogenic or likely pathogenic variants were identified for report generation based on the provided data and HPO terms in this simulation.</p>")
    else:
        for i, variant_doc in enumerate(variants):
            gene_symbol = "N/A"
            myvariant_anno = variant_doc.get("annotation_myvariant_info", {})
            if "ensembl" in myvariant_anno and myvariant_anno["ensembl"] and "gene" in myvariant_anno["ensembl"]:
                gene_info = myvariant_anno["ensembl"]["gene"]
                if isinstance(gene_info, list) and gene_info:
                    gene_symbol = gene_info[0].get("symbol", "N/A")
                elif isinstance(gene_info, dict):
                    gene_symbol = gene_info.get("symbol", "N/A")
            
            report_parts.append(f"<p><strong>Variant {i + 1}: {gene_symbol} ({variant_doc['chrom']}:{variant_doc['pos']} {variant_doc['ref']}>{variant_doc['alt']})</strong></p>")
            report_parts.append("<ul>")
            report_parts.append(f"    <li><strong>Classification:</strong> {variant_doc.get('pathogenicity_assessment', 'N/A')}</li>")
            
            evidence_html = "<ul>"
            for ev in variant_doc.get('pathogenicity_evidence', []):
                evidence_html += f"<li>{ev}</li>"
            evidence_html += "</ul>"
            report_parts.append(f"    <li><strong>Annotation & Evidence:</strong> {evidence_html}</li>")
            
            # Simulated MedGemma insights
            report_parts.append(f"    <li><strong>Potential Implications (Simulated):</strong> Based on the {gene_symbol} gene's role and the nature of this variant, it is considered {variant_doc.get('pathogenicity_assessment', 'N/A').lower()}. This variant may contribute to the patient's phenotype or predisposition to conditions associated with the selected HPO terms. For example, alterations in {gene_symbol} are known to be involved in [simulated disease area, e.g., cancer development, metabolic disorders].</li>")
            report_parts.append(f"    <li><strong>Recommendations (Conceptual):</strong> Consider confirmatory testing if not already performed. Genetic counseling is recommended. Specific therapeutic options targeting pathways involving {gene_symbol} may be relevant (e.g., [simulated therapy type]).</li>")
            report_parts.append("</ul>")

    report_parts.append("<hr>")
    report_parts.append("<h4>Disclaimer:</h4>")
    report_parts.append("<p>This report is generated for demonstration purposes using simulated data and a conceptual MedGemma model via the Python backend. It is not a substitute for professional medical advice, diagnosis, or treatment. All interpretations and clinical decisions should be made by qualified healthcare providers.</p>")
    
    return "\\n".join(report_parts)

# --- 主註釋工作流程 ---
def run_annotation_workflow():
    """
    執行完整的 VCF 註釋和 MongoDB 儲存工作流程。
    """
    collection = get_mongo_collection()
    # 修正: 檢查 Collection 物件是否為 None
    if collection is None:
        logging.error("無法初始化 MongoDB Collection，退出。")
        return

    processed_count = 0
    inserted_count = 0

    for chunk in parse_vcf_in_chunks(VCF_FILE_PATH, VCF_CHUNK_SIZE):
        variants_to_insert = []
        for record in chunk:
            if len(record.ALT) != 1:
                logging.warning(f"跳過多等位基因變異: {record.CHROM}-{record.POS}-{record.REF}-{record.ALT}")
                continue
            
            chrom = str(record.CHROM)
            pos = record.POS
            ref = str(record.REF)
            alt = str(record.ALT[0])  # 獲取第一個變異等位基因

            variant_doc = {
                "chrom": chrom,
                "pos": pos,
                "id": record.ID if record.ID else ".",
                "ref": ref,
                "alt": alt,
                "qual": record.QUAL,
                "filter": record.FILTER,
                "info": dict(record.INFO),
                "samples": []
            }

            for sample in record.samples:
                sample_data = {
                    "sample_id": sample.sample,
                    "genotype": sample['GT']
                }
                variant_doc["samples"].append(sample_data)

            logging.info(f"正在註釋變異: {chrom}-{pos}-{ref}-{alt} (ID: {record.ID})")

            # --- 執行 MyVariant.info 註釋 ---
            myvariant_annotation = annotate_myvariant_info(chrom, pos, ref, alt)
            if myvariant_annotation.get("_id"):  # 確保有找到資訊
                variant_doc["annotation_myvariant_info"] = myvariant_annotation

                # --- 執行致病性評估 ---
                pathogenicity, evidence = assess_pathogenicity(
                    variant_doc, 
                    LUNG_ADENOCARCINOMA_GENE_PANEL, 
                    LUNG_ADENOCARCINOMA_HPO_TERMS
                )
                variant_doc["pathogenicity_assessment"] = pathogenicity
                variant_doc["pathogenicity_evidence"] = evidence

                variants_to_insert.append(variant_doc)
                inserted_count += 1
            else:
                logging.info(f"MyVariant.info 未找到變異資訊: {chrom}-{pos}-{ref}-{alt}")

            processed_count += 1

        # 批量插入當前 chunk 的變異
        if variants_to_insert:
            try:
                collection.insert_many(variants_to_insert)
                logging.info(f"已插入 {len(variants_to_insert)} 個變異到 MongoDB。總計已處理: {processed_count}，已插入: {inserted_count}")
            except Exception as e:
                logging.error(f"批量插入 MongoDB 時發生錯誤: {e}")
            variants_to_insert = []  # 清空以避免重複嘗試失敗的文檔

    logging.info(f"註釋工作流程完成。總計處理了 {processed_count} 個變異，成功插入 {inserted_count} 個變異。")
    logging.info(f"您可以透過 MongoDB Compass 或 Atlas UI 檢查 '{DB_NAME}.{COLLECTION_NAME}' Collection。")

    # --- 產生模擬 MedGemma 報告 ---
    if inserted_count > 0 and collection is not None:
        logging.info("正在從 MongoDB 檢索已註釋的致病性/可能致病性變異以產生報告...")
        try:
            report_variants = list(collection.find({
                "pathogenicity_assessment": {"$in": ["Pathogenic (致病性)", "Likely Pathogenic (可能致病性)"]}
            }))
            
            if report_variants:
                logging.info(f"找到 {len(report_variants)} 個變異用於報告生成。")
                # VCF_FILE_PATH 和 LUNG_ADENOCARCINOMA_HPO_TERMS 是全域變數
                simulated_report = generate_report_with_medgemma(report_variants, VCF_FILE_PATH, LUNG_ADENOCARCINOMA_HPO_TERMS)
                
                logging.info("--- SIMULATED MEDGEMMA REPORT (Python Backend) ---")
                # 為了方便在控制台閱讀，逐行打印（或整體打印，取決於日誌配置）
                # print 函数用於在腳本執行時直接輸出到控制台
                print("\\n" + "="*30 + " Simulated MedGemma Report " + "="*30)
                # simulated_report 包含 HTML 換行符，直接 print 可能不易閱讀，但日誌會保留格式
                # 這裡我們直接 print，讓使用者看到原始的 HTML 結構報告
                print(simulated_report.replace("\\n", "\\n")) # 確保 Python 的換行符被正確處理
                print("="*80 + "\\n")
                logging.info("模擬 MedGemma 報告已生成並記錄到控制台。")
                # 你也可以將 simulated_report 寫入一個 .html 檔案
                # with open("simulated_gene_report.html", "w", encoding="utf-8") as f_report:
                # f_report.write(simulated_report)
                # logging.info("模擬報告也已儲存到 simulated_gene_report.html")

            else:
                logging.info("在 MongoDB 中未找到致病性/可能致病性變異以產生報告。")
        except Exception as e:
            logging.error(f"從 MongoDB 檢索變異或生成報告時發生錯誤: {e}")
    elif collection is None:
        logging.warning("MongoDB collection 未初始化，跳過報告生成。")
    else:
        logging.info("沒有插入任何變異，跳過報告生成。")


if __name__ == "__main__":
    import os # Required for os.path.exists

    # 提醒使用者替換連接字串和 VCF 檔案
    if MONGO_URI == "YOUR_MONGODB_ATLAS_CONNECTION_STRING": # Ensure this default string is different if you have a real one set
        logging.error("請在程式碼中替換為您的 MongoDB Atlas 連接字串。")
    elif not VCF_FILE_PATH:
        logging.error("請提供 VCF 檔案路徑。")
    else:
        try:
            if not os.path.exists(VCF_FILE_PATH):
                logging.info(f"'{VCF_FILE_PATH}' 不存在，正在創建一個簡易的 GZIPPED 測試文件。")
                with gzip.open(VCF_FILE_PATH, 'wt', encoding='utf-8') as f: # 'wt' for text mode gzip
                    f.write("##fileformat=VCFv4.2\n")
                    f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
                    f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
                    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA12878\n")
                    f.write("1\t10000\t.\tA\tG\t100\tPASS\tDP=10\tGT\t0/1\n")
                    f.write("1\t10001\trs123\tT\tC\t90\tPASS\tDP=8\tGT\t1/1\n")
                    f.write("X\t20000\t.\tC\tA\t100\tPASS\tDP=15\tGT\t0/0\n")
                    f.write("7\t55242468\t.\tT\tTG\t100\tPASS\tDP=20\tGT\t0/1\n")
                    f.write("7\t55249071\t.\tC\tT\t100\tPASS\tDP=25\tGT\t0/1\n")
                    f.write("12\t25398284\t.\tG\tT\t100\tPASS\tDP=18\tGT\t0/1\n")
                    f.write("17\t7674903\t.\tC\tA\t100\tPASS\tDP=12\tGT\t0/1\n")
                logging.info(f"已創建一個簡易的 GZIPPED '{VCF_FILE_PATH}' 測試文件。")
            else:
                logging.info(f"'{VCF_FILE_PATH}' 文件已存在，將使用現有文件。")
        except Exception as e:
            logging.error(f"創建/檢查測試 VCF 文件時發生錯誤 ({type(e).__name__}): {e}")
            # Exit or prevent run_annotation_workflow if dummy file creation fails critically
            exit() # Or handle more gracefully
            
        run_annotation_workflow()
