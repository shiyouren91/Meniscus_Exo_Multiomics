"""
Nature期刊投稿文件全面审核脚本
检查所有关键要求和格式规范
"""

import sys
import io
import re
from docx import Document

# 设置UTF-8输出
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def check_title_page(doc):
    """检查标题页"""
    print("\n" + "="*80)
    print("📋 标题页检查 (Title_Page)")
    print("="*80)
    
    issues = []
    warnings = []
    passes = []
    
    # 提取所有文本
    full_text = ' '.join([para.text for para in doc.paragraphs])
    
    # 检查标题
    title_patterns = [
        r'meniscus|injury|repair|exosome',
        r'degeneration|fibrocartilage'
    ]
    title_found = any(re.search(pattern, full_text, re.I) for pattern in title_patterns)
    if title_found:
        passes.append("✅ 标题存在")
    else:
        issues.append("❌ 未找到标题")
    
    # 检查作者（使用英文名或拼音）
    authors_patterns = [
        r'(任师友|Ren\s+Shiyou|Shiyou\s+Ren)',
        r'(关思瑶|Guan\s+Siyao|Siyao\s+Guan)',
        r'(于海阳|Yu\s+Haiyang|Haiyang\s+Yu)',
        r'(张文涛|Zhang\s+Wentao|Wentao\s+Zhang)'
    ]
    authors_found = []
    for pattern in authors_patterns:
        if re.search(pattern, full_text, re.I):
            authors_found.append(pattern)
    
    if len(authors_found) >= 3:
        passes.append(f"✅ 找到 {len(authors_found)} 位作者")
    elif len(authors_found) > 0:
        warnings.append(f"⚠️ 作者可能不完整: 找到 {len(authors_found)} 位")
    else:
        warnings.append("⚠️ 未找到作者信息（可能使用英文名）")
    
    # 检查单位
    units = ['北京积水潭医院', '北京大学']
    units_found = [unit for unit in units if unit in full_text]
    if units_found:
        passes.append(f"✅ 单位信息: {', '.join(units_found)}")
    else:
        warnings.append("⚠️ 单位信息可能不完整")
    
    # 检查通讯作者邮箱
    email_patterns = [
        r'zhangwt2007@sina\.cn',
        r'fy\.yhy@163\.com'
    ]
    emails = re.findall(r'[\w\.-]+@[\w\.-]+\.\w+', full_text)
    if len(emails) >= 2:
        passes.append(f"✅ 通讯作者邮箱: {emails[:2]}")
    else:
        warnings.append(f"⚠️ 邮箱数量不足: {emails}")
    
    # 检查Running title
    running_title_found = 'running' in full_text.lower() or 'short title' in full_text.lower()
    if running_title_found:
        passes.append("✅ Running title存在")
    else:
        warnings.append("⚠️ 未明确标记Running title")
    
    # 检查关键词
    keywords_found = 'keyword' in full_text.lower() or 'keywords' in full_text.lower()
    if keywords_found:
        passes.append("✅ 关键词存在")
    else:
        warnings.append("⚠️ 未明确标记关键词")
    
    # 统计字数
    word_count = len(full_text.replace(' ', '').replace('\n', ''))
    passes.append(f"✅ 标题页字数: {word_count} 字符")
    
    return issues, warnings, passes, full_text

def check_main_text(doc):
    """检查主文稿"""
    print("\n" + "="*80)
    print("📄 主文稿检查 (Main_Text)")
    print("="*80)
    
    issues = []
    warnings = []
    passes = []
    
    # 提取所有文本
    full_text = ' '.join([para.text for para in doc.paragraphs])
    
    # 检查摘要
    abstract_match = re.search(r'abstract[:\s]*([^\n]*?)(?=\n\s*\n|\n\s*(introduction|keywords))', full_text, re.I)
    if abstract_match:
        abstract_text = abstract_match.group(1).strip()
        abstract_char_count = len(abstract_text.replace(' ', ''))
        passes.append(f"✅ 摘要存在 ({abstract_char_count} 字符)")
        
        # Nature建议摘要不超过250字符
        if abstract_char_count > 250:
            warnings.append(f"⚠️ 摘要过长: {abstract_char_count} 字符 (建议≤250字符)")
        else:
            passes.append(f"✅ 摘要长度符合建议 (≤250字符)")
    else:
        issues.append("❌ 未找到摘要")
    
    # 检查关键词
    keywords_match = re.search(r'keyword[s]?[:\s]*([^\n]*?)(?=\n\s*\n|\n\s*(introduction))', full_text, re.I)
    if keywords_match:
        keywords_text = keywords_match.group(1).strip()
        keyword_list = [k.strip() for k in re.split(r'[,;]', keywords_text)]
        passes.append(f"✅ 关键词存在: {len(keyword_list)} 个")
        if len(keyword_list) >= 3:
            passes.append(f"✅ 关键词数量符合要求 (≥3个)")
        else:
            warnings.append(f"⚠️ 关键词数量不足: {len(keyword_list)} 个 (建议≥3个)")
    else:
        warnings.append("⚠️ 未找到关键词")
    
    # 检查主要章节
    required_sections = {
        'introduction': '引言',
        'results': '结果',
        'discussion': '讨论',
        'methods': '方法'
    }
    found_sections = []
    for section_en, section_cn in required_sections.items():
        if re.search(r'\b' + section_en + r'\b', full_text, re.I):
            found_sections.append(section_cn)
    
    if len(found_sections) == 4:
        passes.append(f"✅ 所有主要章节存在: {', '.join(found_sections)}")
    else:
        missing = set(required_sections.values()) - set(found_sections)
        warnings.append(f"⚠️ 章节可能缺失: {missing}")
    
    # 检查参考文献引用
    ref_citations = re.findall(r'\[(\d+)\]', full_text)
    if ref_citations:
        max_ref = max(map(int, ref_citations))
        unique_refs = len(set(map(int, ref_citations)))
        passes.append(f"✅ 参考文献引用: {unique_refs} 篇, 最大编号 {max_ref}")
        if max_ref >= 50:
            passes.append(f"✅ 引用覆盖所有50篇参考文献")
        else:
            warnings.append(f"⚠️ 最大引用编号 {max_ref} 可能不足50篇")
    else:
        warnings.append("⚠️ 未找到参考文献引用格式")
    
    # 检查上标（化学式等）
    superscript_count = len([run for para in doc.paragraphs for run in para.runs if run.font.superscript])
    if superscript_count > 0:
        passes.append(f"✅ 上标格式: {superscript_count} 处")
    else:
        warnings.append("⚠️ 未检测到上标格式（可能手动检查）")
    
    # 检查希腊字母
    greek_letters = ['α', 'β', 'γ', 'δ', 'ε', 'θ', 'μ', 'π', 'σ', 'ω']
    greek_count = sum(full_text.count(letter) for letter in greek_letters)
    if greek_count > 0:
        passes.append(f"✅ 希腊字母: {greek_count} 个")
    else:
        warnings.append("⚠️ 未检测到希腊字母")
    
    # 统计总字数
    word_count = len(full_text.replace(' ', '').replace('\n', ''))
    passes.append(f"✅ 主文稿总字数: {word_count} 字符")
    
    return issues, warnings, passes, full_text

def check_references(doc):
    """检查参考文献"""
    print("\n" + "="*80)
    print("📚 参考文献检查 (References)")
    print("="*80)
    
    issues = []
    warnings = []
    passes = []
    
    # 提取所有文本
    full_text = '\n'.join([para.text for para in doc.paragraphs])
    
    # 检查编号格式
    numbered_refs = re.findall(r'^\s*(\d+)\.\s+', full_text, re.MULTILINE)
    ref_numbers = []
    if numbered_refs:
        ref_numbers = list(map(int, numbered_refs))
        passes.append(f"✅ 找到 {len(ref_numbers)} 篇参考文献（编号格式）")
        
        # 检查编号连续性
        if sorted(ref_numbers) == list(range(1, len(ref_numbers) + 1)):
            passes.append(f"✅ 参考文献编号连续: 1-{len(ref_numbers)}")
        else:
            issues.append(f"❌ 参考文献编号不连续: {ref_numbers}")
        
        # 检查数量
        if len(ref_numbers) == 50:
            passes.append(f"✅ 参考文献数量符合要求: 50篇")
        elif len(ref_numbers) > 50:
            warnings.append(f"⚠️ 参考文献数量超标: {len(ref_numbers)} 篇 (要求≤50篇)")
        else:
            warnings.append(f"⚠️ 参考文献数量不足: {len(ref_numbers)} 篇 (要求≥30篇)")
    else:
        warnings.append("⚠️ 未找到编号格式参考文献")
    
    # 检查DOI链接
    doi_links = re.findall(r'https?://doi\.org/[^,\s\)]+', full_text)
    doi_count = len(doi_links)
    if doi_count > 0:
        passes.append(f"✅ DOI链接: {doi_count} 个")
        if ref_numbers:
            if doi_count >= 40:
                passes.append(f"✅ DOI覆盖率足够: {doi_count}/{len(ref_numbers)}")
            else:
                warnings.append(f"⚠️ DOI覆盖率不足: {doi_count}/{len(ref_numbers)}")
        else:
            passes.append(f"✅ DOI链接存在")
    else:
        warnings.append("⚠️ 未找到DOI链接")
    
    # 检查格式一致性（Vancouver格式）
    vancouver_pattern = r'^\s*\d+\.\s+[A-Z][^,]+,\s+[^,]+,\s+et\s+al\.'  # 基本格式
    vancouver_refs = re.findall(vancouver_pattern, full_text, re.MULTILINE)
    if ref_numbers and len(vancouver_refs) >= len(ref_numbers) * 0.8:  # 80%匹配
        passes.append(f"✅ 格式符合Vancouver标准")
    elif ref_numbers:
        warnings.append(f"⚠️ 部分参考文献格式可能不符合Vancouver标准")
    
    # 检查期刊名称
    journal_match = re.search(r'[A-Z][a-z]{3,20}\s+\d{4}', full_text)
    if journal_match:
        passes.append(f"✅ 包含期刊名称和年份")
    else:
        warnings.append(f"⚠️ 部分参考文献可能缺少期刊信息")
    
    # 统计字数
    word_count = len(full_text.replace(' ', '').replace('\n', ''))
    passes.append(f"✅ 参考文献总字数: {word_count} 字符")
    
    return issues, warnings, passes, full_text

def main():
    print("="*80)
    print("🔍 Nature期刊投稿文件全面审核")
    print("="*80)
    print("审核时间:", "2026-04-04")
    print()
    
    # 文件路径
    files = {
        'Title_Page': r'C:\Users\89367\Desktop\Nature_Submission_Meniscus_20260404\Manuscript\Title_Page_new.docx',
        'Main_Text': r'C:\Users\89367\Desktop\Nature_Submission_Meniscus_20260404\Manuscript\Main_Text_new.docx',
        'References': r'C:\Users\89367\Desktop\Nature_Submission_Meniscus_20260404\Manuscript\References_new.docx'
    }
    
    all_issues = []
    all_warnings = []
    all_passes = []
    
    # 检查每个文件
    for file_name, file_path in files.items():
        try:
            doc = Document(file_path)
            
            if file_name == 'Title_Page':
                issues, warnings, passes, _ = check_title_page(doc)
            elif file_name == 'Main_Text':
                issues, warnings, passes, _ = check_main_text(doc)
            elif file_name == 'References':
                issues, warnings, passes, _ = check_references(doc)
            
            all_issues.extend([f"[{file_name}] {issue}" for issue in issues])
            all_warnings.extend([f"[{file_name}] {warning}" for warning in warnings])
            all_passes.extend([f"[{file_name}] {pass_item}" for pass_item in passes])
            
        except Exception as e:
            all_issues.append(f"[{file_name}] ❌ 文件读取失败: {str(e)}")
    
    # 总结
    print("\n" + "="*80)
    print("📊 审核总结")
    print("="*80)
    
    print(f"\n✅ 通过项 ({len(all_passes)}):")
    for item in all_passes:
        print(f"  {item}")
    
    if all_warnings:
        print(f"\n⚠️ 警告项 ({len(all_warnings)}):")
        for item in all_warnings:
            print(f"  {item}")
    
    if all_issues:
        print(f"\n❌ 问题项 ({len(all_issues)}):")
        for item in all_issues:
            print(f"  {item}")
    
    print("\n" + "="*80)
    print("📈 投稿准备度评估")
    print("="*80)
    
    # 计算准备度
    total_checks = len(all_passes) + len(all_issues) * 5
    pass_score = len(all_passes)
    issue_penalty = len(all_issues) * 5
    
    if total_checks > 0:
        readiness = min(100, max(0, int((pass_score - issue_penalty) / total_checks * 100)))
    else:
        readiness = 0
    
    print(f"\n总体准备度: {readiness}%")
    
    if len(all_issues) == 0:
        print("\n🎉 没有严重问题，可以投稿！")
    elif len(all_issues) <= 2:
        print("\n⚠️ 存在少量问题，建议修复后投稿")
    else:
        print("\n❌ 存在多个问题，必须修复后再投稿")
    
    print("\n" + "="*80)

if __name__ == '__main__':
    main()
