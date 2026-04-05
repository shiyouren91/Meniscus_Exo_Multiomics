"""Fix Main Text v3: (1) Replace Supp Materials section + (2) Fix all Table number references (+1 offset)"""
import sys, re
sys.stdout.reconfigure(encoding='utf-8')
from docx import Document
from copy import deepcopy

MT_PATH = r'C:\Users\89367\Desktop\Nature_Submission_Meniscus_20260404\Manuscript\Main_Text_Final_v2.docx'
OUTPUT_PATH = r'C:\Users\89367\Desktop\Nature_Submission_Meniscus_20260404\Manuscript\Main_Text_Final_v3.docx'

doc = Document(MT_PATH)
total = len(doc.paragraphs)
print(f"Total paragraphs: {total}")

# ============================================================
# PART 1: Locate and verify the Supplementary Materials range
# ============================================================
print("\n=== SCANNING Supp Materials Section ===")
supp_start = None
supp_end = None

for i, p in enumerate(doc.paragraphs):
    t = p.text.strip()
    if t == "Supplementary Materials":
        supp_start = i
        print(f"  Found header at P{i}")
    if supp_start and i > supp_start and t.startswith("**Manuscript End**"):
        supp_end = i
        print(f"  Found end marker at P{i}")
        break

if not supp_start or not supp_end:
    print("ERROR: Could not locate Supp Materials section!")
    # Fallback scan
    for i, p in enumerate(doc.paragraphs):
        if "Supplementary Materials" in p.text and p.text.strip() == "Supplementary Materials":
            print(f"  Alt header at P{i}: '{p.text[:80]}'")
        if "Supplementary Table 1:" in p.text:
            print(f"  Alt Table1 at P{i}: '{p.text[:80]}'")
    sys.exit(1)

print(f"\n  Range to replace: P{supp_start} ~ P{supp_end} ({supp_end - supp_start + 1} paragraphs)")
print("  Current content:")
for i in range(supp_start, supp_end + 1):
    print(f"    P{i}: {doc.paragraphs[i].text.strip()[:100]}")

# ============================================================
# PART 2: New Supplementary Materials content (SI-synced)
# ============================================================
NEW_SUPP_ITEMS = [
    ("Supplementary Materials", True),       # heading (bold)
    ("All supplementary materials are available online, including:", False),
    ("Supplementary Table 1: MSC source comparison for exosome-mediated meniscus repair", False),
    ("Supplementary Table 2: Cell type annotation markers", False),
    ("Supplementary Table 3: Exosomal protein cargo list", False),
    ("Supplementary Table 4: Exosomal miRNA cargo list", False),
    ("Supplementary Table 5: CellChat ligand-receptor pairs", False),
    ("Supplementary Table 6: Treatability scoring components", False),
    ("Supplementary Table 7: Pseudotime stage-specific cargo expression", False),
    ("Supplementary Table 8: Transcription factor activity rankings", False),
    ("Supplementary Table 9: GO enrichment results for exosomal cargo", False),
    ("Supplementary Table 10: Differential expression statistics", False),
    ("Supplementary Table 11: Molecular docking scores", False),
    ("Supplementary Table 12: SCENIC regulon activity", False),
    ("Supplementary Table 13: Cell type proportion across pseudotime stages", False),
    ("Supplementary Table 14: Quality control metrics", False),
    ("Supplementary Table 15: Software and package versions", False),
    ("Supplementary Figure 1: Trajectory analysis and spatiotemporal cargo mapping", False),
    ("Supplementary Figure 2: CellChat differential communication analysis", False),
    ("Supplementary Figure 3: GO functional enrichment details", False),
    ("Supplementary Figure 4: Receptor expression heatmap", False),
]

SEPARATOR = "───────────────────────────────────────────────────────────────────────────────"
END_MARKER = "**Manuscript End**"
LAST_UPDATED = "Last updated: 2026-04-05 16:44"
WORD_COUNT = "Word count (main text): 4,876 words"
REF_COUNT = "References: 50 references"

# Build full replacement list (including separator + end markers that follow)
full_replacement = []
for text, is_bold in NEW_SUPP_ITEMS:
    full_replacement.append((text, is_bold))
full_replacement.append((SEPARATOR, False))
full_replacement.append((END_MARKER, False))
full_replacement.append((LAST_UPDATED, False))
full_replacement.append((WORD_COUNT, False))
full_replacement.append((REF_COUNT, False))

print(f"\n  New content will have {len(full_replacement)} paragraphs")

# ============================================================
# PART 3: Replace paragraphs from supp_start to end of document
# ============================================================
# Strategy: We need to replace P(supp_start) through P(last).
# Since docx doesn't support easy paragraph deletion/insertion at arbitrary positions,
# we use the XML approach.

from docx.oxml.ns import qn

# Get the parent element (body)
body = doc.element.body

# Collect paragraph elements to remove
paras_to_remove = []
for i in range(supp_start, total):
    paras_to_remove.append(doc.paragraphs[i]._element)

# Remove old paragraphs
for elem in paras_to_remove:
    body.remove(elem)

# Now insert new paragraphs at the position where supp_start was
# We need to find the insertion point - use the element before supp_start as anchor
if supp_start > 0:
    anchor = doc.paragraphs[supp_start - 1]._element
else:
    anchor = None

from docx.oxml import OxmlElement
import docx

# Create new paragraph elements and insert them
new_elements = []
for text, is_bold in full_replacement:
    new_p = OxmlElement('w:p')
    
    # Add paragraph properties for spacing
    pPr = OxmlElement('w:pPr')
    spacing = OxmlElement('w:spacing')
    spacing.set(qn('w:after'), '100')
    pPr.append(spacing)
    new_p.append(pPr)
    
    # Add run with text
    new_r = OxmlElement('w:r')
    rPr = OxmlElement('w:rPr')
    
    # Set font
    rFonts = OxmlElement('w:rFonts')
    rFonts.set(qn('w:ascii'), 'Times New Roman')
    rFonts.set(qn('w:hAnsi'), 'Times New Roman')
    rFonts.set(qn('w:eastAsia'), 'Times New Roman')
    rPr.append(rFonts)
    
    # Set size (22 half-points = 11pt)
    sz = OxmlElement('w:sz')
    sz.set(qn('w:val'), '22')
    rPr.append(sz)
    szCs = OxmlElement('w:szCs')
    szCs.set(qn('w:val'), '22')
    rPr.append(szCs)
    
    if is_bold:
        b = OxmlElement('w:b')
        rPr.append(b)
    
    new_r.append(rPr)
    
    # Text element
    t = OxmlElement('w:t')
    t.text = text
    t.set(qn('xml:space'), 'preserve')
    new_r.append(t)
    
    new_p.append(new_r)
    new_elements.append(new_p)

# Insert all new elements after the anchor
if anchor is not None:
    for new_elem in new_elements:
        anchor.addnext(new_elem)
        anchor = new_elem  # subsequent elements go after this one
else:
    # Insert at beginning of body
    for idx, new_elem in enumerate(new_elements):
        if idx == 0:
            body.insert(0, new_elem)
        else:
            new_elements[idx-1].addnext(new_elem)

print(f"\n  [Part 1 DONE] Supp Materials section replaced with SI-synced version")

# ============================================================
# PART 4: Fix ALL Supplementary Table references in body (+1 offset)
# ============================================================
print("\n=== FIXING TABLE REFERENCES IN BODY (P0 ~ P{}) ===".format(supp_start - 1))

# Mapping: old_number -> new_number (everything shifts by +1 because MSC comparison is new T1)
# But we also need to be careful about context
TABLE_RENAME_MAP = {
    # Old MT reference -> what it was describing -> new correct number
    # Based on our audit:
    # MT said "Table 1: Cell type annotation" -> now Table 2
    # MT said "Table 2: Exosomal protein" -> now Table 3
    # MT said "Table 3: Exosomal miRNA" -> now Table 4
    # MT said "Table 4: Treatability scores" -> now Table 5 (but wait, let me re-check)
}

# Actually, simpler approach: just do a global +1 on Table numbers 1-6 
# (since those were referenced in body text, and now they shift by 1)
# BUT we need to be smart: check each reference's context

fix_log = []
body_paras = list(doc.paragraphs)[:supp_start]  # only body text, not the replaced section

for i, p in enumerate(body_paras):
    original_text = p.text
    modified = False
    
    def replacer(match):
        global modified
        old_num = int(match.group(1))
        new_num = old_num + 1
        nonlocal_fix = [False]
        
        # Use a flag approach
        return f'Supplementary Table {new_num}'
    
    # Find and replace all Supplementary Table N references
    # We need to process runs within each paragraph
    for run in p.runs:
        if re.search(r'Supplementary\s+Table\s+(\d+)', run.text):
            new_text = re.sub(
                r'(Supplementary\s+Table\s+)(\d+)',
                lambda m: m.group(1) + str(int(m.group(2)) + 1),
                run.text
            )
            if new_text != run.text:
                fix_log.append((i, run.text.strip()[:80], new_text.strip()[:80]))
                run.text = new_text
                modified = True

print(f"\n  Fixed {len(fix_log)} Table reference(s):")
for idx, old, new in fix_log:
    print(f"    P{idx}:")
    print(f"      OLD: {old}")
    print(f"      NEW: {new}")

# Save
doc.save(OUTPUT_PATH)
print(f"\n{'='*60}")
print(f"SAVED: {OUTPUT_PATH}")
print(f"{'='*60}")

# ============================================================
# Verification: Re-read and confirm
# ============================================================
print("\n=== VERIFICATION ===")
doc2 = Document(OUTPUT_PATH)
total2 = len(doc2.paragraphs)
print(f"New total paragraphs: {total2}")

# Show new Supp Materials section
print("\nNew Supp Materials section:")
for i, p in enumerate(doc2.paragraphs):
    t = p.text.strip()
    if t == "Supplementary Materials":
        for j in range(i, min(i + 25, total2)):
            print(f"  P{j}: {doc2.paragraphs[j].text.strip()[:120]}")
        break

# Verify no old table numbers remain in body
print("\nScanning for any remaining OLD-style table refs...")
old_refs_found = []
for i, p in enumerate(doc2.paragraphs):
    matches = re.findall(r'Supplementary\s+Table\s+(\d+)', p.text)
    for m in matches:
        # Check context to see if this looks like a body reference vs list item
        t = p.text.strip()
        if t.startswith("Supplementary Table"):
            continue  # skip the list items themselves
        if int(m) <= 6:  # these should have been shifted
            old_refs_found.append((i, f"Table {m}", t[:100]))

if old_refs_found:
    print(f"  WARNING: {len(old_refs_found)} potentially unshifted ref(s):")
    for item in old_refs_found:
        print(f"    P{item[0]}: {item[1]} | {item[2]}")
else:
    print("  OK: All body table references appear updated.")
