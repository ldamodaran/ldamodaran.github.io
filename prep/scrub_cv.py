import sys, re, fitz

src, dst = sys.argv[1], sys.argv[2]
doc = fitz.open(src)

# 1. drop the references page (other people's phone numbers and emails)
drop = [i for i in range(doc.page_count)
        if "Professional References" in doc[i].get_text()]
print("references page(s):", [i + 1 for i in drop])

# 2. redact personal gmail and phone number on remaining pages
PATTERNS = [
    re.compile(r"[Ll]ambodhar1@gmail\.com"),
    re.compile(r"Phone:\s*\(?919\)?[\s\-–—]*672[\s\-–—]*6182"),
    re.compile(r"\(?919\)?[\s\-–—]*672[\s\-–—]*6182"),
]
hits = 0
for page in doc:
    for pat in PATTERNS:
        for m in re.finditer(pat, page.get_text()):
            for r in page.search_for(m.group(0)):
                page.add_redact_annot(r, text="", fill=(1, 1, 1))
                hits += 1
    # strip mailto/tel link annotations pointing at the personal address
    for link in page.get_links():
        u = (link.get("uri") or "").lower()
        if "gmail" in u or "6182" in u or u.startswith("tel:"):
            page.delete_link(link)
    page.apply_redactions()
print("redactions applied:", hits)

for i in sorted(drop, reverse=True):
    doc.delete_page(i)

doc.save(dst, garbage=4, deflate=True, clean=True)
doc.close()

# 3. verify the output really is clean
out = fitz.open(dst)
txt = "\n".join(p.get_text() for p in out)
print("pages:", out.page_count)
for probe in ("gmail", "672", "6182", "Professional References",
              "JUSTIN.BAHL", "lhmoncla", "katia.koelle"):
    print(f"  {probe:24s} {'STILL PRESENT' if probe.lower() in txt.lower() else 'removed'}")
print("  emory address kept:", "ldamoda@emory.edu" in txt)
