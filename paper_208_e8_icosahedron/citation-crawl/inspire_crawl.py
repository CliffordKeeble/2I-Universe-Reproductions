"""
Paper 208 §8 item 2 — Inspire-HEP citation crawl (the spine).

Forward-citation extraction for the two Dechant anchors, community tagging,
and the four cross-citation claim tests. Inspire-HEP only; ADS requires a
Bearer token (HTTP 401 without one) and is run separately if a token is
supplied — see crosscite_findings.md for the database-coverage statement.

No interpretation past the data (hard gate 1): bridging papers are reported,
not adjudicated. Run: python inspire_crawl.py
"""
import json, time, csv, os, urllib.parse, urllib.request

OUT = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(OUT, "raw")
os.makedirs(RAW, exist_ok=True)

API = "https://inspirehep.net/api/literature"
FIELDS = ("titles,arxiv_eprints,dois,control_number,publication_info,"
          "earliest_date,authors.full_name,number_of_pages")

def fetch(q, size=200, fields=FIELDS, tag=None):
    url = f"{API}?q={urllib.parse.quote(q)}&size={size}&fields={fields}"
    req = urllib.request.Request(url, headers={"Accept": "application/json",
                                               "User-Agent": "paper208-crawl/1.0"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                data = json.loads(r.read().decode())
            if tag:
                with open(os.path.join(RAW, f"{tag}.json"), "w", encoding="utf-8") as f:
                    json.dump(data, f, indent=1)
            return data
        except Exception as e:
            print(f"  retry {attempt+1} for q={q!r}: {e}")
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"failed: {q}")

# ---- anchors ----------------------------------------------------------------
ANCHORS = {
    "Dechant2016": "arxiv:1602.05985",   # [4]  E8 out of spinors of icosahedron, RSPA 472:20150504
    "Dechant2017": "arxiv:1603.04805",   # [31] E8 geometry from Clifford perspective, AACA 27
    "CHL1995":     "arxiv:hep-th/9506048",  # [22] Chaudhuri-Polchinski (CHL)
    "KN1990":      "Kronheimer Nakajima Yang-Mills instantons ALE",  # [36] no arXiv (1990)
}

def recid_of(data):
    h = data["hits"]["hits"]
    return h[0]["metadata"]["control_number"] if h else None

print("=== resolving anchor recids ===")
recids = {}
for name, q in ANCHORS.items():
    d = fetch(q, size=3, tag=f"anchor_{name}")
    rid = recid_of(d)
    recids[name] = rid
    cc = d["hits"]["hits"][0]["metadata"].get("citation_count") if d["hits"]["hits"] else "?"
    title = d["hits"]["hits"][0]["metadata"]["titles"][0]["title"] if d["hits"]["hits"] else "?"
    print(f"  {name:12s} recid={rid}  cites={cc}  {title[:50]}")
    time.sleep(1)

# ---- community tagging ------------------------------------------------------
# P=polytope/Clifford-spinor  A=ALE/ADHM/instanton  S=heterotic string
# F=field-theory twin/mirror  E=emergent-E8/integrable  U=EJA/octonion
KW = {
    "P": ["h4", "600-cell", "polytope", "coxeter", "clifford", "spinor",
          "icosahed", "root system", "quaternion", "reflection group"],
    "A": ["ale ", "kronheimer", "nakajima", "adhm", "instanton", "orbifold",
          "resolution", "mckay", "du val", "gravitational instanton", "hyperkahler"],
    "S": ["heterotic", "e8 x e8", "e8xe8", "anomaly cancel", "string",
          "compactification", "m-theory", "supergravity"],
    "F": ["twin higgs", "mirror matter", "shadow matter", "mirror sector",
          "twin sector", "dark sector", "hidden sector", "particle pairing"],
    "E": ["emergent", "zamolodchikov", "coldea", "conb2o6", "cobalt niobate",
          "integrable", "ising", "quantum criticality", "e8 spectrum", "golden ratio mass"],
    "U": ["jordan algebra", "exceptional jordan", "octonion", "furey", "gresnigt",
          "boyle", "singh", "division algebra", "standard model from", "trace dynamics"],
}

def tag_communities(title, cats):
    t = (title or "").lower()
    c = " ".join(cats or []).lower()
    blob = t + " || " + c
    tags = set()
    for comm, kws in KW.items():
        if any(k in blob for k in kws):
            tags.add(comm)
    # category-driven hints
    if "math.rt" in c or "math-ph" in c or "math.mp" in c or "math.gr" in c:
        tags.add("P")  # representation/Clifford math community
    return "".join(sorted(tags))

def row_from_hit(h, citing_anchor):
    m = h["metadata"]
    arx = m.get("arxiv_eprints", [{}])
    arxid = arx[0].get("value", "") if arx else ""
    cats = arx[0].get("categories", []) if arx else []
    title = m["titles"][0]["title"] if m.get("titles") else ""
    doi = m.get("dois", [{}])
    doi = doi[0].get("value", "") if doi else ""
    pubs = m.get("publication_info", [])
    venue = ""
    for p in pubs:
        if p.get("journal_title"):
            venue = p["journal_title"]
            break
    year = (m.get("earliest_date", "") or "")[:4]
    auths = [a.get("full_name", "") for a in m.get("authors", [])]
    first = auths[0] if auths else ""
    return {
        "recid": m["control_number"],
        "arxiv": arxid,
        "primary_cat": cats[0] if cats else "",
        "all_cats": "|".join(cats),
        "title": title,
        "first_author": first,
        "venue": venue,
        "year": year,
        "doi": doi,
        "communities": tag_communities(title, cats),
        "cites_anchor": citing_anchor,
    }

# ---- forward citers of the two Dechant anchors ------------------------------
print("\n=== forward citers (the spine) ===")
spine = {}   # recid -> row
for name in ("Dechant2016", "Dechant2017"):
    rid = recids[name]
    d = fetch(f"refersto:recid:{rid}", size=200, tag=f"citers_{name}")
    hits = d["hits"]["hits"]
    print(f"  {name}: {d['hits']['total']} citers")
    for h in hits:
        r = row_from_hit(h, name)
        key = r["recid"]
        if key in spine:
            spine[key]["cites_anchor"] += "+" + name
        else:
            spine[key] = r
    time.sleep(1)

rows = sorted(spine.values(), key=lambda r: (r["year"], r["recid"]))
print(f"  deduped union: {len(rows)} distinct citers")

cols = ["recid","arxiv","primary_cat","all_cats","title","first_author",
        "venue","year","doi","communities","cites_anchor"]
with open(os.path.join(OUT, "inspire_ads_citers_dechant.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    w.writerows(rows)
print(f"  wrote inspire_ads_citers_dechant.csv ({len(rows)} rows)")

# ---- community histogram ----------------------------------------------------
from collections import Counter
hist = Counter()
for r in rows:
    for ch in (r["communities"] or "_"):
        hist[ch] += 1
print("  community spread:", dict(hist))

# ---- cross-citation claim tests --------------------------------------------
# Method: anchor each "side" of a pair to canonical Inspire recid(s), pull the
# full citing set of each side, and intersect. A paper in the intersection
# cites BOTH sides; whether it does so FOR THE QUALIFIED PURPOSE is a Scout
# read, not ours (hard gate 1). Empty intersection = claim holds at this depth.
print("\n=== cross-citation tests ===")

# canonical side-anchors (resolved 12 Jun 2026, see raw/ and brief)
SIDE = {
    "KN_ALE":        [415625],            # A: Kronheimer-Nakajima 1990 (ALE/ADHM)
    "heterotic_E8":  [213013, 207031],    # S: Gross-Harvey-Martinec-Rohm 1985 + Heterotic String
    "twin_Higgs":    [685922, 699848],    # F: Chacko-Goh-Harnik + LR twin Higgs
    "emergent_E8":   [842444],            # E: Coldea et al CoNb2O6 (Ising-chain E8)
    "EJA_octonion":  [2095962, 2667356, 1658866, 555166],  # U: E8xE8/octonion-SM, Furey, EJA-matrix-string
}

def citers_recids(rids, tag=None):
    """union of citing-records over a list of side-anchor recids -> {recid: hit}"""
    out = {}
    for rid in rids:
        if not rid:
            continue
        d = fetch(f"refersto:recid:{rid}", size=1000, tag=(f"{tag}_{rid}" if tag else None))
        for h in d["hits"]["hits"]:
            out[h["metadata"]["control_number"]] = h
    return out

def title_of(h):
    m = h["metadata"]
    return m["titles"][0]["title"] if m.get("titles") else str(m["control_number"])

def short_cats(h):
    ar = h["metadata"].get("arxiv_eprints", [{}])
    return "|".join(ar[0].get("categories", [])) if ar else ""

def intersect_report(left, right, label, cap=40):
    inter = set(left) & set(right)
    print(f"   |{label.split('INTERSECT')[0].strip()}|={len(left)}  "
          f"|other|={len(right)}  intersection={len(inter)}")
    rows_out = []
    for rid in sorted(inter):
        h = right.get(rid) or left.get(rid)
        rows_out.append((rid, title_of(h), short_cats(h)))
    for rid, t, c in rows_out[:cap]:
        print(f"     recid={rid} [{c[:24]}] {t[:60]}")
    if len(rows_out) > cap:
        print(f"     ... +{len(rows_out)-cap} more (see claim_findings.json)")
    if not inter:
        print("   EMPTY")
    return rows_out

# §3.4 P<->A : Dechant-citers INTERSECT Kronheimer-Nakajima-citers
print("\n[§3.4 P<->A]  Dechant-citers INTERSECT ALE/Kronheimer-Nakajima-citers")
print("   (qualified purpose: S3/2I orbit-space / Poincare-homology-sphere reading;")
print("    Allen-Sutcliffe 1302.4664 & Choi-Lee are KNOWN instanton bridges, excluded)")
kn = citers_recids(SIDE["KN_ALE"], tag="side_KN")
findings_34 = intersect_report(spine, kn, "Dechant INTERSECT ALE")

# §4.4 S<->F : heterotic-E8xE8 INTERSECT twin/mirror-matter
print("\n[§4.4 S<->F]  heterotic-E8xE8-citers INTERSECT twin-Higgs/mirror-matter-citers")
print("   (qualified purpose: two-E8-factor EXCHANGE read as particle pairing,")
print("    not lattice/moduli device; purpose-call is Scout's per gate 1)")
het = citers_recids(SIDE["heterotic_E8"], tag="side_het")
twin = citers_recids(SIDE["twin_Higgs"], tag="side_twin")
findings_44 = intersect_report(het, twin, "heterotic INTERSECT twin")
hits44 = [(r[0], r[1]) for r in findings_44]

# §4.6 E<->U : emergent-E8 (Coldea) INTERSECT EJA/octonion
print("\n[§4.6 E<->U]  emergent-E8(Coldea)-citers INTERSECT EJA/octonion-citers")
print("   (qualified purpose: SUBSTANTIVE interaction emergent-E8 <-> E8-as-fundamental,")
print("    not a passing citation; purpose-call is Scout's per gate 1)")
emer = citers_recids(SIDE["emergent_E8"], tag="side_emer")
eja = citers_recids(SIDE["EJA_octonion"], tag="side_eja")
findings_46 = intersect_report(emer, eja, "emergent INTERSECT EJA")
hits46 = [(r[0], r[1]) for r in findings_46]

# §5.2 internal : single source unifying the five §5.2 involutions
print("\n[§5.2 internal]  NOT testable from citation metadata alone")
print("   The five involutions are internal to Paper 208 §5.2; deciding whether one")
print("   source unifies them needs the involution list + a body read (Scout/Dig).")
print("   DEFERRED — out of scope for the citation crawl (declared).")

# §4.8 Dechant-forward : any hep-ph/hep-th citer attaching a physical parameter
print("\n[§4.8]  hep-ph/hep-th Dechant citers (physical-parameter attachers?)")
phys = [r for r in rows if r["primary_cat"] in ("hep-ph", "hep-th")
        or "hep-ph" in r["all_cats"] or "hep-th" in r["all_cats"]]
for r in phys:
    print(f"   recid={r['recid']} cat={r['primary_cat']} [{r['communities']}] {r['title'][:55]}")
print(f"   ({len(phys)} citers carry a hep-ph/hep-th tag; "
      "physical-parameter attachment requires Scout read of each)")

# dump claim findings for the markdown step
findings = {
    "recids": recids,
    "side_anchors": SIDE,
    "n_citers_union": len(rows),
    "community_hist": dict(hist),
    "claim_34_intersection": [list(x) for x in findings_34],
    "claim_44_intersection": [list(x) for x in findings_44],
    "claim_46_intersection": [list(x) for x in findings_46],
    "claim_48_physcat_citers": [(r["recid"], r["arxiv"], r["primary_cat"],
                                 r["title"], r["communities"]) for r in phys],
}
with open(os.path.join(OUT, "claim_findings.json"), "w", encoding="utf-8") as f:
    json.dump(findings, f, indent=2)
print("\nwrote claim_findings.json")
print("DONE")
