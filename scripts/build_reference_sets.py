#!/usr/bin/env python3
"""
Build zebra finch reference sets for ShinyGO and mitch.

Approach developed with guidance from Fiona McCarthy (Oct 2026):
  * New reference genome = bTaeGut7.mat (GCA_048771995.1 / RefSeq GCF_048771995.1)
  * Parse the GFF3s (Ensembl + NCBI) for Ensembl gene ID, NCBI GeneID, symbol, description
  * BioMart is being retired -> rebuild orthology from Ensembl FTP homology files
  * Keep only "strict" orthologs (1:1) for cross-species work

Run from the "Gene Comparisons" folder:   python3 scripts/build_reference_sets.py
Inputs : reference_sources/   (raw public downloads, see README in reference_sets/)
Outputs: reference_sets/
"""
import gzip, hashlib, os, re, sys, urllib.parse
from collections import defaultdict, Counter
import pandas as pd

SRC, OUT = "reference_sources", "reference_sets"
os.makedirs(OUT, exist_ok=True)
GMT_MIN, GMT_MAX = 5, 1000
log_lines = []
def log(*a):
    s = " ".join(str(x) for x in a); print(s); log_lines.append(s)

def op(p):
    p = os.path.join(SRC, p)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def attrs(s):
    d = {}
    for kv in s.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1); d[k] = urllib.parse.unquote(v)
    return d

def clean_desc(d):
    return re.sub(r"\s*\[(Ensembl NN prediction|Source:)[^\]]*\]", "", d or "").strip()

# ---------------------------------------------------------------- Ensembl GFF3
def ens_genes(fn):
    rows = []
    with op(fn) as f:
        for l in f:
            if l.startswith("#"): continue
            c = l.rstrip("\n").split("\t")
            if len(c) < 9 or not c[8].startswith("ID=gene:"): continue
            a = attrs(c[8])
            rows.append(dict(ensembl_gene_id=a["gene_id"], symbol=a.get("Name", ""),
                             biotype=a.get("biotype", ""), description=clean_desc(a.get("description", "")),
                             chrom=c[0], start=int(c[3]), end=int(c[4]), strand=c[6]))
    return pd.DataFrame(rows)

new = ens_genes("ENS2026_bTaeGut7mat_genes.gff3.gz")
old = ens_genes("ENS104_bTaeGut1v1_genes.gff3.gz")
log(f"Ensembl 2026_05 bTaeGut7.mat genes: {len(new)} ({(new.biotype=='protein_coding').sum()} protein-coding, {(new.symbol!='').sum()} with symbol)")
log(f"Ensembl 104 bTaeGut1_v1 genes:      {len(old)} ({(old.biotype=='protein_coding').sum()} protein-coding, {(old.symbol!='').sum()} with symbol)")

# ---------------------------------------------------------------- Ensembl xrefs
def ens_xref(fn):
    hgnc, go = {}, defaultdict(set)
    with op(fn) as f:
        next(f)
        for l in f:
            c = l.rstrip("\n").split("\t")
            if len(c) < 9: continue
            if c[6] == "HGNC Symbol": hgnc.setdefault(c[0], c[4])
            elif c[6] == "GO" and c[3].startswith("GO:"): go[c[0]].add(c[3])
    return hgnc, go
new_hgnc, new_go_ens = ens_xref("ENS2026_bTaeGut7mat_xref.tsv.gz")
old_hgnc, old_go_ens = ens_xref("ENS104_bTaeGut1v1_xref.tsv.gz")
new["hgnc_symbol_xref"] = new.ensembl_gene_id.map(new_hgnc).fillna("")
old["hgnc_symbol_xref"] = old.ensembl_gene_id.map(old_hgnc).fillna("")

# ---------------------------------------------------------------- new-assembly homology (Ensembl 2026_05, BLAST RBBH)
hom = pd.read_csv(os.path.join(SRC, "ENS2026_bTaeGut7mat_homology.tsv.gz"), sep="\t", dtype=str)
hom = hom[hom.ref_species.isin(["homo_sapiens", "mus_musculus", "gallus_gallus"])]
for sp, tag in [("homo_sapiens", "human"), ("mus_musculus", "mouse"), ("gallus_gallus", "chicken")]:
    h = hom[(hom.ref_species == sp) & (hom.homology_type == "homolog_rbbh")]
    # an RBBH should be 1:1; drop any finch gene or ref gene that occurs more than once
    h = h[~h.query_gene_stable_id.duplicated(keep=False) & ~h.ref_gene_stable_id.duplicated(keep=False)]
    m = h.set_index("query_gene_stable_id")
    new[f"{tag}_ensembl_id"] = new.ensembl_gene_id.map(m.ref_gene_stable_id).fillna("")
    new[f"{tag}_symbol"] = new.ensembl_gene_id.map(m.ref_gene_name).fillna("")
    new[f"{tag}_perc_id"] = new.ensembl_gene_id.map(m.query_perc_id).fillna("")
    log(f"  bTaeGut7 strict (RBBH 1:1) {tag} orthologs: {(new[f'{tag}_ensembl_id']!='').sum()}")

# ---------------------------------------------------------------- NCBI RefSeq genes on same assembly
asm = {}
with op("NCBI_GCF_048771995.1_bTaeGut7.mat_assembly_report.txt") as f:
    for l in f:
        if l.startswith("#"): continue
        c = l.rstrip("\n").split("\t")
        asm[c[6]] = c[2]          # RefSeq accn -> assigned molecule (Ensembl seq name)
ncbi = []
with op("NCBI_GCF_048771995.1_bTaeGut7.mat_genomic.gff.gz") as f:
    for l in f:
        if l.startswith("#"): continue
        c = l.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] not in ("gene", "pseudogene"): continue
        a = attrs(c[8])
        gid = re.search(r"GeneID:(\d+)", a.get("Dbxref", ""))
        if not gid: continue
        ncbi.append(dict(ncbi_gene_id=gid.group(1), ncbi_symbol=a.get("Name", ""),
                         ncbi_description=a.get("description", ""), ncbi_biotype=a.get("gene_biotype", ""),
                         chrom=asm.get(c[0], c[0]), start=int(c[3]), end=int(c[4]), strand=c[6]))
ncbi = pd.DataFrame(ncbi).drop_duplicates("ncbi_gene_id")
log(f"NCBI RefSeq GCF_048771995.1 (RS_2025_04) genes: {len(ncbi)}")

# coordinate matching (same assembly, same strand): best reciprocal overlap (Jaccard)
BIN = 100_000
idx = defaultdict(list)
for r in ncbi.itertuples():
    for b in range(r.start // BIN, r.end // BIN + 1):
        idx[(r.chrom, r.strand, b)].append(r)
best = {}
for r in new.itertuples():
    seen, cand = set(), None
    for b in range(r.start // BIN, r.end // BIN + 1):
        for q in idx.get((r.chrom, r.strand, b), []):
            if q.ncbi_gene_id in seen: continue
            seen.add(q.ncbi_gene_id)
            ov = min(r.end, q.end) - max(r.start, q.start) + 1
            if ov <= 0: continue
            j = ov / (max(r.end, q.end) - min(r.start, q.start) + 1)
            sym_ok = bool(r.symbol) and r.symbol.upper() == q.ncbi_symbol.upper()
            score = (sym_ok, j)
            if cand is None or score > cand[0]: cand = (score, q, j, sym_ok)
    if cand and (cand[3] or cand[2] >= 0.3):
        best[r.ensembl_gene_id] = cand
new["ncbi_gene_id"] = new.ensembl_gene_id.map(lambda g: best[g][1].ncbi_gene_id if g in best else "")
new["ncbi_symbol"] = new.ensembl_gene_id.map(lambda g: best[g][1].ncbi_symbol if g in best else "")
new["ncbi_description"] = new.ensembl_gene_id.map(lambda g: best[g][1].ncbi_description if g in best else "")
new["ncbi_overlap_jaccard"] = new.ensembl_gene_id.map(lambda g: round(best[g][2], 3) if g in best else "")
new["ncbi_symbol_agrees"] = new.ensembl_gene_id.map(lambda g: ("yes" if best[g][3] else "no") if g in best else "")
log(f"  Ensembl bTaeGut7 genes matched to an NCBI GeneID by coordinates: {(new.ncbi_gene_id!='').sum()} "
    f"(symbol agrees for {(new.ncbi_symbol_agrees=='yes').sum()})")

# NCBI orthologs (finch GeneID -> human GeneID); human symbols from NCBI gene_info
ortho = pd.read_csv(os.path.join(SRC, "NCBI_gene_orthologs_taxid59729_subset.tsv"), sep="\t", dtype=str)
ortho = ortho[(ortho["#tax_id"] == "9606") & (ortho.Other_tax_id == "59729")]
hs_sym = {}
if os.path.exists(os.path.join(SRC, "NCBI_Homo_sapiens.gene_info.gz")):
    gi = pd.read_csv(os.path.join(SRC, "NCBI_Homo_sapiens.gene_info.gz"), sep="\t", dtype=str, usecols=[1, 2])
    hs_sym = dict(zip(gi.iloc[:, 0], gi.iloc[:, 1]))
f2h = dict(zip(ortho.Other_GeneID, ortho.GeneID))
new["ncbi_human_gene_id"] = new.ncbi_gene_id.map(f2h).fillna("")
new["ncbi_human_symbol"] = new.ncbi_human_gene_id.map(hs_sym).fillna("")

# consensus human ortholog (Ensembl RBBH first, then NCBI ortholog)
def consensus(r):
    if r.human_symbol and r.ncbi_human_symbol:
        return r.human_symbol, ("Ensembl RBBH + NCBI agree" if r.human_symbol.upper() == r.ncbi_human_symbol.upper()
                                else "Ensembl RBBH (NCBI disagrees: %s)" % r.ncbi_human_symbol)
    if r.human_symbol: return r.human_symbol, "Ensembl RBBH only"
    if r.ncbi_human_symbol: return r.ncbi_human_symbol, "NCBI ortholog only"
    return "", ""
cons = new.apply(consensus, axis=1, result_type="expand")
new["human_ortholog_symbol_consensus"], new["human_ortholog_source"] = cons[0], cons[1]
# a consensus symbol used by two finch genes is not 1:1 -> blank it
dup = new.human_ortholog_symbol_consensus.ne("") & new.human_ortholog_symbol_consensus.duplicated(keep=False)
new.loc[dup, "human_ortholog_source"] = new.loc[dup, "human_ortholog_source"] + " [dropped: not 1:1]"
new.loc[dup, "human_ortholog_symbol_consensus"] = ""
log(f"  bTaeGut7 genes with a 1:1 human ortholog symbol (consensus): {(new.human_ortholog_symbol_consensus!='').sum()}")

# best display name for each finch gene
def pick(s, n, h):
    real_n = n if n and not n.startswith("LOC") else ""
    for cand in (s, real_n):              # prefer whichever symbol agrees with the 1:1 human ortholog
        if cand and h and cand.upper() == h.upper(): return cand
    return s or n
new["display_symbol"] = [pick(s, n, h) for s, n, h in zip(new.symbol, new.ncbi_symbol, new.human_ortholog_symbol_consensus)]
new["ensembl_ncbi_symbol_conflict"] = ["yes" if s and n and not n.startswith("LOC") and s.upper() != n.upper() else ""
                                       for s, n in zip(new.symbol, new.ncbi_symbol)]
log(f"  Ensembl vs NCBI symbol conflicts (both named, different): {(new.ensembl_ncbi_symbol_conflict=='yes').sum()}")
new.loc[new.display_symbol == "", "display_symbol"] = new.ensembl_gene_id

# ---------------------------------------------------------------- old assembly (Ensembl 104) strict Compara orthologs
def compara(fn, sp):
    p = os.path.join(SRC, fn)
    if not os.path.exists(p): return None
    c = pd.read_csv(p, sep="\t", dtype=str)
    # file is written from the ref species' side; finch is homology_species
    c = c[(c.homology_species == "taeniopygia_guttata") & (c.homology_type == "ortholog_one2one")]
    return c
# Mouse: prefer Compara 116 when present - same bTaeGut1_v1 finch IDs, current mouse gene set.
MOUSE116 = "ENS116_Compara_mus_musculus_vs_tguttata.tsv"
mouse_fn = MOUSE116 if os.path.exists(os.path.join(SRC, MOUSE116)) else "ENS104_Compara_mus_musculus_vs_tguttata.tsv"
mouse_rel = "116" if mouse_fn == MOUSE116 else "104"
old_ortho_release = {}
for fn, tag in [("ENS104_Compara_homo_sapiens_vs_tguttata.tsv", "human"), (mouse_fn, "mouse")]:
    c = compara(fn, tag)
    if c is None:
        log(f"  (missing {fn}; skipping old-assembly {tag} orthologs)"); continue
    rel = mouse_rel if tag == "mouse" else "104"
    old_ortho_release[tag] = rel
    m = c.drop_duplicates("homology_gene_stable_id").set_index("homology_gene_stable_id")
    old[f"{tag}_ensembl_id"] = old.ensembl_gene_id.map(m.gene_stable_id).fillna("")
    old[f"{tag}_perc_id"] = old.ensembl_gene_id.map(m.homology_identity).fillna("")
    old[f"{tag}_high_confidence"] = old.ensembl_gene_id.map(m.is_high_confidence).fillna("")
    log(f"  bTaeGut1_v1 Compara {rel} ortholog_one2one {tag}: {(old[f'{tag}_ensembl_id']!='').sum()} "
        f"(high-confidence: {(old[f'{tag}_high_confidence']=='1').sum()})")

# human/mouse Ensembl gene ID -> symbol, learned from the 2026 homology file (same IDs, current names)
ens2sym = dict(zip(hom.ref_gene_stable_id, hom.ref_gene_name))
for tag in ("human", "mouse"):
    if f"{tag}_ensembl_id" in old:
        old[f"{tag}_symbol"] = old[f"{tag}_ensembl_id"].map(ens2sym).fillna("")
# mouse symbols: NCBI GRCm39 mouse annotation (Ensembl_GeneID -> Symbol) takes priority when present
MOUSE_ANN = "NCBI_mouse_gene_annotations_merged_FMcCarthy.tsv"
if "mouse_ensembl_id" in old and os.path.exists(os.path.join(SRC, MOUSE_ANN)):
    ma = pd.read_csv(os.path.join(SRC, MOUSE_ANN), sep="\t", dtype=str).dropna(subset=["Ensembl_GeneID"])
    ma = ma.drop_duplicates("Ensembl_GeneID")
    msym, mgid = dict(zip(ma.Ensembl_GeneID, ma.Symbol)), dict(zip(ma.Ensembl_GeneID, ma.GeneID))
    s = old.mouse_ensembl_id.map(msym)
    old["mouse_symbol"] = s.fillna(old.mouse_symbol)
    old["mouse_ncbi_gene_id"] = old.mouse_ensembl_id.map(mgid).fillna("")
    nosym = ((old.mouse_ensembl_id != "") & (old.mouse_symbol == "")).sum()
    log(f"  mouse symbols from {MOUSE_ANN}; mouse partners still without a symbol: {nosym}")
if "human_ensembl_id" in old:
    # fall back to HGNC xref when the ID is not in the 2026 file
    miss = (old.human_ensembl_id != "") & (old.human_symbol == "")
    old.loc[miss, "human_symbol"] = old.loc[miss, "hgnc_symbol_xref"]

# ---------------------------------------------------------------- crosswalk bTaeGut1_v1 (Ens 104) -> bTaeGut7.mat (Ens 2026_05)
def pep_md5(fn):
    d = defaultdict(set); g = None; seq = []
    def flush():
        if g and seq: d[hashlib.md5("".join(seq).encode()).hexdigest()].add(g)
    with op(fn) as f:
        for l in f:
            if l.startswith(">"):
                flush(); seq = []
                g = re.search(r"gene:(\S+?)(\.\d+)?\s", l).group(1)
            else: seq.append(l.strip())
    flush(); return d
pn, po = pep_md5("ENS2026_bTaeGut7mat_pep.fa.gz"), pep_md5("ENS104_bTaeGut1v1_pep.fa.gz")
ev = defaultdict(lambda: defaultdict(set))   # old -> new -> evidence
for h, og in po.items():
    for o in og:
        for n in pn.get(h, ()): ev[o][n].add("identical_protein")
# DIAMOND blastp best hits between the two proteomes (precomputed; see README)
rbh_pid = {}
p = os.path.join(SRC, "crosswalk_protein_besthits_diamond.tsv")
if os.path.exists(p):
    dm = pd.read_csv(p, sep="\t")
    dm = dm[dm.reciprocal & (dm.pident >= 80) & (dm.qcov >= 0.5)]
    for o, n, pid in zip(dm.old_gene, dm.new_gene, dm.pident):
        ev[o][n].add("protein_RBH"); rbh_pid[o] = pid
    log(f"  DIAMOND reciprocal best protein hits (>=80% id, >=50% cov): {len(dm)}")
else:
    log("  (crosswalk_protein_besthits_diamond.tsv missing - protein RBH evidence skipped)")
def uniq_map(df, col):
    s = df[df[col] != ""]; s = s[~s[col].str.upper().duplicated(keep=False)]
    return dict(zip(s[col].str.upper(), s.ensembl_gene_id))
nsym = uniq_map(new, "symbol")
osym = uniq_map(old, "symbol")   # symbol must be unique on BOTH sides (fix 26 Oct 2026)
for r in old[old.symbol != ""].itertuples():
    if osym.get(r.symbol.upper()) != r.ensembl_gene_id: continue
    n = nsym.get(r.symbol.upper())
    if n: ev[r.ensembl_gene_id][n].add("same_symbol")
if "human_ensembl_id" in old:
    nh = {h: g for g, h in zip(new.ensembl_gene_id, new.human_ensembl_id) if h}
    for r in old[old.human_ensembl_id != ""].itertuples():
        n = nh.get(r.human_ensembl_id)
        if n: ev[r.ensembl_gene_id][n].add("same_human_ortholog")
cw = []
newinfo = new.set_index("ensembl_gene_id")
for r in old.itertuples():
    cands = ev.get(r.ensembl_gene_id, {})
    if not cands:
        cw.append(dict(old_ensembl_id=r.ensembl_gene_id, old_symbol=r.symbol, new_ensembl_id="", new_symbol="",
                       evidence="", protein_rbh_pident="", n_candidates=0,
                       confidence="unmapped" + ("" if r.biotype == "protein_coding" else " (non-coding)"))); continue
    ranked = sorted(cands.items(), key=lambda kv: (len(kv[1]), "identical_protein" in kv[1]), reverse=True)
    top = [k for k, v in ranked if len(v) == len(ranked[0][1])]
    n, e = ranked[0]
    conf = ("high" if len(e) >= 2 or "identical_protein" in e else "medium")
    if len(top) > 1: conf = "ambiguous"
    cw.append(dict(old_ensembl_id=r.ensembl_gene_id, old_symbol=r.symbol, new_ensembl_id=";".join(top) if len(top) > 1 else n,
                   new_symbol=newinfo.display_symbol.get(n, ""), evidence="+".join(sorted(e)),
                   protein_rbh_pident=rbh_pid.get(r.ensembl_gene_id, "") if "protein_RBH" in e else "",
                   n_candidates=len(cands), confidence=conf))
cw = pd.DataFrame(cw)
log("Crosswalk bTaeGut1_v1 -> bTaeGut7.mat: " + ", ".join(f"{k}={v}" for k, v in cw.confidence.value_counts().items()))

# ---------------------------------------------------------------- Gene Ontology
def read_obo(fn):
    terms, cur = {}, None
    with op(fn) as f:
        for l in f:
            l = l.rstrip("\n")
            if l == "[Term]": cur = {"parents": set(), "alt": []}; continue
            if l.startswith("[") : cur = None; continue
            if cur is None or ": " not in l: continue
            k, v = l.split(": ", 1)
            if k == "id": cur["id"] = v; terms[v] = cur
            elif k == "name": cur["name"] = v
            elif k == "namespace": cur["ns"] = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}[v]
            elif k == "is_a": cur["parents"].add(v.split()[0])
            elif k == "relationship" and v.startswith("part_of "): cur["parents"].add(v.split()[1])
            elif k == "alt_id": cur["alt"].append(v)
            elif k == "is_obsolete" and v == "true": cur["obs"] = True
    alt = {a: t for t, d in terms.items() for a in d["alt"]}
    return terms, alt
terms, alt = read_obo("GO_go-basic.obo")
anc_cache = {}
def ancestors(t):
    if t in anc_cache: return anc_cache[t]
    out = {t}
    for p in terms.get(t, {}).get("parents", ()): out |= ancestors(p)
    anc_cache[t] = out; return out
ROOTS = {"GO:0008150", "GO:0003674", "GO:0005575"}
def propagate(gene2go):
    sets = defaultdict(set)
    for g, ts in gene2go.items():
        for t in ts:
            t = alt.get(t, t)
            if t not in terms or terms[t].get("obs"): continue
            for a in ancestors(t):
                if a not in ROOTS: sets[a].add(g)
    return sets
def write_gmt(sets, fn, note):
    n = 0
    with open(os.path.join(OUT, fn), "w") as f:
        for t in sorted(sets):
            gs = sets[t]
            if not (GMT_MIN <= len(gs) <= GMT_MAX): continue
            d = terms[t]
            f.write(f"{t} {d['name']} [{d['ns']}]\t{note}\t" + "\t".join(sorted(gs)) + "\n"); n += 1
    log(f"  {fn}: {n} GO sets ({GMT_MIN}-{GMT_MAX} genes), {len(set().union(*sets.values())) if sets else 0} genes annotated")

# human GOA (experimental + curated + IEA; NOT-qualified and ND annotations dropped)
hgo = defaultdict(set)
p = os.path.join(SRC, "GO_goa_human.gaf.gz")
if os.path.exists(p):
    with op("GO_goa_human.gaf.gz") as f:
        for l in f:
            if l.startswith("!"): continue
            c = l.split("\t")
            if "NOT" in c[3] or c[6] == "ND": continue
            hgo[c[2]].add(c[4])
else:
    log("  (GO_goa_human.gaf.gz missing - human-projected sets skipped)")
# NCBI finch GAF (GeneID-keyed) -> Ensembl bTaeGut7 via coordinate match
ngo = defaultdict(set)
with op("NCBI_GCF_048771995.1-RS_2025_04_gene_ontology.gaf.gz") as f:
    for l in f:
        if l.startswith("!"): continue
        c = l.split("\t")
        if "NOT" in c[3] or c[6] == "ND": continue
        ngo[c[1]].add(c[4])

finch_native_new = defaultdict(set)
for g, ts in new_go_ens.items(): finch_native_new[g] |= ts
for g, nid in zip(new.ensembl_gene_id, new.ncbi_gene_id):
    if nid in ngo: finch_native_new[g] |= ngo[nid]
proj_new = {g: hgo[s] for g, s in zip(new.ensembl_gene_id, new.human_ortholog_symbol_consensus) if s and s in hgo}
combo_new = defaultdict(set)
for d in (finch_native_new, proj_new):
    for g, ts in d.items(): combo_new[g] |= ts

log("GO gene-set files (GMT):")
write_gmt(propagate(combo_new), "GO_finch_bTaeGut7_EnsemblID_combined.gmt",
          "finch native GO (Ensembl xref + NCBI RefSeq GAF) + human GOA projected via 1:1 orthologs")
write_gmt(propagate(finch_native_new), "GO_finch_bTaeGut7_EnsemblID_finchNativeOnly.gmt",
          "finch native GO only (Ensembl xref + NCBI RefSeq GAF)")
ngo_sets = propagate(ngo)
write_gmt(ngo_sets, "GO_finch_bTaeGut7_NCBIGeneID_RefSeqGAF.gmt", "NCBI RefSeq RS_2025_04 GAF, NCBI GeneID keyed")

if "human_symbol" in old:
    combo_old = defaultdict(set)
    for g, ts in old_go_ens.items(): combo_old[g] |= ts
    for g, s in zip(old.ensembl_gene_id, old.human_symbol):
        if s and s in hgo: combo_old[g] |= hgo[s]
    write_gmt(propagate(combo_old), "GO_finch_bTaeGut1v1_ENS104_EnsemblID_combined.gmt",
              "OLD IDs (bTaeGut1_v1 / Ensembl 104): Ensembl 104 GO xref + human GOA via Compara one2one")
if hgo:
    write_gmt(propagate(hgo), "GO_human_symbol_goa.gmt", "human GOA, HGNC symbol keyed (for ortholog-converted lists)")

# ---------------------------------------------------------------- write tables
new_cols = ["ensembl_gene_id", "display_symbol", "symbol", "biotype", "description", "chrom", "start", "end", "strand",
            "ncbi_gene_id", "ncbi_symbol", "ncbi_description", "ncbi_overlap_jaccard", "ncbi_symbol_agrees",
            "ensembl_ncbi_symbol_conflict",
            "hgnc_symbol_xref", "human_ortholog_symbol_consensus", "human_ortholog_source",
            "human_ensembl_id", "human_symbol", "human_perc_id", "ncbi_human_gene_id", "ncbi_human_symbol",
            "mouse_ensembl_id", "mouse_symbol", "mouse_perc_id", "chicken_ensembl_id", "chicken_symbol", "chicken_perc_id"]
new[new_cols].to_csv(os.path.join(OUT, "finch_bTaeGut7_gene_annotation.tsv"), sep="\t", index=False)
old_cols = [c for c in ["ensembl_gene_id", "symbol", "biotype", "description", "chrom", "start", "end", "strand", "hgnc_symbol_xref",
            "human_ensembl_id", "human_symbol", "human_perc_id", "human_high_confidence",
            "mouse_ensembl_id", "mouse_symbol", "mouse_ncbi_gene_id", "mouse_perc_id", "mouse_high_confidence"] if c in old]
oldo = old[old_cols].merge(cw[["old_ensembl_id", "new_ensembl_id", "new_symbol", "evidence", "confidence"]]
                          .rename(columns={"old_ensembl_id": "ensembl_gene_id", "new_ensembl_id": "bTaeGut7_ensembl_id",
                                           "new_symbol": "bTaeGut7_symbol", "evidence": "crosswalk_evidence",
                                           "confidence": "crosswalk_confidence"}), on="ensembl_gene_id", how="left")
oldo.to_csv(os.path.join(OUT, "finch_bTaeGut1v1_ENS104_gene_annotation.tsv"), sep="\t", index=False)
cw.to_csv(os.path.join(OUT, "crosswalk_bTaeGut1v1_to_bTaeGut7.tsv"), sep="\t", index=False)

# strict 1:1 ortholog tables for finch-mouse comparisons
st = new[new.mouse_ensembl_id != ""][["ensembl_gene_id", "display_symbol", "mouse_ensembl_id", "mouse_symbol", "mouse_perc_id",
                                      "human_ensembl_id", "human_symbol"]]
st.to_csv(os.path.join(OUT, "orthologs_finch_bTaeGut7_to_mouse_1to1.tsv"), sep="\t", index=False)
if "mouse_ensembl_id" in old:
    cols = [c for c in ["ensembl_gene_id", "symbol", "mouse_ensembl_id", "mouse_symbol", "mouse_ncbi_gene_id", "mouse_perc_id",
                        "mouse_high_confidence", "human_ensembl_id", "human_symbol", "bTaeGut7_ensembl_id"] if c in oldo]
    so = oldo[oldo.mouse_ensembl_id != ""][cols]
    rel = old_ortho_release.get("mouse", "104")
    # primary set = one2one AND high confidence; sensitivity set = all one2one
    so.to_csv(os.path.join(OUT, f"orthologs_finch_bTaeGut1v1_to_mouse_1to1_Compara{rel}_all.tsv"), sep="\t", index=False)
    so[so.mouse_high_confidence == "1"].to_csv(
        os.path.join(OUT, f"orthologs_finch_bTaeGut1v1_to_mouse_1to1_Compara{rel}_highconf.tsv"), sep="\t", index=False)
    log(f"Finch-mouse ortholog files (Compara {rel}): all one2one = {len(so)}, high-confidence = {(so.mouse_high_confidence=='1').sum()}")

# ShinyGO: genome-wide backgrounds (prefer your own expressed-gene list when you have one)
bg = sorted(set(new.human_ortholog_symbol_consensus) - {""})
open(os.path.join(OUT, "ShinyGO_background_bTaeGut7_all_human_symbols.txt"), "w").write("\n".join(bg) + "\n")
if "human_symbol" in old:
    bgo = sorted(set(old.human_symbol) - {""})
    open(os.path.join(OUT, "ShinyGO_background_bTaeGut1v1_all_human_symbols.txt"), "w").write("\n".join(bgo) + "\n")
# mitch geneTable (profile gene ID -> ID used in gene sets)
new[new.human_ortholog_symbol_consensus != ""][["ensembl_gene_id", "human_ortholog_symbol_consensus"]] \
    .to_csv(os.path.join(OUT, "mitch_geneTable_bTaeGut7_to_human_symbol.tsv"), sep="\t", index=False, header=False)
if "human_symbol" in old:
    old[old.human_symbol != ""][["ensembl_gene_id", "human_symbol"]] \
        .to_csv(os.path.join(OUT, "mitch_geneTable_bTaeGut1v1_to_human_symbol.tsv"), sep="\t", index=False, header=False)

open(os.path.join(OUT, "build_log.txt"), "w").write("\n".join(log_lines) + "\n")
print("done")
