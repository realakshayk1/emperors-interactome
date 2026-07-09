import pandas as pd, json, urllib.request, time

df = pd.read_parquet("data/processed/certified.parquet")
hc = df[df["high_conf"]].copy()

# collect all UniProt accessions in high-conf pairs
accs = sorted(set(hc["uniprot_a"]) | set(hc["uniprot_b"]))
print("unique accessions in high-conf:", len(accs))

# fetch lengths from UniProt REST (batch)
lengths = {}
CHUNK = 100
for i in range(0, len(accs), CHUNK):
    batch = accs[i:i + CHUNK]
    q = "+OR+".join(f"accession:{a}" for a in batch)
    url = f"https://rest.uniprot.org/uniprotkb/search?query={q}&fields=accession,length&format=tsv&size=500"
    try:
        txt = urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "u"}), timeout=60).read().decode()
        for line in txt.splitlines()[1:]:
            parts = line.split("\t")
            if len(parts) >= 2:
                lengths[parts[0]] = int(parts[1])
    except Exception as e:
        print("chunk err", i, str(e)[:80])
    time.sleep(0.3)

print("lengths fetched:", len(lengths))
hc["len_a"] = hc["uniprot_a"].map(lengths)
hc["len_b"] = hc["uniprot_b"].map(lengths)
hc["total_len"] = hc["len_a"] + hc["len_b"]
hc = hc.dropna(subset=["total_len"])
hc["total_len"] = hc["total_len"].astype(int)

cert = hc[hc["certified@0.1"]].copy()
drop = hc[~hc["certified@0.1"]].copy()

# Score and certification are nearly separable (dropped all <0.5, certified all >=0.5),
# so a score-matched test is impossible from the map. Refocus: does Boltz-2 (independent
# architecture) corroborate the conformal verdict? Pick smallest-by-length in each group.
print(f"dropped score range [{drop['score'].min():.3f},{drop['score'].max():.3f}]; "
      f"certified [{cert['score'].min():.3f},{cert['score'].max():.3f}] — near-separable")

# smallest by total_len within each group (cheapest folds)
drop_pick = drop.nsmallest(6, "total_len")
cert_pick = cert.nsmallest(6, "total_len")
pilot = pd.concat([drop_pick.assign(label="dropped"), cert_pick.assign(label="certified")])
pilot = pilot[["pair", "gene_a", "gene_b", "uniprot_a", "uniprot_b", "len_a", "len_b", "total_len", "score", "conf_pvalue", "label"]]
pilot.to_csv("data/structures/pilot_pairs.csv", index=False)
print("\nPILOT SET (", len(pilot), "pairs):")
print(pilot.to_string(index=False))
print("\ntotal_len stats: min", int(pilot.total_len.min()), "median", int(pilot.total_len.median()), "max", int(pilot.total_len.max()))
