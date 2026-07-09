import pandas as pd, urllib.request, os, time

pilot = pd.read_csv("data/structures/pilot_pairs.csv")
os.makedirs("data/structures/pilot", exist_ok=True)

def fetch_seq(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    txt = urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "u"}), timeout=60).read().decode()
    return "".join(txt.splitlines()[1:])

seqs = {}
accs = sorted(set(pilot["uniprot_a"]) | set(pilot["uniprot_b"]))
for a in accs:
    seqs[a] = fetch_seq(a); time.sleep(0.2)
print("fetched", len(seqs), "sequences")

# write one Boltz FASTA per pair (chains A|protein / B|protein)
manifest = []
for _, r in pilot.iterrows():
    name = r["pair"]
    sa, sb = seqs[r["uniprot_a"]], seqs[r["uniprot_b"]]
    path = f"data/structures/pilot/{name}.fasta"
    with open(path, "w") as f:
        f.write(f">A|protein\n{sa}\n>B|protein\n{sb}\n")
    manifest.append({"pair": name, "label": r["label"], "total_len": int(r["total_len"]),
                     "afm_score": float(r["score"]), "conf_pvalue": float(r["conf_pvalue"]),
                     "fasta": path})

import json
json.dump(manifest, open("data/structures/pilot/manifest.json", "w"), indent=1)
print("wrote", len(manifest), "boltz FASTA inputs")
for m in manifest:
    print(f"  {m['pair']:18s} {m['label']:9s} len={m['total_len']:4d} afm={m['afm_score']:.3f}")
