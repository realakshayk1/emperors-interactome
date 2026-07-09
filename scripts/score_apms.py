import glob, json, os, re
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from emperor.mimicry import pdockq

fs = sorted(glob.glob("data/mimicry/models/apms/apms_hpidb_afmultimer_ranked0/*.pdb"))

def one(f):
    r = pdockq(f)
    if not r: return None
    name = os.path.basename(f)[:-4]  # strip .pdb
    # name = <pathogen_protein>_<host_uniprot>; host UniProt is the last _-token if it looks like an accession
    m = re.match(r"^(.*)_([A-NR-Z0-9]{6}|[OPQ][0-9][A-Z0-9]{3}[0-9])$", name)
    if m:
        r["pathogen_protein"], r["host_uniprot"] = m.group(1), m.group(2)
    else:
        r["pathogen_protein"], r["host_uniprot"] = name, None
    r["model"] = name
    return r

out=[]
with ProcessPoolExecutor(max_workers=8) as ex:
    for i,r in enumerate(ex.map(one, fs, chunksize=20)):
        if r: out.append(r)
        if (i+1)%1000==0: print(f"  scored {i+1}/{len(fs)}", flush=True)

json.dump(out, open("data/mimicry/apms_pdockq.json","w"))
pds=np.array([r["pdockq"] for r in out])
print("TOTAL scored:", len(out))
print("host_uniprot parsed:", sum(1 for r in out if r["host_uniprot"]))
print("pDockQ pctiles [5,25,50,75,90,95,99]:", [round(float(np.percentile(pds,p)),3) for p in (5,25,50,75,90,95,99)])
print("n>=0.23:", int((pds>=0.23).sum()), "| n>=0.5:", int((pds>=0.5).sum()))
