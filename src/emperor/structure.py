"""Structure-based physical-validity check for the nominated edge (KANSL3-KANSL1).

Run remotely on Modal (A100-80GB) via Boltz-2 v2.2.1 — an architecture independent
of the AlphaFold-Multimer that produced the audited interactome. NOT a local CPU step;
this module documents the exact remote invocation and parses the harvested outputs into
data/processed/ + results/figures/. See DATA.md §5 and DECISIONS.md.

Remote job (host.compute.create('byoc:modal', image=<boltz2>, gpu='A100-80GB')):
    boltz predict kansl3_kansl1_precomputed.yaml --out_dir ./out \
        --output_format pdb --cache /root/.boltz --num_workers 2
where the YAML supplies both chains + precomputed ColabFold MSA (cached from a prior run).

Physical-validity term for the nonconformity score:
    phys_penalty = max(0, PHYS_IPTM_REF - iptm)      # 0 when interface confidence is high
It stays 0 in the shipped audit (single-metric ipTM path); this module makes the nominee's
interface an explicit, reproducible corroboration rather than a stub.
"""
import json
from pathlib import Path
from . import config as C

# Boltz-2 image built on the user's Modal account (see compute_details ledger):
BOLTZ_IMAGE = "im-HUmTG4YGPRXCeHPD6lJUXV"   # proteomics_boltz_gpu@fc3f36390514cc77
BOLTZ_GPU = "A100-80GB"
PHYS_IPTM_REF = 0.5   # reference: edges below this contribute a physical-validity penalty


def phys_penalty(iptm: float) -> float:
    """Physical-validity penalty from predicted interface ipTM (0 when confident)."""
    return max(0.0, PHYS_IPTM_REF - float(iptm))


def load_physical_validity():
    p = C.ROOT / "data" / "structures" / "physical_validity.json"
    return json.loads(p.read_text()) if p.exists() else None


if __name__ == "__main__":
    pv = load_physical_validity()
    if pv is None:
        print("no physical_validity.json — run the Boltz-2 Modal job first (see module docstring)")
    else:
        g = pv["global_confidence"]; i = pv["confident_interface"]
        print(f"nominee interface (Boltz-2, {pv['gpu']}):")
        print(f"  global ipTM={g['iptm']}  pTM={g['ptm']}")
        print(f"  confident interface KANSL3[{i['kansl3_residues']}]-KANSL1[{i['kansl1_residues']}], "
              f"{i['n_contacts']} contacts, iface pLDDT {i['interface_plddt_kansl3']}/{i['interface_plddt_kansl1']}")
        print(f"  phys_penalty(iptm) = {phys_penalty(g['iptm']):.3f}")
