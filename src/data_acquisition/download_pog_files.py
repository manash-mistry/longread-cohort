"""Download the in-scope files from listing.json into data/pog/, preserving
the site's folder structure and verifying each file's size on disk."""
import json, os, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor

DEST = "/Users/mmistry/longread-cohort/data/pog"
ALLOWED = (".tsv", ".csv", ".txt", ".bed", ".vcf",
           ".tsv.gz", ".csv.gz", ".txt.gz", ".bed.gz", ".vcf.gz")
SKIP_EXT = (".bam", ".cram", ".fastq", ".fq", ".fast5", ".pod5", ".tar.gz", ".rmd.gz")

L = json.load(open(sys.argv[1]))
todo = [f for f in L if f["rel"].lower().endswith(ALLOWED) and not f["rel"].lower().endswith(SKIP_EXT)]
skipped = [f for f in L if f not in todo]

def fetch(f):
    path = os.path.join(DEST, f["rel"])
    if os.path.exists(path) and os.path.getsize(path) == f["bytes"]:
        return (f["rel"], "cached")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    req = urllib.request.Request(f["url"], headers={"User-Agent": "curl/8"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=120) as r, open(path, "wb") as out:
                while chunk := r.read(1 << 20):
                    out.write(chunk)
            if os.path.getsize(path) == f["bytes"]:
                return (f["rel"], "ok")
        except Exception as e:
            err = e
    return (f["rel"], f"FAILED after 3 attempts")

with ThreadPoolExecutor(6) as ex:
    results = list(ex.map(fetch, todo))

ok = [r for r in results if r[1] in ("ok", "cached")]
bad = [r for r in results if r not in ok]
print(f"downloaded/verified {len(ok)} of {len(todo)} files; "
      f"{sum(f['bytes'] for f in todo)/1e9:.2f} GB")
for r in bad: print("FAILED:", r[0])
print("\nSkipped (not in allowed types):")
for f in skipped: print(f"  {f['bytes']/1e9:.2f} GB  {f['rel']}")
