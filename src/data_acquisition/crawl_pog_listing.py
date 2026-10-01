"""Recursively crawl the Apache index at bcgsc.ca/downloads/nanopore_pog/ and
record every file with its exact byte size (from a HEAD request)."""
import re, sys, json, urllib.request, urllib.parse
from concurrent.futures import ThreadPoolExecutor

BASE = "https://www.bcgsc.ca/downloads/nanopore_pog/"
HREF = re.compile(r'<a href="([^"?][^"]*)">')

def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "curl/8"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "replace")

def head_size(url):
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "curl/8"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers.get("Content-Length", -1))

files, dirs = [], [BASE]
while dirs:
    d = dirs.pop()
    for href in HREF.findall(get(d)):
        if href.startswith("/") or href.startswith("http"):
            continue  # parent dir / logo
        full = urllib.parse.urljoin(d, href)
        (dirs if href.endswith("/") else files).append(full)

with ThreadPoolExecutor(8) as ex:
    sizes = list(ex.map(head_size, files))

out = [{"url": u, "rel": u[len(BASE):], "bytes": s} for u, s in zip(files, sizes)]
json.dump(out, open(sys.argv[1], "w"), indent=1)
print(len(out), "files,", sum(o["bytes"] for o in out), "bytes total")
