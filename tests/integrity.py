"""Verify preserved source/assets, publication links, and saved examples."""
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]


def main():
    checked = 0
    for name in ("source_manifest.json", "figure_manifest.json"):
        manifest = json.loads((ROOT / "docs" / name).read_text(encoding="utf-8-sig"))
        for entry in manifest["files"]:
            path = ROOT / entry["path"]
            assert path.is_file(), path
            assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"], path
            checked += 1
    links = 0
    for path in ROOT.rglob("*.md"):
        if ".git" in path.parts or "results" in path.parts:
            continue
        for target in re.findall(r"\]\(([^)]+)\)", path.read_text(encoding="utf-8")):
            if urlsplit(target).scheme or target.startswith("#"):
                continue
            relative = unquote(target.split("#", 1)[0])
            assert (path.parent / relative).exists(), (path, target)
            links += 1
    for path in (ROOT / "examples" / "expected").glob("*.csv"):
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        assert rows, path
        assert all(math.isfinite(float(v)) for r in rows for v in r.values()), path
    if "--plots" in sys.argv:
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "trace.png"
            command = [sys.executable, str(ROOT / "examples" / "plot.py"),
                       str(ROOT / "examples" / "expected"), str(output)]
            subprocess.run(command, check=True)
            before = output.read_bytes()
            assert before.startswith(b"\x89PNG\r\n\x1a\n")
            again = subprocess.run(command, capture_output=True)
            assert again.returncode != 0
            assert output.read_bytes() == before
    print(f"Verified {checked} hashes, {links} local links, and saved CSV data.")


if __name__ == "__main__":
    main()
