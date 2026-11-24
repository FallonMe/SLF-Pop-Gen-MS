#!/usr/bin/env python3
"""Validate archived sample metadata without modifying files (Python 3, stdlib)."""
import csv
import math
import sys
from pathlib import Path

def read_metadata(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"SLF_ID", "lat", "lon"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{path.name}: missing required columns")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path.name}: no sample records")
    identifiers = [row["SLF_ID"].strip() for row in rows]
    if any(not key for key in identifiers) or len(set(identifiers)) != len(identifiers):
        raise ValueError(f"{path.name}: blank or duplicate sample IDs")
    for number, row in enumerate(rows, 2):
        lat, lon = float(row["lat"]), float(row["lon"])
        if not (math.isfinite(lat) and math.isfinite(lon) and -90 <= lat <= 90 and -180 <= lon <= 180):
            raise ValueError(f"{path.name}: invalid coordinates at row {number}")
    return set(identifiers), len(rows)

def main():
    directory = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data")
    first, first_n = read_metadata(directory / "Final_SLF_MS_Info.csv")
    second, second_n = read_metadata(directory / "Environmental_PC1s.csv")
    if first != second:
        raise ValueError("Sample-ID sets differ between metadata tables")
    print(f"PASS: {first_n} and {second_n} records; unique matching IDs; valid coordinates.")

if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        sys.exit(1)
