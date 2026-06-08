from __future__ import annotations

import argparse
import csv
import json
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen


PROJECT_ROOT = Path(__file__).resolve().parents[1]
API = "https://gtexportal.org/api/v2/association"
DEFAULT_TISSUES = ["Whole_Blood", "Liver", "Heart_Left_Ventricle", "Muscle_Skeletal", "Brain_Cortex"]


def fetch_json(path: str, params: dict[str, object]) -> dict[str, object]:
    url = f"{API}/{path}?{urlencode(params)}"
    with urlopen(url, timeout=60) as handle:
        return json.loads(handle.read().decode("utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch a small real GTEx v8 sQTL case-study table.")
    parser.add_argument("--out-dir", type=Path, default=PROJECT_ROOT / "data/raw")
    parser.add_argument("--items-per-tissue", type=int, default=5)
    parser.add_argument("--variants-per-gene", type=int, default=10)
    parser.add_argument("--sleep", type=float, default=0.2)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    sgene_rows: list[dict[str, object]] = []
    sqtl_rows: list[dict[str, object]] = []
    event_rows: list[dict[str, object]] = []

    for tissue in DEFAULT_TISSUES:
        payload = fetch_json(
            "sgene",
            {
                "datasetId": "gtex_v8",
                "tissueSiteDetailId": tissue,
                "itemsPerPage": args.items_per_tissue,
                "page": 0,
            },
        )
        genes = payload.get("data", [])
        for gene in genes:
            sgene_rows.append(gene)
            phenotype = str(gene.get("phenotypeId", ""))
            event_rows.append(
                {
                    "event_id": phenotype,
                    "tissue": tissue,
                    "gene_symbol": gene.get("geneSymbol", ""),
                    "gencode_id": gene.get("gencodeId", ""),
                    "p_value": gene.get("pValue", ""),
                    "q_value": gene.get("qValue", ""),
                    "phenotype_id": phenotype,
                    "data_source": "GTEx v8 Portal API /association/sgene",
                }
            )
            sqtl_payload = fetch_json(
                "singleTissueSqtl",
                {
                    "datasetId": "gtex_v8",
                    "tissueSiteDetailId": tissue,
                    "gencodeId": gene.get("gencodeId", ""),
                    "itemsPerPage": args.variants_per_gene,
                    "page": 0,
                },
            )
            for variant in sqtl_payload.get("data", []):
                sqtl_rows.append(variant)
            time.sleep(args.sleep)

    write_tsv(args.out_dir / "gtex_sgenes.tsv", sgene_rows)
    write_tsv(args.out_dir / "gtex_sqtl.tsv", sqtl_rows)
    write_tsv(args.out_dir / "known_splice_events.tsv", event_rows)
    print(f"wrote {len(sgene_rows)} sgenes, {len(sqtl_rows)} sQTL variants, {len(event_rows)} splice events")


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows to write: {path}")
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()

