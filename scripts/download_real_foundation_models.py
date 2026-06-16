from __future__ import annotations

import argparse
import os
from pathlib import Path

from huggingface_hub import snapshot_download


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MODEL_SPECS = {
    "rnafm": "multimolecule/rnafm",
    "rnabert": "multimolecule/rnabert",
}


def load_dotenv(path: Path) -> None:
    if not path.exists():
        return
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        value = value.strip().strip('"').strip("'")
        os.environ.setdefault(key, value)


def get_hf_token() -> str | None:
    load_dotenv(PROJECT_ROOT / ".env")
    load_dotenv(PROJECT_ROOT / ".env.local")
    return os.environ.get("HF_TOKEN") or os.environ.get("HUGGINGFACE_HUB_TOKEN")


def download_model(model_key: str, token: str | None, force: bool) -> Path:
    repo_id = MODEL_SPECS[model_key]
    target = PROJECT_ROOT / "models" / "hf" / model_key
    if target.exists() and not force:
        print(f"{model_key}: already present at {target}")
        return target
    target.mkdir(parents=True, exist_ok=True)
    print(f"{model_key}: downloading {repo_id} -> {target}")
    snapshot_download(
        repo_id=repo_id,
        local_dir=target,
        local_dir_use_symlinks=False,
        token=token,
    )
    return target


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download real RNA foundation model weights for the project.")
    parser.add_argument("--models", nargs="+", choices=sorted(MODEL_SPECS), default=sorted(MODEL_SPECS))
    parser.add_argument("--force", action="store_true", help="Re-download even if the target directory exists.")
    parser.add_argument("--require-token", action="store_true", help="Fail if HF_TOKEN/HUGGINGFACE_HUB_TOKEN is missing.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    token = get_hf_token()
    if args.require_token and not token:
        raise RuntimeError(
            "HF_TOKEN is not configured. Set it in the environment or create .env/.env.local with HF_TOKEN=..."
        )
    for model_key in args.models:
        download_model(model_key, token=token, force=args.force)


if __name__ == "__main__":
    main()
