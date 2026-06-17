from __future__ import annotations

import argparse
import re
from pathlib import Path


BEGIN_RE = re.compile(r"\\begin\{([^}]+)\}")
END_RE = re.compile(r"\\end\{([^}]+)\}")
GRAPHICS_RE = re.compile(r"\\includegraphics(?:\[[^\]]*\])?\{(?:\\detokenize\{)?([^{}]+)(?:\})?\}")
CITE_RE = re.compile(r"\\(?:cite|nocite)\{([^}]+)\}")
BIB_ENTRY_RE = re.compile(r"@\w+\{([^,]+),")


def check_report(tex_path: Path, bib_path: Path) -> list[str]:
    errors: list[str] = []
    text = tex_path.read_text(encoding="utf-8")
    bib = bib_path.read_text(encoding="utf-8") if bib_path.exists() else ""

    begin_stack: list[str] = []
    tokens = []
    for match in BEGIN_RE.finditer(text):
        tokens.append((match.start(), "begin", match.group(1)))
    for match in END_RE.finditer(text):
        tokens.append((match.start(), "end", match.group(1)))
    for _pos, kind, name in sorted(tokens):
        if kind == "begin":
            begin_stack.append(name)
        elif not begin_stack or begin_stack.pop() != name:
            errors.append(f"Mismatched LaTeX environment near \\end{{{name}}}.")
            break
    if begin_stack:
        errors.append(f"Unclosed LaTeX environment(s): {', '.join(begin_stack)}")

    for raw_path in GRAPHICS_RE.findall(text):
        image_path = (tex_path.parent / raw_path).resolve()
        if not image_path.exists():
            errors.append(f"Missing image file: {raw_path}")

    bib_keys = set(BIB_ENTRY_RE.findall(bib))
    for cite_group in CITE_RE.findall(text):
        for key in [item.strip() for item in cite_group.split(",") if item.strip()]:
            if key != "*" and key not in bib_keys:
                errors.append(f"Missing bibliography key: {key}")

    if text.count(r"\begin{document}") != 1 or text.count(r"\end{document}") != 1:
        errors.append("Expected exactly one \\begin{document} and one \\end{document}.")

    return errors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Sanity-check the LaTeX report when no TeX engine is available.")
    parser.add_argument("--tex", type=Path, default=Path("report_letax/njuthesis-sample.tex"))
    parser.add_argument("--bib", type=Path, default=Path("report_letax/njuthesis-sample.bib"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    errors = check_report(args.tex, args.bib)
    if errors:
        print("LaTeX sanity check failed:")
        for error in errors:
            print(f"- {error}")
        raise SystemExit(1)
    print("LaTeX sanity check passed: environments, images, and citations are valid.")


if __name__ == "__main__":
    main()
