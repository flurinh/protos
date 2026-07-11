#!/usr/bin/env python3
"""Verify executable examples and local links in first-party documentation."""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCS = [
    ROOT / "README.md",
    *sorted((ROOT / "docs").rglob("*.md")),
    ROOT / "src/protos/models/protein_design.md",
    ROOT / "src/protos/models/remote.md",
    ROOT / "src/protos/reference_data/grn/README.md",
]

FENCE = re.compile(
    r"^```(?P<language>python|bash|sh)\s*\n(?P<body>.*?)^```\s*$",
    re.MULTILINE | re.DOTALL,
)
MARKDOWN_LINK = re.compile(r"\[[^\]]*\]\(([^)]+)\)")
HTML_IMAGE = re.compile(r'<img\s+[^>]*src="([^"]+)"')


def is_network_example(code: str) -> bool:
    return "download_and_register(\"1ubq\")" in code


def verify_python(path: Path, block_number: int, code: str) -> None:
    with tempfile.TemporaryDirectory(prefix="protos-docs-") as data_root:
        env = os.environ.copy()
        env["PROTOS_DATA_ROOT"] = data_root
        env["PYTHONPATH"] = str(ROOT / "src")
        env["PYTHONDONTWRITEBYTECODE"] = "1"
        result = subprocess.run(
            [sys.executable, "-c", code],
            cwd=ROOT,
            env=env,
            text=True,
            capture_output=True,
            timeout=90,
            check=False,
        )
    if result.returncode:
        raise RuntimeError(
            f"{path.relative_to(ROOT)} Python block {block_number} failed\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )


def verify_shell(path: Path, block_number: int, code: str) -> None:
    result = subprocess.run(
        ["bash", "-n"],
        input=code,
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode:
        raise RuntimeError(
            f"{path.relative_to(ROOT)} shell block {block_number} failed syntax check\n"
            f"{result.stderr}"
        )


def verify_links(path: Path, text: str) -> int:
    checked = 0
    targets = MARKDOWN_LINK.findall(text) + HTML_IMAGE.findall(text)
    for target in targets:
        if target.startswith(("http://", "https://", "mailto:", "#")):
            continue
        checked += 1
        local_target = target.split("#", 1)[0]
        if not (path.parent / local_target).resolve().exists():
            raise FileNotFoundError(
                f"{path.relative_to(ROOT)} points to missing local target: {target}"
            )
    return checked


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--network",
        action="store_true",
        help="also execute examples that download 1ubq from RCSB",
    )
    args = parser.parse_args()

    python_count = 0
    network_skipped = 0
    shell_count = 0
    link_count = 0

    for path in DOCS:
        text = path.read_text(encoding="utf-8")
        link_count += verify_links(path, text)
        language_counts: dict[str, int] = {}
        for match in FENCE.finditer(text):
            language = match.group("language")
            code = match.group("body")
            language_counts[language] = language_counts.get(language, 0) + 1
            block_number = language_counts[language]
            if language == "python":
                if is_network_example(code) and not args.network:
                    network_skipped += 1
                    continue
                verify_python(path, block_number, code)
                python_count += 1
            else:
                verify_shell(path, block_number, code)
                shell_count += 1

    print(
        f"verified {python_count} Python blocks, {shell_count} shell blocks, "
        f"and {link_count} local links/assets"
    )
    if network_skipped:
        print(f"skipped {network_skipped} network blocks; pass --network to execute them")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
