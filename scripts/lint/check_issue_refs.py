#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""check_issue_refs.py -- source comments and reference docs cite no issue or pull-request numbers.

A comment states why the code is as it is, in the present tense; an issue number records when and
why it changed, which is history. History lives in CHANGELOG.md, deferred work in docs/ROADMAP.md
(each item with its issue), and design records in docs/dev_plans/. Those three are the only places
an issue or pull-request number may appear.

Flagged: `#123`, `issue 123`, `PR 123`, `pull request 123`, and GitHub issue or pull URLs. Not
flagged: a Markdown anchor such as `(#4-water-stress)`, an HTML entity, a colour such as `'#333'`.

Usage: check_issue_refs.py REPO_ROOT   (exit 1 and a list when a reference is found)
"""
import re
import subprocess
import sys
from pathlib import Path

# The only files that may cite issue and pull-request numbers (this script quotes the patterns).
ALLOWED = ("CHANGELOG.md", "docs/ROADMAP.md", "scripts/lint/check_issue_refs.py")
ALLOWED_DIRS = ("docs/dev_plans/",)
TEXT = (".f90", ".F90", ".py", ".ipynb", ".toml", ".cmake", ".txt", ".md", ".sh", ".yml", ".yaml", ".cfg",
        ".ini", ".in", ".c", ".h")

REFS = [
    re.compile(r"(?<![\w&'\"/])#\d{1,5}\b(?![-\w])"),                       # #123 (not an anchor)
    re.compile(r"\b(?:issues?|PRs?|pull requests?)\s+#?\d{1,5}\b", re.I),     # issue 123, PR 12
    re.compile(r"github\.com/[\w.-]+/[\w.-]+/(?:issues|pull)/\d+"),           # a GitHub URL
]


def tracked_files(root: Path):
    out = subprocess.run(["git", "-C", str(root), "ls-files"], capture_output=True, text=True)
    if out.returncode == 0 and out.stdout:
        return [root / f for f in out.stdout.split("\n") if f]
    # No git (a source tarball): the source folders and the top-level files, never a build tree.
    folders = ("src", "test", "python", "scripts", "post_proc", "examples", "docs", "cmake", ".claude")
    files = [p for p in root.iterdir() if p.is_file()]
    for name in folders:
        if (root / name).is_dir():
            files += [p for p in (root / name).rglob("*") if p.is_file()]
    return files


def issue_refs(root: Path):
    found = []
    for path in tracked_files(root):
        rel = path.relative_to(root).as_posix()
        if rel in ALLOWED or rel.startswith(ALLOWED_DIRS) or not rel.endswith(TEXT) and path.name != "CMakeLists.txt":
            continue
        try:
            text = path.read_text(errors="replace")
        except (FileNotFoundError, IsADirectoryError):
            continue
        for no, line in enumerate(text.splitlines(), 1):
            for pat in REFS:
                m = pat.search(line)
                if m:
                    found.append(f"{rel}:{no}: {m.group(0)}: {line.strip()[:110]}")
                    break
    return found


if __name__ == "__main__":
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    bad = issue_refs(root)
    for line in bad:
        print(line)
    if bad:
        print(f"\n{len(bad)} issue or pull-request reference(s) outside CHANGELOG.md, docs/ROADMAP.md and "
              "docs/dev_plans/. Say why the code is as it is; put the history in CHANGELOG.md.")
        sys.exit(1)
