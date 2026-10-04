# SPDX-License-Identifier: Apache-2.0
"""The site declaration's settings, checked against and completed from site_reference.toml
(the revision plan MEDS_FAST_CALIBRATION_REVISION_PLAN.md §2).

site_reference.toml is the schema: a key a site file sets that is not there stops the tool (a
misspelling would otherwise silently keep the default), and its values are the defaults -- except
the REQUIRED keys, whose values there are examples, and the EXAMPLE lists, which default to empty.
FREE tables take any keys: model config keys ([overrides], [variants.*], [calibrated]) or the
registry's key names ([priors.*]).
"""
from __future__ import annotations

import copy
from pathlib import Path

from meds.config import load_toml
from residuals import FILTERS

REFERENCE = Path(__file__).resolve().parent / "site_reference.toml"
#: keys a site must set (their reference values are examples)
REQUIRED = {("base", "main"), ("base", "registry"), ("tower", "site")}
#: lists whose reference entries are examples: the default is empty
EXAMPLES = {("windows", "list"), ("windows", "seasonal", "list")}
#: tables that take any keys
FREE = {("overrides",), ("variants",), ("calibrated",), ("priors",)}
#: the keys of one entry of a list of tables
LIST_ENTRY = {("windows", "list"): {"name", "start", "role", "days"},
              ("windows", "seasonal", "list"): {"name", "start", "days"}}
PRIOR_ENTRY = {"centre", "sd", "log_sd", "source", "range"}
#: every target takes the filters (residuals.FILTERS), whether or not its reference table lists them


def reference() -> dict:
    return load_toml(REFERENCE)


def _check(decl: dict, ref: dict, path: tuple, errors: list):
    for k, v in decl.items():
        here = path + (k,)
        if path in FREE:
            if path == ("priors",):
                if not isinstance(v, dict):
                    errors.append(f"[priors.{k}] must be a table")
                else:
                    bad = set(v) - PRIOR_ENTRY
                    if bad:
                        errors.append(f"[priors.{k}]: unknown settings {sorted(bad)}; known: {sorted(PRIOR_ENTRY)}")
            elif path == ("variants",) and not isinstance(v, dict):
                errors.append(f"[variants.{k}] must be a table of model config keys")
            continue
        if len(path) == 2 and path[0] == "targets" and k in FILTERS:
            continue
        if k not in ref:
            where = ".".join(path) or "(top level)"
            errors.append(f"[{where}] has no setting '{k}' (site_reference.toml lists them all)")
            continue
        if isinstance(v, dict) and isinstance(ref[k], dict):
            _check(v, ref[k], here, errors)
        elif here in LIST_ENTRY:
            if not isinstance(v, list):
                errors.append(f"[{'.'.join(path)}].{k} must be a list of tables")
                continue
            for i, e in enumerate(v):
                bad = set(e) - LIST_ENTRY[here] if isinstance(e, dict) else {"(not a table)"}
                if bad:
                    errors.append(f"[{'.'.join(path)}].{k}[{i}]: unknown settings {sorted(bad)}; "
                                  f"known: {sorted(LIST_ENTRY[here])}")


def _merge(defaults: dict, decl: dict, path: tuple) -> dict:
    out = copy.deepcopy(defaults)
    for k, v in decl.items():
        here = path + (k,)
        if isinstance(v, dict) and isinstance(out.get(k), dict) and here not in FREE:
            out[k] = _merge(out[k], v, here)
        else:
            out[k] = copy.deepcopy(v)
    return out


def _defaults(ref: dict) -> dict:
    d = copy.deepcopy(ref)
    for path in REQUIRED:
        node = d
        for k in path[:-1]:
            node = node.get(k, {})
        node.pop(path[-1], None)
    for path in EXAMPLES:
        node = d
        for k in path[:-1]:
            node = node.get(k, {})
        if path[-1] in node:
            node[path[-1]] = []
    return d


def complete(decl: dict) -> dict:
    """The site declaration checked against the reference and completed with its defaults. Raises
    ValueError listing every unknown or missing setting at once."""
    ref = reference()
    errors = []
    _check(decl, ref, (), errors)
    for path in REQUIRED:
        node = decl
        for k in path:
            if not isinstance(node, dict) or k not in node:
                errors.append(f"[{'.'.join(path[:-1])}].{path[-1]} is required")
                break
            node = node[k]
    if errors:
        raise ValueError("the site declaration does not match site_reference.toml:\n  " + "\n  ".join(errors))
    return _merge(_defaults(ref), decl, ())
