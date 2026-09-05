"""Every number quoted in a results table of the docs must match the table a script generated.

The scripts in ``scripts/`` write their result tables to ``figures/*.md``. The docs quote
selections of those tables, sometimes with shorter row labels or fewer columns. This test reads
both and compares each quoted value with the generated one, rounded to the precision the docs
show. It does not rerun any simulation; ``scripts/reproduce.py`` regenerates the tables.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
Table = dict[str, dict[str, str]]  # row label -> column header -> cell text


def read_tables(path: Path) -> list[tuple[list[str], list[list[str]]]]:
    """Return every Markdown table in a file as (header, rows)."""
    lines = path.read_text(encoding="utf-8").splitlines()
    tables, i = [], 0
    while i < len(lines):
        is_header = lines[i].startswith("|") and i + 1 < len(lines)
        if is_header and re.fullmatch(r"\|[-| :]+\|", lines[i + 1].strip()):
            header = [c.strip() for c in lines[i].strip().strip("|").split("|")]
            rows, j = [], i + 2
            while j < len(lines) and lines[j].startswith("|"):
                rows.append([c.strip() for c in lines[j].strip().strip("|").split("|")])
                j += 1
            tables.append((header, rows))
            i = j
        else:
            i += 1
    return tables


def find_table(path: Path, first_column: str, contains: str) -> Table:
    """Return the table whose first header is ``first_column`` and that has a row ``contains``."""
    for header, rows in read_tables(path):
        if header[0] == first_column and any(r[0] == contains for r in rows):
            return {r[0]: dict(zip(header, r, strict=True)) for r in rows}
    raise AssertionError(f"no table in {path.name} with row {contains!r}")


def matches(quoted: str, generated: str) -> bool:
    """Compare a quoted cell with the generated one at the precision of the quoted cell."""
    q, g = quoted.replace(" ", ""), generated.replace(" ", "")
    number = r"-?\d+(?:\.\d+)?"
    if re.fullmatch(number, q) and re.fullmatch(number, g):
        decimals = len(q.split(".")[1]) if "." in q else 0
        return f"{float(g):.{decimals}f}" == q
    return q == g


# Each entry: (docs file, first column of the docs table, generated file, first column of the
# generated table, {docs row: generated row}, {docs column: generated column}).
DOCS = ROOT / "docs"
FIG = ROOT / "figures"
CASES = [
    (
        DOCS / "results.md", "Controller", FIG / "controller_comparison.md", "Controller",
        {
            "On/off thermostat (±0.5 K, full speed)": "On/off thermostat",
            "On/off, fan continuous": "On/off, fan continuous",
            "PI + hysteresis low-load handling (A)": "PI + hysteresis (A)",
            "PI + sigma-delta low-load handling (C)": "PI + sigma-delta (C)",
            "PI (C), fan continuous": "PI (C), fan continuous",
        },
        {c: c for c in ("Energy (kWh)", "Starts", "T rms occ. (K)", "RH mean occ. (%)",
                        "h RH > 60 % occ.")},
    ),
    (
        DOCS / "results.md", "Controller", FIG / "humidity_control.md", "Controller",
        {r: r for r in ("PI (C), full fan", "H1 fan RH-PI", "H2 cascade", "On/off (reference)")},
        {c: c for c in ("Energy (kWh)", "Water removed (kg)", "RH mean occ. (%)",
                        "RH std occ. (%)", "h RH > 60 % occ.")},
    ),
    (
        DOCS / "results.md", "Controller", FIG / "mpc_comparison.md", "Controller",
        {
            "PI, 25 °C constant": "PI, 25 °C constant",
            "PI, set-back, fixed start 6:00": "PI, set-back, start 6:00",
            "MPC": "MPC",
            "MPC, time-of-use tariff": "MPC, time-of-use tariff",
            "MPC, model parameters off by 20 %": "MPC, model off by 20 %",
            "MPC, forecast +2 K and -30 % sun": "MPC, forecast +2 K, -30 % sun",
        },
        {"Energy (kWh)": "Energy (kWh)", "Cost (EUR)": "Cost, ToU (EUR)",
         "K·h outside band": "K·h outside band", "Starts": "Starts"},
    ),
    (
        DOCS / "results.md", "Quantity", FIG / "rated_comparison.md", "Quantity",
        "same", {"Published": "Published", "Model": "Model"},
    ),
    (
        DOCS / "results.md", "Point", FIG / "rated_comparison.md", "Point",
        "same", "same",
    ),
    (
        DOCS / "results.md", "Bins", FIG / "seer.md", "Bins",
        "same", "same",
    ),
    (
        DOCS / "results.md", "City, controller", FIG / "season_comparison.md", "City, controller",
        {f"{city}, {d}": f"{city}, {g}" for city in ("Munich", "Madrid")
         for d, g in (("on/off thermostat", "On/off thermostat"), ("PI (C)", "PI (C)"),
                      ("PI (C) + H2 cascade", "PI (C) + H2 cascade"))},
        "same",
    ),
]

VENT_COLUMNS = {"Total energy (kWh)": "Total (kWh)", "Air (m³)": "Air (m³)",
                "CO₂ max occ. (ppm)": "CO₂ max occ. (ppm)", "h > 1 000 ppm": "h > 1000 ppm"}
# Ventilation: the docs table combines rows of three generated tables (variable occupancy).
VENT_ROWS = {
    "V0 fixed schedule (sized for 2 people)": "V0 fixed schedule",
    "V0 sized for the 3-person peak (exploratory)": "V0 sized for 3 people (70 m³/h)",
    "V1 demand-controlled PI": "V1 demand-controlled PI",
    "V2 staged levels": "V2 staged levels",
    "V3 V1 + night flush at 108 m³/h": "V3 V1 + night flush",
}


def _check(docs_table: Table, gen_table: Table, rows: dict[str, str],
           columns: dict[str, str]) -> list[str]:
    errors = []
    for d_row, g_row in rows.items():
        for d_col, g_col in columns.items():
            quoted, generated = docs_table[d_row][d_col], gen_table[g_row][g_col]
            if not matches(quoted, generated):
                errors.append(f"{d_row} / {d_col}: docs {quoted!r}, generated {generated!r}")
    return errors


@pytest.mark.parametrize("case", CASES, ids=lambda c: f"{c[0].name}->{c[2].name}:{c[3]}")
def test_docs_table_matches_generated(case: tuple) -> None:
    docs_file, d_first, gen_file, g_first, rows, columns = case
    gen_any = next(r for h, rs in read_tables(gen_file) if h[0] == g_first for r in rs)[0]
    gen = find_table(gen_file, g_first, gen_any)
    if rows == "same":
        rows = {r: r for r in gen}
    docs = find_table(docs_file, d_first, next(iter(rows)))
    if columns == "same":
        columns = {c: c for c in next(iter(docs.values())) if c != d_first}
    assert not _check(docs, gen, rows, columns)


def test_ventilation_table_matches_generated() -> None:
    docs = find_table(DOCS / "results.md", "Strategy", "V1 demand-controlled PI")
    tables = [
        {r[0]: dict(zip(h, r, strict=True)) for r in rs}
        for h, rs in read_tables(FIG / "ventilation_comparison.md")
    ]
    variable_occupancy = tables[1]
    combined = {**variable_occupancy, **tables[2]}
    assert not _check(docs, combined, VENT_ROWS, VENT_COLUMNS)
