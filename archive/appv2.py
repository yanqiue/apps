# app.py
import io
import csv
from typing import List, Tuple, Dict, Any, Optional
import streamlit as st
import pandas as pd

#Set page config
st.set_page_config(page_title=None, page_icon="📖", layout="wide", initial_sidebar_state="expanded", menu_items=None)

st.title("Algal CSV Parser")
st.caption("Reads a messy Algal CSV line-by-line.")

# ----------------------------
# Helpers
# ----------------------------
def normalize_cell(s: Optional[str]) -> str:
    if s is None:
        return ""
    return str(s).strip()

def row_contains_all(row: List[str], required: List[str], case_insensitive=True, substring=False) -> bool:
    cells = [normalize_cell(c) for c in row]
    if case_insensitive:
        cells_lookup = [c.lower() for c in cells]
        req = [r.lower() for r in required]
    else:
        cells_lookup = cells
        req = required

    for token in req:
        if substring:
            if not any(token in c for c in cells_lookup):
                return False
        else:
            if token not in cells_lookup:
                return False
    return True

def row_contains_any_substring(row: List[str], tokens: List[str]) -> bool:
    cells = [normalize_cell(c).lower() for c in row]
    return any(any(tok.lower() in c for c in cells) for tok in tokens)

def first_non_empty_row(rows: List[List[str]]) -> Optional[int]:
    for i, r in enumerate(rows):
        if any(normalize_cell(c) for c in r):
            return i
    return None

def parse_metadata_pairs(row: List[str]) -> Dict[str, str]:
    """
    Interpret a single line like: key1, value1, key2, value2, ...
    Any trailing odd item without a pair is ignored.
    """
    cells = [normalize_cell(c) for c in row]
    md = {}
    for i in range(0, len(cells) - 1, 2):
        k = cells[i]
        v = cells[i + 1]
        if k:
            md[k] = v
    return md


def safe_number(x: str) -> Optional[float]:
    x = normalize_cell(x)
    if not x:
        return None
    # allow commas or stray symbols
    x = x.replace(",", "")
    try:
        return float(x)
    except ValueError:
        return None

# ----------------------------
# Line-by-line CSV loader (no pandas yet)
# ----------------------------
def read_csv_lines(file_bytes: bytes, encoding_guess="utf-8") -> List[List[str]]:
    # Try the provided encoding first; fallback to latin-1
    for enc in (encoding_guess, "utf-8-sig", "latin-1"):
        try:
            text = file_bytes.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    reader = csv.reader(io.StringIO(text))
    return [list(row) for row in reader]

# ----------------------------
# Core parsing logic
# ----------------------------
def extract_sections(all_rows: List[List[str]]):
    """
    Returns:
      metadata_dict, metadata_block_text, df1, df2, df3, df2_plus_df3
    """
    # Normalize all rows to trimmed strings
    rows = [[normalize_cell(c) for c in row] for row in all_rows]

    # ---- Find Table 1 header ----
    idx_table1_header = None
    for i, r in enumerate(rows):
        if row_contains_all(r, ["Taxon", "Station"], case_insensitive=True, substring=False):
            idx_table1_header = i
            break

    if idx_table1_header is None:
        raise ValueError("Could not find Table 1 header row containing both 'Taxon' and 'Station'.")

    # ---- Metadata handling: everything above Table 1 header ----
    pre_rows = rows[:idx_table1_header]
    md_row_idx = first_non_empty_row(pre_rows)
    metadata_dict = {}
    if md_row_idx is not None:
        metadata_dict = parse_metadata_pairs(pre_rows[md_row_idx])

    # Also keep a raw text block for all pre-table rows (useful for audit)
    metadata_block_text = "\n".join([",".join(r) for r in pre_rows])

    # ---- Extract Table 1 ----
    header1 = rows[idx_table1_header]
    # Build a mapping of col name -> index for Table 1
    idx_taxon = None
    for j, col in enumerate(header1):
        if col.lower() == "taxon":
            idx_taxon = j
            break
    if idx_taxon is None:
        raise ValueError("Table 1 header located, but no 'Taxon' column was found.")

    # Data for Table 1 goes until BEFORE a line containing both 'Phyto' and 'Diversity:'
    data1: List[List[str]] = []
    i = idx_table1_header + 1
    while i < len(rows):
        r = rows[i]
        if row_contains_all(r, ["Phyto", "Diversity:"], case_insensitive=True, substring=True):
            break  # stop before this row
        data1.append(r)
        i += 1

    # Keep only rows where Taxon value is present
    data1_filtered = [r for r in data1 if idx_taxon < len(r) and normalize_cell(r[idx_taxon])]

    # Convert Table 1 to DataFrame (align widths)
    width1 = max(len(header1), *(len(r) for r in data1_filtered)) if data1_filtered else len(header1)
    header1 = (header1 + [""] * (width1 - len(header1)))[:width1]
    data1_filtered = [(r + [""] * (width1 - len(r)))[:width1] for r in data1_filtered]
    df1 = pd.DataFrame(data1_filtered, columns=header1)

    # ---- Extract Table 2 ----
    # Look for header row like: ["", "mg/m^3", "%", " units or Cells/L", "%"]
    # We'll match loosely: first cell empty-ish, contains mg/m^3, exact "%" at 3rd and 5th positions (loose), one header has "units or Cells/L"
    idx_table2_header = None
    for k in range(i, len(rows)):
        r = rows[k]
        first_empty = (len(r) > 0 and r[0] == "")
        has_mg = any("mg/m^3" == c.lower() for c in r if c)
        has_units_or_cells = any("units or cells/l" in c.lower() for c in r)
        count_percents = sum(1 for c in r if c == "%")
        if first_empty and has_mg and has_units_or_cells and count_percents >= 2:
            idx_table2_header = k
            break

    if idx_table2_header is None:
        raise ValueError("Could not find Table 2 header row (empty, 'mg/m^3', '%', ' units or Cells/L', '%').")

    # Normalize Table 2 headers as required
    # Target headers: ["Alga", "mg/m^3", "mg/m^3_%", "Cells-units/L", "Cells-units/L_%", "ratio"]
    # Note: we will add 'ratio' at the end; the table rows stop before the line containing '===='
    t2_headers_raw = rows[idx_table2_header]
    # Build normalized headers
    t2_headers_norm = []
    # We will scan and map according to the spec
    # Expect roughly 5 first columns of interest
    header_map = []
    percent_seen = 0
    for c in t2_headers_raw:
        lc = c.lower()
        if c == "" and len(t2_headers_norm) == 0:
            t2_headers_norm.append("Alga")
            header_map.append("Alga")
        elif lc == "mg/m^3":
            t2_headers_norm.append("mg/m^3")
            header_map.append("mg/m^3")
        elif c == "%":
            percent_seen += 1
            if percent_seen == 1:
                t2_headers_norm.append("mg/m^3_%")
                header_map.append("mg/m^3_%")
            else:
                t2_headers_norm.append("Cells-units/L_%")
                header_map.append("Cells-units/L_%")
        elif "units or cells/l" in lc:
            t2_headers_norm.append("Cells-units/L")
            header_map.append("Cells-units/L")
        else:
            # Anything else to the right will be dropped later per requirement
            header_map.append(None)

    # Truncate to the five key columns (Alga, mg/m^3, mg/m^3_%, Cells-units/L, Cells-units/L_%)
    keep_headers = ["Alga", "mg/m^3", "mg/m^3_%", "Cells-units/L", "Cells-units/L_%"]
    t2_headers_final = keep_headers + ["ratio"]

    # Collect Table 2 rows until BEFORE a line containing '===='
    data2 = []
    r_idx = idx_table2_header + 1
    while r_idx < len(rows):
        r = rows[r_idx]
        if row_contains_any_substring(r, ["===="]):
            break
        # Build a trimmed row with only the first 5 columns we care about
        trimmed = []
        # Map current row cells by walking raw header_map
        collected = {}
        for raw_idx, col_name in enumerate(header_map):
            if col_name in keep_headers:
                val = rows[r_idx][raw_idx] if raw_idx < len(rows[r_idx]) else ""
                collected[col_name] = val
        trimmed = [collected.get(h, "") for h in keep_headers]
        # Append ratio placeholder as empty
        trimmed.append("")
        data2.append(trimmed)
        r_idx += 1

    df2 = pd.DataFrame(data2, columns=t2_headers_final)

    # ---- Skip lines after '====' until Table 3 header ----
    idx_table3_header = None
    for k in range(r_idx, len(rows)):
        r = rows[k]
        # Match header: ["", "mg/m^3", "Cells-units/L", "ratio"]
        if (len(r) >= 4 and
            r[0] == "" and
            any(c.lower() == "mg/m^3" for c in r) and
            any(c.lower() == "cells-units/l" for c in r) and
            any(c.lower() == "ratio" for c in r)):
            idx_table3_header = k
            break

    if idx_table3_header is None:
        raise ValueError("Could not find Table 3 header row (empty, 'mg/m^3', 'Cells-units/L', 'ratio').")

    # ---- Extract Table 3 ----
    # Target (after insertion): ["Alga", "mg/m^3", "mg/m^3_%", "Cells-units/L", "Cells-units/L_%", "ratio"]
    t3_target_headers = ["Alga", "mg/m^3", "mg/m^3_%", "Cells-units/L", "Cells-units/L_%", "ratio"]
    data3 = []
    k = idx_table3_header + 1
    # We'll read until an all-empty row or end of file
    while k < len(rows):
        r = rows[k]
        if not any(normalize_cell(c) for c in r):
            break
        # Build mapping from the raw position:
        # Expect raw columns in order: ["", "mg/m^3", "Cells-units/L", "ratio"]
        alga = r[0] if len(r) > 0 else ""
        mg = r[1] if len(r) > 1 else ""
        cells_units = r[2] if len(r) > 2 else ""
        ratio = r[3] if len(r) > 3 else ""
        row_out = [alga, mg, "", cells_units, "", ratio]  # Insert the % columns as empty
        data3.append(row_out)
        k += 1

    df3 = pd.DataFrame(data3, columns=t3_target_headers)

    # ---- Concatenate Table 3 under Table 2 ----
    # Ensure same columns/order
    df2_aligned = df2[t3_target_headers] if list(df2.columns) != t3_target_headers else df2
    df2_plus_df3 = pd.concat([df2_aligned, df3], ignore_index=True)

    return metadata_dict, metadata_block_text, df1, df2, df3, df2_plus_df3

# ----------------------------
# Streamlit UI
# ----------------------------
with st.sidebar:
    st.header("Upload CSV")
    f = st.file_uploader("Choose a CSV file", type=["csv"])

if not f:
    st.info("Upload a CSV to begin. The app parses it line-by-line and builds the outputs.")
    st.stop()

# Read file as lines (no pandas yet)
try:
    all_rows = read_csv_lines(f.read())
except Exception as e:
    st.error(f"Failed to read CSV: {e}")
    st.stop()

# Parse into sections
try:
    metadata_dict, metadata_block_text, df1, df2, df3, df2_plus_df3 = extract_sections(all_rows)
except Exception as e:
    st.error(f"Parsing error: {e}")
    st.stop()

# ----------------------------
# Display & Downloads
# ----------------------------
st.subheader("Metadata (key/value pairs from the first non-empty line above Table 1)")
col_md1, col_md2 = st.columns([1,1])
with col_md1:
    st.json(metadata_dict or {"(none)": ""})
with col_md2:
    st.text_area("Raw metadata block (all lines above Table 1)", metadata_block_text, height=180)

# Download metadata
md_json = io.BytesIO(pd.Series(metadata_dict).to_json(indent=2).encode("utf-8"))
st.download_button("Download metadata.json", md_json, file_name="metadata.json", mime="application/json")

raw_md_txt = io.BytesIO(metadata_block_text.encode("utf-8"))
st.download_button("Download metadata_raw.txt", raw_md_txt, file_name="metadata_raw.txt", mime="text/plain")

# Table 1
st.subheader("Table 1 (from 'Taxon'/'Station' header, stops before 'Phyto' + 'Diversity:')")
st.dataframe(df1, use_container_width=True)
csv1 = io.BytesIO(df1.to_csv(index=False).encode("utf-8"))
st.download_button("Download table1.csv", csv1, file_name="table1.csv", mime="text/csv")

# Table 2
st.subheader("Table 2 (normalized headers, truncated after 'Cells-units/L_%', + 'ratio')")
st.caption("Stops before a line containing '===='. Header changes: empty→'Alga', 'units or Cells/L'→'Cells-units/L', '%'→'mg/m^3_%' and 'Cells-units/L_%'.")
st.dataframe(df2, use_container_width=True)
csv2 = io.BytesIO(df2.to_csv(index=False).encode("utf-8"))
st.download_button("Download table2.csv", csv2, file_name="table2.csv", mime="text/csv")

# Table 3
st.subheader("Table 3 (header ['', 'mg/m^3', 'Cells-units/L', 'ratio'] → insert '% columns to match Table 2')")
st.dataframe(df3, use_container_width=True)
csv3 = io.BytesIO(df3.to_csv(index=False).encode("utf-8"))
st.download_button("Download table3.csv", csv3, file_name="table3.csv", mime="text/csv")

# Combined
st.subheader("Combined: Table 2 ⊕ Table 3 (same schema)")
st.dataframe(df2_plus_df3, use_container_width=True)
csv23 = io.BytesIO(df2_plus_df3.to_csv(index=False).encode("utf-8"))
st.download_button("Download table2_plus_table3.csv", csv23, file_name="table2_plus_table3.csv", mime="text/csv")

# Notes & next steps
with st.expander("Notes & adjustments"):
    st.markdown("""
- **Ratio column**: For Table 2, this app adds the column but does **not** compute it (you didn’t specify a formula).  
  If you'd like, we can auto-calc `ratio = (mg/m^3) / (Cells-units/L)` when both are numeric.
- **Header matching is robust but conservative**: it’s case-insensitive, trims spaces, and allows substring checks where appropriate.
- **Stops & starts**:
  - Table 1 stops just before the first row that has both “Phyto” and “Diversity:” anywhere in the row.
  - Table 2 stops before the first row that contains “====”.
  - After “====”, the app skips rows until it finds the Table 3 header row.
- If your real files have slight variations, tell me the exact quirks and I’ll tweak the matchers.
""")
