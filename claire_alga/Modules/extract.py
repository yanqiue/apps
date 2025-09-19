# app.py
import streamlit as st
import pandas as pd

# ----------------------------
# Core parsing logic
# ----------------------------
def extract_sections(all_rows):
    """
    Returns:
      metadata_dict, metadata_block_text, df1, df2, df3, df2_plus_df3
    """

    # ---- Normalize all rows to trimmed strings ---- 
    normalized_rows=[]
    for row in all_rows:
        normalized_strings=[]
        for row_string in row:
            
            if row_string is None:
                normalized_str=""
            else:
                normalized_str = str(row_string).strip()

            normalized_strings.append(normalized_str)
        normalized_rows.append(normalized_strings)


    # ---- Find Table 1 header ------
    idx_table1_header = None
    for index, row in enumerate(normalized_rows):
        if any("taxon" in cell.lower() for cell in row) and any("station" in cell.lower() for cell in row): # Substring matching, not just an exact match.case-insensitive
            idx_table1_header = index
            break
        
    if idx_table1_header is None:
        raise ValueError("Could not find Table 1 header row containing both 'Taxon' and 'Station'.")
    
    else:
        # ---- Metadata handling: everything above Table 1 header ----
        first_non_empty_mdata_row_idx=None
        potentail_metadata_rows = normalized_rows[:idx_table1_header]

        for idx, row in enumerate(potentail_metadata_rows):
            if any(cell.strip() != "" for cell in row): # Row is not empty
                first_non_empty_mdata_row_idx=idx
                break

        # ---- Create a metadata Dictionary ----
        if first_non_empty_mdata_row_idx is not None:
            """
            Interpret a single line like: key1, value1, key2, value2, ...
            Any trailing odd item without a pair is ignored.
            """
            metadata_row=potentail_metadata_rows[first_non_empty_mdata_row_idx]
            cells=[c for c in metadata_row]
            
            metadata_dict = {}
            for i in range(0, len(cells) - 1, 2): #range(start, stop, step) Go up to the second-to-last index (to avoid IndexError when accessing i + 1).Increment by 2 each time
                k = cells[i]
                v = cells[i + 1]
                if k:
                    metadata_dict[k] = v
    

    # ---- Extract Table 1 ----
    header1 = normalized_rows[idx_table1_header]
    # Build a mapping of col name -> index for Table 1
    idx_taxon = None
    for j, col in enumerate(header1):
        if col.lower() == "taxon":
            idx_taxon = j
            break
    if idx_taxon is None:
        raise ValueError("Table 1 header located, but no 'Taxon' column was found.")

    # Data for Table 1 goes until BEFORE a line containing both 'Phyto' and 'Diversity:'
    data_table1 = []
    i = idx_table1_header + 1
    while i < len(normalized_rows):
        table1_row = normalized_rows[i]
        if any("phyto" in cell.lower() for cell in table1_row) and any("diversity" in cell.lower() for cell in table1_row): # Substring matching, not just an exact match.case-insensitive
            break  # stop before this row
        data_table1.append(table1_row)
        i += 1

    # Keep only rows where Taxon value is present
    data_table1_filtered = [r for r in data_table1 if idx_taxon < len(r) and r[idx_taxon]!=""]

    # Convert Table 1 to DataFrame (align widths) _ we want to ensure the same number of columns for each row
    width1 = max(len(header1), *(len(r) for r in data_table1_filtered)) if data_table1_filtered else len(header1) #Determine the widest row (i.e., max number of columns).
    header1 = (header1 + [""] * (width1 - len(header1)))[:width1] #Pad the header with empty strings ("") if it's too short.
    data_table1_filtered = [(r + [""] * (width1 - len(r)))[:width1] for r in data_table1_filtered] #Do the same padding for each data row.
    df1 = pd.DataFrame(data_table1_filtered, columns=header1)

    # st.markdown("##### First data table")
    # st.write(df1)

    # ---- Extract Table 2 ----
    # Look for header row like: ["", "mg/m^3", "%", " units or Cells/L", "%"]
    # We'll match loosely: first cell empty-ish, contains mg/m^3, exact "%" at 3rd and 5th positions (loose), one header has "units or Cells/L"
    idx_table2_header = None
    for k in range(i, len(normalized_rows)):
        table2_row = normalized_rows[k]  #Grab the current row to inspect.
        first_empty = (len(table2_row) > 0 and table2_row[0] == "") # Check if the first cell is empty-ish. This is a loose signal that it might be table 2 header.
        has_mg = any("mg/m^3" == c.lower() for c in table2_row if c) # Look for a cell that exactly matches "mg/m^3", case-insensitively:
        has_units_or_cells = any("units or cells/l" in c.lower() for c in table2_row) # Check if any cell contains the phrase "units or cells/l" (case-insensitive, substring match).
        count_percents = sum(1 for c in table2_row if c == "%") # Count how many cells are exactly equal to "%". You want at least two.
        if first_empty and has_mg and has_units_or_cells and count_percents >= 2: 
            idx_table2_header = k # Save the index of this header row.
            break


    if idx_table2_header is None:
        raise ValueError("Could not find Table 2 header row (empty, 'mg/m^3', '%', ' units or Cells/L', '%').")
    

    # Normalize Table 2 headers as required_terms
    # Target headers: ["Alga", "mg/m^3", "mg/m^3_%", "Cells-units/L", "Cells-units/L_%", "ratio"]
    # Note: we will add 'ratio' at the end; the table rows stop before the line containing '===='
    t2_headers_raw = normalized_rows[idx_table2_header]

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
    while r_idx < len(normalized_rows):
        table2_row = normalized_rows[r_idx]

        if any("====" in cell.lower() for cell in table2_row):
            break

        # Build a trimmed row with only the first 5 columns we care about
        trimmed = []

        # Map current row cells by walking raw header_map
        collected = {}
        for raw_idx, col_name in enumerate(header_map):
            if col_name in keep_headers:
                val = normalized_rows[r_idx][raw_idx] if raw_idx < len(normalized_rows[r_idx]) else ""
                collected[col_name] = val
        trimmed = [collected.get(h, "") for h in keep_headers] #For each header h in keep_headers, look up collected[h] — the value for that column.
                                                              # If h isn’t found (e.g., missing data), it returns "" as a fallback. Ensures The row is ordered according to keep_headers.
        # Append ratio placeholder as empty
        trimmed.append("")
        data2.append(trimmed)
        r_idx += 1

    df2 = pd.DataFrame(data2, columns=t2_headers_final)

    # st.markdown("##### Second data table")
    # st.write(df2)

    # ---- Skip lines after '====' until Table 3 header ----
    idx_table3_header = None
    for k in range(r_idx, len(normalized_rows)):
        table3_rows = normalized_rows[k]
        # Match header: ["", "mg/m^3", "Cells-units/L", "ratio"]
        if (len(table3_rows) >= 4 and
            table3_rows[0] == "" and
            any(c.lower() == "mg/m^3" for c in table3_rows) and
            any(c.lower() == "cells-units/l" for c in table3_rows) and
            any(c.lower() == "ratio" for c in table3_rows)):
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
    while k < len(normalized_rows):
        table3_rows = normalized_rows[k]

        def normalize_cell(s):
            if s is None:
                return ""
            return str(s).strip()
        
        if not any(normalize_cell(c) for c in table3_rows):
            break

        # Build mapping from the raw position:
        # Expect raw columns in order: ["", "mg/m^3", "Cells-units/L", "ratio"]
        alga = table3_rows[0] if len(table3_rows) > 0 else ""
        mg = table3_rows[1] if len(table3_rows) > 1 else ""
        cells_units = table3_rows[2] if len(table3_rows) > 2 else ""
        ratio = table3_rows[3] if len(table3_rows) > 3 else ""
        row_out = [alga, mg, "", cells_units, "", ratio]  # Insert the % columns as empty
        data3.append(row_out)
        k += 1

    df3 = pd.DataFrame(data3, columns=t3_target_headers)
    # st.markdown("##### Third data table")
    # st.write(df3)

    # ---- Concatenate Table 3 under Table 2 ----
    # Ensure same columns/order
    df2_aligned = df2[t3_target_headers] if list(df2.columns) != t3_target_headers else df2
    df2_plus_df3 = pd.concat([df2_aligned, df3], ignore_index=True)

    return metadata_dict, df1, df2, df3, df2_plus_df3
