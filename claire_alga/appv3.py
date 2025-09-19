# app.py
import io
import csv
import streamlit as st
import pandas as pd
import os
import sys

# Set page config
st.set_page_config(page_title=None, page_icon="📖", layout="wide", initial_sidebar_state="expanded", menu_items=None)

#Output Path
path = os.path.dirname(os.path.abspath(__file__))
#Add Modules
sys.path.append(f'{path}/Modules') #adding the Modules directory to Python's search path at runtime.

#Module Imports for the different sections
import extract

st.title("Messy Algae CSV Parser")
st.caption("Reads messy algal data CSV line-by-line and extracts: Metadata, Table 1, Table 2, Table 3 (then Table2 ⊕ Table3).")

# ----------------------------
# Streamlit UI
# ----------------------------
with st.sidebar:
    st.header("Upload CSV")
    f = st.file_uploader("Choose a CSV file", type=["csv"])

if not f:
    st.info("Upload a CSV to begin. The app parses it line-by-line and builds the outputs.")
    st.stop()


# ----------------------------
# Read file as lines (no pandas yet)
# ----------------------------
try:
    raw_csv_as_bytes=f.read()
    # Try the provided encoding first; fallback to latin-1
    for enc in ("utf-8", "utf-8-sig", "latin-1"):
        try:
            text = raw_csv_as_bytes.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    reader = csv.reader(io.StringIO(text))
    list_of_rows= [list(row) for row in reader]

except Exception as e:
    st.error(f"Failed to read CSV: {e}")
    st.stop()


# ----------------------------
# Parse into sections
# ----------------------------
try:
    metadata_dict, df1, df2, df3, df2_plus_df3=extract.extract_sections(list_of_rows)

except Exception as e:
    st.error(f"Parsing error: {e}")
    st.stop()


# ----------------------------
# Display & Downloads
# ----------------------------
st.markdown(" ")
st.markdown(" ##### Metadata (key/value pairs from the first non-empty line above Table 1)")
col_md1, col_md2 = st.columns([1,1])
with col_md1:
    st.json(metadata_dict or {"(none)": ""})

# Download metadata
md_json = io.BytesIO(pd.Series(metadata_dict).to_json(indent=2).encode("utf-8"))
st.download_button("Download metadata.json", md_json, file_name="metadata.json", mime="application/json")

# Table 1
st.markdown(" ")
st.markdown(" ##### Table 1 - Header row has 'Taxon' and 'Station', stops before the row with 'Phyto' and  'Diversity'")
st.dataframe(df1, use_container_width=True)
csv1 = io.BytesIO(df1.to_csv(index=False).encode("utf-8"))
st.download_button("Download table1.csv", csv1, file_name="table1.csv", mime="text/csv")

# Table 2
st.markdown(" ")
st.markdown(" ##### Table 2 -  Normalized headers, trimmed of additional rows to the left")
st.caption("Stops before a line containing '===='. Header changes: empty→'Alga', 'units or Cells/L'→'Cells-units/L', '%'→'mg/m^3_%' and 'Cells-units/L_%'.")
st.dataframe(df2, use_container_width=True)
csv2 = io.BytesIO(df2.to_csv(index=False).encode("utf-8"))
st.download_button("Download table2.csv", csv2, file_name="table2.csv", mime="text/csv")

# Table 3
st.markdown(" ")
st.markdown(" ##### Table 3 -  Headers ['', 'mg/m^3', 'Cells-units/L', 'ratio'] → insert '% columns to match Table 2")
st.dataframe(df3, use_container_width=True)
csv3 = io.BytesIO(df3.to_csv(index=False).encode("utf-8"))
st.download_button("Download table3.csv", csv3, file_name="table3.csv", mime="text/csv")

# Combined
st.markdown(" ")
st.subheader("Combined: Table 2 ⊕ Table 3 (same schema)")
st.dataframe(df2_plus_df3, use_container_width=True)
csv23 = io.BytesIO(df2_plus_df3.to_csv(index=False).encode("utf-8"))
st.download_button("Download table2_plus_table3.csv", csv23, file_name="table2_plus_table3.csv", mime="text/csv")

# Notes & next steps
st.markdown(" ")
with st.expander("Notes & adjustments"):
    st.markdown("""
- **Ratio column**: For Table 2, this app adds the column but does **not** compute it.  
- **Stops & starts**:
  - Table 1 stops just before the first row that has both “Phyto” and “Diversity:” anywhere in the row.
  - Table 2 stops before the first row that contains “====”.
  - After “====”, the app skips rows until it finds the Table 3 header row.
""")
