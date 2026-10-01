import re
import pandas as pd

df_raw = pd.read_csv('data/processed/wq_raw_extracted.csv')

# Drop false positives
valid_records = []
for idx, r in df_raw.iterrows():
    text = str(r['raw_line_text'])
    # Skip false positives
    if any(k in text for k in ['LAKSHDWEEP', 'HERO CYCLE', '1334 WEST BENGAL', '37 23 30 4.7 7.8 7.7 8.6 180 1071']):
        continue
    # Skip 2016 unparsed headers
    if r['year'] == 2016 and len(text.split()) <= 4:
        continue
    valid_records.append((idx, r))

print(f"Retained {len(valid_records)} true raw records out of {len(df_raw)}.")
