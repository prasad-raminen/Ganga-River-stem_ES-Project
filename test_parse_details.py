import re
import pandas as pd
import numpy as np

df_raw = pd.read_csv('data/processed/wq_raw_extracted.csv')

def parse_record(row):
    text = str(row['raw_line_text'])
    year = int(row['year'])
    st_code = str(row['station_code'])
    st_name = str(row['station_name'])
    
    # Filter false positives
    if any(k in text for k in ['LAKSHDWEEP', 'HERO CYCLE', '1334 WEST BENGAL', '37 23 30 4.7 7.8 7.7 8.6 180 1071']):
        return None
    if year == 2016 and len(text.split()) <= 4:
        return None
        
    # Replace non-numeric placeholders like '_', '-', 'BDL' with NaN-safe strings
    tokens = text.split()
    
    # Check 2012-2014 format: 3 values per parameter (Min, Max, Mean)
    if year in [2012, 2013, 2014]:
        # Extract all floats/ints
        nums = re.findall(r'\b\d+\.?\d*\b', text)
        if nums and nums[0] == st_code:
            nums = nums[1:]
        # In 2012-2014: Temp(3), DO(3), pH(3), Cond(3), BOD(3), Nitrate(3), Coliform(2 or 3)
        # e.g. Temp: nums[0:3], DO: nums[3:6], pH: nums[6:9], Cond: nums[9:12], BOD: nums[12:15], Nitrate: nums[15:18]
        if len(nums) >= 15:
            temp_min, temp_max, temp_mean = float(nums[0]), float(nums[1]), float(nums[2])
            do_min, do_max, do_mean = float(nums[3]), float(nums[4]), float(nums[5])
            ph_min, ph_max, ph_mean = float(nums[6]), float(nums[7]), float(nums[8])
            cond_min, cond_max, cond_mean = float(nums[9]), float(nums[10]), float(nums[11])
            bod_min, bod_max, bod_mean = float(nums[12]), float(nums[13]), float(nums[14])
            nitrate_min = float(nums[15]) if len(nums) > 15 else np.nan
            nitrate_max = float(nums[16]) if len(nums) > 16 else np.nan
            nitrate_mean = float(nums[17]) if len(nums) > 17 else np.nan
            fecal_min = float(nums[18]) if len(nums) > 18 else np.nan
            fecal_max = float(nums[19]) if len(nums) > 19 else np.nan
            return {
                'station_code': st_code, 'station_name': st_name, 'year': year,
                'do_min': do_min, 'do_max': do_max, 'do_mean': do_mean,
                'bod_min': bod_min, 'bod_max': bod_max, 'bod_mean': bod_mean,
                'ph_min': ph_min, 'ph_max': ph_max,
                'cond_min': cond_min, 'cond_max': cond_max,
                'nitrate_min': nitrate_min, 'nitrate_max': nitrate_max,
                'fecal_col_min': fecal_min, 'fecal_col_max': fecal_max
            }
    else:
        # 2015+ format: min and max pairs
        # Temp(2), DO(2), pH(2), Cond(2), BOD(2), Nitrate(2), Fecal(2), Total(2)
        # Handle cases where some values are '_' or missing
        # Special case for Haridwar 2017 (row 32): 10148 UTTARAKHAND 15 24 8.6 9.8 7 7.6 1 1 70 170
        # Temp: 15, 24; DO: 8.6, 9.8; pH: 7, 7.6; BOD: 1, 1; Coliform: 70, 170
        if '10148' in text and year == 2017:
            return {
                'station_code': st_code, 'station_name': st_name, 'year': year,
                'do_min': 8.6, 'do_max': 9.8, 'do_mean': 9.2,
                'bod_min': 1.0, 'bod_max': 1.0, 'bod_mean': 1.0,
                'ph_min': 7.0, 'ph_max': 7.6,
                'cond_min': np.nan, 'cond_max': np.nan,
                'nitrate_min': np.nan, 'nitrate_max': np.nan,
                'fecal_col_min': 70.0, 'fecal_col_max': 170.0
            }
        # Haridwar 2018 (row 6): 10148 UTTARAKHAND 16.0 21.0 8.8 10.0 7.4 8.4 _ _ 1.0 1.0 _ _ _ _ 60 170
        if '10148' in text and year == 2018:
            return {
                'station_code': st_code, 'station_name': st_name, 'year': year,
                'do_min': 8.8, 'do_max': 10.0, 'do_mean': 9.4,
                'bod_min': 1.0, 'bod_max': 1.0, 'bod_mean': 1.0,
                'ph_min': 7.4, 'ph_max': 8.4,
                'cond_min': np.nan, 'cond_max': np.nan,
                'nitrate_min': np.nan, 'nitrate_max': np.nan,
                'fecal_col_min': np.nan, 'fecal_col_max': np.nan,
                'total_col_min': 60.0, 'total_col_max': 170.0
            }
        # Haridwar 2019 (row 11): 10148 UTTARAKHAND 17 22 9.0 10.2 7.5 8.4 121 155 1.0 1.0 __ __ 26 50 40 110 - -
        if '10148' in text and year == 2019:
            return {
                'station_code': st_code, 'station_name': st_name, 'year': year,
                'do_min': 9.0, 'do_max': 10.2, 'do_mean': 9.6,
                'bod_min': 1.0, 'bod_max': 1.0, 'bod_mean': 1.0,
                'ph_min': 7.5, 'ph_max': 8.4,
                'cond_min': 121.0, 'cond_max': 155.0,
                'nitrate_min': np.nan, 'nitrate_max': np.nan,
                'fecal_col_min': 26.0, 'fecal_col_max': 50.0
            }
        # Standard 2015-2024 lines with >= 16 numbers
        nums = re.findall(r'\b\d+\.?\d*\b', text)
        if nums and nums[0] == st_code:
            nums = nums[1:]
        if len(nums) >= 16:
            vals = [float(x) for x in nums[-16:]]
            return {
                'station_code': st_code, 'station_name': st_name, 'year': year,
                'do_min': vals[2], 'do_max': vals[3], 'do_mean': (vals[2]+vals[3])/2,
                'bod_min': vals[8], 'bod_max': vals[9], 'bod_mean': (vals[8]+vals[9])/2,
                'ph_min': vals[4], 'ph_max': vals[5],
                'cond_min': vals[6], 'cond_max': vals[7],
                'nitrate_min': vals[10], 'nitrate_max': vals[11],
                'fecal_col_min': vals[12], 'fecal_col_max': vals[13],
                'total_col_min': vals[14], 'total_col_max': vals[15]
            }
        else:
            # Let's print any other rows that didn't match
            print(f"Skipped in 2015+: {year} {st_name}: {text}")
            return None

parsed = []
for idx, r in df_raw.iterrows():
    res = parse_record(r)
    if res:
        parsed.append(res)

df_p = pd.DataFrame(parsed)
print(f"\nParsed {len(df_p)} records successfully!")
print(df_p.groupby('station_name')['year'].count())
print(df_p.head(5))
