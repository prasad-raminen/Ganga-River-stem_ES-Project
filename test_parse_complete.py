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

    # Replace placeholder symbols
    clean_text = re.sub(r'[_—–-]+', ' ', text)
    nums = re.findall(r'\b\d+\.?\d*\b', clean_text)
    if nums and nums[0] == st_code:
        nums = nums[1:]

    # 2012-2014 Format (Triplets: min, max, mean)
    if year in [2012, 2013, 2014]:
        if len(nums) >= 15:
            return {
                'station_code': st_code, 'station_name': st_name, 'year': year,
                'temp_min': float(nums[0]), 'temp_max': float(nums[1]),
                'do_min': float(nums[3]), 'do_max': float(nums[4]), 'do_mean': float(nums[5]),
                'ph_min': float(nums[6]), 'ph_max': float(nums[7]),
                'cond_min': float(nums[9]), 'cond_max': float(nums[10]),
                'bod_min': float(nums[12]), 'bod_max': float(nums[13]), 'bod_mean': float(nums[14]),
                'nitrate_min': float(nums[15]) if len(nums) > 15 else np.nan,
                'nitrate_max': float(nums[16]) if len(nums) > 16 else np.nan,
                'fecal_col_min': float(nums[18]) if len(nums) > 18 else np.nan,
                'fecal_col_max': float(nums[19]) if len(nums) > 19 else np.nan
            }

    # 2017 Haridwar special
    if '10148' in text and year == 2017:
        return {
            'station_code': st_code, 'station_name': st_name, 'year': year,
            'temp_min': 15.0, 'temp_max': 24.0,
            'do_min': 8.6, 'do_max': 9.8, 'do_mean': 9.2,
            'ph_min': 7.0, 'ph_max': 7.6,
            'cond_min': np.nan, 'cond_max': np.nan,
            'bod_min': 1.0, 'bod_max': 1.0, 'bod_mean': 1.0,
            'nitrate_min': np.nan, 'nitrate_max': np.nan,
            'fecal_col_min': 70.0, 'fecal_col_max': 170.0
        }

    # 2017 Gandhi Ghat (missing nitrate)
    if '2552' in text and year == 2017:
        return {
            'station_code': st_code, 'station_name': st_name, 'year': year,
            'temp_min': 15.0, 'temp_max': 33.0,
            'do_min': 6.5, 'do_max': 9.2, 'do_mean': 7.85,
            'ph_min': 7.3, 'ph_max': 8.8,
            'cond_min': 206.0, 'cond_max': 644.0,
            'bod_min': 2.4, 'bod_max': 2.9, 'bod_mean': 2.65,
            'nitrate_min': np.nan, 'nitrate_max': np.nan,
            'fecal_col_min': 1700.0, 'fecal_col_max': 1700000.0
        }

    # 2018/2019 missing nitrate cases (len(nums) == 14: temp2, do2, ph2, cond2, bod2, fecal2, total2)
    if len(nums) == 14:
        vals = [float(x) for x in nums]
        return {
            'station_code': st_code, 'station_name': st_name, 'year': year,
            'temp_min': vals[0], 'temp_max': vals[1],
            'do_min': vals[2], 'do_max': vals[3], 'do_mean': (vals[2]+vals[3])/2,
            'ph_min': vals[4], 'ph_max': vals[5],
            'cond_min': vals[6], 'cond_max': vals[7],
            'bod_min': vals[8], 'bod_max': vals[9], 'bod_mean': (vals[8]+vals[9])/2,
            'nitrate_min': np.nan, 'nitrate_max': np.nan,
            'fecal_col_min': vals[10], 'fecal_col_max': vals[11]
        }

    # Standard 16 numbers
    if len(nums) >= 16:
        vals = [float(x) for x in nums[-16:]]
        return {
            'station_code': st_code, 'station_name': st_name, 'year': year,
            'temp_min': vals[0], 'temp_max': vals[1],
            'do_min': vals[2], 'do_max': vals[3], 'do_mean': (vals[2]+vals[3])/2,
            'ph_min': vals[4], 'ph_max': vals[5],
            'cond_min': vals[6], 'cond_max': vals[7],
            'bod_min': vals[8], 'bod_max': vals[9], 'bod_mean': (vals[8]+vals[9])/2,
            'nitrate_min': vals[10], 'nitrate_max': vals[11],
            'fecal_col_min': vals[12], 'fecal_col_max': vals[13]
        }
    return None

parsed = [parse_record(r) for _, r in df_raw.iterrows() if parse_record(r) is not None]
df_p = pd.DataFrame(parsed)
print(f"Total parsed: {len(df_p)}")
print(df_p.groupby('station_name')['year'].count())
