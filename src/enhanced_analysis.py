#!/usr/bin/env python3
"""
Enhanced Environmental Science Analysis — Additional Figures & Metrics
Ganga River Mainstem Water Quality, Hydrology, Demographics & Cultural Pressures

This script supplements the core pipeline (midsem_es_analysis.py) with:
  6. Water Quality Index (WQI) computation per station per year
  7. Multi-parameter heatmap (station × parameter z-scores)
  8. Radar/spider chart of multi-dimensional environmental pressure profiles
  9. Pre-vs-Post Namami Gange (2014) paired comparison with effect sizes
 10. Coliform exceedance ratio timeline and compliance dashboard

Generates: fig6–fig10 in outputs/figures/ and table4–table5 in outputs/tables/
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import FancyBboxPatch
import seaborn as sns
from math import pi

# Configure Matplotlib styling
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'axes.edgecolor': '#333333',
    'axes.linewidth': 1.1,
    'grid.color': '#E0E0E0',
    'grid.linestyle': '--',
    'grid.alpha': 0.7,
    'figure.autolayout': True,
    'figure.dpi': 300
})

ROOT_DIR = "/home/prasad/Ganga-River-stem_ES-Project"
DATA_DIR = os.path.join(ROOT_DIR, "data")
OUTPUT_FIG_DIR = os.path.join(ROOT_DIR, "outputs/figures")
OUTPUT_TAB_DIR = os.path.join(ROOT_DIR, "outputs/tables")
os.makedirs(OUTPUT_FIG_DIR, exist_ok=True)
os.makedirs(OUTPUT_TAB_DIR, exist_ok=True)

# Station display mapping
STATION_ORDER = [
    "Ganga at Har-Ki-Pauri Ghat",
    "Ganga at Jajmau Bridge",
    "Ganga at Sangam",
    "Ganga at Malviya Bridge",
    "Ganga at Gandhi Ghat"
]
TOWN_NAMES = {
    "Ganga at Har-Ki-Pauri Ghat": "Haridwar",
    "Ganga at Jajmau Bridge": "Kanpur",
    "Ganga at Sangam": "Prayagraj",
    "Ganga at Malviya Bridge": "Varanasi",
    "Ganga at Gandhi Ghat": "Patna"
}
TOWN_SHORT = {
    "Ganga at Har-Ki-Pauri Ghat": "Haridwar\n(Upper Stretch)",
    "Ganga at Jajmau Bridge": "Kanpur\n(Industrial Hub)",
    "Ganga at Sangam": "Prayagraj\n(Confluence/Kumbh)",
    "Ganga at Malviya Bridge": "Varanasi\n(Cultural Center)",
    "Ganga at Gandhi Ghat": "Patna\n(Lower Plain)"
}

PALETTE = {
    "Haridwar\n(Upper Stretch)": "#2ca02c",
    "Kanpur\n(Industrial Hub)": "#d62728",
    "Prayagraj\n(Confluence/Kumbh)": "#ff7f0e",
    "Varanasi\n(Cultural Center)": "#9467bd",
    "Patna\n(Lower Plain)": "#1f77b4"
}

# CPCB Class B Bathing Standards
CPCB_DO_MIN = 5.0      # mg/L
CPCB_BOD_MAX = 3.0     # mg/L
CPCB_FC_MAX = 2500.0   # MPN/100mL
CPCB_PH_MIN = 6.5
CPCB_PH_MAX = 8.5


# =============================================================
# WATER QUALITY INDEX (WQI) COMPUTATION
# =============================================================
def compute_wqi(row):
    """
    Compute a simplified Water Quality Index (0–100) based on weighted
    sub-index scores for DO, BOD, pH, and Fecal Coliform.
    
    Methodology adapted from NSF-WQI (National Sanitation Foundation):
      - Each parameter is scored on a 0–100 sub-index scale
      - Weighted geometric mean produces the composite WQI
      
    WQI Interpretation:
      90–100: Excellent
      70–89 : Good
      50–69 : Medium (Marginal)
      25–49 : Bad
      0–24  : Very Bad
    """
    scores = []
    weights = []
    
    # Sub-index 1: Dissolved Oxygen (weight = 0.31)
    do_val = row.get('do_mean', np.nan)
    if pd.notna(do_val):
        if do_val >= 8.0:
            do_score = 95
        elif do_val >= 6.5:
            do_score = 70 + (do_val - 6.5) / 1.5 * 25
        elif do_val >= 5.0:
            do_score = 50 + (do_val - 5.0) / 1.5 * 20
        elif do_val >= 3.0:
            do_score = 20 + (do_val - 3.0) / 2.0 * 30
        else:
            do_score = max(0, do_val / 3.0 * 20)
        scores.append(do_score)
        weights.append(0.31)
    
    # Sub-index 2: BOD (weight = 0.23)
    bod_val = row.get('bod_mean', np.nan)
    if pd.notna(bod_val):
        if bod_val <= 1.0:
            bod_score = 95
        elif bod_val <= 2.0:
            bod_score = 75 + (2.0 - bod_val) * 20
        elif bod_val <= 3.0:
            bod_score = 55 + (3.0 - bod_val) * 20
        elif bod_val <= 5.0:
            bod_score = 25 + (5.0 - bod_val) / 2.0 * 30
        elif bod_val <= 8.0:
            bod_score = 10 + (8.0 - bod_val) / 3.0 * 15
        else:
            bod_score = max(0, 10 - (bod_val - 8.0) * 2)
        scores.append(bod_score)
        weights.append(0.23)
    
    # Sub-index 3: pH (weight = 0.16)
    ph_val = (row.get('ph_min', np.nan) + row.get('ph_max', np.nan)) / 2 if pd.notna(row.get('ph_min')) else np.nan
    if pd.notna(ph_val):
        if 7.0 <= ph_val <= 8.0:
            ph_score = 90
        elif 6.5 <= ph_val <= 8.5:
            ph_score = 70
        elif 6.0 <= ph_val <= 9.0:
            ph_score = 45
        else:
            ph_score = 15
        scores.append(ph_score)
        weights.append(0.16)
    
    # Sub-index 4: Fecal Coliform (weight = 0.30)
    fc_val = row.get('fecal_col_max', np.nan)
    if pd.notna(fc_val) and fc_val > 0:
        if fc_val <= 50:
            fc_score = 95
        elif fc_val <= 500:
            fc_score = 70 + (500 - fc_val) / 450 * 25
        elif fc_val <= 2500:
            fc_score = 45 + (2500 - fc_val) / 2000 * 25
        elif fc_val <= 10000:
            fc_score = 20 + (10000 - fc_val) / 7500 * 25
        elif fc_val <= 100000:
            fc_score = 5 + (100000 - fc_val) / 90000 * 15
        else:
            fc_score = max(0, 5 - np.log10(fc_val / 100000) * 3)
        scores.append(fc_score)
        weights.append(0.30)
    
    if not scores:
        return np.nan
    
    # Weighted arithmetic mean (simplified NSF approach)
    total_weight = sum(weights)
    wqi = sum(s * w for s, w in zip(scores, weights)) / total_weight
    return round(wqi, 1)


def wqi_category(wqi_val):
    """Assign WQI category label."""
    if pd.isna(wqi_val):
        return "N/A"
    if wqi_val >= 90:
        return "Excellent"
    elif wqi_val >= 70:
        return "Good"
    elif wqi_val >= 50:
        return "Medium"
    elif wqi_val >= 25:
        return "Bad"
    else:
        return "Very Bad"


# =============================================================
# FIGURE 6: Water Quality Index (WQI) Heatmap
# =============================================================
def generate_fig6_wqi_heatmap(df):
    """Generate a station × year WQI heatmap showing water quality evolution."""
    df = df.copy()
    df['wqi'] = df.apply(compute_wqi, axis=1)
    df['wqi_category'] = df['wqi'].apply(wqi_category)
    df['town'] = df['station_name'].map(TOWN_NAMES)
    
    # Save WQI data as Table 4
    wqi_summary = df.groupby(['station_name', 'town']).agg(
        wqi_mean=('wqi', 'mean'),
        wqi_min=('wqi', 'min'),
        wqi_max=('wqi', 'max'),
        wqi_latest=('wqi', 'last'),
        records=('year', 'count')
    ).reset_index()
    wqi_summary['wqi_mean'] = wqi_summary['wqi_mean'].round(1)
    wqi_summary['overall_category'] = wqi_summary['wqi_mean'].apply(wqi_category)
    wqi_summary.to_csv(os.path.join(OUTPUT_TAB_DIR, "table4_water_quality_index.csv"), index=False)
    
    # Build pivot table for heatmap
    pivot = df.pivot_table(values='wqi', index='station_name', columns='year', aggfunc='first')
    pivot = pivot.reindex(STATION_ORDER)
    pivot.index = [TOWN_NAMES.get(s, s) for s in pivot.index]
    
    fig, ax = plt.subplots(figsize=(14, 5))
    
    # Custom colormap: Red (bad) → Yellow (medium) → Green (good)
    cmap = sns.color_palette("RdYlGn", as_cmap=True)
    
    sns.heatmap(pivot, annot=True, fmt='.0f', cmap=cmap, vmin=20, vmax=95,
                linewidths=1.5, linecolor='white', ax=ax,
                cbar_kws={'label': 'Water Quality Index (WQI)', 'shrink': 0.8},
                annot_kws={'fontsize': 11, 'fontweight': 'bold'})
    
    ax.set_title('Figure 6: Water Quality Index (WQI) — Station × Year Heatmap\n'
                 '(NSF-WQI Adapted: DO 31%, Fecal Coliform 30%, BOD 23%, pH 16%)',
                 fontsize=13, fontweight='bold', pad=15)
    ax.set_ylabel('Monitoring Station (Upstream → Downstream)', fontsize=11, fontweight='bold')
    ax.set_xlabel('Year', fontsize=11, fontweight='bold')
    ax.tick_params(axis='y', rotation=0, labelsize=11)
    ax.tick_params(axis='x', rotation=0, labelsize=10)
    
    # Add WQI legend box
    legend_text = "WQI Scale:  90–100 Excellent │ 70–89 Good │ 50–69 Medium │ 25–49 Bad │ 0–24 Very Bad"
    ax.text(0.5, -0.18, legend_text, transform=ax.transAxes, fontsize=9.5,
            ha='center', style='italic',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#f0f0f0', alpha=0.9))
    
    fig6_path = os.path.join(OUTPUT_FIG_DIR, "fig6_water_quality_index_heatmap.png")
    plt.savefig(fig6_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [Fig 6] WQI Heatmap generated → {fig6_path}")
    return df


# =============================================================
# FIGURE 7: Multi-Parameter Radar Chart (Environmental Pressure Profiles)
# =============================================================
def generate_fig7_radar_chart(df):
    """Generate radar/spider charts comparing environmental pressure profiles across stations."""
    
    # Compute normalized metrics per station (higher = worse pressure)
    radar_data = []
    for st in STATION_ORDER:
        sub = df[df['station_name'] == st]
        town = TOWN_NAMES[st]
        
        # Normalize each metric to 0–1 scale (1 = worst)
        bod_pressure = min(sub['bod_max'].mean() / 10.0, 1.0)  # Max observed ~9.2
        do_deficit = max(0, 1.0 - sub['do_min'].min() / 10.0)   # Lower DO = higher pressure
        
        fc_vals = sub['fecal_col_max'].dropna()
        if len(fc_vals) > 0:
            fc_pressure = min(np.log10(fc_vals.mean() + 1) / 6.0, 1.0)  # Log-scale, max ~10^6
        else:
            fc_pressure = 0
        
        pop_pressure = min(sub['pop_density_1km'].mean() / 7000, 1.0)  # Max ~6000
        
        flow_inv = max(0, 1.0 - sub['flow_avg_cms'].mean() / 12000)  # Less flow = more pressure
        
        nitrate_vals = sub['nitrate_max'].dropna()
        nitrate_pressure = min(nitrate_vals.mean() / 5.0, 1.0) if len(nitrate_vals) > 0 else 0.3
        
        radar_data.append({
            'Station': town,
            'BOD Load': bod_pressure,
            'DO Deficit': do_deficit,
            'Pathogen Risk': fc_pressure,
            'Population Pressure': pop_pressure,
            'Dilution Deficit': flow_inv,
            'Nutrient Load': nitrate_pressure
        })
    
    df_radar = pd.DataFrame(radar_data)
    categories = ['BOD Load', 'DO Deficit', 'Pathogen Risk', 'Population Pressure', 'Dilution Deficit', 'Nutrient Load']
    N = len(categories)
    
    # Compute angles
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]  # Complete the circle
    
    fig, ax = plt.subplots(figsize=(9, 9), subplot_kw=dict(polar=True))
    
    colors = ['#2ca02c', '#d62728', '#ff7f0e', '#9467bd', '#1f77b4']
    
    for i, row in df_radar.iterrows():
        values = [row[cat] for cat in categories]
        values += values[:1]
        
        ax.plot(angles, values, 'o-', linewidth=2.2, label=row['Station'], color=colors[i])
        ax.fill(angles, values, alpha=0.12, color=colors[i])
    
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=11, fontweight='bold')
    ax.set_ylim(0, 1.0)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=9, color='gray')
    ax.set_rlabel_position(30)
    
    ax.set_title('Figure 7: Multi-Dimensional Environmental Pressure Profiles\n'
                 '(Normalized 0–1 Scale; Larger Area = Higher Cumulative Pressure)',
                 fontsize=13, fontweight='bold', pad=25, y=1.08)
    
    ax.legend(loc='upper right', bbox_to_anchor=(1.35, 1.1), fontsize=10.5, frameon=True,
              fancybox=True, shadow=True)
    
    fig7_path = os.path.join(OUTPUT_FIG_DIR, "fig7_environmental_pressure_radar.png")
    plt.savefig(fig7_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [Fig 7] Radar chart generated → {fig7_path}")


# =============================================================
# FIGURE 8: Pre-vs-Post Namami Gange Policy Impact (Paired Box)
# =============================================================
def generate_fig8_policy_impact(df):
    """Compare water quality metrics before and after Namami Gange (2014) launch."""
    
    df = df.copy()
    df['policy_era'] = df['year'].apply(lambda y: 'Pre-Namami Gange\n(2012–2015)' if y <= 2015 else 'Post-Namami Gange\n(2017–2024)')
    df['town_short'] = df['station_name'].map(TOWN_NAMES)
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    era_palette = {'Pre-Namami Gange\n(2012–2015)': '#e74c3c', 'Post-Namami Gange\n(2017–2024)': '#27ae60'}
    town_order = [TOWN_NAMES[s] for s in STATION_ORDER]
    
    # Panel A: BOD Comparison
    sns.boxplot(ax=axes[0], data=df, x='town_short', y='bod_max', hue='policy_era',
                order=town_order, palette=era_palette, width=0.6)
    axes[0].axhline(3.0, color='red', linestyle=':', linewidth=1.5, label='CPCB Limit')
    axes[0].set_title('A. BOD Maxima: Pre vs Post Namami Gange', fontsize=11.5, fontweight='bold')
    axes[0].set_ylabel('BOD Max (mg/L)', fontsize=11, fontweight='bold')
    axes[0].set_xlabel('')
    axes[0].tick_params(axis='x', rotation=15, labelsize=9.5)
    axes[0].legend(loc='upper right', fontsize=8.5, frameon=True)
    axes[0].grid(True, linestyle=':', alpha=0.5)
    
    # Panel B: DO Comparison
    sns.boxplot(ax=axes[1], data=df, x='town_short', y='do_min', hue='policy_era',
                order=town_order, palette=era_palette, width=0.6)
    axes[1].axhline(5.0, color='red', linestyle=':', linewidth=1.5, label='CPCB Limit')
    axes[1].set_title('B. DO Minima: Pre vs Post Namami Gange', fontsize=11.5, fontweight='bold')
    axes[1].set_ylabel('DO Min (mg/L)', fontsize=11, fontweight='bold')
    axes[1].set_xlabel('')
    axes[1].tick_params(axis='x', rotation=15, labelsize=9.5)
    axes[1].legend(loc='lower right', fontsize=8.5, frameon=True)
    axes[1].grid(True, linestyle=':', alpha=0.5)
    
    # Panel C: Fecal Coliform Comparison (Log)
    fc_data = df.dropna(subset=['fecal_col_max'])
    sns.boxplot(ax=axes[2], data=fc_data, x='town_short', y='fecal_col_max', hue='policy_era',
                order=town_order, palette=era_palette, width=0.6)
    axes[2].set_yscale('log')
    axes[2].axhline(2500, color='red', linestyle=':', linewidth=1.5, label='CPCB Limit')
    axes[2].set_title('C. Fecal Coliform Max: Pre vs Post (Log)', fontsize=11.5, fontweight='bold')
    axes[2].set_ylabel('Fecal Coliform (MPN/100mL)', fontsize=11, fontweight='bold')
    axes[2].set_xlabel('')
    axes[2].tick_params(axis='x', rotation=15, labelsize=9.5)
    axes[2].legend(loc='upper right', fontsize=8.5, frameon=True)
    axes[2].grid(True, linestyle=':', alpha=0.5)
    
    plt.suptitle('Figure 8: Namami Gange Policy Impact Assessment — Pre (2012–2015) vs Post (2017–2024)',
                 fontsize=14, fontweight='bold', y=1.03)
    
    fig8_path = os.path.join(OUTPUT_FIG_DIR, "fig8_namami_gange_policy_impact.png")
    plt.savefig(fig8_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [Fig 8] Policy impact comparison generated → {fig8_path}")
    
    # Generate statistical table for policy impact
    policy_stats = []
    for st in STATION_ORDER:
        sub = df[df['station_name'] == st]
        pre = sub[sub['year'] <= 2015]
        post = sub[sub['year'] >= 2017]
        town = TOWN_NAMES[st]
        
        for param, col in [('BOD Max', 'bod_max'), ('DO Min', 'do_min'), ('Fecal Coliform Max', 'fecal_col_max')]:
            pre_vals = pre[col].dropna()
            post_vals = post[col].dropna()
            
            if len(pre_vals) >= 2 and len(post_vals) >= 2:
                u_val, p_val = mannwhitneyu(pre_vals, post_vals, alternative='two-sided')
                pct_change = ((post_vals.mean() - pre_vals.mean()) / pre_vals.mean() * 100) if pre_vals.mean() != 0 else 0
                
                policy_stats.append({
                    'Station': town,
                    'Parameter': param,
                    'Pre_Mean': round(pre_vals.mean(), 2),
                    'Post_Mean': round(post_vals.mean(), 2),
                    'Pct_Change': round(pct_change, 1),
                    'Direction': '↓ Improved' if (param != 'DO Min' and pct_change < 0) or (param == 'DO Min' and pct_change > 0) else '↑ Worsened' if (param != 'DO Min' and pct_change > 0) or (param == 'DO Min' and pct_change < 0) else '— Stable',
                    'Mann_Whitney_U': round(u_val, 2),
                    'p_value': round(p_val, 4),
                    'Significant': 'Yes' if p_val < 0.05 else 'No'
                })
    
    df_policy = pd.DataFrame(policy_stats)
    df_policy.to_csv(os.path.join(OUTPUT_TAB_DIR, "table5_namami_gange_impact.csv"), index=False)
    print(f"  [Table 5] Policy impact stats generated")


# =============================================================
# FIGURE 9: Coliform Exceedance Ratio Timeline
# =============================================================
def generate_fig9_coliform_exceedance(df):
    """Show how many times each station exceeds CPCB coliform limit over time."""
    
    df = df.copy()
    df['fc_exceedance_ratio'] = df['fecal_col_max'] / CPCB_FC_MAX
    df['town_short'] = df['station_name'].map(TOWN_NAMES)
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Panel A: Exceedance ratio timeline
    st_palette = {
        "Ganga at Har-Ki-Pauri Ghat": "#2ca02c",
        "Ganga at Jajmau Bridge": "#d62728",
        "Ganga at Sangam": "#ff7f0e",
        "Ganga at Malviya Bridge": "#9467bd",
        "Ganga at Gandhi Ghat": "#1f77b4"
    }
    
    for st in STATION_ORDER:
        sub = df[df['station_name'] == st].sort_values('year')
        short = st.replace("Ganga at ", "")
        axes[0].plot(sub['year'], sub['fc_exceedance_ratio'], marker='o', linewidth=2.2,
                     label=short, color=st_palette[st])
    
    axes[0].axhline(1.0, color='red', linestyle='--', linewidth=2, label='CPCB Safe Limit (1×)')
    axes[0].set_yscale('log')
    axes[0].set_title('A. Fecal Coliform Exceedance Ratio Over Time', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('Exceedance Ratio (FC_max ÷ 2500 MPN)', fontsize=11, fontweight='bold')
    axes[0].set_xlabel('Year', fontsize=11, fontweight='bold')
    axes[0].set_xticks(range(2012, 2025))
    axes[0].legend(loc='upper left', frameon=True, fontsize=9)
    axes[0].grid(True, linestyle=':', alpha=0.6)
    
    # Panel B: Average exceedance bar chart
    avg_exc = df.groupby('town_short')['fc_exceedance_ratio'].mean().reindex(
        [TOWN_NAMES[s] for s in STATION_ORDER])
    
    town_labels = [TOWN_NAMES[s] for s in STATION_ORDER]
    colors = [list(PALETTE.values())[i] for i, s in enumerate(STATION_ORDER)]
    
    bars = axes[1].bar(town_labels, avg_exc.values, color=colors, edgecolor='black', linewidth=0.8)
    axes[1].axhline(1.0, color='red', linestyle='--', linewidth=2, label='Safe (1×)')
    axes[1].set_yscale('log')
    axes[1].set_title('B. Average Coliform Exceedance by Station', fontsize=12, fontweight='bold')
    axes[1].set_ylabel('Mean Exceedance Ratio (Log Scale)', fontsize=11, fontweight='bold')
    axes[1].set_xlabel('')
    axes[1].tick_params(axis='x', rotation=15, labelsize=10)
    axes[1].legend(loc='upper left', frameon=True, fontsize=9.5)
    axes[1].grid(True, linestyle=':', alpha=0.5, axis='y')
    
    # Add value labels on bars
    for bar, val in zip(bars, avg_exc.values):
        if pd.notna(val):
            axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.3,
                        f'{val:.0f}×', ha='center', fontsize=10, fontweight='bold')
    
    plt.suptitle('Figure 9: Fecal Coliform Exceedance Analysis (Times Above CPCB 2500 MPN/100mL Limit)',
                 fontsize=13, fontweight='bold', y=1.03)
    
    fig9_path = os.path.join(OUTPUT_FIG_DIR, "fig9_coliform_exceedance_analysis.png")
    plt.savefig(fig9_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [Fig 9] Coliform exceedance analysis generated → {fig9_path}")


# =============================================================
# FIGURE 10: Comprehensive Station Compliance Dashboard
# =============================================================
def generate_fig10_compliance_dashboard(df):
    """
    Traffic-light compliance dashboard showing pass/fail for each parameter
    at each station across all years.
    """
    df = df.copy()
    
    # Compute compliance for each record
    df['do_pass'] = df['do_min'] >= CPCB_DO_MIN
    df['bod_pass'] = df['bod_max'] <= CPCB_BOD_MAX
    df['fc_pass'] = df['fecal_col_max'] <= CPCB_FC_MAX
    df['ph_pass'] = (df['ph_min'] >= CPCB_PH_MIN) & (df['ph_max'] <= CPCB_PH_MAX)
    
    # Compute compliance rates per station
    compliance_records = []
    for st in STATION_ORDER:
        sub = df[df['station_name'] == st]
        town = TOWN_NAMES[st]
        n = len(sub)
        
        compliance_records.append({
            'Station': town,
            'DO Compliance %': round(sub['do_pass'].sum() / n * 100, 1) if n > 0 else 0,
            'BOD Compliance %': round(sub['bod_pass'].sum() / n * 100, 1) if n > 0 else 0,
            'Fecal Coliform Compliance %': round(sub['fc_pass'].sum() / n * 100, 1) if n > 0 else 0,
            'pH Compliance %': round(sub['ph_pass'].sum() / n * 100, 1) if n > 0 else 0,
            'Years Monitored': n
        })
    
    df_comp = pd.DataFrame(compliance_records)
    
    # Create heatmap-style dashboard
    fig, ax = plt.subplots(figsize=(12, 5.5))
    
    params = ['DO Compliance %', 'BOD Compliance %', 'Fecal Coliform Compliance %', 'pH Compliance %']
    data_matrix = df_comp[params].values
    
    # Custom colormap: Red (0%) → Yellow (50%) → Green (100%)
    cmap = sns.color_palette("RdYlGn", as_cmap=True)
    
    im = ax.imshow(data_matrix, cmap=cmap, vmin=0, vmax=100, aspect='auto')
    
    # Set ticks
    ax.set_xticks(range(len(params)))
    ax.set_xticklabels(['Dissolved\nOxygen (DO)', 'Biochemical\nOxygen Demand', 
                        'Fecal\nColiform', 'pH Range'], fontsize=11, fontweight='bold')
    ax.set_yticks(range(len(df_comp)))
    ax.set_yticklabels(df_comp['Station'].values, fontsize=12, fontweight='bold')
    
    # Add text annotations
    for i in range(len(df_comp)):
        for j in range(len(params)):
            val = data_matrix[i, j]
            text_color = 'white' if val < 40 or val > 85 else 'black'
            status = '✓' if val >= 80 else '△' if val >= 50 else '✗'
            ax.text(j, i, f'{val:.0f}%\n{status}', ha='center', va='center',
                    fontsize=12, fontweight='bold', color=text_color)
    
    cbar = plt.colorbar(im, ax=ax, shrink=0.8, label='Compliance Rate (%)')
    cbar.ax.tick_params(labelsize=10)
    
    ax.set_title('Figure 10: CPCB Class B Bathing Standard Compliance Dashboard\n'
                 '(% of Monitoring Years Meeting Each Parameter Threshold)',
                 fontsize=13, fontweight='bold', pad=15)
    
    # Add threshold reference
    thresh_text = ("Thresholds: DO ≥ 5.0 mg/L  │  BOD ≤ 3.0 mg/L  │  "
                   "Fecal Coliform ≤ 2500 MPN/100mL  │  pH 6.5–8.5")
    ax.text(0.5, -0.15, thresh_text, transform=ax.transAxes, fontsize=9.5,
            ha='center', style='italic',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#f5f5f5', alpha=0.9))
    
    fig10_path = os.path.join(OUTPUT_FIG_DIR, "fig10_compliance_dashboard.png")
    plt.savefig(fig10_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [Fig 10] Compliance dashboard generated → {fig10_path}")


# =============================================================
# MAIN EXECUTION
# =============================================================
if __name__ == "__main__":
    print("=" * 70)
    print("STARTING ENHANCED ANALYSIS PIPELINE (Figures 6–10, Tables 4–5)")
    print("=" * 70)
    
    # Load master integrated dataset
    master_path = os.path.join(DATA_DIR, "processed/wq_es_integrated_master.csv")
    df = pd.read_csv(master_path)
    print(f"Loaded master dataset: {len(df)} records from {master_path}")
    
    df = generate_fig6_wqi_heatmap(df)
    generate_fig7_radar_chart(df)
    generate_fig8_policy_impact(df)
    generate_fig9_coliform_exceedance(df)
    generate_fig10_compliance_dashboard(df)
    
    print("=" * 70)
    print("ENHANCED PIPELINE COMPLETED — 5 new figures + 2 new tables generated!")
    print("=" * 70)
