"""
Plot 2 for infant vaccination scenarios
"""

import pandas as pd
import sciris as sc
from run_scenarios import coverage_array, target_age_list, mac_list
import utils as ut
import seaborn as sns
import numpy as np


def preprocess_data(msim_dict):

    # What to store
    start_year = 2025
    metrics = ['cancers', 'cancer_deaths']
    records = sc.autolist()
    base_label = 'Baseline'
    si = sc.findinds(msim_dict[base_label].year, start_year)[0]

    for target_age in target_age_list:
        age_label = f'Girls {target_age[0]}-{target_age[1]}'

        for cn, cov_val in enumerate(coverage_array):
            scen_label = age_label + f': {np.round(cov_val, decimals=1)} coverage'

            for pn, metric in enumerate(metrics):
                base_vals = msim_dict[base_label][metric].values[si:]
                scen_vals = msim_dict[scen_label][metric].values[si:]
                n_averted = sum(base_vals - scen_vals)
                records += {'age': f'{target_age[0]}', 'coverage': int(round(cov_val, 1)*100), 'metric': f'{metric.replace("_"," ").capitalize()}', 'val': n_averted, 'mac': 'Age window: 1 year'}

    for target_age in target_age_list:
        for mac in mac_list:
            catchup_age = (target_age[1], target_age[1]+mac)
            age_label = f'Girls {target_age[0]}-{target_age[1]} + MAC {catchup_age[0]}-{catchup_age[1]}'

            for cn, cov_val in enumerate(coverage_array):

                scen_label = age_label + f': {np.round(cov_val, decimals=1)} coverage'

                for pn, metric in enumerate(metrics):
                    base_vals = msim_dict[base_label][metric].values[si:]
                    scen_vals = msim_dict[scen_label][metric].values[si:]
                    n_averted = sum(base_vals - scen_vals)
                    records += {'age': f'{target_age[0]}', 'coverage': int(round(cov_val, 1)*100), 'metric': f'{metric.replace("_"," ").capitalize()}', 'val': n_averted, 'mac': 'Age window: '+str(mac) +' years'}

    df = pd.DataFrame.from_dict(records)

    return df


def plot_age_impact(df):

    sns.set_style("whitegrid")
    ut.set_font(18)
    g = sns.catplot(
        data=df,
        kind="bar",
        x="coverage",
        y="val",
        row="metric",
        col="mac",
        hue="age",
        sharey=True,
        # height=12, aspect=1,
        )
    g.set_axis_labels("Coverage (%)", "")
    g.set_titles("{row_name} averted 2025-2100\nrelative to no vaccination\n{col_name}")

    for ax in g.axes.flat:
        sc.SIticks(ax)
    g.legend.set_title("Lower age")

    # fig.tight_layout()
    fig_name = 'figures/fig1_vx_impact.png'
    sc.savefig(fig_name, dpi=100)
    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    msim_dict = sc.loadobj('results/vx_scens.obj')
    df = preprocess_data(msim_dict)

    plot_age_impact(df)
