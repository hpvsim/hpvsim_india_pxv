"""
Plot impact of shifting year and rate of introduction of vaccination
"""

import pandas as pd
import sciris as sc
from run_year_scenarios import start_year_list, n_years_to_scaleup_list
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

    for start_year in start_year_list:
        for n_years_to_scaleup in n_years_to_scaleup_list:
            scen_label = f'Start year: {start_year}, scaleup: {n_years_to_scaleup} years'

            for pn, metric in enumerate(metrics):
                base_vals = msim_dict[base_label][metric].values[si:]
                scen_vals = msim_dict[scen_label][metric].values[si:]
                n_averted = sum(base_vals - scen_vals)
                records += {'Introduction year': start_year, 'Years to scale up': n_years_to_scaleup, 'metric': f'{metric.replace("_"," ").capitalize()}', 'val': n_averted}

    df = pd.DataFrame.from_dict(records)

    return df


def plot_year_impact(df):

    # sns.set_style("whitegrid")
    ut.set_font(18)
    g = sns.catplot(
        data=df,
        kind="bar",
        x="Introduction year",
        y="val",
        row="metric",
        hue="Years to scale up",
        palette="viridis",
        # col="Years to scale up",
        sharey=True,
        # height=12, aspect=1,
        )
    g.set_axis_labels("Introduction year", "")
    g.set_titles("{row_name} averted 2025-2100\nrelative to no vaccination")  #\n{col_name}")

    for ax in g.axes.flat:
        sc.SIticks(ax)
    g.legend.set_title("Years to scale up")

    # fig.tight_layout()
    fig_name = 'figures/fig2_vx_year_impact.png'
    sc.savefig(fig_name, dpi=100)
    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    msim_dict = sc.loadobj('results/vx_year_scens.obj')
    df = preprocess_data(msim_dict)

    plot_year_impact(df)
