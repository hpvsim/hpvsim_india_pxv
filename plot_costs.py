"""
Plot costs for infant vaccination scenarios
"""

import sciris as sc
import utils as ut
import seaborn as sns
import matplotlib.pyplot as plt
from plot_fig2 import preprocess_data


def plot_costs(df):
    sns.set_style("whitegrid")
    ut.set_font(24)
    df2 = df.loc[df['metric'].isin(['DALYs', 'cost'])]
    df2 = df2.groupby(['coverage', 'efficacy', 'metric']).val.first().unstack().reset_index()
    df2['Cost per DALY averted'] = df2['cost'] / df2['DALYs']
    df2['DALYs'] = df2['DALYs']/1e6
    df2['Adolescent coverage'] = df2['coverage']
    df2['Infant efficacy'] = df2['efficacy']
    dfplot = df2.loc[df2['Cost per DALY averted'] >= 0]
    dfplot = dfplot.loc[dfplot['Adolescent coverage'].isin([20,40,60,80])]

    fig, ax = plt.subplots(1, 1, figsize=(11, 10))
    sns.scatterplot(
        ax=ax,
        data=dfplot,
        x="DALYs",
        y="Cost per DALY averted",
        hue="Infant efficacy",
        palette='viridis_r',
        size="Adolescent coverage",
        sizes=(20, 500),
        legend="full"
    )
    ax.legend(frameon=False)
    ax.set_xlabel('DALYs averted (M), 2025-2100')
    ax.set_ylabel('Incremental cost / DALY averted (USD), 2025-2100')
    fig.tight_layout()
    fig_name = 'figures/vx_econ_impact.png'
    sc.savefig(fig_name, dpi=100)

    return df2


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    msim_dict = sc.loadobj('results/vx_scens.obj')
    cost_dict = sc.objdict({
        'Routine vx': 9,
        'Catchup vx': 9,
        'Infant vx': 5,
        'excision': 41.76,
        'ablation': 11.76,
        'radiation': 450
    })
    df = preprocess_data(msim_dict, cost_dict)
    df2 = plot_costs(df)
