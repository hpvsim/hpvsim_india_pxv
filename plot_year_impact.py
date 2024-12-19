"""
Plot impact of shifting year and rate of introduction of vaccination
"""

import pandas as pd
import sciris as sc
from run_year_scenarios import start_year_list, n_years_to_scaleup_list, product_list
import utils as ut
import seaborn as sns
import numpy as np
import pylab as pl


def preprocess_data(msim_dict):

    # What to store
    start_year = 2025
    metrics = ['cancers', 'cancer_deaths']
    records = sc.autolist()
    base_label = 'Baseline'
    si = sc.findinds(msim_dict[base_label].year, start_year)[0]

    for start_year in start_year_list:
        for n_years_to_scaleup in n_years_to_scaleup_list:
            for product in product_list:
                scen_label = f'Start year: {start_year}, scaleup: {n_years_to_scaleup} years, {product}'

                for pn, metric in enumerate(metrics):
                    base_vals = msim_dict[base_label][metric].values[si:]
                    scen_vals = msim_dict[scen_label][metric].values[si:]
                    n_averted = sum(base_vals - scen_vals)
                    records += {'Introduction year': start_year, 'Years to scale up': n_years_to_scaleup, 'Product': product, 'metric': f'{metric.replace("_"," ").capitalize()}', 'val': n_averted}

    df = pd.DataFrame.from_dict(records)

    return df


def plot_single(ax, mres, to_plot, si, ei, color, label=None, smooth=True, al=None):
    years = mres.year[si:ei]
    ts = 0.5  # detection rate / sensitivity

    if to_plot == 'precin_incidence':
        best = mres.n_precin_by_age[3:11, si:ei].sum(axis=0) / mres.n_females_alive_by_age[3:11, si:ei].sum(axis=0) * ts
        low = mres.n_precin_by_age.low[3:11, si:ei].sum(axis=0) / mres.n_females_alive_by_age.low[3:11, si:ei].sum(axis=0) * ts
        high = mres.n_precin_by_age.high[3:11, si:ei].sum(axis=0) / mres.n_females_alive_by_age.high[3:11, si:ei].sum(axis=0) * ts
    else:
        best = mres[to_plot][si:ei]
        low = mres[to_plot].low[si:ei]
        high = mres[to_plot].high[si:ei]

    if smooth:
        best = np.convolve(list(best), np.ones(5), "valid")/5
        low = np.convolve(list(low), np.ones(5), "valid")/5
        high = np.convolve(list(high), np.ones(5), "valid")/5
        years = years[4:]

    ax.plot(years, best, color=color, label=label)
    # ax.fill_between(years, low, high, alpha=0.1, color=color)
    return ax


def plot_year_impact(df):

    sns.set_style("whitegrid")
    ut.set_font(100)
    g = sns.catplot(
        data=df,
        kind="bar",
        x="Introduction year",
        y="val",
        row="metric",
        col="Product",
        hue="Years to scale up",
        palette="viridis",
        sharey=True,
        height=20, aspect=1,
        )
    g.set_axis_labels("Introduction year", "# averted 2025-2100")
    g.set_titles("Cumulative {row_name}\n{col_name} vx")

    for ax in g.axes.flat:
        sc.SIticks(ax)
        ax.tick_params(axis='x', which='both', labelsize=100, rotation=90)
    g.legend.set_title("Years to scale up")

    # fig.tight_layout()
    fig_name = 'figures/fig2_vx_year_impact.png'
    sc.savefig(fig_name, dpi=100)
    return


def plot_year_ts(msim_dict):

    ut.set_font(20)
    plot_start_year_list = np.arange(2025, 2036, 5)
    colors = sc.vectocolor(len(plot_start_year_list), reverse=True)
    # covcolors = sc.vectocolor(len(plot_coverage_arr), reverse=True)
    plot_dict = sc.objdict(
        precin_incidence='Detectable HPV prevalence (F15+)',
        asr_cancer_incidence='ASR cancer incidence'
    )

    fig, axes = pl.subplots(len(plot_dict), len(product_list), figsize=(17, 10))
    # axes = axes.ravel()

    # What to plot
    start_year = 2016
    end_year = 2100
    mbase = msim_dict['Baseline']
    si = sc.findinds(mbase.year, start_year)[0]
    ei = sc.findinds(mbase.year, end_year)[0]
    cn = 0

    for rn, to_plot, plot_label in plot_dict.enumitems():

        for cn, product in enumerate(product_list):

            ax = axes[rn, cn]

            # Plot baseline
            baseline_label = 'Baseline'
            mres = msim_dict[baseline_label]
            ax = plot_single(ax, mres, to_plot, si, ei, 'k', label='No vaccination')

            # Plot with vaccination
            for syn, start_year in enumerate(plot_start_year_list):
                scen_label = f'Start year: {start_year}, scaleup: 0 years, {product}'
                mres = msim_dict[scen_label]
                ax = plot_single(ax, mres, to_plot, si, ei, colors[syn], label=f'90% PxV coverage from {start_year}', al=scen_label)

            ax.set_ylim(bottom=0)
            # ax.set_ylabel(plot_label)
            ax.set_title(plot_label+f'\n{product} vx')
            if (cn == 0) & (rn == 0): ax.legend(frameon=False)
            if to_plot == 'asr_cancer_incidence': ax.axhline(y=4, color='k', ls='--')

    fig.tight_layout()
    fig_name = 'figures/fig3_vx_year_ts.png'
    sc.savefig(fig_name, dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    msim_dict = sc.loadobj('results/vx_year_scens.obj')

    # Cumulative impact
    df = preprocess_data(msim_dict)
    plot_year_impact(df)

    # Time series
    plot_year_ts(msim_dict)

