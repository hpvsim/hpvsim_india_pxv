"""
Plot 1 for infant vaccination scenarios
"""


import pylab as pl
import sciris as sc
from run_scenarios import coverage_arr, efficacy_arr
import utils as ut
import numpy as np


def plot_single(ax, mres, to_plot, si, ei, color, label=None, smooth=True):
    years = mres.year[si:ei]
    best = mres[to_plot][si:ei]
    low = mres[to_plot].low[si:ei]
    high = mres[to_plot].high[si:ei]
    if smooth:
        best = np.convolve(list(best), np.ones(5), "valid")/5
        low = np.convolve(list(low), np.ones(5), "valid")/5
        high = np.convolve(list(high), np.ones(5), "valid")/5
        years = years[4:]

    ax.plot(years, best, color=color, label=label)
    # ax.fill_between(years, low, high, alpha=0.5, color=color)
    return ax


def plot_fig3_old(msim_dict):

    ut.set_font(16)
    colors = sc.vectocolor(len(efficacy_arr), reverse=True)
    plot_coverage_arr = coverage_arr[1::2]  # which ones to plot

    n_plots = len(plot_coverage_arr)
    n_rows, n_cols = sc.get_rows_cols(n_plots)
    fig, axes = pl.subplots(n_rows, n_cols, figsize=(11, 10))
    if n_plots > 1: axes = axes.flatten()

    # What to plot
    start_year = 2016
    end_year = 2060
    to_plot = 'asr_cancer_incidence'
    mbase = msim_dict['Baseline']
    si = sc.findinds(mbase.year, start_year)[0]
    ei = sc.findinds(mbase.year, end_year)[0]

    for pn, cov_val in enumerate(plot_coverage_arr):

        ax = axes[pn] if n_plots > 1 else axes

        # Plot baseline
        baseline_label = 'Baseline'
        mres = msim_dict[baseline_label]
        ax = plot_single(ax, mres, to_plot, si, ei, 'k', label='Baseline')

        # Plot adolescents
        adolescent_label = f'Adolescent: {cov_val} coverage'
        mres = msim_dict[adolescent_label]
        ax = plot_single(ax, mres, to_plot, si, ei, 'k', label='Adolescents')

        for ie, eff_val in enumerate(efficacy_arr):
            infant_label = f'Adolescent: {cov_val} coverage, Infants: {eff_val} efficacy'
            mres = msim_dict[infant_label]
            ax = plot_single(ax, mres, to_plot, si, ei, colors[ie], label=f'Infants, {int(eff_val*100)}% efficacy')

        ax.set_ylim(bottom=0, top=23)
        ax.set_ylabel('ASR cancer incidence')
        ax.set_title(f'{int(cov_val*100)}% vaccine coverage')
        ax.axhline(y=4, color='k', ls='--')
        if pn == 0: ax.legend()

    fig.tight_layout()
    fig_name = 'figures/fig3_vx_scens.png'
    sc.savefig(fig_name, dpi=100)

    return


def plot_fig3(msim_dict):

    ut.set_font(16)
    colors = sc.vectocolor(len(efficacy_arr), reverse=True)
    plot_coverage_arr = coverage_arr[1::2]  # which ones to plot
    covcolors = sc.vectocolor(len(plot_coverage_arr), reverse=True)

    fig, axes = pl.subplots(2, 1, figsize=(11, 10))

    # What to plot
    start_year = 2016
    end_year = 2060
    to_plot = 'asr_cancer_incidence'
    mbase = msim_dict['Baseline']
    si = sc.findinds(mbase.year, start_year)[0]
    ei = sc.findinds(mbase.year, end_year)[0]

    # Plot adolescent scenarios
    pn = 0
    ax = axes[pn]

    # Plot baseline
    baseline_label = 'Baseline'
    mres = msim_dict[baseline_label]
    ax = plot_single(ax, mres, to_plot, si, ei, 'k', label='Baseline')

    # Plot adolescents
    for cn, cov_val in enumerate(plot_coverage_arr):
        adolescent_label = f'Adolescent: {cov_val} coverage'
        mres = msim_dict[adolescent_label]
        ax = plot_single(ax, mres, to_plot, si, ei, covcolors[cn], label=f'Adolescents, {int(np.ceil(cov_val*100))}% coverage')

    ax.set_ylim(bottom=0, top=23)
    ax.set_ylabel('ASR cancer incidence')
    ax.set_title('Adolescent vaccination scenarios')
    ax.axhline(y=4, color='k', ls='--')
    ax.legend()
    pn += 1

    # Plot infants
    ax = axes[pn]

    # Plot baseline
    baseline_label = 'Baseline'
    mres = msim_dict[baseline_label]
    ax = plot_single(ax, mres, to_plot, si, ei, 'k', label='Baseline')

    for ie, eff_val in enumerate(efficacy_arr):
        infant_label = f'Adolescents: {cov_val} coverage, Infants: {eff_val} efficacy'
        mres = msim_dict[infant_label]
        ax = plot_single(ax, mres, to_plot, si, ei, colors[ie], label=f'Infants, {int(np.ceil(eff_val*100))}% efficacy')

    ax.set_ylim(bottom=0, top=23)
    ax.set_ylabel('ASR cancer incidence')
    ax.set_title('Infant vaccination scenarios')
    ax.axhline(y=4, color='k', ls='--')
    ax.legend()

    fig.tight_layout()
    fig_name = 'figures/fig3a_vx_scens.png'
    sc.savefig(fig_name, dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    # Load scenarios and construct figure
    msim_dict = sc.loadobj('results/vx_scens.obj')
    plot_fig3(msim_dict)




