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
    ts = 0.5  # detection rate / sensitivity

    if to_plot == 'precin_incidence':
        best = mres.n_precin_by_age[3:10, si:ei].sum(axis=0) / mres.n_females_alive_by_age[3:10, si:ei].sum(axis=0) * ts
        low = mres.n_precin_by_age.low[3:10, si:ei].sum(axis=0) / mres.n_females_alive_by_age.low[3:10, si:ei].sum(axis=0) * ts
        high = mres.n_precin_by_age.high[3:10, si:ei].sum(axis=0) / mres.n_females_alive_by_age.high[3:10, si:ei].sum(axis=0) * ts
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


def plot_fig3(msim_dict):

    ut.set_font(16)
    colors = sc.vectocolor(len(efficacy_arr), reverse=True)
    plot_coverage_arr = coverage_arr  #[[0,4,8]]  # which ones to plot
    covcolors = sc.vectocolor(len(plot_coverage_arr), reverse=True)
    plot_dict = sc.objdict(
        precin_incidence='Detectable HPV prevalence, females 15-49',
        asr_cancer_incidence='ASR cancer incidence'
    )

    fig, axes = pl.subplots(len(plot_dict), 2, figsize=(11, 10))

    # What to plot
    start_year = 2016
    end_year = 2060
    mbase = msim_dict['Baseline']
    si = sc.findinds(mbase.year, start_year)[0]
    ei = sc.findinds(mbase.year, end_year)[0]

    for rn, to_plot, plot_label in plot_dict.enumitems():

        # Plot adolescent scenarios
        cn = 0
        ax = axes[rn, cn]

        # Plot baseline
        baseline_label = 'Baseline'
        mres = msim_dict[baseline_label]
        ax = plot_single(ax, mres, to_plot, si, ei, 'k', label='Baseline')

        # Plot adolescents
        for cvn, cov_val in enumerate(plot_coverage_arr):
            adolescent_label = f'Adolescent: {cov_val} coverage'
            mres = msim_dict[adolescent_label]
            ax = plot_single(ax, mres, to_plot, si, ei, covcolors[cvn], label=f'Adolescents, {int(np.ceil(cov_val*100))}% coverage')

        ax.set_ylim(bottom=0)  #, top=23)
        ax.set_ylabel(plot_label)
        ax.set_title('Adolescent vaccination scenarios')
        if to_plot == 'asr_cancer_incidence': ax.axhline(y=4, color='k', ls='--')
        ax.legend()
        cn += 1

        # Plot infants
        ax = axes[rn, cn]

        # Plot baseline
        baseline_label = 'Baseline'
        mres = msim_dict[baseline_label]
        ax = plot_single(ax, mres, to_plot, si, ei, 'k', label='Baseline')

        for ie, eff_val in enumerate(efficacy_arr):
            cov_val = eff_val*0.9/0.95
            infant_label = f'Adolescents: {cov_val} coverage, Infants: {eff_val} efficacy'
            mres = msim_dict[infant_label]
            ax = plot_single(ax, mres, to_plot, si, ei, colors[ie], label=f'Infants, {int(np.ceil(eff_val*100))}% efficacy')

        ax.set_ylim(bottom=0)  #, top=23)
        ax.set_ylabel(plot_label)
        ax.set_title('Infant vaccination scenarios')
        if to_plot == 'asr_cancer_incidence': ax.axhline(y=4, color='k', ls='--')
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




