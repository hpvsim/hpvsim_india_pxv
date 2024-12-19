"""
Run HPVsim scenarios varying the year of introduction of prophylactic vaccination
"""


# %% General settings

import os

os.environ.update(
    OMP_NUM_THREADS='1',
    OPENBLAS_NUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1',
    MKL_NUM_THREADS='1',
)

# Standard imports
import numpy as np
import sciris as sc
import hpvsim as hpv

# Imports from this repository
import run_sim as rs

# Settings - used here and imported elsewhere
debug = 0
n_seeds = [20, 1][debug]  # How many seeds to run per cluster
start_year_list = np.arange(2025, 2036)
n_years_to_scaleup_list = [0, 5, 10]
product_list = ['bivalent', 'quadrivalent', 'nonavalent']


# %% Create interventions
def make_vx_scenarios(start_year_list=None, n_years_to_scaleup_list=None, product_list=None, final_coverage=0.9, efficacy=0.98, end=2100):

    vx_scenarios = dict()

    # Baseline
    vx_scenarios['Baseline'] = []

    # Vaccine eligibility
    routine_age = [13, 14]
    catchup_age = [14, 18]
    eligibility = lambda sim: (sim.people.doses == 0)

    for start_year in start_year_list:
        for n_years_to_scaleup in n_years_to_scaleup_list:

            scaleup_years = np.arange(start_year-1, start_year+n_years_to_scaleup+1)
            constant_years = np.arange(start_year+n_years_to_scaleup+1, end+1)
            scaleup_coverage = np.linspace(0, final_coverage, n_years_to_scaleup+2)
            constant_coverage = np.ones(len(constant_years)) * final_coverage
            years = np.concatenate([scaleup_years, constant_years])
            coverage = np.concatenate([scaleup_coverage, constant_coverage])

            for product in product_list:
                prod = hpv.default_vx(prod_name=product)
                prod.imm_init = dict(dist='uniform', par1=efficacy-0.01, par2=efficacy+0.01)
                label = f'Start year: {start_year}, scaleup: {n_years_to_scaleup} years, {product}'

                routine_vx = hpv.routine_vx(
                    prob=coverage,
                    years=years,
                    product=prod,
                    age_range=routine_age,
                    eligibility=eligibility,
                    label='Routine vx'
                )
                catchup_vx = hpv.campaign_vx(
                                prob=coverage,
                                years=years,
                                product=prod,
                                age_range=catchup_age,
                                eligibility=eligibility,
                                label='Catchup vx'
                )

                vx_scenarios[label] = [routine_vx]  #, catchup_vx]

    return vx_scenarios


def make_sims(calib_pars=None, vx_scenarios=None):
    """ Set up scenarios """

    st_intv = []  # make_st()

    all_msims = sc.autolist()
    for name, vx_intv in vx_scenarios.items():
        sims = sc.autolist()
        for seed in range(n_seeds):
            interventions = vx_intv + st_intv
            sim = rs.make_sim(calib_pars=calib_pars, debug=debug, interventions=interventions, end=2100, seed=seed)
            sim.label = name
            sims += sim
        all_msims += hpv.MultiSim(sims)

    msim = hpv.MultiSim.merge(all_msims, base=False)

    return msim


def run_sims(calib_pars=None, vx_scenarios=None, verbose=0.2):
    """ Run the simulations """
    msim = make_sims(calib_pars=calib_pars, vx_scenarios=vx_scenarios)
    msim.run(verbose=verbose)
    return msim


# %% Run as a script
if __name__ == '__main__':

    T = sc.timer()
    do_run = True
    do_save = False
    do_process = True

    # Run scenarios (usually on VMs, runs n_seeds in parallel over M scenarios)
    if do_run:
        calib_pars = sc.loadobj('results/india_pars.obj')
        vx_scenarios = make_vx_scenarios(start_year_list=start_year_list, n_years_to_scaleup_list=n_years_to_scaleup_list, product_list=product_list, end=2100)
        msim = run_sims(vx_scenarios=vx_scenarios)
        if do_save: msim.save('results/vs.msim')

        if do_process:

            metrics = ['year', 'asr_cancer_incidence', 'n_precin_by_age', 'n_females_alive_by_age', 'cancers', 'cancer_deaths']

            # Process results
            vx_scenarios = make_vx_scenarios(start_year_list=start_year_list, n_years_to_scaleup_list=n_years_to_scaleup_list, product_list=product_list, end=2100)
            scen_labels = list(vx_scenarios.keys())
            mlist = msim.split(chunks=len(scen_labels))

            msim_dict = sc.objdict()
            for si, scen_label in enumerate(scen_labels):
                reduced_sim = mlist[si].reduce(output=True)
                mres = sc.objdict({metric: reduced_sim.results[metric] for metric in metrics})
                msim_dict[scen_label] = mres

            sc.saveobj(f'results/vx_year_scens.obj', msim_dict)

    print('Done.')
