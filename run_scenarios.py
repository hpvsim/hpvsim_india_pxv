"""
Run HPVsim scenarios varying the age of prophylactic vaccination
Note: requires an HPC to run with debug=False; with debug=True, should take 5-15 min
to run.
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
target_age_list = [[9, 10], [13, 14]]
mac_list = [4, 8]
coverage_array = np.linspace(0.1, 0.9, 9)


# %% Create interventions
def make_vx_scenarios(target_age_list=None, efficacy=0.98, mac_list=None, coverage_array=None, product='bivalent', start_year=2025):

    prod = hpv.default_vx(prod_name=product)
    prod.imm_init = dict(dist='uniform', par1=efficacy-0.01, par2=efficacy+0.01)
    eligibility = lambda sim: (sim.people.doses == 0)

    vx_scenarios = dict()

    # Baseline
    vx_scenarios['Baseline'] = []

    # Scenarios without catchup
    for age_range in target_age_list:
        routine_age = (age_range[0], age_range[1])

        for cov_val in coverage_array:
            label = f'Girls {age_range[0]}-{age_range[1]}: {np.round(cov_val, decimals=2)} coverage'
            routine_vx = hpv.routine_vx(
                prob=cov_val,
                start_year=start_year,
                product=prod,
                age_range=routine_age,
                eligibility=eligibility,
                label='Routine vx'
            )
            vx_scenarios[label] = [routine_vx]

    # Scenarios with catchup
    for age_range in target_age_list:
        routine_age = (age_range[0], age_range[1])
        for cov_val in coverage_array:
            routine_vx = hpv.routine_vx(
                        prob=cov_val,
                        start_year=start_year,
                        product=prod,
                        age_range=routine_age,
                        eligibility=eligibility,
                        label='Routine vx'
                    )
            for mac_age in mac_list:
                catchup_age = (age_range[1], age_range[1]+mac_age)
                label = f'Girls {age_range[0]}-{age_range[1]} + MAC {catchup_age[0]}-{catchup_age[1]}: {np.round(cov_val, decimals=2)} coverage'
                catchup_vx = hpv.campaign_vx(
                    prob=cov_val,
                    years=start_year,
                    product=prod,
                    age_range=catchup_age,
                    eligibility=eligibility,
                    label='Catchup vx'
                )

            vx_scenarios[label] = [routine_vx, catchup_vx]

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
        vx_scenarios = make_vx_scenarios(target_age_list=target_age_list, mac_list=mac_list, coverage_array=coverage_array)
        msim = run_sims(vx_scenarios=vx_scenarios)

        if do_save: msim.save('results/vs.msim')

        if do_process:

            metrics = ['year', 'asr_cancer_incidence', 'n_precin_by_age', 'n_females_alive_by_age', 'cancers', 'cancer_deaths']

            # Process results
            vx_scenarios = make_vx_scenarios(target_age_list=target_age_list, mac_list=mac_list, coverage_array=coverage_array)
            scen_labels = list(vx_scenarios.keys())
            mlist = msim.split(chunks=len(scen_labels))

            msim_dict = sc.objdict()
            for si, scen_label in enumerate(scen_labels):
                reduced_sim = mlist[si].reduce(output=True)
                mres = sc.objdict({metric: reduced_sim.results[metric] for metric in metrics})
                msim_dict[scen_label] = mres

            sc.saveobj(f'results/vx_scens.obj', msim_dict)

    print('Done.')
