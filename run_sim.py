"""
Define the HPVsim simulation
"""

# Additions to handle numpy multithreading
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

# %% Settings and filepaths
# Debug switch
debug = 0  # Run with smaller population sizes and in serial
do_shrink = True  # Do not keep people when running sims (saves memory)

# Save settings
do_save = True
save_plots = True


# %% Simulation creation functions
def make_sim(calib_pars=None, debug=0, interventions=None, datafile=None, seed=1, end=2020):
    """
    Define parameters, analyzers, and interventions for the simulation
    """

    # Basic parameters
    pars = sc.objdict(
        n_agents=[20e3, 1e3][debug],
        dt=[0.25, 1.0][debug],
        start=[1960, 1980][debug],
        end=end,
        genotypes=[16, 18, 'hi5', 'ohr'],
        location='india',
        init_hpv_dist=dict(hpv16=0.4, hpv18=0.25, hi5=0.25, ohr=.1),
        init_hpv_prev={
            'age_brackets': np.array([12, 17, 24, 34, 44, 64, 80, 150]),
            'm': np.array([0.0, 0.25, 0.6, 0.25, 0.05, 0.01, 0.0005, 0]),
            'f': np.array([0.0, 0.35, 0.7, 0.25, 0.05, 0.01, 0.0005, 0]),
        },
        ms_agent_ratio=100,
        verbose=0.0,
        rand_seed=seed,
    )

    # Sexual behavior parameters
    # Debut: derived by fitting to 2019-21 DHS
    # Women:
    #           Age:   15,   18,   20,   22,   25
    #   Prop_active: 10.3, 38.5, 60.0, 73.9, 85.7
    # Men:
    #           Age:  15,   18,   20,   22,   25
    #   Prop_active: 0.8,  6.2, 16.2, 30.2, 51.5
    # For fitting, see https://doi.org/10.1371/journal.pcbi.1012181
    pars.debut = dict(
        f=dict(dist='lognormal', par1=15, par2=2),
        m=dict(dist='lognormal', par1=20, par2=2),
    )

    # Participation in marital and casual relationships
    # Derived to fit 2019-21 DHS data
    # For fitting, see https://doi.org/10.1371/journal.pcbi.1012181
    pars.layer_probs = dict(
        # Share of people of each age who are married
        m=np.array([
            [0, 5,   10,   15,   20,   25,   30,   35,   40,   45,   50,   55,   60,   65,   70,   75],
            [0, 0, 0.05, 0.25, 0.60, 0.80, 0.95, 0.80, 0.80, 0.65, 0.55, 0.40, 0.40, 0.40, 0.40, 0.40],  # Share f
            [0, 0, 0.01, 0.05, 0.10, 0.70, 0.90, 0.90, 0.90, 0.90, 0.80, 0.60, 0.60, 0.60, 0.60, 0.60]]  # Share m
        ),
        # Share of people of each age in casual partnerships
        c=np.array([
            [0, 5,   10,   15,   20,   25,   30,   35,   40,   45,   50,   55,   60,   65,   70,   75],
            [0, 0, 0.10, 0.50, 0.60, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.50, 0.01, 0.01],  # Share f
            [0, 0, 0.10, 0.20, 0.25, 0.35, 0.40, 0.70, 0.90, 0.90, 0.95, 0.95, 0.70, 0.30, 0.10, 0.10]],  # Share m
        ),
    )

    pars.m_partners = dict(
        m=dict(dist='poisson1', par1=0.01),
        c=dict(dist='poisson1', par1=0.1),
    )
    pars.f_partners = dict(
        m=dict(dist='poisson1', par1=0.01),
        c=dict(dist='neg_binomial', par1=2, par2=0.025),
    )

    if calib_pars is not None:
        pars = sc.mergedicts(pars, calib_pars)

    analyzers = [
        hpv.dalys(start=2000, life_expectancy=88.8),
        hpv.age_pyramid(
            timepoints=['2025', '2050', '2075', '2100'],
            datafile='india_age_pyramid.csv',
            edges=np.linspace(0, 100, 21),
        )
        ]

    # Interventions
    sim = hpv.Sim(pars=pars, interventions=interventions, analyzers=analyzers, datafile=datafile)

    return sim


# %% Simulation running functions
def run_sim(
        location=None, interventions=None, debug=0, seed=1, verbose=0.2,
        do_save=False, end=2020, calib_pars=None, meta=None):

    dflocation = location.replace(' ', '_')

    # Make sim
    sim = make_sim(
        debug=debug,
        interventions=interventions,
        calib_pars=calib_pars,
        end=end,
    )
    sim['rand_seed'] = seed
    sim.label = f'{location}--{seed}'

    # Store metadata
    sim.meta = sc.objdict()
    if meta is not None:
        sim.meta = meta  # Copy over meta info
    else:
        sim.meta = sc.objdict()
    sim.meta.location = location  # Store location in an easy-to-access place

    # Run
    sim['verbose'] = verbose
    sim.run()
    sim.shrink()

    if do_save:
        sim.save(f'results/india.sim')

    return sim


def run_parsets(debug=False, verbose=.1, analyzers=None, save_results=True, **kwargs):
    ''' Run multiple simulations in parallel '''

    parsets = sc.loadobj(f'results/india_pars_all.obj')
    kwargs = sc.mergedicts(dict(debug=debug, end=2040, verbose=verbose, analyzers=analyzers), kwargs)
    simlist = sc.parallelize(run_sim, iterkwargs=dict(calib_pars=parsets), kwargs=kwargs, serial=debug, die=True)
    msim = hpv.MultiSim(simlist)
    msim.reduce()
    if save_results:
        sc.saveobj(f'results/india_msim.obj', msim.results)

    return msim


# %% Run as a script
if __name__ == '__main__':

    T = sc.timer()

    location = 'india'
    sim = run_sim(location=location, end=2100, debug=debug)
    sim.plot()
    a = sim.get_analyzer()
    fig = a.plot(percentages=True)

    T.toc('Done')

