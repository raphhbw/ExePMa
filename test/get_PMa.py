from PMa import circular_PMa
import numpy as np
import json

iparam = False
odir = '/Users/rbendahan/Astronomy/tools/PMa/test/'

if iparam:
    param_path = odir + 'params.json'
    f = open(param_path)
    discdata = json.load(f)
else:
    star = 'HD107146'
    discdata = {star: {"name": star, 'inc': None, 'inc_err':None,
              "PA": None, "PA_err": None}}
    # discdata = {star: {"name": star, 'inc': 66.8, 'inc_err':1.15,
    #           "PA": 100.0, "PA_err": 1.5}}

    # star = 'HD110058'
    # discdata = {star: {"name": star, 'inc': None, 'inc_err':None,
    #           "PA": None, "PA_err": None}}

gaia = 'eDR3'

# MC params (default in PMa.circular_PMa)
Na = 300 # number of semi-major axis to sample
Nrandom = 10000 # number of random PMa per radius
aps = np.logspace(-0.5, 2.5, Na) # semi major axis [au]

# target list
targets=[f for f in discdata.keys()]

for i, targeti in enumerate(targets):
    PMa = circular_PMa.PMa(discdata[targeti], gaia=gaia, iparam=iparam, odir=odir)

    ms_hip = PMa.mass_retrieval(aps, Nrandom, epoch='Hipparcos', subdir='log')
    ms_gaia = PMa.mass_retrieval(aps, Nrandom, epoch=gaia, subdir='log')

    PMa.plot(ms_hip, ms_gaia, figsize=(4,3), subdir='figures', only_save = True)
