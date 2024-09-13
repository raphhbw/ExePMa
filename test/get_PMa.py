from PMa import Database
# import numpy as np

geometry = {"inc": 30.,
            "PA": 20.}

GetData = Database.Database(star='HD92945', 
                            geometry=geometry
                            )
print(GetData.data)

# # MC params (default in PMa.circular_PMa)
# Na = 300 # number of semi-major axis to sample
# Nrandom = 10000 # number of random PMa per radius
# aps = np.logspace(-0.5, 2.5, Na) # semi major axis [au]

# # target list
# targets=[f for f in discdata.keys()]

# for i, targeti in enumerate(targets):
#     PMa = circular_PMa.PMa(discdata[targeti], gaia=gaia, iparam=iparam, odir=odir)

#     ms_hip = PMa.mass_retrieval(aps, Nrandom, epoch='Hipparcos', subdir='log')
#     ms_gaia = PMa.mass_retrieval(aps, Nrandom, epoch=gaia, subdir='log')

#     PMa.plot(ms_hip, ms_gaia, figsize=(4,3), subdir='figures', only_save = True)
