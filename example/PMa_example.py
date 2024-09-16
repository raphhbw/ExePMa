from PMa import Database
from PMa import circular_PMa
import matplotlib.pyplot as plt
# Personal matplotlib style lib
plt.style.use('images_inv')

## System information
star = 'HD92945'

# Geometry
# {"inc": value [default None],
#  "inc_err": value [default None],
#  "PA": value [default None],
#  "PA_err": value [default None]}
geometry = {"inc": 30.,
            "PA": 20.}
# Disc information
# {"r_in": value, Requires value
#  "r_out": value, Requires value
#  "NRhill": value [default 3.0]}
disc_data = {"r_in":45,
             "r_out":120}

## PMa parameters
# {"Na": value [default 300],
#  "aps": array [default np.logspace(-0.5, 2.5, Na)],
#  "Nrandom": value [default 10000],
#  "savelog": Boolean [default True],
#  "subdir": subdir name [default 'log'],
#  "odir": path to repo [default './']}
pma_params = {"savelog": False}

## Plotting parameters
# {"snr": value [default 3.0],
#  "color": array colour [default ['C1', 'C0']],
#  "alpha": value [default 0.1],
#  "xlim": array xmin, xmax [default [0.3,300]],
#  "ylim": array ymin, ymax [default [3.0e-1, 1.0e2]],
#  "disc_color": colour [default 'C1'],
#  "disc_alpha": value [default 0.4]}
plotting_params = {"snr": 3.0}

### Initialise main PMa function
GetPMa = circular_PMa.PMa(star=star,
                          geometry=geometry
                          )
# Retrive database with PMa info + stellar info
# print(GetPMa.data)

### Test mass_retrival function
# gaia_ms = GetPMa.mass_retrieval(epoch='eDR3', pma_params=pma_params)

### Plotting test with PMa curves
fig, ax = plt.subplots(1, figsize=(8,6))
# Generate PMa curve for one epoch [either eDR3 or Hipparcos]
ax, aps = GetPMa.plot_pma(ax, epoch='eDR3', pma_params=pma_params, plotting_params=plotting_params)
ax, aps = GetPMa.plot_pma(ax, epoch='Hipparcos', pma_params=pma_params, plotting_params=plotting_params)
# Add axis information
ax = GetPMa.add_axis_info(ax)
# Add disc inner and outer edges information + planet truncation argument
ax = GetPMa.plot_disc_extent(ax=ax, aps=aps, discdata=disc_data)
plt.show()