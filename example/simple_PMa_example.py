from ExePMa import Database
from ExePMa import circular_PMa
import matplotlib.pyplot as plt

## System information
star = 'HIP544'

## PMa parameters
pma_params = {"savelog": False}

## Plotting parameters
plotting_params = {"snr": 3.0} # significance in sigma

### Initialise main PMa function
GetPMa = circular_PMa.PMa(star=star)

### Plotting test with PMa curves
fig, ax = plt.subplots(1, figsize=(4,3))
# Generate PMa curve for one epoch [either eDR3 or Hipparcos]
ax, aps = GetPMa.plot_pma(ax, epoch='eDR3', pma_params=pma_params, plotting_params=plotting_params)
ax, aps = GetPMa.plot_pma(ax, epoch='Hipparcos', pma_params=pma_params, plotting_params=plotting_params)
# Add axis information
ax = GetPMa.add_axis_info(ax)
ax.text(30, 0.6, star, color='k', fontsize=10)
plt.tight_layout()
plt.show()