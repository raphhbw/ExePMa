# Proper Motion Anomaly (PMa)

Code to plot the Proper Motion anomaly (PMa) of any star in the [Kervella+2022 Vizier catalog](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/A%2bA/657/A7/tablea1). The code follows the method detailed in [Kervella+2019](https://ui.adsabs.harvard.edu/abs/2019A%26A...623A..72K/abstract).

An example of how to use the code is given in `example/PMa_example.py` with the potential kwargs that can be given to the different functions.

##### Important
- `self.data`: DataFrame with queried PMa information from Gaia eDR3 and Hipparcos + any geometry information given for the system  
- `plot_pma(ax, epoch=['eDR3' or 'Hipparcos'])`: Returns ax where the specifc PMa curve is plotted  
- `plot_disc_extent(ax, discdata)`: Returns ax where the disc extent is plotted as well as a Hill radius argument for the stability of the inner edge location.
