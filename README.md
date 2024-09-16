# ExePMa

Code to plot the Proper Motion anomaly (PMa) of any star in the <a href="https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/A%2bA/657/A7/tablea1" target="_blank">Kervella+2022 Vizier catalog</a>. The code follows the method detailed in <a href="https://ui.adsabs.harvard.edu/abs/2019A%26A...623A..72K/abstract" target="_blank"> Kervella+2019</a>.

An example of how to use the code is given in `example/PMa_example.py` with the potential kwargs that can be given to the different functions.

##### Important
- [_self.data_](./ExePMa/circular_PMa.py): DataFrame with queried PMa information from Gaia eDR3 and Hipparcos + any geometry information given for the system  
- [_plot_pma(ax, epoch=['eDR3' or 'Hipparcos'])_](./ExePMa/plotting.py): Returns ax where the specifc PMa curve is plotted  
- [_plot_disc_extent(ax, discdata)_](./ExePMa/plotting.py): Returns ax where the disc extent is plotted as well as a Hill radius argument for the stability of the inner edge location.

##### Requirements
Requirements are listed in [_requirements.txt_](requirements.txt).