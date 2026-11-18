# ExePMa

Constrain the mass and separation of an unseen companion (planet or brown dwarf) from the **Hipparcos–Gaia proper motion anomaly (PMa)** of its host star.

ExePMa takes the PMa of any star in the [Kervella et al. (2022) catalogue](https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=J/A%2bA/657/A7/tablea1) and converts it into companion mass as a function of orbital separation, following the circular-orbit method of [Kervella et al. (2019)](https://scixplorer.org/abs/2019A%26A...623A..72K/abstract). It can also overlay complementary constraints on the same plot: debris-disc extent and Hill stability, inner-edge sculpting, diffusion timescale and the Gaia RUWE limit.

## Installation

```bash
git clone https://github.com/raphhbw/ExePMa.git
cd ExePMa
pip install -r requirements.txt
```

ExePMa is not yet a pip-installable package. Run your scripts from the repository root, or add the root to your `PYTHONPATH`. Querying a star needs an internet connection (SIMBAD and Vizier).

## Quickstart

```python
import matplotlib.pyplot as plt
from ExePMa import circular_PMa

star = circular_PMa.PMa(star='HD 107146')   # queries SIMBAD + Vizier, prints star.data

fig, ax = plt.subplots(figsize=(4, 3))
ax, aps = star.plot_pma(ax, epoch='eDR3', pma_params={'savelog': False})
ax, aps = star.plot_pma(ax, epoch='Hipparcos', pma_params={'savelog': False})
ax = star.add_axis_info(ax)
plt.show()
```

A detection curve (median mass with 1, 2 and 3σ bands) is only drawn when the PMa signal-to-noise is at least `snr` (default 3). Below that, set `upper_limit=True` in `plotting_params` to draw the 3σ upper limit instead.

A fuller example with geometry, disc constraints and the RUWE cut is in [`example/PMa_example.py`](example/PMa_example.py).

## Selecting the star

`PMa(...)` takes the following keyword arguments:

| Keyword | Default | Description |
|---|---|---|
| `star` | – | Star name, resolved to a HIP number via SIMBAD (e.g. `'HD 107146'`) |
| `hip` | – | Hipparcos identifier (e.g. `'HIP 60074'`); skips the SIMBAD query |
| `geometry` | `{}` | Orientation of the system, used to deproject the PMa (see below) |
| `mstar` | catalogue `M1` | Stellar mass [M<sub>☉</sub>] overriding the catalogue value |
| `gaia` | `'eDR3'` | Gaia release (only eDR3 is implemented) |
| `local_catalog` | `None` | Local copy of the Kervella+22 table (astropy `Table`), used instead of querying Vizier |

You must give either `star` or `hip`.

**Geometry.** If the orientation of the orbital plane is known (for example from a resolved debris disc), pass it to deproject the PMa:

```python
star = circular_PMa.PMa(star='HD 92945', geometry={'inc': 30., 'inc_err': 2., 'PA': 20., 'PA_err': 1.})
```

| Key | Default | Description |
|---|---|---|
| `inc`, `PA` | `None` (unknown) | Inclination and position angle [deg]. When unknown, a random orientation is assumed. |
| `inc_err`, `PA_err` | 1.0 | Uncertainties [deg] |

## The data table (`star.data`)

`star.data` is a one-row pandas DataFrame indexed by HIP identifier, with three column groups:

| Group | Columns | Units |
|---|---|---|
| `'Hipparcos'`, `'eDR3'` | `S_N` (PMa signal-to-noise), `PMa_ra`, `PMa_ra_err`, `PMa_dec`, `PMa_dec_err` | mas/yr |
| `'params'` | `parallax` | mas |
| | `dpc` | pc |
| | `ruwe` | – |
| | `mstar` | M<sub>☉</sub> |
| | `PA`, `PA_err`, `inc`, `inc_err` | deg |

Access single values with `.iloc[0]`, e.g. `star.data['params']['dpc'].iloc[0]`.

## Plotting functions

Every function draws on a matplotlib `Axes` (separation [au] on x, mass [M<sub>Jup</sub>] on y) and returns it. These plotting elements are used, and the constraints they show are described, in [Bendahan-West et al. (2026)](https://scixplorer.org/abs/2026MNRAS.546f2255B/abstract), Section 3.3.

| Function | What it plots |
|---|---|
| `plot_pma(ax, epoch, ...)` | Companion mass needed to explain the PMa at `epoch` (`'eDR3'` or `'Hipparcos'`). Also returns the separations `aps`. |
| `add_axis_info(ax, ...)` | Log axes, limits, labels, legend, and a top axis in arcsec |
| `plot_disc_extent(ax, aps, discdata, mstar)` | Disc edges, and the region where a planet's `NRhill` Hill radii would overlap the disc |
| `single_planet_inner_edge_sculptor(ax, aps, discdata, mstar, min_mp)` | Mass–separation relation of a single planet whose Hill radii reach the disc inner edge |
| `plot_diffusion(ax, aps, discdata, mstar, age)` | Minimum mass to sculpt the inner edge within the system age [Myr] ([Pearce et al. 2022](https://scixplorer.org/abs/2022A%26A...659A.135P/abstract), Eq. 7) |
| `ruwe_cutoff(ax, dpc, mstar, ruwe)` | Region excluded by a low Gaia RUWE (< 1.4) ([Limbach et al. 2024](https://scixplorer.org/abs/2024ApJ...973L..11L/abstract); [Kiefer et al. 2025](https://scixplorer.org/abs/2025A%26A...702A..76K/abstract)) |

### Main options

Options are passed as dictionaries. The most commonly used ones are:

| Dictionary | Key | Default | Description |
|---|---|---|---|
| `pma_params` | `aps` | `np.logspace(-0.5, 2.5, 300)` | Separations sampled [au] |
| | `Nrandom` | 10000 | Monte Carlo samples per separation |
| | `savelog` | `True` | Save the mass percentiles to `./log/` |
| `plotting_params` | `snr` | 3.0 | Minimum PMa S/N for a detection curve |
| | `upper_limit` | `False` | Draw the 3σ upper limit when S/N < `snr` |
| `discdata` | `r_in`, `r_out` | required | Disc inner and outer edges [au] |
| | `NRhill` | 3.0 | Number of Hill radii cleared by a planet |
| | `gaps` | 0 | Number of gaps in the disc |
| | `gaps_in`, `gaps_out` | `None` | Inner and outer edges of each gap [au], as lists of length `gaps` (required if `gaps > 0`) |

For example, a disc from 40 to 150 au with one gap between 70 and 90 au is:

```python
discdata = {'r_in': 40, 'r_out': 150, 'gaps': 1, 'gaps_in': [70], 'gaps_out': [90]}
```

Every function has more options (colours, labels, sigma bands, output paths, ...). These are listed in the function's docstring, e.g. `help(star.plot_pma)`.

## Requirements

Listed in [`requirements.txt`](requirements.txt): numpy, pandas, matplotlib, astropy, astroquery.

## References

If you use ExePMa, please take a look at [Bendahan-West et al. 2026, MNRAS, 546](https://scixplorer.org/abs/2026MNRAS.546f2255B/abstract), where the combined PMa and disc constraints are presented (Section 3.3), together with the papers its methods are based on:

- [Kervella et al. 2019, A&A, 623, A72](https://scixplorer.org/abs/2019A%26A...623A..72K/abstract): PMa method
- [Kervella et al. 2022, A&A, 657, A7](https://scixplorer.org/abs/2022A%26A...657A...7K/abstract): Hipparcos–Gaia eDR3 PMa catalogue
- [Pearce et al. 2022, A&A, 659, A135](https://scixplorer.org/abs/2022A%26A...659A.135P/abstract): diffusion and inner-edge sculpting
- [Limbach et al. 2024, ApJ, 973, L11](https://scixplorer.org/abs/2024ApJ...973L..11L/abstract): RUWE constraints
- [Kiefer et al. 2025, A&A, 702, A76](https://scixplorer.org/abs/2025A%26A...702A..76K/abstract): RUWE constraints

## License

MIT, see [LICENSE](LICENSE).
