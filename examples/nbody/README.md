# N-body Retrievals 

EXOTIC provides an [N-body API](../../exotic/api/nbody.py) for simulation and retrieval development.

## Prior estimate

The Lomb-Scargle periodogram from our ephemeris fitting code can help estimate a prior for an N-body retrieval. Transit timing variations (TTV) alone can admit multiple solutions with the same amplitude and periodicity. Radial velocity measurements can help distinguish those solutions.

The API's `estimate_prior()` function reads a precomputed simulation grid from a FITS file. Its default path is `exotic/api/ttv_grid.fits`, but that grid is not included in this repository. Supply an existing grid with `estimate_prior(..., file_name="/path/to/ttv_grid.fits")`; the function does not generate it.

The [historical demo script](../../older_for_review/examples/nbody/example.py) and its [grid illustration](../../older_for_review/examples/nbody/grid_500_150.png) are preserved in `older_for_review/`. The demo depends on the missing grid, so its complete workflow cannot run from the bundled files alone. With the grid available, it saves a prior-estimate figure as `ttv_prior.png`.

The demo was adapted from [Nbody-AI](https://github.com/pearsonkyle/Nbody-ai) and remains development material.

If you make use of this code please cite:

https://ui.adsabs.harvard.edu/abs/2019AJ....158..243P/abstract
