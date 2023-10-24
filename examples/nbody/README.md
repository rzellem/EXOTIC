# N-body Retrievals 

Coming soon to EXOTIC...

## Prior estimate

The lomb-scargle periodogram from our ephemeris fitting code can be used to estimate a prior for an N-body retrieval. Transit timing variations (TTV) alone are not sufficient to uniquely constrain a single solution, instead multiple modes can exist which give the same amplitude/periodicity of perturbation. The distinguishing factor will rely on radial velocity measurements but none the less this approach can help limit the search space for potential planets. Below is a figure highlighting results from a grid of N-body simulations where we construct a mask based on the amplitude and periodicity found in the O-C diagram. The mask is used to constrain the N-body retrieval (i.e. search space) which significantly reduces the number of simulations needed to find a solution.

![](ttv_prior.png)

The example code has recently been refactored from [Nbody-AI](https://github.com/pearsonkyle/Nbody-ai) and is still in development. Some retrievals may not work entirely... please be paitent and report any issues to our slack.

If you make use of this code please cite:

https://ui.adsabs.harvard.edu/abs/2019AJ....158..243P/abstract
