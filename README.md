# The Music of the Sphere

A speculative project to map the entire universe over time with low spatial resolution using CMB, 21 cm and local surveys. The simple time-evolution of the perturbations in the gravitational potential makes this a natural quantity to map; it (and its gradient) should be well constrained in a thin shell at z=1000 by the all-sky CMB data, and also by large scale structure surveys at low redshift. The longest wavelength fourier modes of the 3D potential will be best constrained; a question we are interested in is how high the resolution of a truly global evolving map of the universe can be made.

![](https://github.com/rogerblandford/Music/blob/master/doc/figures/hirezrealunivslice.png) ![](https://github.com/rogerblandford/Music/blob/master/Demos/Demo_figures/opac_phi3Ddomain.png)

Applications of such a global map could include: predicting very large scale structures present in the reionization era, visible with 21cm surveys; cleaning the CMB and low redshift survey data by insisting on physical consistency; exploring the universe pre-recombination; improved tests of non-Gaussianity; we might anticipate many more.

3D maps constrained by CMB data have been investigated by Wandelt et al, at the epoch of last scattering. We are interested in pushing this idea to its limits, to make a map that is either constrained by, or (equivalently) capable of predicting, all the cosmological data we have and to understand what are the ultimate limits with future data.

## People

* Roger Blandford (KIPAC)
* Phil Marshall (KIPAC)
* Laurence Perreault Levasseur (KIPAC)
* Hans Kristian Eriksen (Oslo)
* Ingunn Kathrine Wehus (Oslo)
* Ryan Keisler (KIPAC)
* Michael Schneider (LLNL)
* Roland de Putter (Caltech)

## Activities

* Project description, plus some theoretical explorations (Blandford) [here](https://github.com/rogerblandford/Music/blob/master/doc/music-allegro.pdf), specifically, connections to and consequences for inflation (Perreault Levasseur) [here](https://github.com/rogerblandford/Music/blob/master/doc/music-inflation.pdf).
* Simple universe-in-a-box investigation. IPython notebook demo [here](http://nbviewer.ipython.org/github/rogerblandford/Music/blob/master/Potential_Demo.ipynb).
* Application to Planck CMB temperature anisotropy maps. Work in progress here, [extracting the covariance matrix of the Gaussian approximation to the Planck posterior over maps at various $l_{max}$](https://github.com/rogerblandford/Music/blob/master/Demos/PlanckSamples_Demo.ipynb) and then [reconstructing the maximum posterior 3D map](https://github.com/rogerblandford/Music/blob/master/Demos/Reconstruction_Demo.ipynb).

Note: Prior to executing the IPython notebooks, please ensure to include the Music folder to your python path.

## Credits, License, etc

This is work in progress! If you use any of the ideas, code etc you find here, please cite us as (Blandford et al in preparation), and do please get in touch! We welcome feedback via the [issues](https://github.com/rogerblandford/Music/issues).

All content Copyright 2014, 2015 the authors. The code is distributed under the MIT license.
