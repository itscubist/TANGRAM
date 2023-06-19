# TANGRAM
An unbinned fitter (of 3 observables or less) I developed for my thesis: First measurement of atmospheric flux weighted electron neutrino - oxygen nuclei cross section below 125 MeV.

Since data is not public and study is not complete, I am not filling the input and ouput
directories. This is just the code and some example card files.

For more details of the project see the pdf or jpg file I presented in NEUTRINO2022 Conference:
BaranBodur_neutrino2022poster.pdf
BaranBodur_neutrino2022poster.jpg

The fits and results in the middle row of the poster is done with this fitter. However fitter
is improved a lot since then. Now it can deal with many nuissance parameters (systematics),
evalute their effects and importance to the result etc...


Requires CERN ROOT 6 (It has libraries for weighted KDEs which are essential)

This is unbinned likelihood where templates are obtained from MC events by kernel density
estimation. It can handle systematic effects changing templates by using event weights, where
each event will reweighted as the fitter progresses.

The fitter uses ntag,evis and gtag (gamma tagging) information for SK1-5. Fits everything
simultaneously.

The fit has 4 free variables for the relative fraction of 4 types of interactions, nueoxygen,
IBD and nuebar oxygen, atmpd backgrounds, and nueCC>125MeV

Systematics are restricted with pull terms based on their sigmas.

The fitter can be run in various modes, based on the card file and executables in the main dir:
It can do mock data studies to see fit variance.
It can make profile fits (though each step should typically be given separately for speed)
It can scan the effect of systematics on the MC templates at +-1 sigma
It can estimate the impact of systematics to nue-16O result, by generating and reweighting mock
data and observing how much % change occured on nue-16O result.
See some example cards and executables in the directory...

Requires root 6 (not root5), so use makeMainFit.sh and runMainFit.sh (or queue)

	**cards: Card files to determine and input every aspect of fitting 
	**inputs: Data and MC trees
	**src: Source files for the fitter
	**include: Include files for the fitter
	**outputs: Some test outputs
	**withQueue: Queue scripts and outputs (log, error, job and output files) of queued jobs







