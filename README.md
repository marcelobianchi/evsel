evsel
=====

Event Selection Tool using SHELL script.

The main tool is evsel.bshm. The package is a BASH set of function that can
fetch seismological earthquake catalog file from major sources making the
information available for future filtering and plotting.

For general usage of evsel.bshm try:

% source evsel.bshm
% evhelp

With the main package we disponibilize a tool for plotting the output using
the pgplot5 package (evgraph). The tool can make

	* Map plots (longitude, latitude)
	* Section plots (coordinate vs Depth)
	* General plots
	* Magnitude Histograms

The program also can plot log-log and semilog plots, and perform a line-fit
into the plotted data.  In the Map mode, user has the ability to overlay the
earthquakes points with the continents boundaries (GMT) or the plate
tectonics as mapped from Bird et al., 2012.

Bianchi,
24-05-2014