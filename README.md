# dc2_satellite_census

Estimating the LSST survey selection function for resolved ultra-faint dwarf galaxies using DC2 and the Rubin Science Platform.

The simulation inputs (found in `data/`) were generated on the Rubin Science Platform (RSP).

Mock satellites are simulated using [ugali](https://github.com/DarkEnergySurvey/ugali) and injected into DC2 at the catalog level.

Satellite search is run using [simple](https://github.com/sidneymau/simple_adl/tree/kb)

Predictions for the number of observable galaxies are performed using [subhalo_satellite_connection](https://github.com/eonadler/subhalo_satellite_connection)
