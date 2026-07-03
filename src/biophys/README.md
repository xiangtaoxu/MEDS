# biophys (placeholder)

**Fast, sub-daily, mostly-stateless physical flux calculators**: canopy radiative transfer
(reads canopy structure from `src/state/`, writes a diagnostic absorbed-PAR field consumed by
`plant/leaf`), leaf/canopy energy balance, and soil hydrology. Device-eligible and netCDF-free
(forcing enters via a passed-in boundary struct, never a direct `use netcdf`). Empty until
implemented.
