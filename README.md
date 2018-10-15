# HighResMIP-extreme_indices
Use ICCLIM python package to calculate various climate extreme indices
ICCLIM documentation at https://icclim.readthedocs.io/en/latest/index.html

Standard indices are simple to calculate.
User defined indices are also straightforward.

The script simply loops over models, years and indices required, and uses icclim.
The metadata from the source file is added back into the final indices files.
The files are also converted to netcdf4 for extra compression options.
