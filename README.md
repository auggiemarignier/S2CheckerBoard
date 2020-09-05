# S2CheckerBoard

A very simple script to build checkerboards on the sphere

## Usage
To install all dependencies using [poetry](https://python-poetry.org/)
```
poetry install
source .venv/bin/install
```
Otherwise the main dependency you need is [healpy](https://healpy.readthedocs.io/en/latest/index.html) which you can install using pip.

To make checkerboards, run
```
python main.py <options>
```

## Options
`--checkersize` - Size of the checkers in degrees e.g. 15 is 15X15 degrees.

`--nside` - Healpix resolution parameter, must be a power of 2.  Max/Default is 128.

`--L` - Bandlimit.  All spherical harmonic coeffs for l>=L are 0

`--save_map` - Saves map in a .fits file.

`--save_alm` - Saves harmonic coeffs in a .fits file.  Argument --L must be given.

`--save_img` - Saves a .png image of the map.

If any of the `save` options are specified, an `outputs` directory will be created in the directory where the script is run.