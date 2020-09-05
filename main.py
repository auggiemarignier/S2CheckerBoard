from scipy.stats import multivariate_normal
import numpy as np
import healpy as hp
from matplotlib import cm
import matplotlib.pyplot as plt
import argparse
import os


def pixelise(signal, Nside, longs, lats):
    Npix = hp.nside2npix(Nside)
    pixnum = hp.ang2pix(Nside, longs, lats, lonlat=True)
    amap = np.zeros(Npix)
    count = np.zeros(Npix)
    nsample = len(signal)
    for i in range(nsample):
        pix = pixnum[i]
        amap[pix] += signal[i]
        count[pix] += 1.0
    for i in range(Npix):
        if count[i] > 0:
            amap[i] = amap[i] / count[i]
        else:
            amap[i] = hp.UNSEEN
    return amap


def healpy_lm(el, em, L):
    return int(em * (2 * L - 1 - em) / 2 + el)


def lm_hp2lm(flm_hp, L):
    f_lm = np.zeros([L * L], dtype=complex)
    for el in range(L):
        for em in range(el + 1):
            f_lm[el * el + el - em] = (
                pow(-1.0, -em) * (flm_hp[healpy_lm(el, em, L)]).conjugate()
            )
            f_lm[el * el + el + em] = flm_hp[healpy_lm(el, em, L)]
    return f_lm


parser = argparse.ArgumentParser()
parser.add_argument(
    "--checkersize",
    type=int,
    default=30,
    help="size of the checkers in degrees e.g. 15 is 15X15 degrees",
)
parser.add_argument(
    "--nside",
    type=int,
    default=128,
    help="Healpix resolution parameter, must be a power of 2.  Max/Default is 128",
)
parser.add_argument(
    "--L",
    type=int,
    default=32,
    help="Bandlimit.  All spherical harmonic coeffs for l>=L are 0",
)
parser.add_argument(
    "--save_map", action="store_true", help="Saves map in a .fits file."
)
parser.add_argument(
    "--save_alm",
    action="store_true",
    help="Saves harmonic coeffs in a .fits file.  Argument --L must be given",
)
parser.add_argument(
    "--save_img", action="store_true", help="Saves a .png image of the map."
)

args = parser.parse_args()

x, y = np.mgrid[-1:1:0.01, -1:1:0.01]
pos = np.empty(x.shape + (2,))
pos[:, :, 0] = x
pos[:, :, 1] = y
rv = multivariate_normal([0.0, 0.0], np.eye(2) * 0.15)
pdf = rv.pdf(pos)

n_in_row = 360 // (2 * args.checkersize)
n_in_col = 180 // (2 * args.checkersize)

row = np.hstack([pdf, -pdf] * n_in_row)
col = np.vstack([row, -row] * n_in_col)

longs = np.linspace(-180, 180, col.shape[1])
lats = np.linspace(-90, 90, col.shape[0])
longs, lats = np.meshgrid(longs, lats)
longs, lats = longs.flatten(), lats.flatten()
cboard = pixelise(col.flatten(), args.nside, longs, lats)

if args.save_map or args.save_alm or args.save_img:
    if not os.path.exists("outputs"):
        os.makedirs("outputs")

if args.save_map:
    hp.write_map(f"outputs/checkerboard{args.checkersize}.fits", cboard, overwrite=True)
if args.save_alm:
    try:
        alm_hp = hp.map2alm(cboard, args.L - 1)
        alm = lm_hp2lm(alm_hp, args.L)
        np.savetxt(f"outputs/checkerboardL{args.L}.txt", alm)
    except NameError:
        raise NameError("Please specify a bandlimit with the --L option.")
if args.save_img:
    cbar_end = max([abs(min(cboard)), max(cboard)])
    hp.mollview(
        cboard,
        cmap=cm.seismic_r,
        flip="geo",
        title=f"{args.checkersize}x{args.checkersize} Checkerboard",
        min=-cbar_end,
        max=cbar_end,
    )
    hp.graticule(args.checkersize)
    plt.savefig(f"outputs/checkerboard{args.checkersize}.png")
