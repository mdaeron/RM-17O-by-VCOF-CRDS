#! /usr/bin/env python3

import numpy as np
from scipy.optimize import fsolve


def ratios_to_deltas(
	R45,
	R46,
	R13_VPDB = 0.01118,  # (Chang & Li, 1990)
	R18_VSMOW = 0.0020052,  # (Baertschi, 1976)
	R17_VSMOW = 0.00038475,  # (Assonov & Brenninkmeijer, 2003, rescaled to R13_VPDB)
	LAMBDA_17=0.528,  # (Barkan & Luz, 2005)
	D17O = 0,  # in permil
):
	R45 = np.asarray(R45)
	R46 = np.asarray(R46)
	if R45.shape != R46.shape:
		raise ValueError(
			'R45 and R46 must both be floats or both be arrays of the same shape.'
		)

	def f(R18):
		K = np.exp(D17O / 1e3) * R17_VSMOW / R18_VSMOW**LAMBDA_17
		return (
			-3 * K**2 * R18 ** (2 * LAMBDA_17)
			+ 2 * K * R45 * R18**LAMBDA_17
			+ 2 * R18
			- R46
		)

	R18 = fsolve(f, R46 / R18_VSMOW / 2, xtol = 1e-16)
	R17 = np.exp(D17O / 1e3) * R17_VSMOW * (R18 / R18_VSMOW) ** LAMBDA_17
	R13 = R45 - 2 * R17

	d13C_VPDB = (R13 / R13_VPDB - 1) * 1e3
	d18O_VSMOW = (R18 / R18_VSMOW - 1) * 1e3

	if d13C_VPDB.size == 1:
		return (float(d13C_VPDB), float(d18O_VSMOW))

	return (d13C_VPDB, d18O_VSMOW)


def deltas_to_ratios(
	d13C_VPDB,
	d18O_VSMOW,
	R13_VPDB = 0.01118,  # (Chang & Li, 1990)
	R18_VSMOW = 0.0020052,  # (Baertschi, 1976)
	R17_VSMOW = 0.00038475,  # (Assonov & Brenninkmeijer, 2003, rescaled to R13_VPDB)
	LAMBDA_17=0.528,  # (Barkan & Luz, 2005)
	D17O = 0,  # in permil
):
	d13C_VPDB = np.asarray(d13C_VPDB)
	d18O_VSMOW = np.asarray(d18O_VSMOW)
	if d13C_VPDB.shape != d18O_VSMOW.shape:
		raise ValueError(
			'd13C_VPDB and d18O_VSMOW must both be floats or both be arrays of the same shape.'
		)

	R13 = R13_VPDB * (1 + d13C_VPDB / 1e3)
	R18 = R18_VSMOW * (1 + d18O_VSMOW / 1e3)
	R17 = np.exp(D17O / 1e3) * R17_VSMOW * (1 + d18O_VSMOW / 1e3) ** LAMBDA_17

	R45 = 2 * R17 + R13
	R46 = 2 * R18 + 2 * R17 * R13 + R17**2

	return (R45, R46)


if __name__ == '__main__':
	x1, y1 = deltas_to_ratios(-10, -20)
	x2, y2 = deltas_to_ratios(10, 20)
	print(ratios_to_deltas(np.linspace(x1, x2, 4), np.linspace(y1, y2, 4)))
