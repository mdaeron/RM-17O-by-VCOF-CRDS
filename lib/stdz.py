#! /usr/bin/env python3

# from pylab import *

import scipy.linalg, lmfit
from numpy import log, array, exp, cov, ix_, sqrt
from scipy.stats import t as tstudent
from scipy.stats import f as ssf
from sklearn.covariance import MinCovDet


def estimate_covariance(X, robust = True):
	# X.shape must be (N,2)
	if robust:
		return MinCovDet().fit(X).covariance_
	else:
		return cov(X.T)


class IsotopicParameters:
	def __init__(self):
		self.R13_VPDB = 0.01118  # (Chang & Li, 1990)
		self.R18_VSMOW = 0.0020052  # (Baertschi, 1976)
		self.R17_VSMOW = (
			0.00038475  # (Assonov & Brenninkmeijer, 2003, rescaled to R13_VPDB)
		)
		self.LAMBDA_17 = 0.528  # (Barkan & Luz, 2005)
		self.d18OVSMOW_VPDB = 30.92  # (Coplen et al., 1983))

	R18_VPDB = property(
		fget = lambda self: self.R18_VSMOW * (1 + self.d18OVSMOW_VPDB / 1e3)
	)
	R17_VPDB = property(
		fget = lambda self: self.R17_VSMOW
		* (1 + self.d18OVSMOW_VPDB / 1e3) ** self.LAMBDA_17
	)


ISOPARAMETERS = IsotopicParameters()

RM_NBS_19 = {
	'name': 'NBS-19',
	'd13C_VPDB': {'value': +1.95, 'sigma': 0.0},
	'd18O_VPDB': {'value': -2.2, 'sigma': 0.0},
	'D17O_VSMOW': {
		'value': -0.102,
		'sigma': 0.010,
	},  # Wostbrock et al. (2020) 10.1016/j.chemgeo.2019.119432
}

RM_NBS_19_CO2_25C = {
	'name': 'NBS-19-CO2-25C',
	'd13C_VPDB': {'value': +1.95, 'sigma': 0.0},
	'd18O_VSMOW': {'value': (1000 - 2.2) * 1.03092 * 1.01025 - 1000, 'sigma': 0.0},
	'D17O_VSMOW': {
		'value': -0.182,
		'sigma': 0.010,
	},  # Barkan et al. (2019) 10.1002/rcm.8391
	'D17O_VSMOW': {
		'value': -0.155,
		'sigma': 0.010,
	},  # Wostbrock et al. (2020) 10.1016/j.chemgeo.2019.119432
}

RM_NBS_19_CO2_90C = {
	'name': 'NBS-19-CO2-90C',
	'd13C_VPDB': {'value': +1.95, 'sigma': 0.0},
	'd18O_VSMOW': {'value': (1000 - 2.2) * 1.03092 * 1.008146 - 1000, 'sigma': 0.0},
	'D17O_VSMOW': {
		'value': -0.135,
		'sigma': 0.010,
	},  # Passey et al. (2014) 10.1016/j.gca.2014.06.006
}

RM_NBS_18 = {
	'name': 'NBS-18',
	'd13C_VPDB': {'value': -5.014, 'sigma': 0.035},
	'd18O_VPDB': {'value': -23.2, 'sigma': 0.1},
	'D17O_VSMOW': {
		'value': -0.048,
		'sigma': 0.010,
	},  # Wostbrock et al. (2020) 10.1016/j.chemgeo.2019.119432
}

RM_NBS_18_CO2_25C = {
	'name': 'NBS-18-CO2-25C',
	'd13C_VPDB': {'value': -5.014, 'sigma': 0.035},
	'd18O_VSMOW': {'value': (1000 - 23.2) * 1.03092 * 1.01025 - 1000, 'sigma': 0.1},
	'D17O_VSMOW': {
		'value': -0.163,
		'sigma': 0.010,
	},  # Barkan et al. (2019) 10.1002/rcm.8391
}

RM_NBS_18_CO2_90C = {
	'name': 'NBS-18-CO2-90C',
	'd13C_VPDB': {'value': -5.014, 'sigma': 0.0},
	'd18O_VSMOW': {'value': (1000 - 23.2) * 1.03092 * 1.008146 - 1000, 'sigma': 0.0},
	'D17O_VSMOW': {
		'value': -0.098,
		'sigma': 0.010,
	},  # Passey et al. (2014) 10.1016/j.gca.2014.06.006
}

RM_IAEA_603 = {
	'name': 'IAEA-603',
	'd13C_VPDB': {'value': +2.46, 'sigma': 0.01},
	'd18O_VPDB': {'value': -2.37, 'sigma': 0.04},
	'D17O_VSMOW': {
		'value': -0.100,
		'sigma': 0.010,
	},  # Wostbrock et al. (2020) 10.1016/j.chemgeo.2019.119432
}

RM_IAEA_603_CO2_25C = {
	'name': 'IAEA-603-CO2-25C',
	'd13C_VPDB': {'value': +2.46, 'sigma': 0.01},
	'd18O_VSMOW': {'value': (1000 - 2.37) * 1.03092 * 1.01025 - 1000, 'sigma': 0.04},
	'D17O_VSMOW': {
		'value': -0.194,
		'sigma': 0.010,
	},  # Barkan et al. (2019) 10.1002/rcm.8391
	'D17O_VSMOW': {
		'value': -0.147,
		'sigma': 0.010,
	},  # Wostbrock et al. (2020) 10.1016/j.chemgeo.2019.119432
}

RM_IAEA_610 = {
	'name': 'IAEA-610',
	'd13C_VPDB': {'value': -9.109, 'sigma': 0.012},
	'd18O_VPDB': {'value': -18.834, 'sigma': 0.045},
}

RM_IAEA_610_CO2_90C = {
	'name': 'IAEA-610-CO2-90C',
	'd13C_VPDB': {'value': -9.109, 'sigma': 0.012},
	'd18O_VSMOW': {
		'value': (1000 - 18.834) * 1.03092 * 1.008146 - 1000,
		'sigma': 0.045,
	},
}

RM_IAEA_610_CO2_25C = {
	'name': 'IAEA-610-CO2-25C',
	'd13C_VPDB': {'value': -9.109, 'sigma': 0.012},
	'd18O_VSMOW': {'value': (1000 - 18.834) * 1.03092 * 1.01025 - 1000, 'sigma': 0.045},
}

RM_IAEA_611 = {
	'name': 'IAEA-611',
	'd13C_VPDB': {'value': -30.795, 'sigma': 0.013},
	'd18O_VPDB': {'value': -4.224, 'sigma': 0.046},
}

RM_IAEA_611_CO2_90C = {
	'name': 'IAEA-611-CO2-90C',
	'd13C_VPDB': {'value': -30.795, 'sigma': 0.013},
	'd18O_VSMOW': {'value': (1000 - 4.224) * 1.03092 * 1.008146 - 1000, 'sigma': 0.046},
}

RM_IAEA_611_CO2_25C = {
	'name': 'IAEA-611-CO2-25C',
	'd13C_VPDB': {'value': -30.795, 'sigma': 0.013},
	'd18O_VSMOW': {'value': (1000 - 4.224) * 1.03092 * 1.01025 - 1000, 'sigma': 0.046},
}

RM_IAEA_612 = {
	'name': 'IAEA-612',
	'd13C_VPDB': {'value': -36.722, 'sigma': 0.015},
	'd18O_VPDB': {'value': -12.079, 'sigma': 0.062},
}

RM_IAEA_612_CO2_90C = {
	'name': 'IAEA-612-CO2-90C',
	'd13C_VPDB': {'value': -36.722, 'sigma': 0.015},
	'd18O_VSMOW': {
		'value': (1000 - 12.079) * 1.03092 * 1.008146 - 1000,
		'sigma': 0.062,
	},  # 26.764
}

RM_IAEA_612_CO2_25C = {
	'name': 'IAEA-612-CO2-25C',
	'd13C_VPDB': {'value': -36.722, 'sigma': 0.015},
	'd18O_VSMOW': {'value': (1000 - 12.079) * 1.03092 * 1.01025 - 1000, 'sigma': 0.062},
}


def sanitize(x):
	return x.replace('-', '_').replace('.', '_')


def standardize(
	data,
	anchors,
	key_in,
	key_out = 'd13C_VPDB',
	constraints = {},
	method = 'correct_observations',  # 'predict_observations'
	verbose = False,
):
	if verbose:
		if method == 'correct_observations':
			print('Using "correct_observations" method.')
		if method == 'predict_observations':
			print('Using "predict_observations" method.')

	out = {}

	sessions = sorted({r['Session'] for r in data})
	samples = sorted({r['Sample'] for r in data})
	usamples = [s for s in samples if s not in anchors]
	asamples = [s for s in samples if s in anchors]

	out['data'] = [r.copy() for r in data]

	fitparams = lmfit.Parameters()

	for s in sessions:
		fitparams.add('delta_scaling_' + sanitize(s), value = 1.0)
		fitparams.add('delta_of_wg_' + sanitize(s), value = 0.0)
	for s in usamples:
		fitparams.add('x_' + sanitize(s), value = 0.0)

	for p in fitparams:
		if p in constraints:
			fitparams[p].expr = constraints[p]

	def observations(r):
		return r[key_in]

	def truevalues(p, r):
		if r['Sample'] in anchors:
			return anchors[r['Sample']][key_out]
		else:
			return p[f"x_{sanitize(r['Sample'])}"]

	def correct_observations(p, r):
		delta_obs = observations(r)
		delta_corrected = 1e3 * (
			(1 + delta_obs / 1e3 / p[f"delta_scaling_{sanitize(r['Session'])}"])
			* (1 + p[f"delta_of_wg_{sanitize(r['Session'])}"] / 1e3)
			- 1
		)
		return delta_corrected

	def predict_observations(p, r):
		delta_true = truevalues(p, r)
		delta_predicted = (
			(
				(1 + delta_true / 1e3)
				/ (1 + p[f"d18O_of_wg_{sanitize(r['Session'])}"] / 1e3)
				- 1
			)
			* 1e3
			* p[f"d18O_scaling_{sanitize(r['Session'])}"]
		)
		return delta_predicted

	def prediction_residuals(p):
		return array([observations(r) - predict_observations(p, r) for r in data])

	def correction_residuals(p):
		return array([correct_observations(p, r) - truevalues(p, r) for r in data])

	residuals = {
		'correct_observations': correction_residuals,
		'predict_observations': prediction_residuals,
	}[method]

	fitresult = lmfit.minimize(
		residuals, fitparams, method = 'least_squares', scale_covar = True
	)

	out['bestfit'] = fitresult
	out['fitreport'] = lmfit.fit_report(fitresult)

	out['Nf'] = fitresult.nfree
	out['t95'] = tstudent.ppf(1 - 0.05 / 2, out['Nf'])

	if verbose:
		print(f"Student's t-factor is {out['t95']:.2f} (p = 0.95, Nf = {out['Nf']})")

	p = {_: fitresult.params[_].value for _ in fitresult.params}
	for r in out['data']:
		r[f'{key_out}_corrected'] = correct_observations(p, r)
		r[f'{key_out}_true'] = truevalues(p, r)

		r[f'{key_out}_residual'] = r[f'{key_out}_corrected'] - r[f'{key_out}_true']

	out['sessions'] = {}
	for s in sessions:
		_s = sanitize(s)
		out['sessions'][s] = {}

		out['sessions'][s]['N'] = len([r for r in data if r['Session'] == s])
		out['sessions'][s]['Na'] = len(
			[r for r in data if r['Session'] == s and r['Sample'] in anchors]
		)
		out['sessions'][s]['Nu'] = len(
			[r for r in data if r['Session'] == s and r['Sample'] not in anchors]
		)
		out['sessions'][s]['data'] = [r for r in out['data'] if r['Session'] == s]

		out['sessions'][s][f'{key_out}_scaling'] = fitresult.params[
			f'delta_scaling_' + _s
		].value
		out['sessions'][s][f'{key_out}_of_wg'] = fitresult.params[
			'delta_of_wg_' + _s
		].value

		out['sessions'][s][f'SE_{key_out}_scaling'] = fitresult.params[
			'delta_scaling_' + _s
		].stderr
		out['sessions'][s][f'SE_{key_out}_of_wg'] = fitresult.params[
			'delta_of_wg_' + _s
		].stderr

		out['sessions'][s][f'RMSE_{key_out}'] = (
			array([r[f'{key_out}_residual'] for r in out['sessions'][s]['data']]) ** 2
		).mean() ** 0.5

	out['samples'] = {}
	for s in samples:
		out['samples'][s] = {}

		out['samples'][s]['N'] = len([r for r in data if r['Sample'] == s])
		out['samples'][s]['data'] = [r for r in out['data'] if r['Sample'] == s]

	for s in asamples:
		out['samples'][s][key_out] = anchors[s][key_out]
		out['samples'][s]['SD_' + key_out] = array(
			[r[f'{key_out}_residual'] for r in out['samples'][s]['data']]
		).std(ddof = 0)

	for s in usamples:
		_s = sanitize(s)
		out['samples'][s][key_out] = fitresult.params['x_' + _s].value
		out['samples'][s]['SD_' + key_out] = array(
			[r[f'{key_out}_residual'] for r in out['samples'][s]['data']]
		).std(ddof = 1)
		out['samples'][s]['SE_' + key_out] = fitresult.params['x_' + _s].stderr
		out['samples'][s]['95CL_' + key_out] = (
			out['samples'][s]['SE_' + key_out] * out['t95']
		)

	out['anchors'] = asamples[:]
	out['unknowns'] = usamples[:]

	csv = f'Session,N,Na,Nu,{key_out}_scaling,SE_{key_out}_scaling,{key_out}_of_wg,SE_{key_out}_of_wg,RMSE_{key_out}'
	for s in sessions:
		csv += f'\n{s}'
		csv += f',{out["sessions"][s]["N"]}'
		csv += f',{out["sessions"][s]["Na"]}'
		csv += f',{out["sessions"][s]["Nu"]}'
		csv += f',{out["sessions"][s][key_out+"_scaling"]:.4f}'
		csv += f',{out["sessions"][s]["SE_"+key_out+"_scaling"]:.4f}'
		csv += f',{out["sessions"][s][key_out+"_of_wg"]:.3f}'
		csv += f',{out["sessions"][s]["SE_"+key_out+"_of_wg"]:.3f}'
		csv += f',{out["sessions"][s]["RMSE_"+key_out]:.3f}'
	out['csv_sessions'] = csv

	return out


def standardize_D17O(
	data,
	anchors,
	d18O_key_in,
	d17O_key_in,
	d18O_key_out = 'd18O_VSMOW',
	d17O_key_out = 'd17O_VSMOW',
	constraints = {},
	method = 'correct_observations',  # 'predict_observations'
	residuals_17='D17',  # 'D17'
	robust_cov_estimator = True,
	verbose = False,
):
	"""
	Standardize δ18O, δ17O, and Δ17O values using a pooled LS regression approach.

	Arguments:
	----------
	* `data`: a list of dicts, with each list item corresponding to one analysis
	* `anchors`: a dict with the following structure:

	```
	anchors = {
		'FIRST_ANCHOR':
		}
	```
	"""

	D17O_key_out = d17O_key_out.replace('d17', 'D17')

	if verbose:
		if residuals_17 == 'D17':
			print('Using D17 residuals.')
		elif residuals_17 == 'd17':
			print('Using d17 residuals.')
		if method == 'correct_observations':
			print('Using "correct_observations" method.')
		if method == 'predict_observations':
			print('Using "predict_observations" method.')

	out = {}

	sessions = sorted({r['Session'] for r in data})
	samples = sorted({r['Sample'] for r in data})
	usamples = [s for s in samples if s not in anchors]
	asamples = [s for s in samples if s in anchors]

	out['data'] = [r.copy() for r in data]

	fitparams = lmfit.Parameters()

	for s in sessions:
		fitparams.add('d18O_scaling_' + sanitize(s), value = 1.0)
		fitparams.add('d18O_of_wg_' + sanitize(s), value = 0.0)
		fitparams.add('d17O_scaling_' + sanitize(s), value = 1.0)
		fitparams.add('d17O_of_wg_' + sanitize(s), value = 0.0)

	for s in usamples:
		fitparams.add('d18O_' + sanitize(s), value = 0.0)
		if residuals_17 == 'D17':
			fitparams.add(f'D17O_' + sanitize(s), value = 0.0)
		elif residuals_17 == 'd17':
			fitparams.add(f'd17O_' + sanitize(s), value = 0.0)

	for p in fitparams:
		if p in constraints:
			fitparams[p].expr = constraints[p]

	def observations(r):
		d18obs = r[d18O_key_in]
		d17obs = r[d17O_key_in]
		D17obs = 1e3 * (
			log(1 + d17obs / 1e3) - ISOPARAMETERS.LAMBDA_17 * log(1 + d18obs / 1e3)
		)

		return (d18obs, d17obs, D17obs)

	def truevalues(p, r):
		if r['Sample'] in anchors:
			d18true = anchors[r['Sample']][d18O_key_out]
			d17true = anchors[r['Sample']][d17O_key_out]
			D17true = 1e3 * (
				log(1 + d17true / 1e3)
				- ISOPARAMETERS.LAMBDA_17 * log(1 + d18true / 1e3)
			)
		else:
			d18true = p[f"d18O_{sanitize(r['Sample'])}"]
			if residuals_17 == 'D17':
				D17true = p[f"D17O_{sanitize(r['Sample'])}"]
				d17true = 1e3 * (
					exp(D17true / 1e3) * (1 + d18true / 1e3) ** ISOPARAMETERS.LAMBDA_17
					- 1
				)
			elif residuals_17 == 'd17':
				d17true = p[f"d17O_{sanitize(r['Sample'])}"]
				D17true = 1e3 * (
					log(1 + d17true / 1e3)
					- ISOPARAMETERS.LAMBDA_17 * log(1 + d18true / 1e3)
				)

		return (d18true, d17true, D17true)

	def correct_observations(p, r):
		d18obs, d17obs, D17obs = observations(r)

		d18corrected = 1e3 * (
			(1 + d18obs / 1e3 / p[f"d18O_scaling_{sanitize(r['Session'])}"])
			* (1 + p[f"d18O_of_wg_{sanitize(r['Session'])}"] / 1e3)
			- 1
		)
		d17corrected = 1e3 * (
			(1 + d17obs / 1e3 / p[f"d17O_scaling_{sanitize(r['Session'])}"])
			* (1 + p[f"d17O_of_wg_{sanitize(r['Session'])}"] / 1e3)
			- 1
		)
		D17corrected = 1e3 * (
			log(1 + d17corrected / 1e3)
			- ISOPARAMETERS.LAMBDA_17 * log(1 + d18corrected / 1e3)
		)

		return (d18corrected, d17corrected, D17corrected)

	def predict_observations(p, r):
		d18true, d17true, D17true = truevalues(p, r)

		d18predicted = (
			(
				(1 + d18true / 1e3)
				/ (1 + p[f"d18O_of_wg_{sanitize(r['Session'])}"] / 1e3)
				- 1
			)
			* 1e3
			* p[f"d18O_scaling_{sanitize(r['Session'])}"]
		)
		d17predicted = (
			(
				(1 + d17true / 1e3)
				/ (1 + p[f"d17O_of_wg_{sanitize(r['Session'])}"] / 1e3)
				- 1
			)
			* 1e3
			* p[f"d17O_scaling_{sanitize(r['Session'])}"]
		)
		D17predicted = 1e3 * (
			log(1 + d17predicted / 1e3)
			- ISOPARAMETERS.LAMBDA_17 * log(1 + d18predicted / 1e3)
		)

		return (d18predicted, d17predicted, D17predicted)

	def prediction_residuals(p, residuals_17=residuals_17, ChM = None):
		R = []
		for r in data:
			d18predicted, d17predicted, D17predicted = predict_observations(p, r)
			d18obs, d17obs, D17obs = observations(r)

			R.append(d18obs - d18predicted)
			if residuals_17 == 'D17':
				R.append(D17obs - D17predicted)
			elif residuals_17 == 'd17':
				R.append(d17obs - d17predicted)

		R = array(R)

		if ChM is None:
			return R

		R = R.reshape((R.size, 1))

		return ChM @ R

	def correction_residuals(p, residuals_17=residuals_17, ChM = None):
		R = []
		for r in data:
			d18true, d17true, D17true = truevalues(p, r)
			d18corrected, d17corrected, D17corrected = correct_observations(p, r)

			R.append(d18corrected - d18true)
			if residuals_17 == 'D17':
				R.append(D17corrected - D17true)
			elif residuals_17 == 'd17':
				R.append(d17corrected - d17true)

		R = array(R)

		if ChM is None:
			return R

		R = R.reshape((R.size, 1))

		return ChM @ R

	residuals = {
		'correct_observations': correction_residuals,
		'predict_observations': prediction_residuals,
	}[method]

	fitresult = lmfit.minimize(
		residuals, fitparams, method = 'least_squares', scale_covar = True
	)

	for _ in range(2):  # overkill, one pass suffices
		CM = estimate_covariance(
			residuals(fitresult.params).reshape((len(data), 2)),
			robust = robust_cov_estimator,
		)

		# 		if verbose:
		# 			print(f'Residuals: SE = {1e3*CM[0,0]**.5:.1f} ppm, {1e3*CM[1,1]**.5:.1f} ppm, correl = {CM[0,1]/(CM[0,0]*CM[1,1])**.5:.4f}')

		ChM = scipy.linalg.block_diag(
			*[scipy.linalg.cholesky(scipy.linalg.inv(CM))] * len(data)
		)
		fitresult = lmfit.minimize(
			residuals,
			fitparams,
			method = 'least_squares',
			scale_covar = False,
			kws = {'ChM': ChM},
		)

	if verbose:
		CM = estimate_covariance(
			residuals(fitresult.params).reshape((len(data), 2)),
			robust = robust_cov_estimator,
		)
		print(
			f'Residuals: SE = {1e3*CM[0,0]**.5:.1f} ppm, {1e3*CM[1,1]**.5:.1f} ppm, correl = {CM[0,1]/(CM[0,0]*CM[1,1])**.5:.4f}'
		)

	out['bestfit'] = fitresult
	out['fitreport'] = lmfit.fit_report(fitresult)

	out['Nf'] = fitresult.nfree // 2
	out['t95'] = tstudent.ppf(1 - 0.05 / 2, out['Nf'])

	if verbose:
		print(f"Student's t-factor is {out['t95']:.2f} (p = 0.95, Nf = {out['Nf']})")

	p = {_: fitresult.params[_].value for _ in fitresult.params}
	for r in out['data']:
		r['d18corrected'], r['d17corrected'], r['D17corrected'] = correct_observations(
			p, r
		)
		r['d18true'], r['d17true'], r['D17true'] = truevalues(p, r)

		r['d18residual'] = r['d18corrected'] - r['d18true']
		r['d17residual'] = r['d17corrected'] - r['d17true']
		r['D17residual'] = r['D17corrected'] - r['D17true']

	out['sessions'] = {}
	for s in sessions:
		_s = sanitize(s)
		out['sessions'][s] = {}

		out['sessions'][s]['N'] = len([r for r in data if r['Session'] == s])
		out['sessions'][s]['Na'] = len(
			[r for r in data if r['Session'] == s and r['Sample'] in anchors]
		)
		out['sessions'][s]['Nu'] = len(
			[r for r in data if r['Session'] == s and r['Sample'] not in anchors]
		)
		out['sessions'][s]['data'] = [r for r in out['data'] if r['Session'] == s]

		out['sessions'][s]['d18O_scaling'] = fitresult.params[
			'd18O_scaling_' + _s
		].value
		out['sessions'][s]['d18O_of_wg'] = fitresult.params['d18O_of_wg_' + _s].value
		out['sessions'][s]['d17O_scaling'] = fitresult.params[
			'd17O_scaling_' + _s
		].value
		out['sessions'][s]['d17O_of_wg'] = fitresult.params['d17O_of_wg_' + _s].value

		out['sessions'][s]['SE_d18O_scaling'] = fitresult.params[
			'd18O_scaling_' + _s
		].stderr
		out['sessions'][s]['SE_d18O_of_wg'] = fitresult.params[
			'd18O_of_wg_' + _s
		].stderr
		out['sessions'][s]['SE_d17O_scaling'] = fitresult.params[
			'd17O_scaling_' + _s
		].stderr
		out['sessions'][s]['SE_d17O_of_wg'] = fitresult.params[
			'd17O_of_wg_' + _s
		].stderr

		out['sessions'][s]['RMSE_d18O'] = (
			array([r['d18residual'] for r in out['sessions'][s]['data']]) ** 2
		).mean() ** 0.5
		out['sessions'][s]['RMSE_d17O'] = (
			array([r['d17residual'] for r in out['sessions'][s]['data']]) ** 2
		).mean() ** 0.5
		out['sessions'][s]['RMSE_D17O'] = (
			array([r['D17residual'] for r in out['sessions'][s]['data']]) ** 2
		).mean() ** 0.5

	out['samples'] = {}
	for s in samples:
		out['samples'][s] = {}

		out['samples'][s]['N'] = len([r for r in data if r['Sample'] == s])
		out['samples'][s]['data'] = [r for r in out['data'] if r['Sample'] == s]

	for s in asamples:
		out['samples'][s][d18O_key_out] = anchors[s][d18O_key_out]
		out['samples'][s][d17O_key_out] = anchors[s][d17O_key_out]
		out['samples'][s][D17O_key_out] = 1e3 * (
			log(1 + out['samples'][s][d17O_key_out] / 1e3)
			- ISOPARAMETERS.LAMBDA_17 * log(1 + out['samples'][s][d18O_key_out] / 1e3)
		)

	for s in usamples:
		_s = sanitize(s)
		out['samples'][s][d18O_key_out] = fitresult.params['d18O_' + _s].value
		if residuals_17 == 'D17':
			out['samples'][s][D17O_key_out] = fitresult.params['D17O_' + _s].value
			out['samples'][s][d17O_key_out] = 1e3 * (
				exp(out['samples'][s][D17O_key_out] / 1e3)
				* (1 + out['samples'][s][d18O_key_out] / 1e3) ** ISOPARAMETERS.LAMBDA_17
				- 1
			)
		elif residuals_17 == 'd17':
			out['samples'][s][d17O_key_out] = fitresult.params['d17O_' + _s].value
			out['samples'][s][D17O_key_out] = 1e3 * (
				log(1 + out['samples'][s][d17O_key_out] / 1e3)
				- ISOPARAMETERS.LAMBDA_17
				* log(1 + out['samples'][s][d18O_key_out] / 1e3)
			)

		out['samples'][s]['SE_' + d18O_key_out] = fitresult.params['d18O_' + _s].stderr
		out['samples'][s]['95CL_' + d18O_key_out] = (
			out['samples'][s]['SE_' + d18O_key_out] * out['t95']
		)

		if residuals_17 == 'D17':
			out['samples'][s]['SE_' + D17O_key_out] = fitresult.params[
				'D17O_' + _s
			].stderr
			out['samples'][s]['95CL_' + D17O_key_out] = (
				out['samples'][s]['SE_' + D17O_key_out] * out['t95']
			)

			out['samples'][s]['SE_' + d17O_key_out] = 0.0
			out['samples'][s]['95CL_' + d17O_key_out] = 0.0

			try:
				i = fitresult.var_names.index('d18O_' + _s)
				j = fitresult.var_names.index('D17O_' + _s)
				CM = fitresult.covar[ix_([i, j], [i, j])]

				J = array(
					[
						[
							exp(out['samples'][s][D17O_key_out] / 1e3)
							* ISOPARAMETERS.LAMBDA_17
							* (1 + out['samples'][s][d18O_key_out] / 1e3)
							** (ISOPARAMETERS.LAMBDA_17 - 1)
						],
						[
							exp(out['samples'][s][D17O_key_out] / 1e3)
							* (1 + out['samples'][s][d18O_key_out] / 1e3)
							** ISOPARAMETERS.LAMBDA_17
						],
					]
				)

				out['samples'][s]['SE_' + d17O_key_out] = float(J.T @ CM @ J) ** 0.5
				out['samples'][s]['95CL_' + d17O_key_out] = (
					out['samples'][s]['SE_' + d17O_key_out] * out['t95']
				)
			except ValueError:
				pass

		elif residuals_17 == 'd17':
			out['samples'][s]['SE_' + d17O_key_out] = fitresult.params[
				'd17O_' + _s
			].stderr
			out['samples'][s]['95CL_' + d17O_key_out] = (
				out['samples'][s]['SE_' + d17O_key_out] * out['t95']
			)

			out['samples'][s]['SE_' + D17O_key_out] = 0.0
			out['samples'][s]['95CL_' + D17O_key_out] = 0.0

			try:
				i = fitresult.var_names.index('d18O_' + _s)
				j = fitresult.var_names.index('d17O_' + _s)
				CM = fitresult.covar[ix_([i, j], [i, j])]
				J = array(
					[
						[
							-1
							/ (1 + out['samples'][s][d18O_key_out] / 1e3)
							* ISOPARAMETERS.LAMBDA_17
						],
						[1 / (1 + out['samples'][s][d17O_key_out] / 1e3)],
					]
				)

				out['samples'][s]['SE_' + D17O_key_out] = float(J.T @ CM @ J) ** 0.5
				out['samples'][s]['95CL_' + D17O_key_out] = (
					out['samples'][s]['SE_' + D17O_key_out] * out['t95']
				)
			except ValueError:
				pass

	for s in samples:
		out['samples'][s]['SD_' + d18O_key_out] = array(
			[r[f'd18residual'] for r in out['samples'][s]['data']]
		).std(ddof = 0 if s in asamples else 1)
		out['samples'][s]['SD_' + d17O_key_out] = array(
			[r[f'd17residual'] for r in out['samples'][s]['data']]
		).std(ddof = 0 if s in asamples else 1)
		out['samples'][s]['SD_' + D17O_key_out] = array(
			[r[f'D17residual'] for r in out['samples'][s]['data']]
		).std(ddof = 0 if s in asamples else 1)

	out['anchors'] = asamples[:]
	out['unknowns'] = usamples[:]

	csv = f'Session,N,Na,Nu,{d18O_key_out}_scaling,SE_{d18O_key_out}_scaling,{d17O_key_out}_scaling,SE_{d17O_key_out}_scaling,{d18O_key_out}_of_wg,SE_{d18O_key_out}_of_wg,{d17O_key_out}_of_wg,SE_{d17O_key_out}_of_wg,RMSE_{d18O_key_out},RMSE_{d17O_key_out},RMSE_{D17O_key_out}'
	for s in sessions:
		csv += f'\n{s}'
		csv += f',{out["sessions"][s]["N"]}'
		csv += f',{out["sessions"][s]["Na"]}'
		csv += f',{out["sessions"][s]["Nu"]}'
		csv += f',{out["sessions"][s]["d18O_scaling"]:.4f}'
		csv += f',{out["sessions"][s]["SE_d18O_scaling"]:.4f}'
		csv += f',{out["sessions"][s]["d17O_scaling"]:.4f}'
		csv += f',{out["sessions"][s]["SE_d17O_scaling"]:.4f}'
		csv += f',{out["sessions"][s]["d18O_of_wg"]:.3f}'
		csv += f',{out["sessions"][s]["SE_d18O_of_wg"]:.3f}'
		csv += f',{out["sessions"][s]["d17O_of_wg"]:.3f}'
		csv += f',{out["sessions"][s]["SE_d17O_of_wg"]:.3f}'
		csv += f',{out["sessions"][s]["RMSE_d18O"]:.3f}'
		csv += f',{out["sessions"][s]["RMSE_d17O"]:.3f}'
		csv += f',{out["sessions"][s]["RMSE_D17O"]:.4f}'
	out['csv_sessions'] = csv

	return out
