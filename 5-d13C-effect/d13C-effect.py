import sys

sys.path.append('../lib')

from mylib import *

RAWDATA_FILENAME = 'input/d13C-effect-rawdata.csv'

logger.remove(0)
logger.add('output/d13C-effect.log', mode = 'w', format = '[{line:3.0f}]    {message}')

ppl.style.use('../lib/mydefault.mplstyle')


def H2O_eq_CO2_predictions(starting_CO2, eq_waters):
	"""PREDICT COMPOSITIONS OF EQUILIBRATED CO2 SAMPLES"""

	predictions = {}
	for w in eq_waters:
		d18, D17 = isofunctions.D17_CO2_H2O_eq(
			eq_waters[w]['x'],
			starting_CO2['d18O_VSMOW'],
			starting_CO2['D17O_VSMOW'],
			eq_waters[w]['d18O_VSMOW'],
			eq_waters[w]['D17O_VSMOW'],
		)[:2]
		predictions[w] = {
			'd18O_VSMOW': d18,
			'd17O_VSMOW': (np.exp(D17 / 1e3) * (1 + d18 / 1e3) ** LAMBDA_17 - 1) * 1e3,
			'D17O_VSMOW': D17,
		}
	return predictions


def H2O_eq_CO2_D17_corrections(starting_CO2, eq_waters):
	"""COMPUTE CORRECTIONS DUE TO FINITE H2O/CO2 RATIO"""
	D17_corrections = {
		w: isofunctions.D17_CO2_H2O_eq(
			eq_waters[w]['x'],
			starting_CO2['d18O_VSMOW'],
			starting_CO2['D17O_VSMOW'],
			eq_waters[w]['d18O_VSMOW'],
			eq_waters[w]['D17O_VSMOW'],
		)[1]
		- eq_waters[w]['D17O_VSMOW']
		- 1e3 * (np.log(alpha17_CO2g_H2O) - LAMBDA_17 * np.log(alpha18_CO2g_H2O))
		for w in eq_waters
	}

	for w in D17_corrections:
		logger.info(
			f'Δ17O correction is {D17_corrections[w]*1e3:+.1f} ppm for CO2 equilibrated with {w}.'
		)

	d18_corrections = {
		w: isofunctions.D17_CO2_H2O_eq(
			eq_waters[w]['x'],
			starting_CO2['d18O_VSMOW'],
			starting_CO2['D17O_VSMOW'],
			eq_waters[w]['d18O_VSMOW'],
			eq_waters[w]['D17O_VSMOW'],
		)[0]
		- ((1000 + eq_waters[w]['d18O_VSMOW']) * alpha18_CO2g_H2O - 1000)
		for w in eq_waters
	}

	return d18_corrections, D17_corrections


def standardize_triple_oxygen(data, anchors, constraints = {}):
	stdz_triple_oxygen = stdz.standardize_D17O(
		data, anchors, 'd628', 'd627', constraints = constraints
	)

	for S, key, fn, ss in [
		(stdz_triple_oxygen, 'd18O_VSMOW', 'residuals_d18O', ''),
		(stdz_triple_oxygen, 'D17O_VSMOW', 'residuals_D17O', 'sessions_triple_oxygen'),
	]:
		isofunctions.plot_residuals(
			S['data'],
			key = key,
			filename = f'output/{fn}',
			anchors = [_ for _ in S['anchors']],
		)
		if ss:
			with open(f'output/{ss}.csv', 'w') as fid:
				fid.write(S['csv_sessions'])

	return stdz_triple_oxygen


def waters(S, group, starting_CO2, filename = 'output/water_samples.csv'):
	co2_samples = [
		{
			'Sample': s,
			'd18O_VSMOW': S['samples'][s]['d18O_VSMOW'],
			'D17O_VSMOW': S['samples'][s]['D17O_VSMOW'],
			'H2O/CO2': float(s.split('_x')[-1]),
		}
		for s in S[group]
	]

	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		isofunctions.correct_CO2eqH2O_back_to_initial_H2O(
			co2_samples, starting_CO2['d18O_VSMOW'], starting_CO2['D17O_VSMOW']
		)

	with open(filename, 'w') as fid:
		fid.write('Sample,water_d18O_VSMOW,water_D17O_VSMOW')
		for s in co2_samples:
			fid.write(
				f'\n{s["Sample"]},{s["water_d18O_VSMOW"]:.3f},{s["water_D17O_VSMOW"]:.6f}'
			)


def plot_d13C_effect():
	X = [ref_waters[s]['d18O_VSMOW'] for s in ['VSMOW2', 'SLAP2']]
	with open('output/water_samples_dCO2.csv') as fid:
		data = [l.strip().split(',') for l in fid.readlines()]
	data = [{k: v for k, v in zip(data[0], l)} for l in data[1:]]
	data = {r['Sample']: r for r in data}
	Y = np.array(
		[
			float(data[s]['water_D17O_VSMOW']) * 1000
			for s in ['C13D-VSMOW2_x387', 'C13D-SLAP2_x340']
		]
	)
	eY = np.array(
		[
			stdz_triple_oxygen['samples'][s]['95CL_D17O_VSMOW'] * 1000
			for s in ['C13D-VSMOW2_x387', 'C13D-SLAP2_x340']
		]
	)

	fig = ppl.figure(figsize = (3.15, 2))
	fig.subplots_adjust(0.2, 0.05, 0.97, 0.95)
	ax = ppl.subplot(111)
	kwe = dict(
		ls = 'None', marker = 'None', ecolor = 'k', elinewidth = 1, capsize = 3, capthick = 1
	)
	kwm = dict(ls = 'None', marker = 'o', mfc = 'w', mec = 'k', mew = 1, ms = 6)
	ax.axhline(0, color = 'k', lw = 1.5, dashes = (6, 3), alpha = 0.25)
	ax.errorbar(X, Y, eY, **kwe)
	ax.plot(X, Y, **kwm)
	for x, y, ey, t in zip(X, Y, eY, ['VSMOW2', 'SLAP2']):
		ax.text(
			x + np.sign(x + 25) * 20,
			y,
			f'{t}\n{y:.1f} ± {ey:.1f} ppm',
			va = 'bottom' if y > 0 else 'top',
			ha = 'right' if x < 0 else 'left',
			size = 9,
			weight = 'bold',
		)

	ax.set_yticks([-8, -4, 0, 4, 8])
	ax.set_xticks([])
	ax.axis([-255, 200, -9, 9])
	ax.set_ylabel('Δ$^{17}Ο_{VSMOW}$ (ppm)')
	fig.savefig('output/d13C_effect')
	ppl.close(fig)


if __name__ == '__main__':
	"""COMPOSITION OF EQUILIBRATED WATERS"""
	ref_waters = {
		'VSMOW2': {'d18O_VSMOW': 0.0, 'D17O_VSMOW': 0.000},
		'SLAP2': {'d18O_VSMOW': -55.5, 'D17O_VSMOW': 0.000},
	}

	"""LOAD VCOF-CRDS ANALYSES"""
	data = read_rawdata([RAWDATA_FILENAME])
	data = [r for r in data if r['UID'] not in ['833D54CF']]  # problematic outlier

	"""DEFINE WATER-EQUILIBRATED STANDARDS"""
	unique_combinations = {
		(r['Sample'], r['H2O_to_CO2'])
		for r in data
		if r['Sample'] in ['VSMOW2', 'SLAP2', 'C13D-VSMOW2', 'C13D-SLAP2']
	}
	eq_waters = {
		f'{w}_x{x:.0f}': {
			'd18O_VSMOW': ref_waters[w.split('C13D-')[-1]]['d18O_VSMOW'],
			'D17O_VSMOW': ref_waters[w.split('C13D-')[-1]]['D17O_VSMOW'],
			'x': x,
		}
		for w, x in unique_combinations
	}
	anchors = {w: eq_waters[w] for w in eq_waters if not w.startswith('C13D')}
	c13d = {w: eq_waters[w] for w in eq_waters if w.startswith('C13D')}
	for r in data:
		if r['Sample'] not in ['dCO2']:
			r['Sample'] = r['Sample'] + f'_x{r["H2O_to_CO2"]:.0f}'

	POSTULATED_SLAP_D17 = 0  # Fool around by turning this knob
	logger.info(f'Postulated Δ17O_VSMOW value of SLAP = {POSTULATED_SLAP_D17} ‰')

	if POSTULATED_SLAP_D17:  # Adjust D17O values in the case of non-zero SLAP value
		for w in eq_waters:
			eq_waters[w]['D17O_VSMOW'] += (
				log(1 + eq_waters[w]['d18O_VSMOW'] / 1e3)
				/ log(1 - 55.5e-3)
				* POSTULATED_SLAP_D17
			)

	"""COMPOSITION CO2 USED WITH EQUILIBRATED WATERS"""
	match CO2g_H2O_FRACTIONATION:
		case 'Guo & Zhou (2019)':
			starting_CO2 = {
				'eCO2': {'d18O_VSMOW': 24.71, 'D17O_VSMOW': -0.0835},
				'dCO2': {'d18O_VSMOW': 10.71, 'D17O_VSMOW': -0.1381},
			}
		case 'Barkan & Luz (2012)':
			starting_CO2 = {
				'eCO2': {'d18O_VSMOW': 24.25, 'D17O_VSMOW': -0.150},
				'dCO2': {'d18O_VSMOW': 10.30, 'D17O_VSMOW': -0.2053},
			}

	for co2 in ['eCO2', 'dCO2']:
		isofunctions.H2O_eq_CO2_figure(
			starting_CO2[co2],
			ref_waters,
			40e-6,
			40e-6 * 312,
			H2O_to_CO2_ratio_half_width = 1 / 3,
			filename = f'output/H2O_eq_CO2_D17_vs_ratio_{co2}',
		)

	"""H2O_eq_CO2 CORRECTIONS"""
	d18_corrections, D17_corrections = H2O_eq_CO2_D17_corrections(
		starting_CO2['eCO2'], anchors
	)
	c13d_d18_corrections, c13d_D17_corrections = H2O_eq_CO2_D17_corrections(
		starting_CO2['dCO2'], c13d
	)
	d18_corrections = d18_corrections | c13d_d18_corrections
	D17_corrections = D17_corrections | c13d_D17_corrections
	predictions = H2O_eq_CO2_predictions(starting_CO2['eCO2'], anchors)

	"""FLAG d626b OUTLIERS"""
	# first flag them without excluding them:
	data, p_grubbs = isofunctions.d626b_outliers(data, rename = False)
	# standardize the data with these outliers still included:
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')

		stdz_triple_oxygen = standardize_triple_oxygen(
			data,
			predictions,
			constraints = {
				'D17O_C13D_SLAP2_x396': f'D17O_C13D_SLAP2_x340 - {D17_corrections["C13D-SLAP2_x340"]} + {D17_corrections["C13D-SLAP2_x396"]}',
				'D17O_C13D_SLAP2_x438': f'D17O_C13D_SLAP2_x340 - {D17_corrections["C13D-SLAP2_x340"]} + {D17_corrections["C13D-SLAP2_x438"]}',
				'D17O_C13D_VSMOW2_x438': f'D17O_C13D_VSMOW2_x387 - {D17_corrections["C13D-VSMOW2_x387"]} + {D17_corrections["C13D-VSMOW2_x438"]}',
				'D17O_C13D_VSMOW2_x520': f'D17O_C13D_VSMOW2_x387 - {D17_corrections["C13D-VSMOW2_x387"]} + {D17_corrections["C13D-VSMOW2_x520"]}',
			},
		)
	# generate outlier plot:
	isofunctions.d626b_plot(stdz_triple_oxygen, p_grubbs)
	isofunctions.plot_kde(stdz_triple_oxygen, filename = 'output/kde')
	isofunctions.save_sample_results(None, None, stdz_triple_oxygen)

	stdz_triple_oxygen['unknowns'] = [
		_ for _ in stdz_triple_oxygen['unknowns'] if _.startswith('C13D')
	]
	for co2, group in [('eCO2', 'anchors'), ('dCO2', 'unknowns')]:
		waters(
			stdz_triple_oxygen,
			group,
			starting_CO2[co2],
			filename = f'output/water_samples_{co2}.csv',
		)

	plot_d13C_effect()
