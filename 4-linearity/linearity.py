import sys

sys.path.append('../lib')
from mylib import *

RAWDATA_FILENAME = 'input/linearity-rawdata.csv'

logger.remove(0)
logger.add('output/linearity.log', mode = 'w', format = '[{line:3.0f}]    {message}')

ppl.style.use('../lib/mydefault.mplstyle')


def H2O_eq_CO2_predictions(starting_CO2, eq_waters, N_starting_CO2, N_eq_waters):
	"""PREDICT COMPOSITIONS OF EQUILIBRATED CO2 SAMPLES"""

	predictions = {}
	for w in eq_waters:
		d18, D17 = isofunctions.D17_CO2_H2O_eq(
			N_eq_waters / N_starting_CO2,
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


def standardize_triple_oxygen(data, anchors):
	stdz_triple_oxygen = stdz.standardize_D17O(data, anchors, 'd628', 'd627')

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

	isofunctions.save_sample_results(None, None, stdz_triple_oxygen)

	return stdz_triple_oxygen


def H2O_eq_CO2_D17_corrections(starting_CO2, eq_waters, N_starting_CO2, N_eq_waters):
	"""COMPUTE CORRECTIONS DUE TO FINITE H2O/CO2 RATIO"""
	D17_corrections = {
		w: isofunctions.D17_CO2_H2O_eq(
			N_eq_waters / N_starting_CO2,
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

	return D17_corrections


def plot_waters(eq_waters, filename = 'output/waters', metadata = {}):
	fig = ppl.figure(figsize = (3.15, 3.5))
	fig.subplots_adjust(0.23, 0.13, 0.97, 0.78)
	ax = ppl.subplot(111)

	X = [eq_waters[w]['d18O_VSMOW'] for w in ['OC4', 'NEEM', 'HAWAI']]
	Y = [eq_waters[w]['D17O_VSMOW'] for w in ['OC4', 'NEEM', 'HAWAI']]
	ax.plot(X, Y, 'ws', mec = 'k', mew = 1, ms = 6, label = 'mixing end members')

	label = 'mixing lines'
	for a, b in [
		('OC4', 'HAWAI'),
		('NEEM', 'HAWAI'),
		('MIX_OH', 'MIX_NH'),
	]:
		X, Y = zip(
			*[
				isofunctions.mix_waters(
					eq_waters[a]['d18O_VSMOW'],
					eq_waters[b]['d18O_VSMOW'],
					eq_waters[a]['D17O_VSMOW'],
					eq_waters[b]['D17O_VSMOW'],
					x,
				)[:2]
				for x in np.linspace(0, 1)
			]
		)
		ax.plot(
			X, Y, 'k-', dashes = (6, 2, 2, 2), lw = 1, alpha = 0.5, label = label, zorder = -100
		)
		label = None

	X = [eq_waters[w]['d18O_VSMOW'] for w in ['MIX_OH', 'MIX_ONH', 'MIX_NH']]
	Y = [eq_waters[w]['D17O_VSMOW'] for w in ['MIX_OH', 'MIX_ONH', 'MIX_NH']]
	ax.plot(X, Y, 'wo', mec = 'k', mew = 1, ms = 6, label = 'mixing products')

	for w in eq_waters:
		_str = {
			'NEEM': f'{w}   ',
			'HAWAI': f'{w}  \n',
			'OC4': f'{w}\n',
			'MIX_OH': f'\n\n\n{w}',
			'MIX_NH': f'      {w}\n',
			'MIX_ONH': f'{w} \n\n',
		}[w]
		va = {
			'NEEM': 'center',
			'HAWAI': 'bottom',
			'OC4': 'bottom',
			'MIX_OH': 'center',
			'MIX_NH': 'bottom',
			'MIX_ONH': 'center',
		}[w]
		ha = {
			'NEEM': 'right',
			'HAWAI': 'center',
			'OC4': 'center',
			'MIX_OH': 'center',
			'MIX_NH': 'center',
			'MIX_ONH': 'right',
		}[w]
		ax.text(
			eq_waters[w]['d18O_VSMOW'],
			eq_waters[w]['D17O_VSMOW'],
			_str.replace('_', '-'),
			ha = ha,
			va = va,
			linespacing = 0.7,
			size = 9,
		)

	ax.set_xlabel('Water δ$^{18}$O$_{VSMOW}$ (‰)')
	ax.set_ylabel('Water Δ$^{17}$O$_{VSMOW}$ (‰)')

	ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.03))

	ax.legend(
		loc = 'lower center',
		bbox_to_anchor = (0.5, 1.01),
		fontsize = 9,
	)

	ax.axis(
		[
			eq_waters['OC4']['d18O_VSMOW'] - 6,
			eq_waters['HAWAI']['d18O_VSMOW'] + 7,
			eq_waters['MIX_OH']['D17O_VSMOW'] - 0.025,
			eq_waters['NEEM']['D17O_VSMOW'] + 0.015,
		]
	)

	fig.savefig(filename, metadata = metadata)
	ppl.close(fig)


def triple_oxygen_prediction_plot(
	S,
	predictions,
	filename = 'output/predictions_D17O_vs_d18O',
	metadata = {},
):
	fig = ppl.figure(figsize = (3.15, 3.15))
	fig.subplots_adjust(0.23, 0.15, 0.97, 0.84)

	ax = ppl.subplot(111)

	X = [r['d18corrected'] for r in S['data'] if not r['outlier']]
	Y = [r['D17corrected'] for r in S['data'] if not r['outlier']]
	(_obs,) = ax.plot(X, Y, 'kx', ms = 5, mew = 0.75)

	kw = dict(
		ha = 'center',
		va = 'center',
		size = 9,
		alpha = 1,
		color = (0, 0.5, 1),
	)
	for s in predictions:
		ax.text(
			predictions[s]['d18O_VSMOW'],
			predictions[s]['D17O_VSMOW'] + 0.016,
			s.replace('_', '-'),
			**kw,
		)

	# 	ax.text(S['samples']['LINDE']['d18O_VSMOW'], S['samples']['LINDE']['D17O_VSMOW'] - 0.016, 'LINDE', **kw)

	kw = dict(
		ms = 20,
		mfc = 'none',
		mew = 1,
		alpha = 1,
		zorder = -200,
		mec = kw['color'],
	)

	ax.plot(
		[predictions[s]['d18O_VSMOW'] for s in predictions],
		[predictions[s]['D17O_VSMOW'] for s in predictions],
		's',
		**kw,
	)

	ax.plot(
		[predictions[s]['d18O_VSMOW'] for s in predictions],
		[predictions[s]['D17O_VSMOW'] for s in predictions],
		'w+',
		ms = 22,
		mew = 8,
		zorder = -100,
	)

	kw['ms'] = 10

	(_square,) = ax.plot([], [], 's', **kw)
	(_plus,) = ax.plot([], [], 'w+', mew = 4, ms = 11)

	ax.set_xlabel('δ$^{18}$O$_{VSMOW}$ of CO$_2$ (‰)')
	ax.set_ylabel('Δ$^{17}$O$_{VSMOW}$ of CO$_2$ (‰)')

	ax.legend(
		[(_square, _plus), _obs],
		['predictions', 'measurements'],
		loc = 'lower center',
		labelspacing = 0.5,
		borderpad = 0,
		bbox_to_anchor = (0.45, 1.01),
		frameon = False,
	)

	x1, x2, y1, y2 = ax.axis()
	x0 = (x1 + x2) / 2
	y0 = (y1 + y2) / 2
	ax.set_xlim((x0 - 36, x0 + 36))
	ax.set_ylim((y0 - 0.08, y0 + 0.09))

	ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))

	fig.savefig(filename, metadata = metadata)
	ppl.close(fig)


def D17O_prediction_plot(
	S, predictions, filename = 'output/predictions_D17O', metadata = {}
):
	residuals = []
	fig = ppl.figure(figsize = (3.15, 4.37))
	fig.subplots_adjust(0.23, 0.12, 0.97, 0.84)

	moz = fig.subplot_mosaic("""A;A;A;B""")
	ax, ax2 = moz['A'], moz['B']
	ax2.sharex(ax)

	for G, kw in [
		(
			['OC4', 'HAWAI'],
			dict(
				ls = 'None',
				marker = 's',
				ms = 5,
				mew = 1,
				mfc = 'w',
				mec = 'k',
				label = 'Anchors (unfalsifiable)',
			),
		),
		(
			['NEEM'],
			dict(
				ls = 'None',
				marker = 'D',
				ms = 4.8,
				mew = 1,
				mfc = 'w',
				mec = 'k',
				label = 'Predictions from IRMS on H$_2$O',
			),
		),
		(
			['MIX_NH', 'MIX_OH', 'MIX_ONH'],
			dict(
				ls = 'None',
				marker = 'o',
				ms = 5,
				mew = 1,
				mfc = 'w',
				mec = 'k',
				label = 'Predictions from mixing effects',
			),
		),
	]:
		X = np.array([predictions[s]['D17O_VSMOW'] for s in G])
		Y = np.array([S['samples'][s]['D17O_VSMOW'] for s in G])
		if 'Anchors' not in kw['label']:
			eY = np.array([S['samples'][s]['95CL_D17O_VSMOW'] for s in G])
			ax.errorbar(
				X,
				Y,
				eY,
				elinewidth = 1,
				capsize = kw['ms'] / 2 + 1,
				capthick = 1,
				marker = 'None',
				ls = 'None',
				ecolor = kw['mec'],
			)
			ax2.errorbar(
				X,
				(Y - X) * 1e3,
				eY * 1e3,
				elinewidth = 1,
				capsize = kw['ms'] / 2 + 1,
				capthick = 1,
				marker = 'None',
				ls = 'None',
				ecolor = kw['mec'],
			)
			residuals += [y - x for x, y in zip(X, Y)]
			ax2.plot(X, (Y - X) * 1e3, **kw)
		ax.plot(X, Y, **kw)

	residuals = np.array(residuals)

	logger.info(
		f'RMSE (measured vs predicted) = {(residuals**2).mean()**.5 * 1000:.1f} ppm.'
	)

	x1, x2, y1, y2 = ax.axis()
	xmin = min(x1, y1)
	xmax = max(x2, y2)
	ax.plot(
		[xmin, xmax],
		[xmin, xmax],
		'-',
		color = 'k',
		lw = 0.8,
		dashes = (8, 2, 2, 2),
		zorder = -100,
	)
	ax.axis([xmin, xmax, xmin, xmax])
	ax.set_ylabel('Observed Δ$^{17}$O (‰) of $CO_2$')

	ax.legend(
		labelspacing = 0.2,
		loc = 'lower center',
		bbox_to_anchor = (0.5, 1.01),
		fontsize = 9,
		handlelength = 1,
		borderpad = 0.6,
	)

	ax.xaxis.set_major_locator(ticker.MultipleLocator(0.05, offset = 0.02))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05, offset = 0.02))
	ppl.setp(ax.get_xticklabels(), visible = False)

	ax2.set_xlabel('Predicted Δ$^{17}$O (‰) of $CO_2$')
	ax2.set_ylabel('Residuals\n(ppm)')
	ax2.margins(y = 0.25)
	ax2.yaxis.set_major_locator(ticker.MultipleLocator(5))
	ax2.set_ylim(-9, 7)
	ax2.axhline(0, lw = 1, alpha = 0.25, color = 'k')

	fig.savefig(filename, metadata = metadata)
	ppl.close(fig)


if __name__ == '__main__':
	"""H2O/CO2 RATIO"""
	N_starting_CO2 = 40e-6  # mol
	V_eq_waters = 300e-6  # L
	N_eq_waters = V_eq_waters / 18.015e-3  # mol
	H2O_to_CO2_ratio = N_eq_waters / N_starting_CO2

	logger.info(f'H2O/CO2 ratio is {H2O_to_CO2_ratio:.0f}.')

	"""COMPOSITION OF EQUILIBRATED WATERS"""
	eq_waters = {
		'HAWAI': {'d18O_VSMOW': 0.54, 'D17O_VSMOW': 0},
		'NEEM': {'d18O_VSMOW': -32.87, 'D17O_VSMOW': 0.038},
		'OC4': {'d18O_VSMOW': -53.93, 'D17O_VSMOW': 0.009},
	}

	POSTULATED_SLAP_D17 = 0  # Fool around by turning this knob
	logger.info(f'Postulated Δ17O_VSMOW value of SLAP = {POSTULATED_SLAP_D17} ‰')

	if POSTULATED_SLAP_D17:  # Adjust D17O values in the case of non-zero SLAP value
		for w in eq_waters:
			eq_waters[w]['D17O_VSMOW'] += (
				log(1 + eq_waters[w]['d18O_VSMOW'] / 1e3)
				/ log(1 - 55.5e-3)
				* POSTULATED_SLAP_D17
			)

	"""MIX WATERS"""
	d18, D17, txt = isofunctions.mix_waters(
		eq_waters['OC4']['d18O_VSMOW'],
		eq_waters['HAWAI']['d18O_VSMOW'],
		eq_waters['OC4']['D17O_VSMOW'],
		eq_waters['HAWAI']['D17O_VSMOW'],
		0.5,
	)
	eq_waters['MIX_OH'] = {'d18O_VSMOW': d18, 'D17O_VSMOW': D17}

	d18, D17, txt = isofunctions.mix_waters(
		eq_waters['NEEM']['d18O_VSMOW'],
		eq_waters['HAWAI']['d18O_VSMOW'],
		eq_waters['NEEM']['D17O_VSMOW'],
		eq_waters['HAWAI']['D17O_VSMOW'],
		0.5,
	)
	eq_waters['MIX_NH'] = {'d18O_VSMOW': d18, 'D17O_VSMOW': D17}

	d18, D17, txt = isofunctions.mix_waters(
		eq_waters['MIX_OH']['d18O_VSMOW'],
		eq_waters['MIX_NH']['d18O_VSMOW'],
		eq_waters['MIX_OH']['D17O_VSMOW'],
		eq_waters['MIX_NH']['D17O_VSMOW'],
		0.5,
	)
	eq_waters['MIX_ONH'] = {'d18O_VSMOW': d18, 'D17O_VSMOW': D17}

	"""WATER MIXING PLOT"""
	plot_waters(
		eq_waters,
		metadata = {
			'Author': AUTHOR,
			'Title': f"""
				Triple-oxygen mixing plot showing the water standards analyzed to test the
				linearity of our Δ17O measurements. The δ18O and Δ17O values of each of
				the mixing end-members were indepedently determined by IRMS. Because of
				well-understood nonlinear mixing effects, the Δ17O values of these six
				waters are predicted to range from -0.09 to +0.4 ‰.
				""".replace('\t', '').replace('\n', ' ')[1:-1],
		},
	)

	"""WATER COMPOSITION CSV"""
	with open('output/mixwaters.csv', 'w') as fid:
		fid.write(
			'Water_sample,OC4_fraction,NEEM_fraction,HAWAI_fraction,d18O_VSMOW,D17O_VSMOW'
		)
		for k, (x, y, z) in [
			('HAWAI', (0, 0, 1)),
			(
				'OC4',
				(
					1,
					0,
					0,
				),
			),
			('NEEM', (0, 1, 0)),
			('MIX_NH', (0, 0.5, 0.5)),
			('MIX_OH', (0.5, 0, 0.5)),
			('MIX_ONH', (0.25, 0.5, 0.25)),
		]:
			fid.write(
				f'\n{k},{x},{y},{z},{eq_waters[k]["d18O_VSMOW"]:.2f},{eq_waters[k]["D17O_VSMOW"]:.4f}'
			)

	"""COMPOSITION CO2 USED WITH EQUILIBRATED WATERS"""
	match CO2g_H2O_FRACTIONATION:
		case 'Guo & Zhou (2019)':
			starting_CO2 = {'d18O_VSMOW': 24.71, 'D17O_VSMOW': -0.0835}
		case 'Barkan & Luz (2012)':
			starting_CO2 = {'d18O_VSMOW': 24.25, 'D17O_VSMOW': -0.150}

	"""H2O_eq_CO2 PLOT"""
	isofunctions.H2O_eq_CO2_figure(
		starting_CO2,
		eq_waters,
		N_starting_CO2,
		N_eq_waters,
		metadata = {
			'Author': AUTHOR,
			'Title': f"""
				Predicted 17O anomaly of water-equilibrated CO2 as a function of the molecular ratio
				H2O/CO2. Due to nonlinear mixing effects, final Δ17O of equilibrated CO2 depends on
				H2O/CO2 ratio, initial water Δ17O, and relative δ18O values of CO2 and H2O.
				Initial Δ17O of the CO2 tank used in this study is {starting_CO2['D17O_VSMOW']:.3f} ‰.
				Solid and dashed black lines corespond to the different waters used in this study.
				The mixing ratio used in our experiments (~{H2O_to_CO2_ratio/100:.0f}00) is shown as a vertical grey line.
				""".replace('\t', '').replace('\n', ' ')[1:-1],
		},
	)

	"""H2O_eq_CO2 CORRECTIONS"""
	D17_corrections = H2O_eq_CO2_D17_corrections(
		starting_CO2, eq_waters, N_starting_CO2, N_eq_waters
	)
	predictions = H2O_eq_CO2_predictions(
		starting_CO2, eq_waters, N_starting_CO2, N_eq_waters
	)

	"""LOAD VCOF-CRDS ANALYSES"""
	data = read_rawdata([RAWDATA_FILENAME])
	data = [r for r in data if r['Sample'] != 'LINDE']
	for k, r in enumerate(data):
		r['UID'] = f'{r["Run"]}-{k+1:02}'
		r['Session'] = '2022-10'

	"""FLAG d626b OUTLIERS"""
	# first flag them without excluding them:
	data, p_grubbs = isofunctions.d626b_outliers(data, rename = False)
	# standardize the data with these outliers still included:
	stdz_triple_oxygen = standardize_triple_oxygen(
		data, {_: predictions[_] for _ in ['HAWAI', 'OC4']}
	)
	# generate outlier plot:
	isofunctions.d626b_plot(stdz_triple_oxygen, p_grubbs)

	"""STANDARDIZATION"""
	# now exclude the outliers:
	# 	data, p_grubbs = d626b_outliers(data, rename = False)
	data = [_ for _ in data if not _['outlier']]
	# standardize the data once more, with outliers excluded this time:
	stdz_triple_oxygen = standardize_triple_oxygen(
		data, {_: predictions[_] for _ in ['HAWAI', 'OC4']}
	)
	with open('output/residuals.csv', 'w') as fid:
		fid.write('d18,d17,D17')
		for r in stdz_triple_oxygen['data']:
			fid.write(
				f"\n{r['d18residual']:.6f},{r['d17residual']:.6f},{r['D17residual']:.6f}"
			)

	isofunctions.plot_kde(stdz_triple_oxygen, filename = 'output/kde')

	"""COMPARE OBSERVATIONS TO PREDICTIONS"""
	triple_oxygen_prediction_plot(stdz_triple_oxygen, predictions)
	D17O_prediction_plot(
		stdz_triple_oxygen,
		predictions,
		metadata = {
			'Author': AUTHOR,
			'Title': f"""
				Data-prediction comparison.
				""".replace('\t', '').replace('\n', ' ')[1:-1],
		},
	)

	with open('output/results.csv', 'w') as fid:
		fid.write('Sample')
		fid.write(',initial_d18O_VSMOW_water')
		fid.write(',initial_D17O_VSMOW_water')
		fid.write(',initial_d18O_VSMOW_CO2')
		fid.write(',initial_D17O_VSMOW_CO2')
		fid.write(',H2O/CO2')
		fid.write(',final_d18O_VSMOW_CO2_predicted')
		fid.write(',final_D17O_VSMOW_CO2_predicted')
		fid.write(',final_d18O_VSMOW_CO2_measured')
		fid.write(',final_D17O_VSMOW_CO2_measured')
		fid.write(',final_D17O_VSMOW_CO2_measured_95CL')
		fid.write(',D17O_residual')
		for w in eq_waters:
			fid.write(f"\n{w.replace('_', '-')}")
			fid.write(f",{eq_waters[w]['d18O_VSMOW']:.2f}")
			fid.write(f",{eq_waters[w]['D17O_VSMOW']:.4f}")
			fid.write(f",{starting_CO2['d18O_VSMOW']:.2f}")
			fid.write(f",{starting_CO2['D17O_VSMOW']:.2f}")
			fid.write(f',{H2O_to_CO2_ratio:.0f}')
			fid.write(f",{predictions[w]['d18O_VSMOW']:.2f}")
			fid.write(f",{predictions[w]['D17O_VSMOW']:.4f}")
			fid.write(f",{stdz_triple_oxygen['samples'][w]['d18O_VSMOW']:.2f}")
			fid.write(f",{stdz_triple_oxygen['samples'][w]['D17O_VSMOW']:.4f}")
			fid.write(
				',-'
				if w in ['HAWAI', 'OC4']
				else f",{stdz_triple_oxygen['samples'][w]['95CL_D17O_VSMOW']:.4f}"
			)
			fid.write(
				f",{stdz_triple_oxygen['samples'][w]['D17O_VSMOW']-predictions[w]['D17O_VSMOW']:.4f}"
			)
