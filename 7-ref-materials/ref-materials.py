import sys

sys.path.append('../lib')

from mylib import *

RAWDATA_FILENAME = 'input/ref-materials-rawdata.csv'

logger.remove(0)
logger.add('output/ref-materials.log', mode = 'w', format = '[{line:3.0f}]    {message}')

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

	return stdz_triple_oxygen


def standardize_carb_d18O(data, anchors, **kwargs):
	S = stdz.standardize(data, anchors, 'd628', 'd18O_VPDB', **kwargs)

	isofunctions.plot_residuals(
		S['data'],
		key = 'd18O_VPDB',
		filename = f'output/residuals_carb_d18O',
		anchors = [_ for _ in S['anchors']],
	)
	with open(f'output/sessions_carb_d18O.csv', 'w') as fid:
		fid.write(S['csv_sessions'])

	return S


def standardize_d13C(data, anchors, **kwargs):
	S = stdz.standardize(data, anchors, 'd636', 'd13C_VPDB', **kwargs)

	isofunctions.plot_residuals(
		S['data'],
		key = 'd13C_VPDB',
		filename = f'output/residuals_d13C',
		anchors = [_ for _ in S['anchors']],
	)
	with open(f'output/sessions_d13C.csv', 'w') as fid:
		fid.write(S['csv_sessions'])

	return S


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
			f'Δ’17O correction is {D17_corrections[w]*1e3:+.1f} ppm for CO2 equilibrated with {w}.'
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


def triple_oxygen_prediction_plot(
	S, predictions, filename = 'output/predictions_D17O_vs_d18O'
):
	fig = figure(figsize = (3.15, 3.15))
	subplots_adjust(0.23, 0.15, 0.97, 0.84)

	ax = subplot(111)

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
	ax.set_ylabel('Δ’$^{17}$O$_{VSMOW}$ of CO$_2$ (‰)')

	ax.legend(
		[(_square, _plus), _obs],
		['predictions', 'measurements'],
		loc = 'lower center',
		labelspacing = 0.5,
		borderpad = 0,
		bbox_to_anchor = (0.45, 1.01),
		frameon = False,
	)

	ax.set_xlim((-23, 49))
	ax.set_ylim((-0.31, -0.14))

	ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))

	fig.savefig(filename)
	close(fig)


def D17O_prediction_plot(S, predictions, filename = 'output/predictions_D17O'):
	residuals = []
	fig = figure(figsize = (3.15, 4.37))
	fig.subplots_adjust(0.23, 0.12, 0.97, 0.84)

	moz = fig.subplot_mosaic("""A;A;A;B""")
	ax, ax2 = moz['A'], moz['B']
	ax2.sharex(ax)

	for G, kw in [
		(
			['OC4', 'HAWAI'],
			dict(
				ls = 'None',
				marker = 'o',
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
				marker = 's',
				ms = 5,
				mew = 1,
				mfc = 'w',
				mec = 'k',
				label = 'Predictions from mixing effects',
			),
		),
	]:
		X = array([predictions[s]['D17O_VSMOW'] for s in G])
		Y = array([S['samples'][s]['D17O_VSMOW'] for s in G])
		if 'Anchors' not in kw['label']:
			eY = array([S['samples'][s]['95CL_D17O_VSMOW'] for s in G])
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

	residuals = array(residuals)

	logger.info(
		f'RMSE (measured vs predicted) = {(residuals**2).sum()**.5 * 1000:.1f} ppm.'
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
	ax.set_ylabel(f'Observed Δ’$^{{17}}$O (‰) of $CO_2$')

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
	setp(ax.get_xticklabels(), visible = False)

	ax2.set_xlabel('Predicted Δ’$^{17}$O (‰) of $CO_2$')
	ax2.set_ylabel('Residuals\n(ppm)')
	ax2.margins(y = 0.25)
	ax2.yaxis.set_major_locator(ticker.MultipleLocator(5))
	ax2.set_ylim(-9, 7)
	ax2.axhline(0, lw = 1, alpha = 0.25, color = 'k')

	savefig(filename)
	close(fig)


def compare_standardization_methods(data, anchors):
	methods = ['correct_observations', 'predict_observations']
	residuals = ['D17', 'd17']
	samples = sorted({r['Sample'] for r in data if r['Sample'] not in anchors})

	out = {}

	for method in methods:
		out[method] = {}
		for residuals_17 in residuals:
			S = stdz_triple_oxygen = stdz.standardize_D17O(
				data,
				anchors,
				'd628',
				'd627',
				method = method,
				residuals_17=residuals_17,
			)

			out[method][residuals_17] = {
				s: f"{S['samples'][s]['D17O_VSMOW']:.4f}" for s in samples
			}

	with open('output/compare_stdz.csv', 'w') as fid:
		fid.write('method,residuals,')
		fid.write(','.join(samples))
		for method in methods:
			for residuals_17 in residuals:
				fid.write(f'\n{method},{residuals_17}')
				for s in samples:
					fid.write(f',{out[method][residuals_17][s]}')


def revisit_refmaterials(
	refmaterials, stdz_triple_oxygen, filename = 'output/refmaterials.csv'
):
	with open(filename, 'w') as fid:
		fid.write(
			'Reference Material,Nominal d13C_VPDB,Nominal d18O_VPDB,R45,R46,D17O_VSMOW,Recalculated d13C_VPDB,d13C Shift'
		)
		for s in ['NBS18', 'NBS19', 'IAEA603', 'IAEA610', 'IAEA611', 'IAEA612']:
			refmaterials[s]['D17O_VSMOW'] = round(
				stdz_triple_oxygen['samples'][s]['D17O_VSMOW'], 4
			)
			refmaterials[s]['R45'], refmaterials[s]['R46'] = co2irms.deltas_to_ratios(
				refmaterials[s]['d13C_VPDB'],
				(1000 + refmaterials[s]['d18O_VPDB']) * 1.03092 * 1.01025 - 1000,
			)
			(refmaterials[s]['d13C_VPDB_new'], refmaterials[s]['d18O_VSMOW_new']) = (
				co2irms.ratios_to_deltas(
					refmaterials[s]['R45'],
					refmaterials[s]['R46'],
					D17O = refmaterials[s]['D17O_VSMOW'],
				)
			)
			refmaterials[s]['d18O_VPDB_new'] = (
				1000 + refmaterials[s]['d18O_VSMOW_new']
			) / 1.03092 / 1.01025 - 1000
			fid.write(f'\n{s}')
			fid.write(f",{refmaterials[s]['d13C_VPDB']},{refmaterials[s]['d18O_VPDB']}")
			fid.write(f",{refmaterials[s]['R45']:.9f},{refmaterials[s]['R46']:.12f}")
			fid.write(f",{refmaterials[s]['D17O_VSMOW']:.4f}")
			fid.write(f",{refmaterials[s]['d13C_VPDB_new']:.3f}")
			fid.write(
				f",{refmaterials[s]['d13C_VPDB_new'] - refmaterials[s]['d13C_VPDB']:.3f}"
			)

	return refmaterials


def check_nominal_values(stdz_carbon, stdz_oxygen_carb, refmaterials):
	fig = ppl.figure(figsize = (3.15, 5))
	fig.subplots_adjust(0.25, 0.1, 0.95, 0.97, 0.3, 0.3)

	ax = ppl.subplot(211)

	sigma = (
		np.array([r[f'd13C_VPDB_residual'] for r in stdz_carbon['data']]) ** 2
	).mean() ** 0.5
	X = np.array([refmaterials[s]['d13C_VPDB'] for s in refmaterials])
	Y = np.array(
		[
			np.mean(
				[r['d13C_VPDB_corrected'] for r in stdz_carbon['samples'][s]['data']]
			)
			for s in refmaterials
		]
	)
	eY = np.array([sigma / stdz_carbon['samples'][s]['N'] ** 0.5 for s in refmaterials])
	ax.errorbar(
		X,
		Y - X,
		eY,
		ls = 'None',
		marker = 'None',
		ecolor = 'k',
		elinewidth = 1,
		capsize = 2,
		capthick = 1,
	)
	ax.plot(X, Y - X, 'ws', mec = 'k', mew = 1)
	ax.set_xlabel('Nominal δ$^{13}C_{VPDB}$ (‰)')
	ax.set_ylabel('Residual δ$^{13}C_{VPDB}$ (‰)')

	ax = ppl.subplot(212)

	sigma = (
		np.array([r[f'd18O_VPDB_residual'] for r in stdz_oxygen_carb['data']]) ** 2
	).mean() ** 0.5
	X = np.array([refmaterials[s]['d18O_VPDB'] for s in refmaterials])
	Y = np.array(
		[
			np.mean(
				[
					r['d18O_VPDB_corrected']
					for r in stdz_oxygen_carb['samples'][s]['data']
				]
			)
			for s in refmaterials
		]
	)
	eY = np.array(
		[sigma / stdz_oxygen_carb['samples'][s]['N'] ** 0.5 for s in refmaterials]
	)
	ax.errorbar(
		X,
		Y - X,
		eY,
		ls = 'None',
		marker = 'None',
		ecolor = 'k',
		elinewidth = 1,
		capsize = 2,
		capthick = 1,
	)
	ax.plot(X, Y - X, 'ws', mec = 'k', mew = 1)
	ax.set_xlabel('Nominal δ$^{18}O_{VPDB}$ (‰)')
	ax.set_ylabel('Residual δ$^{18}O_{VPDB}$ (‰)')

	fig.savefig('output/check_nominal_values')


def gresp(data, starting_CO2, filename = 'output/water_samples.csv'):
	xGRESP = [r for r in data if r['Sample'] == 'GRESP'][0]['H2O_to_CO2']
	co2_samples = [
		{
			'Sample': s,
			'd18O_VSMOW': stdz_triple_oxygen['samples'][s]['d18O_VSMOW'],
			'D17O_VSMOW': stdz_triple_oxygen['samples'][s]['D17O_VSMOW'],
			'H2O/CO2': xGRESP if s == 'GRESP' else float(s.split('_x')[-1]),
		}
		for s in ['GRESP'] + [_ for _ in stdz_triple_oxygen['anchors']]
	]

	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		isofunctions.correct_CO2eqH2O_back_to_initial_H2O(
			co2_samples, starting_CO2['d18O_VSMOW'], starting_CO2['D17O_VSMOW']
		)

	GRESP = co2_samples[0]
	logger.info(
		f'GRESP yields\n\t\td18O_VSMOW = {GRESP["water_d18O_VSMOW"]:.3f} ‰\n\t\tD17O_VSMOW = {GRESP["water_D17O_VSMOW"]:.4f} ‰'
	)

	with open(filename, 'w') as fid:
		fid.write('Sample,water_d18O_VSMOW,water_D17O_VSMOW')
		for s in co2_samples:
			fid.write(
				f'\n{s["Sample"]},{s["water_d18O_VSMOW"]:.3f},{s["water_D17O_VSMOW"]:.4f}'
			)


def triple_oxygen_plot(X, filename = 'output/triple_oxygen_plot'):
	fig = ppl.figure(figsize = (3.15, 3.7))
	ax = ppl.subplot(111)
	fig.subplots_adjust(0.22, 0.13, 0.98, 0.78)

	kw = dict(
		ls = 'None',
		capthick = 1,
		capsize = 2,
		elinewidth = 1,
		mew = 1,
		mfc = 'w',
		ms = 6,
		ecolor = 'k',
		mec = 'k',
	)

	ax.errorbar(
		[],
		[],
		marker = (4, 0, 45),
		label = '$CO_2$ from $CaCO_3$ (90$\\,$°C acid reaction)',
		**kw,
	)
	ax.errorbar(
		[],
		[],
		marker = (4, 0, 0),
		label = '$CO_2$ equilibrated with $H_2O$ at 25$\\,$°C',
		**kw,
	)

	done = []
	samples = X['samples']
	for s in samples:
		if '95CL_D17O_VSMOW' in samples[s]:
			x = samples[s]['d18O_VSMOW']
			y = samples[s]['D17O_VSMOW']
			if s in D17_corrections:
				y -= D17_corrections[s]
			ey = samples[s]['95CL_D17O_VSMOW']
			ax.errorbar(
				x, y, ey, marker = (4, 0, 0) if s == 'GRESP' else (4, 0, 45), **kw
			)
			if s == 'GRESP':
				ax.text(x - 2.5, y - 0.0003, s, va = 'center', ha = 'right', size = 8)
				logger.info(f'GRESP-CO2: D17O_VSMOW = {y:.4f} ‰')
				logger.info(
					f'GRESP-H2O: D17O_VSMOW = {y - 1000*(np.log(alpha17_CO2g_H2O) - LAMBDA_17*np.log(alpha18_CO2g_H2O)):.4f} ‰'
				)
			elif x in [
				samples[_]['d18O_VSMOW']
				for _ in ['NBS19', 'NBS18', 'IAEA612', 'IAEA611']
			]:
				ax.text(
					x + 2,
					y - 0.0003,
					chop_name(s),
					va = 'center',
					ha = 'left',
					size = 8,
					linespacing = 1,
				)
			else:
				ax.text(
					x - 2.5,
					y - 0.0003,
					chop_name(s),
					va = 'center',
					ha = 'right',
					size = 8,
					linespacing = 1,
				)
		else:
			shortname = s.split('_x')[0]
			if shortname not in done:
				done.append(shortname)
				x = samples[s]['d18O_VSMOW']
				y = samples[s]['D17O_VSMOW']
				x -= d18_corrections[s]
				y -= D17_corrections[s]

				ax.plot(x, y, 'w', marker = (4, 0, 0), mew = 1, mec = 'k', ms = 5)
				ax.text(
					x + 1 - np.sign(x),
					y - 0.0045,
					shortname,
					va = 'top',
					ha = 'center',
					size = 8,
				)

	ax.legend(
		fontsize = 8,
		loc = 'lower center',
		bbox_to_anchor = (0.5, 1.01),
		labelspacing = 1.0,
		handlelength = 1,
		borderpad = 1.4,
	)
	ax.set_xlabel('δ$^{18}O_{VSMOW}$ (‰)')
	ax.set_ylabel('Δ’$^{17}O_{VSMOW}$ (‰)')
	x1, x2, y1, y2 = ax.axis()
	x0 = (x1 + x2) / 2
	y0 = (y1 + y2) / 2
	ax.set_xlim((x0 - 33, x0 + 37))
	ax.set_ylim((y0 - 0.048, y0 + 0.043))
	fig.savefig(filename)
	ppl.close(fig)


def chop_name(x):
	for k, l in enumerate(x):
		if l in '0123456789':
			break
	return f'{x[:k]}\n{x[k:]}'


def plot_IAEA603_NBS18_residuals(S):
	XY1 = [
		(k + 1, r['D17residual'] * 1000)
		for k, r in enumerate(stdz_triple_oxygen['data'])
		if r['Sample'] == 'IAEA603'
	]
	XY2 = [
		(k + 1, r['D17residual'] * 1000)
		for k, r in enumerate(stdz_triple_oxygen['data'])
		if r['Sample'] == 'NBS18'
	]
	XY3 = [
		(k + 1, r['D17corrected'])
		for k, r in enumerate(stdz_triple_oxygen['data'])
		if r['Sample'] == 'IAEA603'
	]
	XY4 = [
		(k + 1, r['D17corrected'])
		for k, r in enumerate(stdz_triple_oxygen['data'])
		if r['Sample'] == 'NBS18'
	]

	fig = ppl.figure(figsize = (3.15, 3.8))
	fig.subplots_adjust(0.23, 0.08, 0.98, 0.93, 0, -0.25)

	ax = ppl.subplot(414)
	ax.plot(*zip(*XY1), 'o', mfc = (0.5, 0.5, 0.5), mec = 'k', mew = 0.8, label = 'IAEA603')
	ax.plot(*zip(*XY2), 'o', mfc = 'w', mec = 'k', mew = 0.8, label = 'NBS18')
	ax.set_xticks([])
	ax.set_yticks([-20, 0, 20])
	ax.axis([None, None, -24, 24])
	ax.set_xlabel('sequence of analyses')
	ax.set_ylabel('Δ’$^{17}O$ residuals (ppm)')

	ax2 = ppl.subplot(211, sharex = ax)
	ax2.plot(*zip(*XY3), 'o', mfc = (0.5, 0.5, 0.5), mec = 'k', mew = 0.8, label = 'IAEA603')
	ax2.plot(*zip(*XY4), 'o', mfc = 'w', mec = 'k', mew = 0.8, label = 'NBS18')
	ax2.set_xticks([])
	ax2.set_xlabel('sequence of analyses')
	ax2.set_ylabel('Δ’$^{17}O_{VSMOW}$ (‰)')
	ax2.legend(
		fontsize = 8,
		bbox_to_anchor = (0.5, 1.01),
		loc = 'lower center',
		ncols = 2,
		frameon = False,
		borderpad = 0,
		handlelength = 0.5,
	)

	fig.savefig('output/residuals')
	ppl.close(fig)


def plot_seasonal_data(sdata, data):
	# 	offset = data['samples']['NBS18']['D17O_VSMOW'] - data['samples']['IAEA603']['D17O_VSMOW']

	seasons = sorted(
		{s.split('_')[-1] for s in sdata['samples'] if s.startswith('IAEA603')}
	)

	summary = {
		s: {
			'X': sdata['samples'][f'IAEA603_{s}']['D17O_VSMOW'],
			'Y': sdata['samples'][f'NBS18_{s}']['D17O_VSMOW'],
			'sX': sdata['samples'][f'IAEA603_{s}']['SE_D17O_VSMOW'],
			'sY': sdata['samples'][f'NBS18_{s}']['SE_D17O_VSMOW'],
			'rXY': sdata['bestfit']
			.params[f'D17O_IAEA603_{s}']
			.correl[f'D17O_NBS18_{s}'],
		}
		for s in seasons
	}

	fig = ppl.figure(figsize = (3.15, 2.95))
	fig.subplots_adjust(0.22, 0.16, 0.97, 0.97)
	ax = ppl.subplot(111, aspect = 1)
	ax.xaxis.set_major_locator(ticker.MultipleLocator(0.01))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.01))

	for s in seasons:
		covXY = summary[s]['sX'] * summary[s]['sY'] * summary[s]['rXY']
		w, h, r = isofunctions.cov_ellipse(
			np.array([[summary[s]['sX'] ** 2, covXY], [covXY, summary[s]['sY'] ** 2]])
		)
		ax.add_patch(
			Ellipse(
				xy = (summary[s]['X'], summary[s]['Y']),
				width = w,
				height = h,
				angle = r,
				ec = 'k',
				fc = 'None',
			)
		)
		if 'early' in s:
			ax.annotate(
				s.replace('2023', '\n2023'),
				xy = (summary[s]['X'] - 0.0025, summary[s]['Y'] - summary[s]['sY'] * 2.4),
				xytext = (
					summary[s]['X'] - 0.0025,
					summary[s]['Y'] - summary[s]['sY'] * 3.6,
				),
				va = 'top',
				ha = 'center',
				size = 9,
				arrowprops = dict(
					arrowstyle = '->',
					lw = 0.8,
					ec = 'k',
				),
				color = 'k',
			)
		elif 'late' in s:
			ax.annotate(
				s.replace('2023', '\n2023'),
				xy = (summary[s]['X'] - summary[s]['sX'] * 2.4, summary[s]['Y'] - 0.002),
				xytext = (
					summary[s]['X'] - summary[s]['sX'] * 5.9,
					summary[s]['Y'] - 0.002,
				),
				va = 'center',
				ha = 'center',
				size = 9,
				arrowprops = dict(
					arrowstyle = '->',
					lw = 0.8,
					ec = 'k',
				),
				color = 'k',
			)

	x0 = -0.130
	y0 = -0.105
	xmin = x0 - 0.019
	xmax = x0 + 0.019
	ymin = y0 - 0.019
	ymax = y0 + 0.019
	ax.plot([xmin, xmax], [xmin, xmax], 'k-', lw = 0.75, dashes = (3, 2), alpha = 0.4)
	ax.text(
		x0 + 0.012,
		x0 + 0.012,
		'1:1\n',
		size = 9,
		alpha = 0.5,
		va = 'center',
		ha = 'center',
		linespacing = 2,
		rotation = 45,
	)

	for o in [0.025, 0.050]:
		ax.plot(
			[xmin, xmax],
			[xmin + o, xmax + o],
			'k-',
			lw = 0.75,
			dashes = (3, 2),
			alpha = 1 if o == 0.025 else 0.4,
		)
		ax.text(
			x0 - 0.012,
			x0 - 0.012 + o,
			f'{o:+.3f} ‰\n',
			size = 9,
			rotation = 45,
			va = 'center',
			ha = 'center',
			linespacing = 2,
			alpha = 1 if o == 0.025 else 0.5,
		)

	ax.axis([xmin, xmax, ymin, ymax])
	ax.set_xlabel('$Δ’^{17}O_{VSMOW}$ of $CO_2$ from IAEA603 (‰)')
	ax.set_ylabel('$Δ’^{17}O_{VSMOW}$ of $CO_2$ from NBS18 (‰)')

	fig.savefig('output/compare_seasons')


if __name__ == '__main__':
	"""COMPOSITION OF EQUILIBRATED WATERS"""
	ref_waters = {
		'VSMOW2': {'d18O_VSMOW': 0.0, 'D17O_VSMOW': 0.000},
		'SLAP2': {'d18O_VSMOW': -55.5, 'D17O_VSMOW': 0.000},
		# 		'GRESP' : {'d18O_VSMOW': -34.0, 'D17O_VSMOW': 0.041},
	}

	"""LOAD VCOF-CRDS ANALYSES"""
	data = read_rawdata([RAWDATA_FILENAME])
	for r in data:
		r['Sample'] = r['Sample'].replace('IAEA-', 'IAEA').replace('NBS-', 'NBS')
	data = [r for r in data if r['Sample'] != 'LINDE']
	data = [r for r in data if r['UID'] != 'DFC25593']
	# 	data = [r for r in data if r['Session'] not in ['2023-06', '2023-08']]
	# 	data = [r for r in data if r['Session'] not in ['2023-10']]
	# 	data = [r for r in data if r['Session'] not in ['2023-12']]

	"""DEFINE WATER-EQUILIBRATED STANDARDS"""
	unique_combinations = {
		(r['Sample'], r['H2O_to_CO2'])
		for r in data
		if r['Sample'] in ['VSMOW2', 'SLAP2']
	}
	eq_waters = {
		f'{w}_x{x:.0f}': {
			'd18O_VSMOW': ref_waters[w]['d18O_VSMOW'],
			'D17O_VSMOW': ref_waters[w]['D17O_VSMOW'],
			'x': x,
		}
		for w, x in unique_combinations
	}
	for r in data:
		if r['Sample'] in ['VSMOW2', 'SLAP2']:
			r['Sample'] = r['Sample'] + f'_x{r["H2O_to_CO2"]:.0f}'

	POSTULATED_SLAP_D17 = 0  # Fool around by turning this knob
	logger.info(f'Postulated Δ’17O_VSMOW value of SLAP = {POSTULATED_SLAP_D17} ‰')

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
			starting_CO2 = {'d18O_VSMOW': 24.71, 'D17O_VSMOW': -0.0835}
		case 'Barkan & Luz (2012)':
			starting_CO2 = {'d18O_VSMOW': 24.25, 'D17O_VSMOW': -0.150}

	"""H2O_eq_CO2 PLOT"""
	isofunctions.H2O_eq_CO2_figure(
		starting_CO2, ref_waters, 40e-6, 40e-6 * 312, H2O_to_CO2_ratio_half_width = 1 / 3
	)

	"""H2O_eq_CO2 CORRECTIONS"""
	d18_corrections, D17_corrections = H2O_eq_CO2_D17_corrections(
		starting_CO2, eq_waters
	)
	predictions = H2O_eq_CO2_predictions(starting_CO2, eq_waters)

	"""FLAG d626b OUTLIERS"""
	# first flag them without excluding them:
	data, p_grubbs = isofunctions.d626b_outliers(data, rename = False)
	# standardize the data with these outliers still included:
	stdz_triple_oxygen = standardize_triple_oxygen(data, predictions)
	# generate outlier plot:
	isofunctions.d626b_plot(stdz_triple_oxygen, p_grubbs)

	"""STANDARDIZATION"""
	# now exclude the outliers:
	# 	data, p_grubbs = d626b_outliers(data, rename = False)
	data = [_ for _ in data if not _['outlier']]
	# standardize the data once more, with outliers excluded this time:

	for r in data:
		if r['Sample'] in ['IAEA603', 'NBS18']:
			if r['Session'].startswith('2023-03'):
				r['Sample'] = r['Sample'] + '_early2023'
			# 			elif r['Session'].startswith('2023-0'):
			# 				r['Sample'] = r['Sample'] + '_summer'
			# 			elif r['Session'].startswith('2023-1'):
			# 				r['Sample'] = r['Sample'] + '_fall'
			else:
				r['Sample'] = r['Sample'] + '_late2023'

	seasonal_stdz_triple_oxygen = standardize_triple_oxygen(data, predictions)
	isofunctions.save_sample_results(None, None, seasonal_stdz_triple_oxygen)

	for s in ['IAEA603', 'NBS18']:
		for r in data:
			if r['Sample'].startswith(s):
				r['Sample'] = s

	stdz_triple_oxygen = standardize_triple_oxygen(data, predictions)

	plot_seasonal_data(seasonal_stdz_triple_oxygen, stdz_triple_oxygen)

	isofunctions.plot_kde(stdz_triple_oxygen, filename = 'output/kde')
	triple_oxygen_plot(stdz_triple_oxygen, filename = 'output/triple_oxygen_plot')
	plot_IAEA603_NBS18_residuals(stdz_triple_oxygen)

	with open('output/analyses.csv', 'w') as fid:
		fid.write(
			'UID,Sample,Session,P,d636,d628,d627,d18O_VSMOW,D17O_VSMOW,d18O_residual,D17O_residual'
		)
		for r in stdz_triple_oxygen['data']:
			fid.write(
				(
					f"\n{r['UID']}"
					+ f",{r['Sample']}"
					+ f",{r['Session']}"
					+ f",{r['P']:.4f}"
					+ f",{r['d636']:.4f}"
					+ f",{r['d628']:.4f}"
					+ f",{r['d627']:.4f}"
					+ f",{r['d18corrected']:.4f}"
					+ f",{r['D17corrected']:.4f}"
					+ f",{r['d18residual']:.4f}"
					+ f",{r['D17residual']:.4f}"
				)
			)

	"""CHECK d18O VALUE OF NBS18"""
	ratio_IAEA603_to_NBS18 = (
		1000 + stdz_triple_oxygen['samples']['IAEA603']['d18O_VSMOW']
	) / (1000 + stdz_triple_oxygen['samples']['NBS18']['d18O_VSMOW'])
	predicted_d18O_VPDB_of_NBS18 = (1000 - 2.37) / ratio_IAEA603_to_NBS18 - 1000
	error_on_predicted_d18O_VPDB_of_NBS18 = (
		stdz_triple_oxygen['samples']['IAEA603']['SE_d18O_VSMOW'] ** 2
		+ stdz_triple_oxygen['samples']['NBS18']['SE_d18O_VSMOW'] ** 2
	) ** 0.5
	logger.info(
		f'Predicted d18O_VPDB of NBS18 is {predicted_d18O_VPDB_of_NBS18:.2f} ± {error_on_predicted_d18O_VPDB_of_NBS18:.2f} ‰ (1SE)'
	)

	"""COMPARE STANDARDIZATION METHODS"""
	compare_standardization_methods(data, predictions)

	with warnings.catch_warnings():
		warnings.simplefilter('ignore')

		"""d18O STANDARDIZATION"""
		refmaterials_d18O = {
			_: isofunctions.REFMATERIALS[_] for _ in ['NBS18', 'NBS19', 'IAEA603']
		}
		stdz_oxygen_carb = standardize_carb_d18O(data, refmaterials_d18O)

		"""REVISIT 13C REF MATERIALS"""
		refmaterials = {
			_: isofunctions.REFMATERIALS[_]
			for _ in ['NBS18', 'NBS19', 'IAEA603', 'IAEA610', 'IAEA611', 'IAEA612']
		}
		refmaterials = revisit_refmaterials(refmaterials, stdz_triple_oxygen)

		"""d13C STANDARDIZATION"""
		refmaterials_d13C = {
			_: refmaterials[_]
			for _ in ['NBS19', 'IAEA603', 'IAEA610', 'IAEA611', 'IAEA612']
		}
		data_d13C = [
			r for r in data if '2023-03' in r['Session']
		]  # only sessions with enough constraints
		stdz_carbon = standardize_d13C(data_d13C, refmaterials_d13C)

	isofunctions.plot_kde(stdz_triple_oxygen, filename = 'output/kde')

	isofunctions.save_sample_results(stdz_carbon, stdz_oxygen_carb, stdz_triple_oxygen)

	"""COMPARISON WITH NOMINAL VALUES"""
	check_nominal_values(stdz_carbon, stdz_oxygen_carb, refmaterials)

	"""GRESP"""
	gresp(data, starting_CO2)
