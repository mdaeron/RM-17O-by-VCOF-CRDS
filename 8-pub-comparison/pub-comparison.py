import sys

sys.path.append('../lib')
from mylib import *

logger.remove(0)
logger.add('output/ref-materials.log', mode = 'w', format = '[{line:3.0f}]    {message}')

ppl.style.use('../lib/mydefault.mplstyle')


def chop_name(x):
	for k, l in enumerate(x):
		if l in '0123456789':
			break
	return f'{x[:k]}\n{x[k:]}'


colorz = {
	'red': (0.920, 0.000, 0.000),
	'orange': (1.000, 0.600, 0.000),
	'yellow': (1.000, 0.800, 0.200),
	'green': (0.013, 0.660, 0.000),
	'blue': (0.000, 0.384, 0.960),
	'purple': (0.538, 0.320, 1.000),
	'pink': (0.990, 0.500, 1.000),
	'grey': (0.670, 0.670, 0.670),
	'dark': (0.330, 0.330, 0.330),
	'black': (0.000, 0.000, 0.000),
}


def standardize_triple_oxygen(data, anchors):
	stdz_triple_oxygen = stdz.standardize_D17O(
		data,
		anchors,
		'd628',
		'd627',
		d18O_key_out = 'd18O_VPDB',
		d17O_key_out = 'd17O_VPDB',
		method = 'correct_observations',  # 'predict_observations'
		residuals_17='d17',  # 'D17'
	)

	for S, key, fn, ss in [
		(stdz_triple_oxygen, 'd18O_VPDB', 'residuals_d18O', ''),
		(stdz_triple_oxygen, 'D17O_VPDB', 'residuals_D17O', 'sessions_triple_oxygen'),
	]:
		plot_residuals(
			S['data'],
			key = key,
			filename = f'output/{fn}',
			anchors = [_ for _ in S['anchors']],
		)
		if ss:
			with open(f'output/{ss}.csv', 'w') as fid:
				fid.write(S['csv_sessions'])

	return stdz_triple_oxygen


def pub_comparison_plot(data, filename = 'output/pub_comparison'):
	fig = ppl.figure(figsize = (6.62, 3))
	fig.subplots_adjust(0.11, 0.05, 0.98, 0.95, 1.65)

	colors = {
		'Passey et al. (2014)': colorz['red'],
		'Barkan et al. (2019)': colorz['orange'],
		'Passey & Ji (2019)': colorz['yellow'],
		'Fosu et al. (2020)': colorz['green'],
		'Sha et al. (2020)': colorz['blue'],
		'Wostbrock et al. (2020)': colorz['purple'],
		# 		"Perdue et al. (2022)*": colorz["grey"],
		# 		"Hare et al. (2022)*": colorz["dark"],
		'Ellis & Passey (2023)': colorz['pink'],
		'this study': colorz['black'],
	}

	markers = {
		None: (4, 0, 45),
		'90': (3, 0, 0),
		'70': (4, 0, 0),
		'25': (3, 0, 180),
	}

	refs = {_: [r for r in data if r['Ref'] == _] for _ in colors}

	X = {
		'NBS19': -1,
		'NBS18': 0,
		'IAEA603': 1,
	}

	kwm = dict(ls = 'None', ms = 6, mew = 1, mfc = 'w', zorder = 200)
	kwe = dict(
		ls = 'None',
		marker = 'None',
		elinewidth = 1,
		zorder = 100,
		capsize = 3,
		capthick = 1,
		alpha = 0.5,
	)
	kwl = dict(ls = '-', marker = 'None', lw = 1, zorder = 100)

	ax = ppl.subplot(121)

	for T in markers:
		ax.plot(
			[],
			[],
			marker = markers[T],
			mec = 'k',
			label = f'$CO_2$ ({T}$\\,$°C)' if T else '$CaCO_3$',
			**kwm,
		)

	ax.plot([], [], ls = 'None', marker = 'None', label = '\n\n\n')

	dx = 0.04 * 0
	for k, ref in enumerate(refs):
		for r in data:
			if r['Ref'] == ref and r['Tacid']:
				if 'eD17O_VSMOW' in r:
					ax.errorbar(
						X[r['Sample']] + k * dx,
						r['D17O_VSMOW'],
						r['eD17O_VSMOW'],
						ecolor = colors[r['Ref']],
						**kwe,
					)
				ax.plot(
					X[r['Sample']] + k * dx,
					r['D17O_VSMOW'],
					marker = markers[r['Tacid']],
					mec = colors[r['Ref']],
					**kwm,
				)

		x = [X[_['Sample']] + k * dx for _ in refs[ref] if _['Tacid']]
		y = [_['D17O_VSMOW'] for _ in refs[ref] if _['Tacid']]
		ax.plot(x, y, color = colors[ref], label = ref, **kwl)

	k += 1

	for ref in [
		'Wostbrock et al. (2020)',
		'Ellis & Passey (2023)',
	]:
		x = [
			X[_['Sample']] + k * dx
			for _ in data
			if _['Tacid'] is None and _['Ref'] == ref
		]
		y = [_['D17O_VSMOW'] for _ in data if _['Tacid'] is None and _['Ref'] == ref]
		ey = [_['eD17O_VSMOW'] for _ in data if _['Tacid'] is None and _['Ref'] == ref]
		ax.errorbar(x, y, ey, ecolor = colors[ref], **kwe)
		ax.plot(x, y, marker = markers[None], mec = colors[ref], **kwm)
		ax.plot(x, y, color = colors[ref], **kwl)

	if dx == 0:
		xmin = -0.2
		xmax = 0.2
	else:
		xmin = -0.1
		xmax = 0.47
	ymin = -0.28
	ymax = -0.03

	for s in X:
		ax.plot(
			[X[s] + xmax, X[s] + xmax, X[s] + xmin, X[s] + xmin, X[s] + xmax],
			[ymax, ymin, ymin, ymax, ymax],
			'k-',
			lw = 0.7,
			dashes = (3, 2),
		)
		ax.text(
			X[s] + (xmin + xmax) / 2,
			ymin - 0.007,
			chop_name(s),
			va = 'top',
			ha = 'center',
			weight = 'bold',
			color = 'k',
			size = 9,
			linespacing = 1,
		)

	ax.set_ylabel('Δ$^{17}O_{VSMOW}$ (‰)')

	ax.set_xticks([])
	ax.set_xlim((-1 + xmin - 0.5, 1 + xmax + 0.5))
	ax.set_ylim((ymin - 0.035, ymax + 0.01))

	ax.legend(
		loc = 'center left',
		bbox_to_anchor = (1.07, 0.48),
		fontsize = 9,
		labelspacing = 0.25,
		frameon = False,
	)

	y1, y2 = ax.get_ylim()
	h = y2 - y1

	ax = ppl.subplot(122)

	dx = 0.00 * 1
	for k, ref in enumerate(refs):
		for r in data:
			if r['Sample'] == 'IAEA603' or 'D17O_603' not in r:
				continue
			if r['Ref'] == ref and r['Tacid']:
				ax.errorbar(
					X[r['Sample']] + k * dx,
					r['D17O_603'],
					r['eD17O_603'],
					ecolor = colors[r['Ref']],
					**kwe,
				)
				ax.plot(
					X[r['Sample']] + k * dx,
					r['D17O_603'],
					marker = markers[r['Tacid']],
					mec = colors[r['Ref']],
					**kwm,
				)

	k += 1

	kwm['zorder'] -= 10
	kwe['zorder'] -= 10

	for ref in [
		'Wostbrock et al. (2020)',
		'Ellis & Passey (2023)',
	]:
		x = [
			X[_['Sample']] + k * dx
			for _ in data
			if _['Tacid'] is None and _['Ref'] == ref and _['Sample'] != 'IAEA603'
		]
		y = [
			_['D17O_603']
			for _ in data
			if _['Tacid'] is None and _['Ref'] == ref and _['Sample'] != 'IAEA603'
		]
		ey = [
			_['eD17O_603']
			for _ in data
			if _['Tacid'] is None and _['Ref'] == ref and _['Sample'] != 'IAEA603'
		]
		ax.errorbar(x, y, ey, ecolor = colors[ref], **kwe)
		ax.plot(x, y, marker = markers[None], mec = colors[ref], **kwm)

	x1, x2, y1, y2 = ax.axis()
	y1, y2 = (y1 + y2 - h) / 2, (y1 + y2 + h) / 2
	ax.axis([x1 - 0.5, x2 + 0.5, y1, y2])

	ymin = -0.085
	ymax = 0.160

	for s in X:
		if s == 'IAEA603':
			continue
		ax.plot(
			[X[s] + xmax, X[s] + xmax, X[s] + xmin, X[s] + xmin, X[s] + xmax],
			[ymax, ymin, ymin, ymax, ymax],
			'k-',
			lw = 0.7,
			dashes = (3, 2),
		)
		ax.text(
			X[s] + (xmin + xmax) / 2,
			ymin - 0.007,
			chop_name(s),
			va = 'top',
			ha = 'center',
			weight = 'bold',
			color = 'k',
			size = 9,
			linespacing = 1,
		)

	ax.set_xticks([])

	ax.set_ylabel('Δ$^{17}O_{603}$ (‰)')

	fig.savefig(filename)
	ppl.close(fig)


if __name__ == '__main__':
	"""IMPORT PUBLISHED DATA"""
	with open('input/published.csv') as fid:
		pubdata = [l.strip().split('\t') for l in fid.readlines() if len(l) > 1]
	pubdata = [{k: v for k, v in zip(pubdata[0], l) if v} for l in pubdata[1:]]
	d18Onominal = {'IAEA603': -2.37, 'NBS18': -23.01, 'NBS19': -2.2}
	acidfrac = {None: 1.0, '25': 1.01025, '70': 1.00871, '90': 1.00813}
	for r in pubdata:
		if 'Tacid' not in r:
			r['Tacid'] = None
		for k in r:
			if 'D17' in k:
				r[k] = float(r[k])
		r['d18O_VSMOW'] = (1000 + d18Onominal[r['Sample']]) * acidfrac[
			r['Tacid']
		] * alpha18_CO2g_H2O - 1000
		r['D17O_VSMOW'] = r['D17O'] - r['D17O_of_SLAP'] * np.log(
			1 + r['d18O_VSMOW'] / 1000
		) / np.log(1 - 55.5e-3)
		if 'eD17O' in r:
			r['eD17O_VSMOW'] = r['eD17O']

	# 	pubdata = [r for r in pubdata if r['Ref'] in [
	# 		"Passey et al. (2014)",
	# 		"Barkan et al. (2019)",
	# 		"Passey & Ji (2019)",
	# 		"Fosu et al. (2020)",
	# 		"Sha et al. (2020)",
	# 		"Wostbrock et al. (2020)",
	# 		"Ellis & Passey (2023)",
	# 	]]

	"""IMPORT OUR DATA"""
	with open('../7-ref-materials/output/samples.csv') as fid:
		ourdata = [l.strip().split(',') for l in fid.readlines()]
	ourdata = [{k: v for k, v in zip(ourdata[0], l) if v} for l in ourdata[1:]]
	ourdata = {r['Sample']: r for r in ourdata}

	for s in ['NBS19', 'NBS18', 'IAEA603']:
		pubdata.append(
			{
				'Sample': s,
				'Tacid': '90',
				'D17O': float(ourdata[s]['D17O_VSMOW']),
				'D17O_VSMOW': float(ourdata[s]['D17O_VSMOW']),
				'eD17O': float(ourdata[s]['95CL_D17O_VSMOW']),
				'eD17O_VSMOW': float(ourdata[s]['95CL_D17O_VSMOW']),
				'Ref': 'this study',
			}
		)

	"""COMPUTE D17O_603"""
	for r in pubdata:
		try:
			ref_IAEA603 = [
				_
				for _ in pubdata
				if _['Ref'] == r['Ref']
				and _['Tacid'] == r['Tacid']
				and _['Sample'] == 'IAEA603'
			][0]
		except IndexError:
			continue
		r['D17O_603'] = r['D17O_VSMOW'] - ref_IAEA603['D17O_VSMOW']
		if 'eD17O_VSMOW' not in r:
			r['eD17O_VSMOW'] = 0
		r['eD17O_603'] = (
			r['eD17O_VSMOW'] ** 2 + ref_IAEA603['eD17O_VSMOW'] ** 2
		) ** 0.5
		if r['Sample'] == 'NBS18':
			logger.info(
				f"{r['Ref']:<24s}: Δ17O_603 of NBS18 = {r['D17O_603']*1e3:.1f} ± {r['eD17O_603']*1e3:.1f}"
			)

	"""SAVE RESULTS"""
	with open('output/compilation.csv', 'w') as fid:
		fid.write(
			f'Sample,Ref,Tacid,D17O_reported,D17O_of_SLAP,D17O_VSMOW,eD17O_VSMOW,D17O_603,eD17O_603'
		)
		for r in pubdata:
			if 'D17O_of_SLAP' not in r:
				r['D17O_of_SLAP'] = 0
			if 'eD17O' not in r:
				stre = ''
			else:
				stre = f"{r['eD17O']:.4f}"
			if 'D17O_603' not in r or r['Sample'] == 'IAEA603':
				str603 = ','
			else:
				str603 = f"{r['D17O_603']:.4f},{r['eD17O_603']:.4f}"

			fid.write(
				f"\n{r['Sample']},{r['Ref']},{r['Tacid'] if r['Tacid'] else ''},{r['D17O']:.4f},{r['D17O_of_SLAP']:.3f},{r['D17O_VSMOW']:.3f},{stre},{str603}"
			)

	"""PLOT EVERYTHING"""
	pub_comparison_plot(pubdata)
