import sys

sys.path.append('../lib')
from mylib import *
from numpy.random import normal, seed
from matplotlib.patches import Rectangle
from matplotlib import transforms

ppl.style.use('../lib/mydefault.mplstyle')

seed(17)


def d13C_example(
	anchors = {
		'IAEA603': {'d13C_VPDB': 2.46, 'N': 3},
		'IAEA612': {'d13C_VPDB': -36.722, 'N': 3},
	},
	unknowns = {
		'X1': {'d13C_VPDB': -5.0, 'N': 3},
		'X2': {'d13C_VPDB': -45.0, 'N': 3},
	},
	sessions = {
		'A': {'d13C_VPDB_wg': -25.0, 'f': 1.1},
		'B': {'d13C_VPDB_wg': -25.0, 'f': 0.9},
	},
	sigma = 2,
	Nr = 3,
	sessioncolor = (0, 0.35, 0.7),
	samplecolor = (0.9, 0, 0),
	anchorcolor = (0, 0, 0),
):
	data = []
	samples = anchors | unknowns
	for s in sessions:
		for r in samples:
			for _ in range(samples[r]['N']):
				data.append(
					{
						'Sample': r,
						'Session': s,
						'd636': (
							1000
							* sessions[s]['f']
							* (
								(1000 + samples[r]['d13C_VPDB'])
								/ (1000 + sessions[s]['d13C_VPDB_wg'])
								- 1
							)
							+ sigma * normal()
						),
					}
				)

	# 	data[0]['d636'] -= 10
	# 	data = [r for r in data if not(r['Session'] == 'B' and r['Sample'] in anchors)]

	S = stdz.standardize(data, anchors, key_in = 'd636', key_out = 'd13C_VPDB')

	fig = ppl.figure(figsize = (3.15, 6.5))
	_top = 0.76
	_bottom = 0.07
	fig.subplots_adjust(0.165, _bottom, 0.98, _top, -0.77, -0.31)

	ax = {
		'A': ppl.subplot(322),
		'B': ppl.subplot(326),
		'C': ppl.subplot(15, 2, 16),
		'a': ppl.subplot(3, 30, 1),
		'b': ppl.subplot(3, 30, 61),
	}

	xmin, xmax, ymin, ymax = 0, 0, 0, 0
	sessionline = {}

	for s in sessions:
		# slope and intercept (_a, _b)
		_a = (
			1000
			* S['sessions'][s]['d13C_VPDB_scaling']
			/ (1000 + S['sessions'][s]['d13C_VPDB_of_wg'])
		)
		_b = (
			1000
			* S['sessions'][s]['d13C_VPDB_scaling']
			* (1000 / (1000 + S['sessions'][s]['d13C_VPDB_of_wg']) - 1)
		)

		# Jacobian matrix of (_a, _b) vs (d13C_VPDB_scaling, d13C_VPDB_of_wg)
		_J = np.array(
			[
				[
					# d(_a)/d(d13C_VPDB_scaling)
					1000 / (1000 + S['sessions'][s]['d13C_VPDB_of_wg']),
					# d(_b)/d(d13C_VPDB_scaling)
					1000 * (1000 / (1000 + S['sessions'][s]['d13C_VPDB_of_wg']) - 1),
				],
				[
					# d(_a)/d(d13C_VPDB_of_wg)
					-1000
					* S['sessions'][s]['d13C_VPDB_scaling']
					/ (1000 + S['sessions'][s]['d13C_VPDB_of_wg']) ** 2,
					# d(_b)/d(d13C_VPDB_of_wg)
					-1000
					* S['sessions'][s]['d13C_VPDB_scaling']
					* 1000
					/ (1000 + S['sessions'][s]['d13C_VPDB_of_wg']) ** 2,
				],
			]
		)

		# Covariance matrix of (_a, _b)
		if s == 'A':
			_i, _j = 0, 1
		else:
			_i, _j = 2, 3
		_CMab = _J.T @ S['bestfit'].covar[[_i, _j], :][:, [_i, _j]] @ _J

		sessionline[s] = (_a, _b, _CMab)

		G = [r for r in S['data'] if r['Session'] == s and r['Sample'] in anchors]
		ax[s].plot(
			[samples[r['Sample']]['d13C_VPDB'] for r in G],
			[r['d636'] for r in G],
			'x',
			mec = anchorcolor,
			mew = 1,
			ms = 6,
			alpha = 0.75,
			zorder = 100,
		)

		(leg_obs_std,) = ax[s.lower()].plot(
			[0 for r in G],
			[r['d636'] for r in G],
			'x',
			mec = anchorcolor,
			mew = 1,
			ms = 6,
			alpha = 0.8,
			zorder = 100,
			label = 'single analysis (standard)',
		)

		G = [r for r in S['data'] if r['Session'] == s and r['Sample'] in unknowns]
		ax[s].plot(
			[(r['d636'] - _b) / _a for r in G],
			[r['d636'] for r in G],
			'x',
			mec = samplecolor,
			mew = 1,
			ms = 6,
			alpha = 0.8,
			zorder = 100,
			# 			label = 'single analysis (unknown sample)',
		)

		(leg_obs_uk,) = ax[s.lower()].plot(
			[0 for r in G],
			[r['d636'] for r in G],
			'x',
			mec = samplecolor,
			mew = 1,
			ms = 6,
			alpha = 0.8,
			zorder = 100,
			label = 'single analysis (unknown sample)',
		)

		x1, x2, y1, y2 = ax[s].axis()
		xmin = np.min((xmin, x1))
		xmax = np.max((xmax, x2))
		ymin = np.min((ymin, y1))
		ymax = np.max((ymax, y2))

	xmin -= S['samples']['X2']['95CL_d13C_VPDB']
	# 	xmax += 4

	for s in sessions:
		xi = np.array([xmin, xmax])
		true_yi = (
			1000
			* sessions[s]['f']
			* ((1000 + xi) / (1000 + sessions[s]['d13C_VPDB_wg']) - 1)
		)
		(leg_true_slope,) = ax[s].plot(
			xi,
			true_yi,
			'-',
			color = sessioncolor,
			lw = 0.75,
			dashes = (6, 2, 2, 2),
			label = 'unknown true session slope',
		)

		xi = np.linspace(xmin, xmax)
		_a, _b, _CMab = sessionline[s]

		est_yi = _a * xi + _b

		est_yi_se = (_CMab[0, 0] * xi**2 + _CMab[1, 1] + 2 * _CMab[0, 1] * xi) ** 0.5

		ax[s].fill_between(
			xi,
			est_yi + S['t95'] * est_yi_se,
			est_yi - S['t95'] * est_yi_se,
			color = sessioncolor,
			lw = 0,
			alpha = 0.2,
		)

		(leg_est_slope,) = ax[s].plot(
			[],
			[],
			'-',
			color = sessioncolor,
			lw = 4,
			alpha = 0.2,
			label = 'estimated session slope',
		)

		# 		ax[s].text(
		# 			0.04, 0.96,
		# 			f"""$δ^{{13}}C_{{VPDB}}^{{WG}}$ = {S['sessions'][s]['d13C_VPDB_of_wg']:.2f} ± {S['sessions'][s]['SE_d13C_VPDB_of_wg']:.2f} ‰
		# 			$δ_{{636}}$ scaling = {S['sessions'][s]['d13C_VPDB_scaling']:.2f} ± {S['sessions'][s]['SE_d13C_VPDB_scaling']:.2f}""".replace('\t', ''),
		# 			va = 'top',
		# 			ha = 'left',
		# 			transform = ax[s].transAxes,
		# 			size = 9,
		# 			color = sessioncolor,
		# 			)

		ax[s].text(
			0.45,
			0.9,
			f'Session {s}',
			va = 'top',
			ha = 'center',
			transform = ax[s].transAxes,
			size = 9,
			color = sessioncolor,
			weight = 'bold',
		)

		ax[s].errorbar(
			S['sessions'][s]['d13C_VPDB_of_wg'],
			0,
			None,
			S['t95'] * S['sessions'][s]['SE_d13C_VPDB_of_wg'],
			ecolor = sessioncolor,
			elinewidth = 1,
			capthick = 1,
			capsize = 3.5,
		)
		(leg_est_wg,) = ax[s].plot(
			S['sessions'][s]['d13C_VPDB_of_wg'],
			0,
			'wo',
			mec = sessioncolor,
			mew = 1,
			ms = 5,
			label = 'estimated WG composition',
		)

		for x in unknowns:
			x0 = S['samples'][x]['d13C_VPDB']
			xe = S['samples'][x]['95CL_d13C_VPDB']
			ax[s].axvspan(
				x0 - xe,
				x0 + xe,
				color = samplecolor,
				lw = 0,
				alpha = 0.1,
			)
		# 			ax[s].errorbar(
		# 				S['samples'][x]['d13C_VPDB'],
		# 				1000 * S['sessions'][s]['d13C_VPDB_scaling'] * (
		# 					(1000 + S['samples'][x]['d13C_VPDB']) / (1000 + S['sessions'][s]['d13C_VPDB_of_wg']) - 1
		# 				),
		# 				None, S['samples'][x]['95CL_d13C_VPDB'],
		# 				ecolor = samplecolor, elinewidth = 1, capthick = 1, capsize = 2)
		# 			ax[s].plot(
		# 				S['samples'][x]['d13C_VPDB'],
		# 				1000 * S['sessions'][s]['d13C_VPDB_scaling'] * (
		# 					(1000 + S['samples'][x]['d13C_VPDB']) / (1000 + S['sessions'][s]['d13C_VPDB_of_wg']) - 1
		# 				),
		# 				'wo', mec = samplecolor, mew = 1, ms = 4)

		ax[s].axis([xmin, xmax, ymin, ymax])

	ax['B'].set_xlabel('True $δ^{13}C_{VPDB}$ (‰)')
	ax['a'].set_ylabel('Observed $δ_{636}$ (‰)', labelpad = -1)
	ax['b'].set_ylabel('Observed $δ_{636}$ (‰)', labelpad = -1)

	ax['B'].sharex(ax['A'])
	ax['C'].sharex(ax['A'])

	ax['a'].sharey(ax['A'])
	ax['b'].sharey(ax['B'])

	ax['A'].tick_params(
		labelbottom = False,
		labelleft = False,
	)

	ax['B'].tick_params(
		labelleft = False,
	)

	ax['C'].tick_params(
		left = False,
		labelleft = False,
		labelbottom = False,
	)

	ax['a'].tick_params(
		bottom = False,
		labelbottom = False,
	)

	ax['b'].tick_params(
		bottom = False,
		labelbottom = False,
	)

	# 	ax['B'].tick_params(
	# 		labeltop = False,
	# 		labelbottom = True,
	# 		top = False,
	# 		bottom = True,
	# 	)

	# 	txt = '\n'
	# 	for x in unknowns:
	# 		txt += f"  $δ^{{13}}C_{{VPDB}}^{{{x}}}$ = {S['samples'][x]['d13C_VPDB']:.2f} ± {S['samples'][x]['SE_d13C_VPDB']:.2f}  \n"

	# 	ax['B'].text(
	# 		0.52, 0.6,
	# 		txt[1:-1],
	# 		va = 'center',
	# 		ha = 'center',
	# 		size = 9,
	# 		color = 'r',
	# 		bbox = dict(facecolor = 'w', edgecolor = 'k', pad = 0.45, boxstyle = 'round'),
	# 		transform = fig.transFigure,
	# 		linespacing = 2,
	# 		)

	for x in anchors:
		ax['C'].plot(
			samples[x]['d13C_VPDB'],
			0,
			'D',
			mfc = (*anchorcolor, 0.5),
			mec = anchorcolor,
			mew = 1,
			ms = 5,
		)
		(leg_true_std,) = ax['A'].plot(
			[],
			[],
			'D',
			mfc = (*anchorcolor, 0.5),
			mec = anchorcolor,
			mew = 1,
			ms = 5,
			label = 'known true composition of standard',
		)

	# 	ax['C'].text(
	# 		np.mean([S['samples'][_]['d13C_VPDB'] for _ in anchors]), 0,
	# 		'known true values',
	# 		va = 'center',
	# 		ha = 'center',
	# 		size = 8,
	# 		color = anchorcolor,
	# 		style = 'italic',
	# # 		bbox = dict(fc = 'w', lw = 0, alpha = 1)
	# 	)

	(leg_true_wg,) = ax['C'].plot(
		sessions['A']['d13C_VPDB_wg'],
		0,
		'D',
		mfc = (*sessioncolor, 0.5),
		mec = sessioncolor,
		mew = 1,
		ms = 5,
		label = 'unknown true composition of WG',
	)

	for x in unknowns:
		ax['C'].plot(
			samples[x]['d13C_VPDB'],
			-1,
			'D',
			mfc = (*samplecolor, 0.5),
			mec = samplecolor,
			mew = 1,
			ms = 5,
		)
		(leg_true_uk,) = ax['A'].plot(
			[],
			[],
			'D',
			mfc = (*samplecolor, 0.5),
			mec = samplecolor,
			mew = 1,
			ms = 5,
			label = 'unknown true composition of sample',
		)

		ax['C'].errorbar(
			S['samples'][x]['d13C_VPDB'],
			1,
			None,
			S['samples'][x]['95CL_d13C_VPDB'],
			ecolor = samplecolor,
			elinewidth = 1,
			capthick = 1,
			capsize = 3.5,
		)
		ax['C'].plot(
			S['samples'][x]['d13C_VPDB'],
			1,
			'wo',
			mec = samplecolor,
			mew = 1,
			ms = 5,
		)
		(leg_est_uk,) = ax['A'].plot(
			[],
			[],
			'wo',
			mec = samplecolor,
			mew = 1,
			ms = 5,
			label = 'estimated composition of sample',
		)

		x0 = S['samples'][x]['d13C_VPDB']
		xe = S['samples'][x]['95CL_d13C_VPDB']

		trans = transforms.blended_transform_factory(ax['C'].transData, fig.transFigure)

		x0 = S['samples'][x]['d13C_VPDB']
		xe = S['samples'][x]['95CL_d13C_VPDB']
		ax['C'].axvspan(
			x0 - xe,
			x0 + xe,
			color = samplecolor,
			lw = 0,
			alpha = 0.1,
		)
	#
	# 		ax['C'].add_patch(
	# 			Rectangle(
	# 				(x0-xe, _bottom),
	# 				2*xe, _top - _bottom,
	# 				fc = samplecolor,
	# 				alpha = 0.1,
	# 				lw = 0,
	# 				clip_on = False,
	# 				transform = trans,
	# 			)
	# 		)

	# 		ax['C'].text(
	# 			S['samples'][x]['d13C_VPDB'], 1, f'{x}\n',
	# 			va = 'center',
	# 			ha = 'center',
	# 			size = 9,
	# 			linespacing = 2,
	# 			color = samplecolor,
	# 		)

	# 	ax['C'].text(
	# 		(S['samples']['X1']['d13C_VPDB'] + S['samples']['X2']['d13C_VPDB'])/2, 1,
	# 		'final estimates',
	# 		va = 'center',
	# 		ha = 'center',
	# 		size = 8,
	# 		color = samplecolor,
	# 		style = 'italic',
	# # 		bbox = dict(fc = 'w', lw = 0, alpha = 1)
	# 	)
	# 	ax['C'].text(
	# 		(S['samples']['X1']['d13C_VPDB'] + S['samples']['X2']['d13C_VPDB'])/2, -1,
	# 		'unknown true values',
	# 		va = 'center',
	# 		ha = 'center',
	# 		size = 8,
	# 		color = samplecolor,
	# 		style = 'italic',
	# # 		bbox = dict(fc = 'w', lw = 0, alpha = 1)
	# 	)

	legs = [
		leg_obs_std,
		leg_obs_uk,
		leg_true_std,
		leg_true_uk,
		leg_true_wg,
		leg_est_uk,
		leg_est_wg,
		leg_true_slope,
		leg_est_slope,
	]
	ax['A'].legend(
		legs,
		[l.get_label() for l in legs],
		loc = 'lower center',
		bbox_to_anchor = (0.304, 1.05),
		borderpad = 0.0,
		fontsize = 8,
		labelspacing = 0.4,
		handlelength = 2.7,
		# 		bbox_to_anchor = (.3851, .96),
		# 		borderpad = 1.7,
		frameon = False,
	)

	ax['C'].axis([None, None, -3, 3])

	fig.savefig('output/d13C_stdz_example')
	ppl.close(fig)


def kde_example(plot_margins = 0.5, robust_cov_estimator = True):
	with open('../4-linearity/output/residuals.csv') as fid:
		data = [l.split(',') for l in fid.readlines()[1:]]
		d18, d17, D17 = np.array(data, dtype = float).T

	fig = ppl.figure(figsize = (3.15, 5.71))
	fig.subplots_adjust(0.23, 0.08, 0.98, 0.99, 0.2, 0.15)

	for x, y, ax, delta in [
		(d18, d17, ppl.subplot(211), 'δ'),
		(d18, D17, ppl.subplot(212), 'Δ’'),
	]:
		ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
		ax.yaxis.set_major_locator(
			ticker.MultipleLocator(0.01 if delta == 'Δ’' else 0.5)
		)
		if delta == 'Δ’':
			ax.yaxis.set_major_formatter(
				ticker.FuncFormatter(lambda x, pos: f'${x:+.2f}$' if x else '$0$')
			)
		else:
			ax.yaxis.set_major_formatter(
				ticker.FuncFormatter(lambda x, pos: f'${x:+.1f}$' if x else '$0$')
			)
		ax.xaxis.set_major_formatter(
			ticker.FuncFormatter(lambda x, pos: f'${x:+.1f}$' if x else '$0$')
		)

		xmean = x.mean()
		xwidth = np.abs(x - xmean).max() * (1 + plot_margins)
		xmin, xmax = (
			xmean - xwidth,
			xmean + xwidth,
		)

		ymean = y.mean()
		ywidth = np.abs(y - ymean).max() * (1 + plot_margins)
		ymin, ymax = (
			ymean - ywidth,
			ymean + ywidth,
		)

		ax.set_xlim(xmin, xmax)
		if delta == 'Δ’':
			ax.set_ylim(ymin, ymax)
		else:
			ax.set_ylim(xmin, xmax)

		ax.plot(x, y, 'k+', mew = 0.7, ms = 9, alpha = 1)

		ax.set_ylabel('$' + delta + '{}^{17}$O$_{VSMOW}$ residuals (‰)')

		CM = isofunctions.estimate_covariance(np.array([x, y]).T)
		kw = dict(fc = 'None', ec = 'k', lw = 1, alpha = 0.25)
		w, h, r = isofunctions.cov_ellipse(CM, r = 1)
		for f in [1, 2, 3, 4] if delta == 'Δ’' else [4]:
			ax.add_patch(Ellipse(xy = (0, 0), width = w * f, height = h * f, angle = r, **kw))

	# 		kw = dict(color = kw["fc"], lw = 6, ls = '-', alpha = kw['alpha'])
	# 		ax.plot([], [], label = "95 % confidence from covariance estimate", **kw)

	# 		ax.legend(fontsize = 7, handlelength = 2.5)

	ax.set_xlabel('δ$^{18}$O$_{VSMOW}$ residuals (‰)')

	fig.savefig('output/kde')
	ppl.close(fig)


if __name__ == '__main__':
	d13C_example()
	kde_example()
