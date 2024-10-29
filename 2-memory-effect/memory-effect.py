import sys

sys.path.append('../lib')

from mylib import *

RAWDATA_FILENAME = 'input/memory-effect-rawdata.csv'

logger.remove(0)
logger.add('output/memory-effect.log', mode = 'w', format = '[{line:3.0f}]    {message}')

ppl.style.use('../lib/mydefault.mplstyle')

if __name__ == '__main__':
	data = read_rawdata([RAWDATA_FILENAME])
	sigma = {}

	fig = ppl.figure(figsize = (4.5, 4.3))
	fig.subplots_adjust(0.15, 0.11, 0.98, 0.94, 0.5, 0.362)

	for field, axa, axb in [
		('d636', ppl.subplot(221), ppl.subplot(222)),
		('d628', ppl.subplot(223), ppl.subplot(224)),
	]:
		data = data[: len(data) // 2 * 2]  # truncate data to an even number of analyses

		# R = d6xx RESIDUALS
		R = np.array(
			[
				r[field]
				- np.mean([_[field] for _ in data if _['Sample'] == r['Sample']])
				for r in data
			]
		)

		# 95% CL for one analysis
		_p = 0.95
		sigma[field] = R.std(ddof = 4)
		logger.info(
			f'Pooled repeatability for {field} measurements is {sigma[field]*1000:.1f} ppm.'
		)
		CL = sigma[field] * stats.t.ppf(1 - (1 - _p) / 2, len(data) - 4)
		logger.info(f'95 % CL for one {field} analysis is {CL*1000:.1f} ppm.')

		AX, AY, BX, BY = [], [], [], []
		for k, r in enumerate(data):
			gas, rank = r['Sample']

			if rank == '1':
				X = r[field]
				if gas == 'A':
					AX.append(X)
				else:
					BX.append(X)
			else:
				Y = r[field]
				if gas == 'A':
					AY.append(Y)
				else:
					BY.append(Y)

			if gas == 'A':
				ax = axa
			else:
				ax = axb

			if rank == '2':
				ax.errorbar(
					X,
					Y,
					CL,
					CL,
					elinewidth = 0.7,
					capsize = 0,
					capthick = 0.7,
					marker = 'None',
					ls = 'None',
					ecolor = [0.6] * 3,
				)
				ax.plot(X, Y, 's', mec = [0.5] * 3, mew = 0.7, mfc = [1] * 3, ms = 4)

		large_arrow = np.array(
			[
				(-0.5, 0.5),
				(0.5, 1.5),
				(1, 1),
				(2, 2),
				(2, -2),
				(-2, -2),
				(-1, -1),
				(-1.5, -0.5),
				(-0.5, 0.5),
			]
		)

		for ax, X, Y, tank, Xarrow, Yarrow, Xarrowtext, Yarrowtext, arrow_direction in [
			(axa, AX, AY, 'A', 0.75, 0.25, 0.8, 0.2, 1),
			(axb, BX, BY, 'B', 0.25, 0.75, 0.23, 0.77, -1),
		]:
			ax.add_patch(
				Ellipse(
					xy = (np.mean(X), np.mean(Y)),
					width = 2
					* CL
					/ (len(X) * stats.chi2.ppf(0.95, 1) / stats.chi2.ppf(0.95, 2))
					** 0.5,
					height = 2
					* CL
					/ (len(X) * stats.chi2.ppf(0.95, 1) / stats.chi2.ppf(0.95, 2))
					** 0.5,
					fc = 'none',
					ec = 'k',
					lw = 1,
					# 					ls = (0, (3, 1)),
					zorder = 200,
				)
			)
			ax.autoscale_view()

			ax.set_xlabel(f'First δ$_{{{field[1:]}}}$ value (‰)')
			ax.set_ylabel(f'Second δ$_{{{field[1:]}}}$ value (‰)')
			if field == 'd636':
				ax.set_title(f'Tank {tank}')

			x1, x2, y1, y2 = ax.axis()
			xmin = min(x1, y1) - 0.022
			xmax = max(x2, y2) + 0.022
			ax.plot(
				[xmin, xmax],
				[xmin, xmax],
				'-',
				color = 'k',
				lw = 0.7,
				dashes = (6, 2, 2, 2),
				alpha = 0.5,
				zorder = -100,
			)
			ax.axis([xmin, xmax, xmin, xmax])

			ax.xaxis.set_major_locator(
				ticker.MultipleLocator(
					3e-2, offset = 0.0 if tank == 'A' or field == 'd628' else 0.01
				)
			)
			ax.yaxis.set_major_locator(
				ticker.MultipleLocator(
					3e-2, offset = 0.0 if tank == 'A' or field == 'd628' else 0.01
				)
			)

			ax.text(
				0.15 if tank == 'B' else 0.86,
				0.495,
				'?',
				transform = ax.transAxes,
				ha = 'center',
				va = 'center',
				size = 10,
				weight = 'bold',
				color = [0.3] * 3,
				linespacing = 0.9,
			)
			ax.arrow(
				*(
					(
						0.24,
						0.5,
						-0.11,
						0,
					)
					if tank == 'B'
					else (
						0.76,
						0.5,
						+0.11,
						0,
					)
				),
				width = 0.18,
				head_length = 0.1,
				head_width = 0.28,
				transform = ax.transAxes,
				lw = 0.75,
				ec = [0.5] * 3,
				fc = [0.9] * 3,
			)

	# 			ax.plot(
	# 				Xarrow,
	# 				Yarrow,
	# 				ls = "None",
	# 				marker = large_arrow * arrow_direction,
	# 				mew = 2,
	# 				mec = [0.75] * 3,
	# 				mfc = [0.85] * 3,
	# 				ms = 45,
	# 				transform = ax.transAxes,
	# 				zorder = -100,
	# 			)
	# 			ax.text(
	# 				Xarrowtext + (.02 if tank == 'A' else 0),
	# 				Yarrowtext - (.02 if tank == 'A' else 0),
	# 				"memory\neffect",
	# 				va = "center",
	# 				ha = "center",
	# 				color = [0.6] * 3,
	# 				size = 8,
	# 				weight = "bold",
	# 				transform = ax.transAxes,
	# 				rotation = 45,
	# 			)

	fig.savefig(
		'output/memory_effect',
		metadata = {
			'Author': AUTHOR,
			'Title': f"""
				Absence of memory effects. When aliquots of the same tank are analyzed consecutively,
				potential memory effects should bias the first analysis as indicated by grey arrows.
				Individual analyses are shown as squares with 95 % error bars. Black ellipses
				correspond to joint 95 % confidence limits for the mean of the first and second
				measured values, based on analytical repeatabilities of {sigma['d636']*1000:.1f} ppm
				and {sigma['d628']*1000:.1f} ppm for δ636 and δ628, respectively.
				""".replace('\t', '').replace('\n', ' ')[1:-1],
		},
	)
	ppl.close(fig)
