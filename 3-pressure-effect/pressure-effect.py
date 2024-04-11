import sys

sys.path.append('../lib')
from mylib import *

RAWDATA_FILENAME = 'input/pressure-effect-rawdata.csv'
INSTRUMENTAL_REPEATABILITY = 4

logger.remove(0)
logger.add('output/pressure-effect.log', mode = 'w', format = '[{line:3.0f}]    {message}')

ppl.style.use('../lib/mydefault.mplstyle')

if __name__ == '__main__':
	data = read_rawdata([RAWDATA_FILENAME])

	fig = ppl.figure(figsize = (3.15, 2.5))
	fig.subplots_adjust(0.185, 0.175, 0.98, 0.98)
	ax = ppl.subplot(111)

	P = np.array([_['P'] for _ in data])
	d627 = np.array([_['d627'] for _ in data])
	d628 = np.array([_['d628'] for _ in data])
	D17O = 1e6 * (np.log(1 + d627 / 1e3) - 0.528 * np.log(1 + d628 / 1e3))
	_p = 0.95
	CL = INSTRUMENTAL_REPEATABILITY * 1.96

	ax.errorbar(
		P,
		(D17O - D17O.mean()),
		CL,
		elinewidth = 1,
		capsize = 4,
		capthick = 1,
		marker = 'None',
		ls = 'None',
		ecolor = 'k',
	)
	ax.plot(P, (D17O - D17O.mean()), 'wo', mec = 'k', mew = 1, ms = 6)
	ax.set_xlabel('Analyte pressure (mbar)')
	ax.set_ylabel('Δ$^{17}$O$_{WG}$ residuals (ppm)')

	ax.axis([4.89, 5.11, -17, 16])
	ax.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
	ax.grid(False, axis = 'x')

	Pmean = 4.9966
	total_P_range = np.abs((P - Pmean)).max()
	Phwidth = 0.002
	pcolor = [0.5, 0.75, 1]
	ax.axvspan(Pmean - Phwidth, Pmean + Phwidth, color = pcolor, lw = 0, alpha = 0.5)
	ax.annotate(
		'P range during\nnormal operation\n',
		xy = (Pmean - Phwidth * 1.5, -15),
		xytext = (4.935, -15),
		va = 'center',
		ha = 'center',
		size = 9,
		arrowprops = dict(
			arrowstyle = '->',
			lw = 0.8,
			ec = [_ * 0.75 for _ in pcolor],
		),
		color = [_ * 0.75 for _ in pcolor],
	)

	logger.info(
		f'SD of Δ17O_WG measurements is {D17O.std(ddof = 1):.4} ‰ (N = {D17O.size})'
	)
	Zvalues = (D17O - D17O.mean()) / INSTRUMENTAL_REPEATABILITY
	pvalue = stats.kstest(Zvalues, 'norm', (0, 1)).pvalue
	logger.info(
		f'Kolmogorov-Smirnoff test that the Z-values are normally distributed: p-value = {pvalue:.2f}'
	)

	fig.savefig(
		'output/P_effects',
		metadata = {
			'Author': AUTHOR,
			'Title': f"""
				Absence of pressure effects. Repeated analyses of the same gas, with pressures
				varying by ±{total_P_range:.3f} mbar, yield statistically indistinguishable
				Δ17O values (SD = {D17O.std(ddof = 1):.1f} ppm). By comparison,
				analyte pressure during routine measurements remain within ±{Phwidth} mbar.
				""".replace('\t', '').replace('\n', ' ')[1:-1],
		},
	)
	ppl.close(fig)
