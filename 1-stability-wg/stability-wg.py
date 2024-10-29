import sys

sys.path.append('../lib')
from mylib import *

RAWDATA_FILENAME = 'input/stability-wg-rawdata.csv'
DRIFTING_RAWDATA_FILENAME = 'input/stability-wg-rawdata_drifting.csv'

ppl.style.use('../lib/mydefault.mplstyle')

data = read_rawdata([DRIFTING_RAWDATA_FILENAME])
for r in data:
	r['t'] /= 60
	r['D627'] = (
		np.log(1 + r['d627'] / 1000)
		- isofunctions.LAMBDA_17 * np.log(1 + r['d628'] / 1000)
	) * 1e6

Xwg = [r['t'] for r in data if r['Sample'] == 'WG']
Ywg = [r['D627'] for r in data if r['Sample'] == 'WG']
Xlinde = [r['t'] for r in data if r['Sample'] == 'LINDE']
Ylinde = [r['D627'] for r in data if r['Sample'] == 'LINDE']

didata = read_rawdata([RAWDATA_FILENAME])
for r in didata:
	r['t'] /= 60
	r['D627'] = (
		np.log(1 + r['d627'] / 1000)
		- isofunctions.LAMBDA_17 * np.log(1 + r['d628'] / 1000)
	) * 1e6

Xdi = [r['t'] for r in didata]
Ydi = [r['D627'] for r in didata]

fig = ppl.figure(figsize = (3.15, 4))
fig.subplots_adjust(0.24, 0.11, 0.99, 0.93, 0.5, 0.05)

ax = ppl.subplot(211)
ax.plot(
	Xlinde,
	Ylinde,
	's',
	ms = 4.2,
	mec = 'k',
	mfc = ISO_COLORS['627'],
	mew = 0.7,
	label = 'Tank #2',
	zorder = 100,
)
ax.plot(
	Xwg,
	Ywg,
	'D',
	ms = 4,
	mec = 'k',
	mfc = [0.7] * 3,
	mew = 0.7,
	label = 'Working reference gas',
	zorder = 100,
)

ax.set_ylabel('Δ’$^{17}Ο_{WG}$ (ppm) relative to\nlong-term WG average')
kw_legend = dict(
	loc = 'lower right',
	bbox_to_anchor = (1, 0),
	borderpad = 0.5,
	framealpha = 1,
	fontsize = 9,
	handlelength = 0.5,
	markerfirst = False,
	frameon = True,
	labelspacing = 0.4,
)
ax.legend(**kw_legend)

ymin, ymax = ax.axis()[-2:]
ymin = -40

ax.axis([None, None, ymin, ymax])

ax2 = ppl.subplot(212, sharex = ax)

ax2.plot(
	Xdi,
	Ydi,
	'o',
	ms = 4.2,
	mec = 'k',
	mfc = ISO_COLORS['627'],
	mew = 0.7,
	label = 'Tank #2, WG-corrected',
	zorder = 100,
)

ax2.legend(**kw_legend)

ax2.set_xlabel('time (h)')
ax2.set_ylabel('Δ’$^{17}Ο_{WG}$ (ppm) relative to\nshort-term WG average')
ax2.axis([None, None, ymin, ymax])

ax.tick_params(
	labelbottom = False,
	labeltop = True,
	bottom = False,
	top = True,
)
ax2.tick_params(
	labelbottom = True,
	labeltop = False,
	bottom = True,
	top = False,
)

ax.xaxis.set_major_locator(ticker.MultipleLocator(6))
ax.yaxis.set_major_locator(ticker.MultipleLocator(15))
ax2.yaxis.set_major_locator(ticker.MultipleLocator(15))

fig.savefig(
	'output/instrumental-stability',
	metadata = {
		'Author': AUTHOR,
		'Title': """
			Instrumental stability over a continuous period of 27 hours.
			Upper panel: uncorrected Δ’17Ο values of repeated aliquots from two CO2 tanks,
			relative to the overall average composition of one of the tanks (“working
			reference gas”).
			Lower panel: Δ’17Ο values of the second tank relative to the preceding
			and subsequent working-gas measurements.
			""".replace('\t', '').replace('\n', ' ')[1:-1],
	},
)
ppl.close(fig)


fig = ppl.figure(figsize = (6.62, 3.0))

moz = fig.subplot_mosaic("""
AAAAAAAAAAA....BBBBBBBBBBBBCC
AAAAAAAAAAA....DDDDDDDDDDDDEE
AAAAAAAAAAA....FFFFFFFFFFFFGG
""")
ax, ax2, ax3, ax4 = moz['D'], moz['C'], moz['B'], moz['E']

fig.subplots_adjust(0.075, 0.15, 1, 0.91, 0.25, 0.1)

ax_D627ts = moz['B']
ax_d628ts = moz['D']
ax_d636ts = moz['F']
ax_D627hist = moz['C']
ax_d628hist = moz['E']
ax_d636hist = moz['G']
ax_avar = moz['A']

ax_d628ts.sharex(ax_D627ts)
ax_d636ts.sharex(ax_D627ts)

ax_D627ts.sharey(ax_D627hist)
ax_d628ts.sharey(ax_d628hist)
ax_d636ts.sharey(ax_d636hist)

ax_D627ts.xaxis.set_major_locator(ticker.MultipleLocator(6))
ax_D627ts.yaxis.set_major_locator(ticker.MultipleLocator(20))
ax_d628ts.yaxis.set_major_locator(ticker.MultipleLocator(0.02, offset = 0.01))
ax_d636ts.yaxis.set_major_locator(ticker.MultipleLocator(0.02))

ax_D627ts.tick_params(
	labelbottom = False,
	labeltop = True,
	bottom = False,
	top = True,
)
ax_d628ts.tick_params(
	labelbottom = False,
	labeltop = False,
	bottom = False,
	top = False,
)
ax_d636ts.tick_params(
	labelbottom = True,
	labeltop = False,
	bottom = True,
	top = False,
)

ax_D627hist.tick_params(
	labelleft = False,
	left = False,
)
ax_d628hist.tick_params(
	labelleft = False,
	left = False,
)
ax_d636hist.tick_params(
	labelleft = False,
	left = False,
)

t = [r['t'] for r in didata]

ax = ax_D627ts
Y = [r['D627'] for r in didata]

ax.plot(t, Y, '-', color = ISO_COLORS['627'], lw = 1.5)
# ppl.setp(ax.get_xticklabels(), visible = False)

ax.text(
	0.5,
	0.05,
	f'Δ’$^{{\\mathbf{{17}}}}$O$_{{WG}}$ = {np.mean(Y):.1f} ppm (SD = {stdev(Y):.1f} ppm)',
	transform = ax.transAxes,
	ha = 'center',
	va = 'bottom',
	color = [_ * 0.75 for _ in ISO_COLORS['627']],
	weight = 'bold',
	size = 10,
)
ax.set_ylabel('Δ’$^{17}O_{WG}$ (ppm)')

ax = ax_D627hist

ax.hist(
	Y,
	bins = np.linspace(-18, 18, 12 + 1) + np.mean(Y),
	histtype = 'stepfilled',
	lw = 0,
	color = ISO_COLORS['627'],
	orientation = 'horizontal',
	alpha = 0.25,
)

ax.hist(
	Y,
	bins = np.linspace(-18, 18, 12 + 1) + np.mean(Y),
	histtype = 'step',
	lw = 1,
	color = ISO_COLORS['627'],
	orientation = 'horizontal',
)

ppl.setp(ax.get_yticklabels(), visible = False)
ax.grid(False)
ax.set_xticks([])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ymin, ymax = ax.axis()[-2:]
yspan = ymax - ymin
yspan *= 1.2
y0 = (ymin + ymax) / 2

ax.axis([None, None, -yspan / 2, +yspan / 2])

ax = ax_d628ts
Y = [r['d628'] for r in didata]
ax.plot(t, Y, '-', color = ISO_COLORS['628'], lw = 1.5)
ppl.setp(ax.get_xticklabels(), visible = False)
ax.text(
	0.5,
	0.05,
	f'δ$_{{\\mathbf{{628}}}}$ = {np.mean(Y):.3f} ‰ (SD = {stdev(Y)*1000:.1f} ppm)',
	transform = ax.transAxes,
	ha = 'center',
	va = 'bottom',
	color = [_ * 0.75 for _ in ISO_COLORS['628']],
	weight = 'bold',
	size = 10,
)
ymin, ymax = ax.axis()[-2:]
y0 = (ymin + ymax) / 2 - 0.006
ax.axis([None, None, y0 - yspan / 2000, y0 + yspan / 2000])
ax.set_ylabel('δ$_{628}$ (‰)')

ax = ax_d628hist
ax.hist(
	Y,
	bins = np.linspace(-18e-3, 18e-3, 12 + 1) + np.mean(Y),
	histtype = 'stepfilled',
	lw = 0,
	color = ISO_COLORS['628'],
	orientation = 'horizontal',
	alpha = 0.25,
)

ax.hist(
	Y,
	bins = np.linspace(-18e-3, 18e-3, 12 + 1) + np.mean(Y),
	histtype = 'step',
	lw = 1,
	color = ISO_COLORS['628'],
	orientation = 'horizontal',
)

ppl.setp(ax.get_yticklabels(), visible = False)
ax.grid(False)
ax.set_xticks([])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ax = ax_d636ts
Y = [r['d636'] for r in didata]
ax.plot(t, Y, '-', color = ISO_COLORS['636'], lw = 1.5)
ax.set_xlabel('time (h)')
ax.text(
	0.5,
	0.05,
	f'δ$_{{\\mathbf{{636}}}}$ = {np.mean(Y):.3f} ‰ (SD = {stdev(Y)*1000:.1f} ppm)',
	transform = ax.transAxes,
	ha = 'center',
	va = 'bottom',
	color = [_ * 0.75 for _ in ISO_COLORS['636']],
	weight = 'bold',
	size = 10,
)
ax.set_ylabel('δ$_{636}$ (‰)')
ymin, ymax = ax.axis()[-2:]
y0 = (ymin + ymax) / 2 - 0.006
ax.axis([None, None, y0 - yspan / 2000, y0 + yspan / 2000])

ax = ax_d636hist
ax.hist(
	Y,
	bins = np.linspace(-18e-3, 18e-3, 12 + 1) + np.mean(Y),
	histtype = 'stepfilled',
	lw = 0,
	color = ISO_COLORS['636'],
	orientation = 'horizontal',
	alpha = 0.25,
)
ax.hist(
	Y,
	bins = np.linspace(-18e-3, 18e-3, 12 + 1) + np.mean(Y),
	histtype = 'step',
	lw = 1,
	color = ISO_COLORS['636'],
	orientation = 'horizontal',
)

ppl.setp(ax.get_yticklabels(), visible = False)
ax.grid(False)
ax.set_xticks([])
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ax = ax_avar

Y = [r['d636'] * 1000 for r in didata]
(tau, ad, ade, ns) = adev(Y, data_type = 'freq', taus = 'all')
ax.plot(tau, ad, '-', color = ISO_COLORS['636'], lw = 1, label = 'δ$_{636}$')

Y = [r['d628'] * 1000 for r in didata]
(tau, ad, ade, ns) = adev(Y, data_type = 'freq', taus = 'all')
ax.plot(tau, ad, '-', color = ISO_COLORS['628'], lw = 1, label = 'δ$_{628}$')

Y = [r['D627'] for r in didata]
(tau, ad, ade, ns) = adev(Y, data_type = 'freq', taus = 'all')
ax.plot(tau, ad, '-', color = ISO_COLORS['627'], lw = 1.5, label = 'Δ’$^{17}$O$_{WG}$')

ax.plot(
	[tau[0], tau[-1]],
	[0.95 * ad[0], 0.95 * ad[0] * (tau[0] / tau[-1]) ** 0.5],
	'k-',
	lw = 1,
	dashes = (6, 3, 2, 3),
	label = 'white noise',
)

ax.set_xscale('log')
ax.set_yscale('log')

ax.legend(
	loc = 'lower left',
	labelspacing = 0.3,
)
ax.set_xlabel('Number of aliquots')
ax.set_ylabel('Allan deviation (ppm)', labelpad = 10)

ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.0f}'))
ax.xaxis.set_minor_formatter(ticker.StrMethodFormatter(''))
ax.yaxis.set_minor_formatter(ticker.StrMethodFormatter(''))
ax.set_xticks([1, 2, 4, 8, 16])
ax.set_yticks([1, 2, 3, 4, 5, 6])
ax.set_ylim(0.4, None)

fig.savefig(
	'output/allan-variance',
	metadata = {
		'Author': AUTHOR,
		'Title': """
			Allan plot (left) and Δ’17O, δ628, δ636, time series corresponding
			to the working-gas measurements of figure 5. Allan deviation is
			computed with no overlap. Analytical scatter of repeated aliquots
			behaves as expected for white noise.
			""".replace('\t', '').replace('\n', ' ')[1:-1],
	},
)
ppl.close(fig)
