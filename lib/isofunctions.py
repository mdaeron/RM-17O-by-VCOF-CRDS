from scipy.optimize import fsolve
from loguru import logger
from scipy.stats import gaussian_kde, chi2
from scipy.linalg import eigh
from pylab import *
from sklearn.covariance import MinCovDet
from matplotlib.patches import Ellipse

R13_VPDB = 0.01118
R18_VSMOW = 0.00200520
R17_VSMOW = 0.00038475
LAMBDA_17 = 0.528

CO2g_H2O_FRACTIONATION = 'Guo & Zhou (2019)'
# CO2g_H2O_FRACTIONATION = 'Barkan & Luz (2012)'

CO2g_H2O_fractionations = {
	'Guo & Zhou (2019)': dict(
		alpha18_CO2g_H2O = 1.041461382,
		theta17_CO2g_H2O = 0.524606996,
	),
	'Barkan & Luz (2012)': dict(
		alpha18_CO2g_H2O = 1.041036,
		theta17_CO2g_H2O = 0.5229,
	),
}
_CO2g_H2O_fractionation = CO2g_H2O_fractionations[CO2g_H2O_FRACTIONATION]

alpha18_CO2g_H2O = _CO2g_H2O_fractionation['alpha18_CO2g_H2O']
theta17_CO2g_H2O = _CO2g_H2O_fractionation['theta17_CO2g_H2O']
alpha17_CO2g_H2O = exp(theta17_CO2g_H2O * log(alpha18_CO2g_H2O))

REFMATERIALS = {
	'NBS18': {'d13C_VPDB': -5.014, 'd18O_VPDB': -23.01},
	'NBS19': {'d13C_VPDB': 1.95, 'd18O_VPDB': -2.2},
	'IAEA603': {'d13C_VPDB': 2.46, 'd18O_VPDB': -2.37},
	'IAEA610': {'d13C_VPDB': -9.109, 'd18O_VPDB': -18.834},
	'IAEA611': {'d13C_VPDB': -30.795, 'd18O_VPDB': -4.224},
	'IAEA612': {'d13C_VPDB': -36.722, 'd18O_VPDB': -12.079},
}


def estimate_covariance(X, robust = True):
	# X.shape must be (N,2)
	if robust:
		return MinCovDet().fit(X).covariance_
	else:
		return cov(X.T)


def cov_ellipse(CM, p = None, r = None):
	"""
	Parameters
	----------
	CM : (2, 2) array
		Covariance matrix.
	p : float
		Confidence level, should be in (0, 1)

	Returns
	-------
	width, height, rotation :
		 The lengths of two axes and the rotation angle in degree
	for the ellipse.
	"""

	assert r is None or p is None

	if r is None and p is None:
		p = 0.95
		r2 = chi2.ppf(p, 2)
	if r is not None and p is None:
		r2 = r**2
	if p is not None and r is None:
		r2 = chi2.ppf(p, 2)

	val, vec = eigh(CM)
	width, height = 2 * sqrt(val[:, None] * r2)
	rotation = degrees(arctan2(*vec[::-1, 0]))

	return width, height, rotation


def mix_waters(
	d18a = 0.0,
	d18b = -55.5,
	D17a = 0.0,
	D17b = 0.0,
	x = 0.5,
	R18_VSMOW = R18_VSMOW,
	R17_VSMOW = R17_VSMOW,
	LAMBDA_17=LAMBDA_17,
):
	R18a = R18_VSMOW * (1 + d18a / 1000)
	R17a = exp(D17a / 1000) * R17_VSMOW * (1 + d18a / 1000) ** LAMBDA_17
	R18b = R18_VSMOW * (1 + d18b / 1000)
	R17b = exp(D17b / 1000) * R17_VSMOW * (1 + d18b / 1000) ** LAMBDA_17

	C16a = 1 / (1 + R17a + R18a)
	C17a = R17a * C16a
	C18a = R18a * C16a
	C16b = 1 / (1 + R17b + R18b)
	C17b = R17b * C16b
	C18b = R18b * C16b

	C16 = x * C16a + (1 - x) * C16b
	C17 = x * C17a + (1 - x) * C17b
	C18 = x * C18a + (1 - x) * C18b

	R17 = C17 / C16
	R18 = C18 / C16

	d18 = 1000 * (R18 / R18_VSMOW - 1)
	D17 = 1000 * (log(R17 / R17_VSMOW) - LAMBDA_17 * log(R18 / R18_VSMOW))

	txt = f'  {x:6} [δ18O = {d18a:6}, Δ17O = {D17a:6}]\n+ {1-x:6} [δ18O = {d18b:6}, Δ17O = {D17b:6}]\n = {" ":8}[δ18O = {d18:6.2f}, Δ17O = {D17:6.6f}]'

	return d18, D17, txt


def D17_CO2_H2O_eq_approx(
	h2o_to_co2_molecular_ratio = 1.0,
	d18i_co2=0.0,
	D17i_co2=0.0,
	d18i_h2o = 0.0,
	D17i_h2o = 0.0,
	alpha18_CO2g_H2O = alpha18_CO2g_H2O,
	alpha17_CO2g_H2O = alpha17_CO2g_H2O,
	R18_VSMOW = R18_VSMOW,
	R17_VSMOW = R17_VSMOW,
	LAMBDA_17=LAMBDA_17,
	xtol = 1e-12,
):
	x_co2 = 1 / (1 + h2o_to_co2_molecular_ratio)
	x_h2o = 1 - x_co2

	r18i_co2 = 2 * R18_VSMOW * (1 + d18i_co2 / 1e3)  # initial [626]/[626]
	r17i_co2 = (
		2 * exp(D17i_co2 / 1e3) * R17_VSMOW * (1 + d18i_co2 / 1e3) ** LAMBDA_17
	)  # initial [627]/[626]
	x16i_co2 = x_co2 / (1 + r18i_co2 + r17i_co2)  # initial [626]/([626]+[627]+[628])
	x17i_co2 = (
		x_co2 * r17i_co2 / (1 + r18i_co2 + r17i_co2)
	)  # initial [627]/([626]+[627]+[628])
	x18i_co2 = (
		x_co2 * r18i_co2 / (1 + r18i_co2 + r17i_co2)
	)  # initial [628]/([626]+[627]+[628])

	r18i_h2o = R18_VSMOW * (1 + d18i_h2o / 1e3)  # initial [181]/[161]
	r17i_h2o = (
		exp(D17i_h2o / 1e3) * R17_VSMOW * (1 + d18i_h2o / 1e3) ** LAMBDA_17
	)  # initial [171]/[161]
	x16i_h2o = x_h2o / (1 + r18i_h2o + r17i_h2o)  # initial [161]/([161]+[171]+[181])
	x17i_h2o = (
		x_h2o * r17i_h2o / (1 + r18i_h2o + r17i_h2o)
	)  # initial [171]/([161]+[171]+[181])
	x18i_h2o = (
		x_h2o * r18i_h2o / (1 + r18i_h2o + r17i_h2o)
	)  # initial [181]/([161]+[171]+[181])

	def equations(p):
		x16_co2, x17_co2, x18_co2, x16_h2o, x17_h2o, x18_h2o = exp(p)
		return (
			x16_co2 + x16_h2o - x16i_co2 - x16i_h2o,  # conservation of N(626) + N(161)
			x17_co2 + x17_h2o - x17i_co2 - x17i_h2o,  # conservation of N(627) + N(171)
			x18_co2 + x18_h2o - x18i_co2 - x18i_h2o,  # conservation of N(628) + N(181)
			alpha18_CO2g_H2O * x18_h2o / x16_h2o
			- x18_co2 / 2 / x16_co2,  # alpha18 = ([628]/[626]/2) / ([181]/[616])
			alpha17_CO2g_H2O * x17_h2o / x16_h2o
			- x17_co2 / 2 / x16_co2,  # alpha17 = ([627]/[626]/2) / ([171]/[616])
			x16_co2 + x17_co2 + x18_co2 - x_co2,  # conservation of CO2
		)

	x16_co2, x17_co2, x18_co2, x16_h2o, x17_h2o, x18_h2o = exp(
		fsolve(
			equations,
			log(array([x16i_co2, x17i_co2, x18i_co2, x16i_h2o, x17i_h2o, x18i_h2o])),
			xtol = 1e-12,
		)
	)

	d18e_h2o = 1e3 * ((x18_h2o / x16_h2o / R18_VSMOW) - 1)
	D17e_h2o = 1e3 * (
		log(x17_h2o / x16_h2o / R17_VSMOW)
		- LAMBDA_17 * log(x18_h2o / x16_h2o / R18_VSMOW)
	)
	d18e_co2 = 1e3 * ((x18_co2 / x16_co2 / 2 / R18_VSMOW) - 1)
	D17e_co2 = 1e3 * (
		log(x17_co2 / x16_co2 / 2 / R17_VSMOW)
		- LAMBDA_17 * log(x18_co2 / x16_co2 / 2 / R18_VSMOW)
	)

	return d18e_co2, D17e_co2, d18e_h2o, D17e_h2o


def D17_CO2_H2O_eq(
	h2o_to_co2_molecular_ratio = 1.0,
	d18i_co2=0.0,
	D17i_co2=0.0,
	d18i_h2o = 0.0,
	D17i_h2o = 0.0,
	alpha18_CO2g_H2O = alpha18_CO2g_H2O,
	alpha17_CO2g_H2O = alpha17_CO2g_H2O,
	R18_VSMOW = R18_VSMOW,
	R17_VSMOW = R17_VSMOW,
	LAMBDA_17=LAMBDA_17,
	D727i = 0.1043,
	D728i = 0.1951,
	D828i = 0.3743,
	D637i = 0.5101,
	D638i = 0.9714,
	D737i = 1.1188,
	D738i = 1.6662,
	D838i = 2.3046,
	D727=0.1043,
	D728=0.1951,
	D828=0.3743,
	D637=0.5101,
	D638=0.9714,
	D737=1.1188,
	D738=1.6662,
	D838=2.3046,
	stochastic = False,
	d13_co2=0.0,
	R13_VPDB = R13_VPDB,
	xtol = 1e-12,
):
	if stochastic:
		rK637, rK637i = 1, 1
		rK638, rK638i = 1, 1
		rK727, rK727i = 1, 1
		rK728, rK728i = 1, 1
		rK737, rK737i = 1, 1
		rK738, rK738i = 1, 1
		rK828, rK828i = 1, 1
		rK838, rK838i = 1, 1
	else:
		rK637, rK637i = 1 + D637 / 1e3, 1 + D637i / 1e3
		rK638, rK638i = 1 + D638 / 1e3, 1 + D638i / 1e3
		rK727, rK727i = 1 + D727 / 1e3, 1 + D727i / 1e3
		rK728, rK728i = 1 + D728 / 1e3, 1 + D728i / 1e3
		rK737, rK737i = 1 + D737 / 1e3, 1 + D737i / 1e3
		rK738, rK738i = 1 + D738 / 1e3, 1 + D738i / 1e3
		rK828, rK828i = 1 + D828 / 1e3, 1 + D828i / 1e3
		rK838, rK838i = 1 + D838 / 1e3, 1 + D838i / 1e3

	x_co2 = 1 / (1 + h2o_to_co2_molecular_ratio)
	x_h2o = 1 - x_co2

	r13i_co2 = R13_VPDB * (1 + d13_co2 / 1e3)  # initial [13]/[12]
	r18i_co2 = R18_VSMOW * (1 + d18i_co2 / 1e3)  # initial [18]/[16]
	r17i_co2 = (
		exp(D17i_co2 / 1e3) * R17_VSMOW * (1 + d18i_co2 / 1e3) ** LAMBDA_17
	)  # initial [17]/[16]

	r628i_co2 = 2 * r18i_co2  # initial [628]/[626]
	r627i_co2 = 2 * r17i_co2  # initial [627]/[626]
	r828i_co2 = rK828i * r18i_co2**2  # initial [828]/[626]
	r727i_co2 = rK727i * r17i_co2**2  # initial [727]/[626]
	r728i_co2 = rK728i * 2 * r17i_co2 * r18i_co2  # initial [728]/[626]
	r636i_co2 = r13i_co2  # initial [628]/[626]
	r638i_co2 = rK638i * 2 * r18i_co2 * r13i_co2  # initial [628]/[626]
	r637i_co2 = rK637i * 2 * r17i_co2 * r13i_co2  # initial [627]/[626]
	r838i_co2 = rK838i * r18i_co2**2 * r13i_co2  # initial [828]/[626]
	r737i_co2 = rK737i * r17i_co2**2 * r13i_co2  # initial [727]/[626]
	r738i_co2 = (
		rK738i * 2 * r17i_co2 * r18i_co2 * r13i_co2
	)  # initial [727]/[626]                                  # initial [728]/[626]

	x626i_co2 = x_co2 / (
		1
		+ r628i_co2
		+ r627i_co2
		+ r828i_co2
		+ r727i_co2
		+ r728i_co2
		+ r636i_co2
		+ r638i_co2
		+ r637i_co2
		+ r838i_co2
		+ r737i_co2
		+ r738i_co2
	)  # initial fraction of [626] over all molecules
	x628i_co2 = x626i_co2 * r628i_co2  # initial fraction of [628] over all molecules
	x627i_co2 = x626i_co2 * r627i_co2  # initial fraction of [627] over all molecules
	x828i_co2 = x626i_co2 * r828i_co2  # initial fraction of [828] over all molecules
	x727i_co2 = x626i_co2 * r727i_co2  # initial fraction of [727] over all molecules
	x728i_co2 = x626i_co2 * r728i_co2  # initial fraction of [728] over all molecules
	x636i_co2 = x626i_co2 * r636i_co2  # initial fraction of [628] over all molecules
	x638i_co2 = x626i_co2 * r638i_co2  # initial fraction of [628] over all molecules
	x637i_co2 = x626i_co2 * r637i_co2  # initial fraction of [627] over all molecules
	x838i_co2 = x626i_co2 * r838i_co2  # initial fraction of [828] over all molecules
	x737i_co2 = x626i_co2 * r737i_co2  # initial fraction of [727] over all molecules
	x738i_co2 = x626i_co2 * r738i_co2  # initial fraction of [728] over all molecules

	r18i_h2o = R18_VSMOW * (1 + d18i_h2o / 1e3)  # initial [181]/[161]
	r17i_h2o = (
		exp(D17i_h2o / 1e3) * R17_VSMOW * (1 + d18i_h2o / 1e3) ** LAMBDA_17
	)  # initial [171]/[161]
	x16i_h2o = x_h2o / (
		1 + r18i_h2o + r17i_h2o
	)  # initial fraction of [161] over all molecules
	x17i_h2o = x16i_h2o * r17i_h2o
	x18i_h2o = x16i_h2o * r18i_h2o

	def equations(p):
		(
			x626_co2,
			x627_co2,
			x628_co2,
			x727_co2,
			x828_co2,
			x728_co2,
			x636_co2,
			x637_co2,
			x638_co2,
			x737_co2,
			x838_co2,
			x738_co2,
			x16_h2o,
			x17_h2o,
			x18_h2o,
		) = exp(p)
		return (
			# conservation of 13C
			(
				x636_co2
				+ x637_co2
				+ x638_co2
				+ x737_co2
				+ x738_co2
				+ x838_co2
				- x636i_co2
				- x637i_co2
				- x638i_co2
				- x737i_co2
				- x738i_co2
				- x838i_co2
			),
			# conservation of 16O
			(
				x16_h2o
				+ 2 * x626_co2
				+ x627_co2
				+ x628_co2
				+ 2 * x636_co2
				+ x637_co2
				+ x638_co2
				- x16i_h2o
				- 2 * x626i_co2
				- x627i_co2
				- x628i_co2
				- 2 * x636i_co2
				- x637i_co2
				- x638i_co2
			),
			# conservation of 17O
			(
				x17_h2o
				+ 2 * x727_co2
				+ x627_co2
				+ x728_co2
				+ 2 * x737_co2
				+ x637_co2
				+ x738_co2
				- x17i_h2o
				- 2 * x727i_co2
				- x627i_co2
				- x728i_co2
				- 2 * x737i_co2
				- x637i_co2
				- x738i_co2
			),
			# conservation of 18O
			(
				x18_h2o
				+ 2 * x828_co2
				+ x628_co2
				+ x728_co2
				+ 2 * x838_co2
				+ x638_co2
				+ x738_co2
				- x18i_h2o
				- 2 * x828i_co2
				- x628i_co2
				- x728i_co2
				- 2 * x838i_co2
				- x638i_co2
				- x738i_co2
			),
			# alpha18 = ([628]/[626]/2) / ([181]/[616])
			alpha18_CO2g_H2O * x18_h2o / x16_h2o - x628_co2 / 2 / x626_co2,
			# alpha17 = ([627]/[626]/2) / ([171]/[616])
			alpha17_CO2g_H2O * x17_h2o / x16_h2o - x627_co2 / 2 / x626_co2,
			# conservation of H2O
			x16_h2o + x17_h2o + x18_h2o - x_h2o,
			# clumped isotopes
			x828_co2 / x626_co2 - rK828 * (x628_co2 / x626_co2 / 2) ** 2,
			x727_co2 / x626_co2 - rK727 * (x627_co2 / x626_co2 / 2) ** 2,
			x728_co2 / x626_co2 - rK728 * x627_co2 / x626_co2 * x628_co2 / x626_co2 / 2,
			x638_co2 / x626_co2 - rK638 * (x628_co2 / x626_co2) * x636_co2 / x626_co2,
			x637_co2 / x626_co2 - rK637 * (x627_co2 / x626_co2) * x636_co2 / x626_co2,
			x838_co2 / x626_co2
			- rK838 * (x628_co2 / x626_co2 / 2) ** 2 * x636_co2 / x626_co2,
			x737_co2 / x626_co2
			- rK737 * (x627_co2 / x626_co2 / 2) ** 2 * x636_co2 / x626_co2,
			x738_co2 / x626_co2
			- rK738
			* x636_co2
			/ x626_co2
			* x627_co2
			/ x626_co2
			* x628_co2
			/ x626_co2
			/ 2,
		)

	(
		x626_co2,
		x627_co2,
		x628_co2,
		x727_co2,
		x828_co2,
		x728_co2,
		x636_co2,
		x637_co2,
		x638_co2,
		x737_co2,
		x838_co2,
		x738_co2,
		x16_h2o,
		x17_h2o,
		x18_h2o,
	) = exp(
		fsolve(
			equations,
			log(
				array(
					[
						x626i_co2,
						x627i_co2,
						x628i_co2,
						x727i_co2,
						x828i_co2,
						x728i_co2,
						x636i_co2,
						x637i_co2,
						x638i_co2,
						x737i_co2,
						x838i_co2,
						x738i_co2,
						x16i_h2o,
						x17i_h2o,
						x18i_h2o,
					]
				)
			),
			xtol = xtol,
		)
	)

	d18e_h2o = 1e3 * ((x18_h2o / x16_h2o / R18_VSMOW) - 1)
	D17e_h2o = 1e3 * (
		log(x17_h2o / x16_h2o / R17_VSMOW)
		- LAMBDA_17 * log(x18_h2o / x16_h2o / R18_VSMOW)
	)
	d18e_co2 = 1e3 * ((x628_co2 / x626_co2 / 2 / R18_VSMOW) - 1)
	D17e_co2 = 1e3 * (
		log(x627_co2 / x626_co2 / 2 / R17_VSMOW)
		- LAMBDA_17 * log(x628_co2 / x626_co2 / 2 / R18_VSMOW)
	)

	return d18e_co2, D17e_co2, d18e_h2o, D17e_h2o


def H2O_eq_CO2_figure(
	starting_CO2,
	eq_waters,
	N_starting_CO2,
	N_eq_waters,
	H2O_to_CO2_ratio_half_width = 0.1,
	filename = 'output/H2O_eq_CO2_D17_vs_ratio',
	metadata = {},
):
	"""GENERATE FIGURE H2O_eq_CO2_D17_vs_ratio"""

	fig = figure(figsize = (3.15, 2.7))
	subplots_adjust(0.27, 0.17, 0.98, 0.98)

	_ = 3
	h2o_to_co2_molar_ratio = geomspace(1e-4 / _, 1e4 * _, 500)

	dashes = [
		(None, None),
		(8, 1),
		(5, 1),
		(8, 1, 2, 1),
		(6, 2, 2, 2, 2, 2),
		(2, 1, 2, 1),
	]

	for k, w in enumerate(sorted(eq_waters, key = lambda _: -eq_waters[_]['D17O_VSMOW'])):
		y = [
			D17_CO2_H2O_eq(
				x,
				starting_CO2['d18O_VSMOW'],
				starting_CO2['D17O_VSMOW'],
				eq_waters[w]['d18O_VSMOW'],
				eq_waters[w]['D17O_VSMOW'],
			)[-1]
			for x in h2o_to_co2_molar_ratio
		]

		semilogx(
			h2o_to_co2_molar_ratio,
			y,
			ls = '-',
			lw = 1.0,
			color = 'k',
			dashes = dashes[k],
			zorder = -100,
			solid_capstyle = 'butt',
			label = w.replace('_', '-'),
		)

	_mycolor = (1, 0, 0)
	_alpha = 0.5
	axvspan(
		N_eq_waters / N_starting_CO2 * (1 - H2O_to_CO2_ratio_half_width),
		N_eq_waters / N_starting_CO2 * (1 + H2O_to_CO2_ratio_half_width),
		color = _mycolor,
		lw = 0,
		alpha = _alpha,
		zorder = -101,
	)
	plot(
		[],
		[],
		'-',
		lw = 2,
		color = _mycolor,
		alpha = _alpha,
		label = 'H$_2$O/CO$_2$ (this study)',
	)

	legend(
		framealpha = 0,
		fontsize = 8,
		handlelength = 2.6,
		loc = 'lower left',
		bbox_to_anchor = (0, 0),
		borderpad = 0.1,
		labelspacing = 0.6,
	)
	y1 = axis()[2]
	axis([h2o_to_co2_molar_ratio[0], h2o_to_co2_molar_ratio[-1], y1 - 0.025, None])
	xlabel('Molar ratio of H$_2$O / CO$_2$')
	ylabel(
		'Δ’$^{17}O_{VSMOW}$ (‰) of $CO_2$ equilibrated\n at 25$\\,$°C with different waters'
	)

	minorticks_off()
	xticks([1e-4, 1e-2, 1, 1e2, 1e4])

	fig.savefig(filename, metadata = metadata)
	close(fig)


def grubbs_two_sided_test(data, alpha):
	from scipy.stats import t

	X = asarray(data)
	mask = array([True for _ in X])

	while True:
		Xmean = X[mask].mean()
		Xsd = X[mask].std(ddof = 1)
		N = X[mask].size
		Z = abs(X - Xmean) / Xsd
		Z[~mask] = 0
		G_calculated = Z.max()
		t_critical = t.isf(alpha / 2 / (2 * N), N - 2)
		G_critical = ((N - 1) / sqrt(N)) * sqrt(t_critical**2 / (N - 2 + t_critical**2))

		if G_calculated > G_critical:
			mask[Z.argmax()] = False
		else:
			break

	return [k for k, v in enumerate(mask) if not v]


def d626b_outliers(data, p_grubbs = 0.95, rename = True, flag = True):
	out = [_.copy() for _ in data]
	# 	outliers = grubbs.two_sided_test_indices([r['d626b'] for r in out], alpha = 1-p_grubbs)
	outliers = grubbs_two_sided_test([r['d626b'] for r in out], alpha = 1 - p_grubbs)
	for k, r in enumerate(out):
		if flag:
			r['outlier'] = k in outliers
		if r['outlier']:
			logger.info(
				f"Analysis {r['UID']} ({r['Sample']}) is a d626b outlier (d626b = {r['d626b']:.1f} ‰)"
			)
			if rename:
				r['Sample'] = r['Sample'] + f'_outlier_{k}'
	return out, p_grubbs


def d626b_plot(S, p_grubbs = 0.95, filename = 'output/d626b_outliers'):
	fig = figure(figsize = (3.15, 2.5))
	subplots_adjust(0.2, 0.2, 0.97, 0.8)

	outliers = [_['outlier'] for _ in S['data']]

	X = array([r['d626b'] for r in S['data']])
	Y = array([r['D17residual'] for r in S['data']]) * 1e3
	kw = dict(mfc = 'w', ms = 5, mew = 1)
	(_all,) = plot(X, Y, 'ko', **kw)

	(_out,) = plot(
		X[outliers],
		Y[outliers],
		'x',
		mew = 1,
		ms = 8,
		mfc = 'none',
		mec = 'k',
	)

	if True in outliers:
		legend(
			[(_all, _out)],
			[f'Grubbs outliers (p = {p_grubbs})\nbased only on δ$_{{626b}}$ values'],
			fontsize = 9,
			handlelength = 1.5,
			loc = 'lower center',
			bbox_to_anchor = (0.5, 1.01),
		)

	axhline(0, color = 'k', lw = 0.75, dashes = (4, 2))
	axvline(0, color = 'k', lw = 0.75, dashes = (4, 2))
	xlabel('δ$_{626b}$ (‰)')
	ylabel('Δ’$^{17}$O residuals (ppm)')
	margins(0.2)

	savefig(filename)
	close(fig)


def plot_residuals(
	data,
	key = 'D17O_VSMOW',
	filename = 'residuals',
	anchors = ['SMOW', 'SLAP'],
	anchorcolor = (1, 0, 0),
	uid_labels = True,
	connecting_lines = False,
):
	fig = figure(figsize = (10, 6))
	subplots_adjust(0.1, 0.1, 0.95, 0.95, 0.05, 0.05)
	moz = fig.subplot_mosaic('AAAAB;CCCC.')

	ax1, ax2, ax3 = moz['C'], moz['A'], moz['B']

	UID = [r['UID'] for r in data]
	S = [r['Sample'] for r in data]
	newsession = [
		k for k in range(len(data) - 1) if data[k]['Session'] != data[k + 1]['Session']
	]
	sessions = []
	for r in data:
		if r['Session'] not in sessions:
			sessions.append(r['Session'])
	samples = sorted(set(S))[::-1]
	Y_sample = {s: k + 1 for k, s in enumerate(samples)}
	sample_colors = {s: anchorcolor if s in anchors else 'k' for s in samples}
	t = array([k + 1 for k, r in enumerate(data)])
	try:
		Y = [r[key + '_residual'] * 1e3 for r in data]
	except KeyError:
		Y = [r[key[:3] + 'residual'] * 1e3 for r in data]
	rmse = (array(Y) ** 2).mean() ** 0.5
	sd = ((array(Y) ** 2).sum() / (len(Y) - len(samples))) ** 0.5

	# 	ax2.axhline(0, color = 'k', lw = 0.5)
	for k in newsession:
		ax1.axvline(k + 1.5, color = 'k', lw = 0.5)
		ax2.axvline(k + 1.5, color = 'k', lw = 0.5)

	for k in range(len(t)):
		_x, _y1, _y2 = t[k], Y_sample[S[k]], Y[k]
		ax1.plot(_x, _y1, 'o', mew = 0, ms = 5, mfc = sample_colors[S[k]])
		ax2.plot(_x, _y2, 'o', mew = 0, ms = 5, mfc = sample_colors[S[k]])
		if uid_labels:
			ax1.text(
				_x,
				_y1,
				'   ' + UID[k],
				rotation = 90,
				va = 'bottom',
				ha = 'center',
				size = 6,
				alpha = 0.5,
				color = sample_colors[S[k]],
			)
		if connecting_lines:
			ax2.add_artist(
				ConnectionPatch(
					xyA = (t[k], Y_sample[S[k]] - 0.3),
					xyB = (t[k], Y[k] + 1),
					coordsA = 'data',
					coordsB = 'data',
					axesA = ax1,
					axesB = ax2,
					color = 'k',
					alpha = 0.25,
					lw = 0.5,
				)
			)

	xmin, xmax = ax2.get_xlim()
	sessionlimts = [xmin] + [_ + 1.5 for _ in newsession] + [xmax]
	for k, s in enumerate(sessions):
		ax1.text(
			(sessionlimts[k] + sessionlimts[k + 1]) / 2,
			0,
			f'\n{s}',
			ha = 'center',
			va = 'center',
		)

	ymin, ymax = ax2.get_ylim()
	bins = linspace(ymin, ymax, 8)
	ax3.hist(
		Y, bins, orientation = 'horizontal', histtype = 'stepfilled', alpha = 0.25, color = 'k'
	)
	ax3.text(
		0, (ymin + ymax) / 2, f'  SD = {sd:.1f} ppm', va = 'center', ha = 'left', size = 10
	)
	ax3.set_ylim((ymin, ymax))
	ax3.axis('off')

	# 	ax1.spines[['right', 'left', 'top', 'bottom']].set_visible(False)
	ax1.set_ylim((None, len(samples) + 4))
	ax1.tick_params(length = 0, labelleft = False, labelright = True, labelsize = 8)
	ax1.set_xticks([])
	ax1.set_yticks([k + 1 for k in range(len(samples))])
	ax1.set_yticklabels(samples)

	ax2.set_xticks([])
	if key == 'D17O_VSMOW':
		ax2.set_ylabel('Δ’$^{17}Ο$ residuals (ppm)')
	if key == 'd18O_VSMOW':
		ax2.set_ylabel('δ$^{18}Ο$ residuals (ppm)')
	if key == 'd13C_VPDB':
		ax2.set_ylabel('δ$^{13}C$ residuals (ppm)')

	savefig(filename)
	# 	show()

	close(fig)


def save_sample_results(
	stdz_carbon,
	stdz_oxygen_carb,
	stdz_triple_oxygen,
	filename = 'output/samples.csv',
	d13key = 'd13C_VPDB',
	d18keycarb = 'd18O_VPDB',
	d18key = 'd18O_VSMOW',
	d17key = 'd17O_VSMOW',
	D17key = 'D17O_VSMOW',
):
	include_d13C = stdz_carbon is not None
	include_d18Ocarb = stdz_oxygen_carb is not None

	samples = stdz_triple_oxygen['samples']
	if include_d13C:
		samples_C = stdz_carbon['samples']
	if include_d18Ocarb:
		samples_O = stdz_oxygen_carb['samples']

	csv = 'Sample'
	csv += ',N'
	if include_d13C:
		csv += ',' + d13key
		csv += ',SD_' + d13key
		csv += ',SE_' + d13key
		csv += ',95CL_' + d13key
	if include_d18Ocarb:
		csv += ',' + d18keycarb
		csv += ',SD_' + d18keycarb
		csv += ',SE_' + d18keycarb
		csv += ',95CL_' + d18keycarb
	csv += ',' + d18key
	csv += ',SD_' + d18key
	csv += ',SE_' + d18key
	csv += ',95CL_' + d18key
	csv += ',' + d17key
	csv += ',SD_' + d17key
	csv += ',SE_' + d17key
	csv += ',95CL_' + d17key
	csv += ',' + D17key
	csv += ',SD_' + D17key
	csv += ',SE_' + D17key
	csv += ',95CL_' + D17key

	for s in samples:
		csv += f'\n{s}'
		csv += f',{samples[s]["N"]}'
		if include_d13C:
			csv += f',{samples_C[s][d13key]:.3f}' if s in samples_C else ',-'
			csv += f',{samples_C[s]["SD_"+d13key]:.3f}' if s in samples_C else ',-'
			csv += (
				f',{samples_C[s]["SE_"+d13key]:.3f}'
				if s in stdz_carbon['unknowns']
				else ',-'
			)
			csv += (
				f',{samples_C[s]["95CL_"+d13key]:.3f}'
				if s in stdz_carbon['unknowns']
				else ',-'
			)
		if include_d18Ocarb:
			if s in stdz_oxygen_carb['samples']:
				csv += f',{samples_O[s][d18keycarb]:.3f}'
				csv += f',{samples_O[s]["SD_"+d18keycarb]:.3f}'
				csv += (
					f',{samples_O[s]["SE_"+d18keycarb]:.3f}'
					if s in stdz_oxygen_carb['unknowns']
					else ',-'
				)
				csv += (
					f',{samples_O[s]["95CL_"+d18keycarb]:.3f}'
					if s in stdz_oxygen_carb['unknowns']
					else ',-'
				)
			else:
				csv += ',-,-,-,-'
		csv += f',{samples[s][d18key]:.4f}'
		csv += f',{samples[s]["SD_"+d18key]:.4f}'
		csv += (
			f',{samples[s]["SE_"+d18key]:.4f}'
			if s in stdz_triple_oxygen['unknowns']
			else ',-'
		)
		csv += (
			f',{samples[s]["95CL_"+d18key]:.4f}'
			if s in stdz_triple_oxygen['unknowns']
			else ',-'
		)
		csv += f',{samples[s][d17key]:.4f}'
		csv += f',{samples[s]["SD_"+d17key]:.4f}'
		csv += (
			f',{samples[s]["SE_"+d17key]:.4f}'
			if s in stdz_triple_oxygen['unknowns']
			else ',-'
		)
		csv += (
			f',{samples[s]["95CL_"+d17key]:.4f}'
			if s in stdz_triple_oxygen['unknowns']
			else ',-'
		)
		csv += f',{samples[s][D17key]:.4f}'
		csv += f',{samples[s]["SD_"+D17key]:.4f}'
		csv += (
			f',{samples[s]["SE_"+D17key]:.4f}'
			if s in stdz_triple_oxygen['unknowns']
			else ',-'
		)
		csv += (
			f',{samples[s]["95CL_"+D17key]:.4f}'
			if s in stdz_triple_oxygen['unknowns']
			else ',-'
		)

	with open(filename, 'w') as fid:
		fid.write(csv)


def plot_kde(S, filename = 'kde', res = 'D17', plot_margins = 0.3, robust_cov_estimator = True):
	data = S['data']

	x = array([_['d18residual'] for _ in data])
	y = array([_[res + 'residual'] for _ in data]) * 1e3

	xmin, xmax = x.min(), x.max()
	xmin, xmax = (
		xmin - (xmax - xmin) * plot_margins,
		xmax + (xmax - xmin) * plot_margins,
	)

	ymin, ymax = y.min(), y.max()
	ymin, ymax = (
		ymin - (ymax - ymin) * plot_margins,
		ymax + (ymax - ymin) * plot_margins,
	)

	xx, yy = mgrid[xmin:xmax:1000j, ymin:ymax:1000j]
	positions = vstack([xx.ravel(), yy.ravel()])
	values = vstack([x, y])
	kernel = gaussian_kde(values)

	f = np.reshape(kernel(positions).T, xx.shape)
	f /= f.max()
	f = 1 - f

	fig = figure(figsize = (5, 5))
	subplots_adjust(0.17, 0.17, 0.97, 0.97)

	ax = subplot(111)
	ax.set_xlim(xmin, xmax)
	ax.set_ylim(ymin, ymax)

	cbins = [-1, 0, 0.25, 0.5, 0.75, 0.95]
	cfset = ax.contourf(xx, yy, f, levels = cbins, cmap = 'Blues_r', alpha = 1)
	## Or kernel density estimate plot instead of the contourf plot
	# ax.imshow(np.rot90(f), cmap = 'Blues', extent = [xmin, xmax, ymin, ymax])
	# Contour plot
	cset = ax.contour(xx, yy, f, levels = cbins, colors = 'k')
	plot(x, y, 'k+', mew = 0.5, ms = 4, alpha = 0.5)

	ax.clabel(cset, inline = 1, fontsize = 10)
	ax.set_xlabel('δ$^{18}$O$_{VSMOW}$ (‰)')
	ax.set_ylabel('Δ’$^{17}$O$_{VSMOW}$ (ppm)')
	# 	ax.grid(False)

	CM = estimate_covariance(array([x, y]).T)
	w, h, r = cov_ellipse(CM)
	kw = dict(ec = 'k', lw = 1, ls = (0, (6, 2, 2, 2)))
	ax.add_patch(Ellipse(xy = (0, 0), width = w, height = h, angle = r, fc = 'None', **kw))
	kw = dict(color = kw['ec'], lw = kw['lw'], ls = kw['ls'])
	plot([], [], label = '95 % confidence from\ncovariance estimate', **kw)

	legend(fontsize = 7, handlelength = 2.5)

	fig.savefig(filename.split('.')[0] + '_' + res + 'O')
	close(fig)


def correct_CO2eqH2O_back_to_initial_H2O(co2_samples, d18Oco2_i, D17Oco2_i):
	"""
	co2_samples = [
		{'Sample': ..., 'd18O_VSMOW' = ..., 'D17O_VSMOW' = ..., 'H2O/CO2' = ...},
		{'Sample': ..., 'd18O_VSMOW' = ..., 'D17O_VSMOW' = ..., 'H2O/CO2' = ...}, ...
		]
	"""

	from scipy.optimize import fsolve

	D17O_change_from_CO2eqH2O = 1e3 * (
		log(alpha17_CO2g_H2O) - LAMBDA_17 * log(alpha18_CO2g_H2O)
	)

	def fun(x, d18Oco2, D17Oco2, d18Oco2_i, D17Oco2_i, wcratio):
		_d18Oco2, _D17Oco2 = D17_CO2_H2O_eq(
			wcratio,
			d18Oco2_i,
			D17Oco2_i,
			x[0],
			x[1],
		)[:2]
		return array([d18Oco2 - _d18Oco2, D17Oco2 - _D17Oco2])

	for s in co2_samples:
		d18Ow_firstguess = (1000 + s['d18O_VSMOW']) / alpha18_CO2g_H2O - 1000
		D17Ow_firstguess = s['D17O_VSMOW'] - D17O_change_from_CO2eqH2O

		d18Ow, D17Ow = fsolve(
			fun,
			array([d18Ow_firstguess, D17Ow_firstguess]),
			args = (
				s['d18O_VSMOW'],
				s['D17O_VSMOW'],
				d18Oco2_i,
				D17Oco2_i,
				s['H2O/CO2'],
			),
			xtol = 1e-12,
		)
		s['water_d18O_VSMOW'] = d18Ow
		s['water_D17O_VSMOW'] = D17Ow


if __name__ == '__main__':
	from pylab import *

	x = logspace(-4, 4, 201)
	A = array(
		[
			D17_CO2_H2O_eq_approx(
				h2o_to_co2_molecular_ratio = _,
				d18i_co2=0.0,
				D17i_co2=0.0,
			)[1]
			for _ in x
		]
	)
	B = array(
		[
			D17_CO2_H2O_eq(
				h2o_to_co2_molecular_ratio = _,
				d18i_co2=0.0,
				D17i_co2=0.0,
				stochastic = True,
			)[1]
			for _ in x
		]
	)
	C = array(
		[
			D17_CO2_H2O_eq(
				h2o_to_co2_molecular_ratio = _,
				d18i_co2=0.0,
				D17i_co2=0.0,
				stochastic = False,
			)[1]
			for _ in x
		]
	)

	ax = subplot(211)
	semilogx(x, A, 'b')
	semilogx(x, B, 'r')

	subplot(212, sharex = ax)
	semilogx(x, B - A, 'g')
	semilogx(x, C - B, 'c')

	show()
