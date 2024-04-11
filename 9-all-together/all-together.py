import sys

sys.path.append('../lib')
from mylib import *

logger.remove(0)
logger.add('output/all-together.log', mode = 'w', format = '[{line:3.0f}]    {message}')

ppl.style.use('../lib/mydefault.mplstyle')

if __name__ == '__main__':
	"""H2O/CO2 RATIO"""
	N_starting_CO2 = 40e-6  # mol
	V_eq_waters = 300e-6  # L
	N_eq_waters = V_eq_waters / 18.015e-3  # mol
	H2O_to_CO2_ratio = N_eq_waters / N_starting_CO2

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

	d18_OHmix, D17_OHmix = zip(
		*[
			isofunctions.mix_waters(
				eq_waters['OC4']['d18O_VSMOW'],
				eq_waters['HAWAI']['d18O_VSMOW'],
				eq_waters['OC4']['D17O_VSMOW'],
				eq_waters['HAWAI']['D17O_VSMOW'],
				x,
			)[:2]
			for x in np.linspace(0, 1)
		]
	)

	d18_NHmix, D17_NHmix = zip(
		*[
			isofunctions.mix_waters(
				eq_waters['NEEM']['d18O_VSMOW'],
				eq_waters['HAWAI']['d18O_VSMOW'],
				eq_waters['NEEM']['D17O_VSMOW'],
				eq_waters['HAWAI']['D17O_VSMOW'],
				x,
			)[:2]
			for x in np.linspace(0, 1)
		]
	)

	d18_ONHmix, D17_ONHmix = zip(
		*[
			isofunctions.mix_waters(
				eq_waters['MIX_OH']['d18O_VSMOW'],
				eq_waters['MIX_NH']['d18O_VSMOW'],
				eq_waters['MIX_OH']['D17O_VSMOW'],
				eq_waters['MIX_NH']['D17O_VSMOW'],
				x,
			)[:2]
			for x in np.linspace(0, 1)
		]
	)

	"""COMPOSITION CO2 USED WITH EQUILIBRATED WATERS"""
	match CO2g_H2O_FRACTIONATION:
		case 'Guo & Zhou (2019)':
			starting_CO2 = {'d18O_VSMOW': 24.71, 'D17O_VSMOW': -0.0835}
		case 'Barkan & Luz (2012)':
			starting_CO2 = {'d18O_VSMOW': 24.25, 'D17O_VSMOW': -0.150}

	d18_OHmix, D17_OHmix = zip(
		*[
			isofunctions.D17_CO2_H2O_eq(
				N_eq_waters / N_starting_CO2,
				starting_CO2['d18O_VSMOW'],
				starting_CO2['D17O_VSMOW'],
				d18,
				D17,
			)[:2]
			for d18, D17 in zip(d18_OHmix, D17_OHmix)
		]
	)

	d18_NHmix, D17_NHmix = zip(
		*[
			isofunctions.D17_CO2_H2O_eq(
				N_eq_waters / N_starting_CO2,
				starting_CO2['d18O_VSMOW'],
				starting_CO2['D17O_VSMOW'],
				d18,
				D17,
			)[:2]
			for d18, D17 in zip(d18_NHmix, D17_NHmix)
		]
	)

	d18_ONHmix, D17_ONHmix = zip(
		*[
			isofunctions.D17_CO2_H2O_eq(
				N_eq_waters / N_starting_CO2,
				starting_CO2['d18O_VSMOW'],
				starting_CO2['D17O_VSMOW'],
				d18,
				D17,
			)[:2]
			for d18, D17 in zip(d18_ONHmix, D17_ONHmix)
		]
	)

	"""IMPORT MIXWATER DATA"""
	with open('../4-linearity/output/samples.csv') as fid:
		mixwaterdata = [l.strip().split(',') for l in fid.readlines()]
	mixwaterdata = [
		{
			k: float(v) if k != 'Sample' else v.replace('_', '-')
			for k, v in zip(mixwaterdata[0], l)
			if v and v != '-'
		}
		for l in mixwaterdata[1:]
	]

	"""IMPORT RM DATA"""
	with open('../7-ref-materials/output/samples.csv') as fid:
		rmdata = [l.strip().split(',') for l in fid.readlines()]
	rmdata = [
		{
			k: float(v) if k != 'Sample' else v
			for k, v in zip(rmdata[0], l)
			if v and v != '-'
		}
		for l in rmdata[1:]
	]

	data = {r['Sample'].split('_')[0]: r for r in mixwaterdata + rmdata}

	kw_plot_water = dict(
		ls = 'None', marker = 'o', mec = 'k', mfc = 'w', ms = 4, zorder = 10, mew = 0.8
	)
	kw_plot_carb = dict(
		ls = 'None', marker = 's', mec = 'k', mfc = 'w', ms = 4, zorder = 10, mew = 0.8
	)
	kw_errorbar = dict(
		ls = 'None',
		marker = 'None',
		ecolor = 'k',
		elinewidth = 0.8,
		capthick = 0.8,
		zorder = 10,
		capsize = 3,
	)
	kw_mixlines = dict(
		ls = '-', marker = 'None', color = 'k', lw = 0.8, dashes = (8, 2, 2, 2), zorder = 5
	)

	anchors = ['VSMOW2', 'SLAP2', 'HAWAI', 'OC4']
	uwaters = ['GRESP', 'MIX-OH', 'MIX-NH', 'MIX-ONH', 'NEEM']
	waters = anchors + uwaters
	carbs = ['NBS18', 'NBS19', 'IAEA603', 'IAEA610', 'IAEA611', 'IAEA612']

	fig = ppl.figure(figsize = (3.15, 4))
	fig.subplots_adjust(0.23, 0.12, 0.98, 0.84)
	ax = ppl.subplot(111)

	for k, s in enumerate(waters):
		if s in uwaters:
			ax.errorbar(
				data[s]['d18O_VSMOW'],
				data[s]['D17O_VSMOW'],
				data[s]['95CL_D17O_VSMOW'],
				**kw_errorbar,
			)
		ax.plot(
			data[s]['d18O_VSMOW'],
			data[s]['D17O_VSMOW'],
			label = None if k else '$CO_2$ equilibrated with $H_2O$ at 25$\\,$°C',
			**kw_plot_water,
		)

	for k, s in enumerate(carbs):
		ax.errorbar(
			data[s]['d18O_VSMOW'],
			data[s]['D17O_VSMOW'],
			data[s]['95CL_D17O_VSMOW'],
			**kw_errorbar,
		)
		ax.plot(
			data[s]['d18O_VSMOW'],
			data[s]['D17O_VSMOW'],
			label = None if k else '$CO_2$ from $CaCO_3$ (90$\\,$°C acid reaction)',
			**kw_plot_carb,
		)

	kw = {
		'VSMOW2': dict(va = 'bottom', ha = 'left', dx = -2, dy = +0.004),
		'SLAP2': dict(va = 'top', ha = 'right', dx = +3, dy = -0.005),
		'HAWAI': dict(va = 'top', ha = 'left', dx = 0, dy = -0.004),
		'OC4': dict(va = 'bottom', ha = 'left', dx = -4, dy = +0.004),
		'GRESP': dict(va = 'top', ha = 'right', dx = -2, dy = 0),
		'MIX-OH': dict(va = 'top', ha = 'center', dx = 0, dy = -0.006),
		'MIX-NH': dict(va = 'center', ha = 'right', dx = -3, dy = -0.001),
		'MIX-ONH': dict(va = 'bottom', ha = 'right', dx = -1, dy = +0.005),
		'NEEM': dict(va = 'bottom', ha = 'center', dx = -3, dy = +0.006),
		'NBS18': dict(va = 'top', ha = 'center', dx = +1, dy = -0.006),
		'NBS19': dict(va = 'center', ha = 'right', dx = -3, dy = -0.003),
		'IAEA603': dict(va = 'center', ha = 'right', dx = -3, dy = +0.001),
		'IAEA610': dict(va = 'center', ha = 'right', dx = -3, dy = -0.001),
		'IAEA611': dict(va = 'top', ha = 'center', dx = 0, dy = -0.010),
		'IAEA612': dict(va = 'center', ha = 'left', dx = +3, dy = -0.001),
	}
	for s in carbs + waters:
		ax.text(
			data[s]['d18O_VSMOW'] + kw[s]['dx'],
			data[s]['D17O_VSMOW'] + kw[s]['dy'],
			s,
			# 			weight = 'bold',
			size = 8,
			color = [0.6] * 3,
			**{k: kw[s][k] for k in ['va', 'ha']},
		)

	ax.plot(d18_OHmix, D17_OHmix, label = 'mixing lines', **kw_mixlines)
	ax.plot(d18_NHmix, D17_NHmix, **kw_mixlines)
	ax.plot(d18_ONHmix, D17_ONHmix, **kw_mixlines)

	ax.set_xlabel('δ$^{18}O_{VSMOW}$ (‰)')
	ax.set_ylabel('Δ$^{17}O_{VSMOW}$ (‰)')
	ax.legend(
		loc = 'lower center',
		bbox_to_anchor = (0.5, 1.01),
		fontsize = 8.5,
		handlelength = 2,
		borderpad = 0,
		frameon = False,
	)

	x1, x2, y1, y2 = ax.axis()
	x0 = (x1 + x2) / 2
	y0 = (y1 + y2) / 2
	ax.axis([x0 - 40, x0 + 45, y0 - 0.098, y0 + 0.090])

	ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05, offset = 0.03))

	fig.savefig('output/all-together')
	ppl.close(fig)

	fig = ppl.figure(figsize = (6.62, 4))
	fig.subplots_adjust(0.11, 0.12, 0.6, 0.97)
	ax = ppl.subplot(111)

	for k, s in enumerate(waters):
		if s in uwaters:
			ax.errorbar(
				data[s]['d18O_VSMOW'],
				data[s]['D17O_VSMOW'],
				data[s]['95CL_D17O_VSMOW'],
				**kw_errorbar,
			)
		ax.plot(
			data[s]['d18O_VSMOW'],
			data[s]['D17O_VSMOW'],
			label = None if k else '$CO_2$ equilibrated with $H_2O$ at 25$\\,$°C',
			**kw_plot_water,
		)

	for k, s in enumerate(carbs):
		ax.errorbar(
			data[s]['d18O_VSMOW'],
			data[s]['D17O_VSMOW'],
			data[s]['95CL_D17O_VSMOW'],
			**kw_errorbar,
		)
		ax.plot(
			data[s]['d18O_VSMOW'],
			data[s]['D17O_VSMOW'],
			label = None if k else '$CO_2$ from $CaCO_3$ (90$\\,$°C acid reaction)',
			**kw_plot_carb,
		)

	kw = {
		'VSMOW2': dict(va = 'bottom', ha = 'left', dx = -2, dy = +0.004),
		'SLAP2': dict(va = 'top', ha = 'right', dx = +3, dy = -0.005),
		'HAWAI': dict(va = 'top', ha = 'left', dx = 0, dy = -0.004),
		'OC4': dict(va = 'bottom', ha = 'left', dx = -4, dy = +0.004),
		'GRESP': dict(va = 'top', ha = 'right', dx = -2, dy = 0),
		'MIX-OH': dict(va = 'top', ha = 'center', dx = 0, dy = -0.006),
		'MIX-NH': dict(va = 'center', ha = 'right', dx = -3, dy = -0.001),
		'MIX-ONH': dict(va = 'bottom', ha = 'right', dx = -1, dy = +0.005),
		'NEEM': dict(va = 'bottom', ha = 'center', dx = -3, dy = +0.006),
		'NBS18': dict(va = 'top', ha = 'center', dx = +1, dy = -0.006),
		'NBS19': dict(va = 'center', ha = 'right', dx = -3, dy = -0.003),
		'IAEA603': dict(va = 'center', ha = 'right', dx = -3, dy = +0.001),
		'IAEA610': dict(va = 'center', ha = 'right', dx = -3, dy = -0.001),
		'IAEA611': dict(va = 'top', ha = 'center', dx = 0, dy = -0.010),
		'IAEA612': dict(va = 'center', ha = 'left', dx = +3, dy = -0.001),
	}
	for s in carbs + waters:
		ax.text(
			data[s]['d18O_VSMOW'] + kw[s]['dx'],
			data[s]['D17O_VSMOW'] + kw[s]['dy'],
			s,
			# 			weight = 'bold',
			size = 8,
			color = [0.6] * 3,
			**{k: kw[s][k] for k in ['va', 'ha']},
		)

	ax.plot(d18_OHmix, D17_OHmix, label = 'mixing lines', **kw_mixlines)
	ax.plot(d18_NHmix, D17_NHmix, **kw_mixlines)
	ax.plot(d18_ONHmix, D17_ONHmix, **kw_mixlines)

	waters = {
		'VSMOW2': {'d18O_VSMOW': 0.0, 'D17O_VSMOW': 0.0},
		'SLAP2': {'d18O_VSMOW': -55.5, 'D17O_VSMOW': 0.0},
	}
	kw_plot_water['mfc'] = [0.6] * 3
	for k, s in enumerate(waters):
		ax.plot(
			waters[s]['d18O_VSMOW'],
			waters[s]['D17O_VSMOW'],
			label = None if k else '$H_2O$',
			**kw_plot_water,
		)

	carbs = {
		'IAEA603': {'d18O_VSMOW': (1000 - 2.37) * 1.03092 - 1000, 'D17O_VSMOW': -0.082},
		'NBS18': {'d18O_VSMOW': (1000 - 23.01) * 1.03092 - 1000, 'D17O_VSMOW': -0.033},
	}
	kw_plot_carb['mfc'] = [0.6] * 3
	for k, s in enumerate(carbs):
		ax.plot(
			carbs[s]['d18O_VSMOW'],
			carbs[s]['D17O_VSMOW'],
			label = None if k else '$CaCO_3$',
			**kw_plot_carb,
		)

	ax.set_xlabel('δ$^{18}O_{VSMOW}$ (‰)')
	ax.set_ylabel('Δ$^{17}O_{VSMOW}$ (‰)')
	ax.legend(
		loc = 'center left',
		bbox_to_anchor = (1.01, 0.5),
		fontsize = 8.5,
		handlelength = 2,
		borderpad = 0,
		frameon = False,
	)

	x1, x2, y1, y2 = ax.axis()
	x0 = (x1 + x2) / 2
	y0 = (y1 + y2) / 2
	ax.axis([x0 - 55, x0 + 65, y0 - 0.14, y0 + 0.14])

	ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05, offset = 0.03))

	fig.savefig('output/all-together-wider')
	ppl.close(fig)

	"""IMPORT PUBLISHED DATA"""
	with open('../8-pub-comparison/input/published.csv') as fid:
		pubdata = [l.strip().split('\t') for l in fid.readlines()]
	pubdata = [
		{k: v for k, v in zip(pubdata[0], l) if v} for l in pubdata[1:] if len(l) > 1
	]
	for r in pubdata:
		if 'Tacid' not in r:
			r['Tacid'] = None
		for k in r:
			if 'D17' in k:
				r[k] = float(r[k])

	unmco2data = {
		_['Sample']: _
		for _ in pubdata
		if _['Ref'] == 'Wostbrock et al. (2020)' and _['Tacid'] is not None
	}
	for s in unmco2data:
		r = unmco2data[s]
		r['d18O_VSMOW'] = (
			1000 + {'IAEA603': -2.37, 'NBS18': -23.01, 'NBS19': -2.2}[r['Sample']]
		) * 1.01025 * alpha18_CO2g_H2O - 1000
		r['D17O_VSMOW'] = r['D17O'] + 0.015 * np.log(
			1 + r['d18O_VSMOW'] / 1000
		) / np.log(1 - 55.5e-3)

	fig = ppl.figure(figsize = (3.15, 4.5))
	fig.subplots_adjust(0.23, 0.11, 0.98, 0.66)
	ax = ppl.subplot(111)

	vsmowslap = np.array([[0, -55.5], [0, 0]])
	co2vsmowslap = np.array(
		[
			(1000 + vsmowslap[0]) * alpha18_CO2g_H2O - 1000,
			(vsmowslap[1] + (theta17_CO2g_H2O - LAMBDA_17) * np.log(alpha18_CO2g_H2O))
			* 1000,
		]
	)

	unmcolor = (1, 0, 0)  # '#ff7f00'
	lscecolor = (0, 0.7, 0.8)
	o2color = [(2 + _) / 3 for _ in unmcolor]
	co2color = [(2 + _) / 3 for _ in lscecolor]

	kw_plot = dict(
		ls = 'None', marker = 'D', ms = 4, mfc = o2color, mec = 'k', mew = 0.8, zorder = 2000
	)
	ax.plot(*vsmowslap, label = '$O_2$ from $H_2O$ (Wostbrock et al., 2020)', **kw_plot)

	kw_plot['marker'] = 's'
	kw_plot['mfc'] = o2color
	for k, s in enumerate(['NBS18', 'IAEA603']):
		unmco2data[s]['d18O_VSMOW'] = (
			1000 + data[s]['d18O_VSMOW']
		) * 1.01025 / 1.00813 - 1000
		if 'eD17O' in unmco2data[s]:
			ax.errorbar(
				unmco2data[s]['d18O_VSMOW'],
				unmco2data[s]['D17O_VSMOW'],
				unmco2data[s]['eD17O'],
				**kw_errorbar,
			)
		ax.plot(
			unmco2data[s]['d18O_VSMOW'],
			unmco2data[s]['D17O_VSMOW'],
			label = None
			if k
			else '$O_2$ from $CO_2$ from $CaCO_3$ (Wostbrock et al., 2020)',
			**kw_plot,
		)

	kw_plot['marker'] = 'D'
	kw_plot['mfc'] = co2color
	ax.plot(*co2vsmowslap, label = '$H_2O$-equilibrated $CO_2$ (this study)', **kw_plot)

	kw_plot['marker'] = 's'
	for k, s in enumerate(['NBS18', 'IAEA603']):
		ax.errorbar(
			data[s]['d18O_VSMOW'],
			data[s]['D17O_VSMOW'],
			data[s]['95CL_D17O_VSMOW'],
			**kw_errorbar,
		)
		ax.plot(
			data[s]['d18O_VSMOW'],
			data[s]['D17O_VSMOW'],
			label = None if k else '$CO_2$ from $CaCO_3$ (this study)',
			**kw_plot,
		)

	X0 = vsmowslap[0].mean()
	X1, Y1 = data['NBS18']['d18O_VSMOW'], data['NBS18']['D17O_VSMOW']
	X2, Y2 = data['IAEA603']['d18O_VSMOW'], data['IAEA603']['D17O_VSMOW']
	X3, Y3 = unmco2data['NBS18']['d18O_VSMOW'], unmco2data['NBS18']['D17O_VSMOW']
	X4, Y4 = unmco2data['IAEA603']['d18O_VSMOW'], unmco2data['IAEA603']['D17O_VSMOW']

	ARBITRARY_FACTOR = 1
	# 	ARBITRARY_FACTOR = 0.6 # this would be enough to reconcile UNM and LSCE observations
	a = ((Y4 - Y3) - (Y2 - Y1)) / ((X3 - X0) ** 2 - (X4 - X0) ** 2) * ARBITRARY_FACTOR
	f = lambda x: -a * (x - X0) ** 2 + a * X0**2

	# 	ax.annotate(
	# 		'',
	# 		xy = (X0, 0),
	# 		xytext = (X0, a * X0**2),
	# 		arrowprops = dict(
	# 			arrowstyle = "->",
	# 			color = unmcolor,
	# 			shrinkA = 0,
	# 			shrinkB = 0,
	# 			patchA = None,
	# 			patchB = None,
	# 			lw = .8,
	# 			mutation_scale = 7,
	# 			),
	# 		)
	# 	ax.plot(X0, a * X0**2, 'o', mew = 0, mfc = unmcolor, ms = 3)

	xi = np.linspace(-60, 45, 111)
	yi = f(xi)
	ax.plot(
		xi,
		yi,
		'-',
		color = unmcolor,
		lw = 0.8,
		dashes = (6, 2, 2, 2),
		label = '\nHypothetical quadratic correction reducing\nthe difference between NBS18 and IAEA603',
	)
	logger.info(
		f'Hypothetical quadratic correction for UNM data would reach {f(X0)*1000:+.1f} ppm between VSMOW2 and SLAP2'
	)
	logger.info(
		f"This would change IAEA603 by {-f(unmco2data['IAEA603']['d18O_VSMOW'])*1000:.1f} ppm and NBS18 by {-f(unmco2data['NBS18']['d18O_VSMOW'])*1000:.1f} ppm."
	)

	kw_plot['mec'] = unmcolor
	kw_plot['mfc'] = 'w'
	kw_errorbar['ecolor'] = unmcolor
	for k, s in enumerate(['NBS18', 'IAEA603']):
		if 'eD17O' in unmco2data[s]:
			ax.errorbar(
				unmco2data[s]['d18O_VSMOW'],
				unmco2data[s]['D17O_VSMOW'] - f(unmco2data[s]['d18O_VSMOW']),
				unmco2data[s]['eD17O'],
				**kw_errorbar,
			)
		ax.plot(
			unmco2data[s]['d18O_VSMOW'],
			unmco2data[s]['D17O_VSMOW'] - f(unmco2data[s]['d18O_VSMOW']),
			**kw_plot,
		)

	for s, aA, aB, shA, shB in [
		('NBS18', 15, -15, 4, 4),
		('IAEA603', 15, -50, 4, 4),
	]:
		ax.annotate(
			'',
			xy = (
				unmco2data[s]['d18O_VSMOW'],
				unmco2data[s]['D17O_VSMOW'] - f(unmco2data[s]['d18O_VSMOW']),
			),
			xytext = (unmco2data[s]['d18O_VSMOW'], unmco2data[s]['D17O_VSMOW']),
			arrowprops = dict(
				arrowstyle = '->',
				color = unmcolor,
				shrinkA = shA,
				shrinkB = shB,
				patchA = None,
				patchB = None,
				lw = 0.8,
				connectionstyle = f'angle3, angleA = {aA:.0f}, angleB = {aB:.0f}',
			),
		)
		ax.annotate(
			'',
			xy = (unmco2data[s]['d18O_VSMOW'], 0),
			xytext = (unmco2data[s]['d18O_VSMOW'], f(unmco2data[s]['d18O_VSMOW'])),
			arrowprops = dict(
				arrowstyle = '->',
				color = unmcolor,
				shrinkA = 0,
				shrinkB = 0,
				patchA = None,
				patchB = None,
				lw = 0.8,
				mutation_scale = 7,
			),
		)
		ax.plot(
			unmco2data[s]['d18O_VSMOW'],
			f(unmco2data[s]['d18O_VSMOW']),
			'o',
			mew = 0,
			mfc = unmcolor,
			ms = 3,
		)

	X0 = co2vsmowslap[0].mean()

	a_lsce = -((Y4 - Y3) - (Y2 - Y1)) / ((X1 - X0) ** 2 - (X2 - X0) ** 2)
	f = lambda x: -a_lsce * (x - X0) ** 2 + a_lsce * (co2vsmowslap[0][0] - X0) ** 2

	logger.info(
		f'Hypothetical quadratic correction for LSCE data would reach {f(X0)*1000:+.1f} ppm between VSMOW2 and SLAP2'
	)
	logger.info(
		f"This would change IAEA603 by {-f(data['IAEA603']['d18O_VSMOW'])*1000:.1f} ppm and NBS18 by {-f(data['NBS18']['d18O_VSMOW'])*1000:.1f} ppm."
	)

	xi = np.linspace(-20, 45)
	yi = f(xi)
	ax.plot(
		xi,
		yi + co2vsmowslap[1].mean(),
		'-',
		color = lscecolor,
		lw = 0.8,
		dashes = (6, 2, 2, 2),
		label = '\nHypothetical quadratic correction increasing\nthe difference between NBS18 and IAEA603',
	)

	kw_plot['mec'] = lscecolor
	kw_plot['mfc'] = 'w'
	kw_errorbar['ecolor'] = lscecolor
	for k, s in enumerate(['NBS18', 'IAEA603']):
		ax.errorbar(
			data[s]['d18O_VSMOW'],
			data[s]['D17O_VSMOW'] - f(data[s]['d18O_VSMOW']),
			data[s]['95CL_D17O_VSMOW'],
			**kw_errorbar,
		)
		ax.plot(
			data[s]['d18O_VSMOW'],
			data[s]['D17O_VSMOW'] - f(data[s]['d18O_VSMOW']),
			**kw_plot,
		)

	for s, aA, aB, shA, shB in [
		('NBS18', 0, 50, 4, 6),
		('IAEA603', -10, 10, 4, 3),
	]:
		ax.annotate(
			'',
			xy = (
				data[s]['d18O_VSMOW'],
				data[s]['D17O_VSMOW'] - f(data[s]['d18O_VSMOW']),
			),
			xytext = (data[s]['d18O_VSMOW'], data[s]['D17O_VSMOW']),
			arrowprops = dict(
				arrowstyle = '->',
				color = lscecolor,
				shrinkA = shA,
				shrinkB = shB,
				patchA = None,
				patchB = None,
				lw = 0.8,
				connectionstyle = f'angle3, angleA = {aA:.0f}, angleB = {aB:.0f}',
			),
		)
		ax.annotate(
			'',
			xy = (data[s]['d18O_VSMOW'], co2vsmowslap[1].mean()),
			xytext = (
				data[s]['d18O_VSMOW'],
				co2vsmowslap[1].mean() + f(data[s]['d18O_VSMOW']),
			),
			arrowprops = dict(
				arrowstyle = '->',
				color = lscecolor,
				shrinkA = 0,
				shrinkB = 0,
				patchA = None,
				patchB = None,
				lw = 0.8,
				mutation_scale = 7,
			),
		)
		ax.plot(
			data[s]['d18O_VSMOW'],
			co2vsmowslap[1].mean() + f(data[s]['d18O_VSMOW']),
			'o',
			mew = 0,
			mfc = lscecolor,
			ms = 3,
		)

	ax.set_xlabel('δ$^{18}O_{VSMOW}$ (‰)')
	ax.set_ylabel('Δ$^{17}O_{VSMOW}$ (‰)')
	ax.legend(
		loc = 'lower center',
		bbox_to_anchor = (0.4, 1.02),
		fontsize = 8,
		handlelength = 1.7,
		borderpad = 0,
		frameon = False,
		labelspacing = 0.4,
	)

	ax.axis([-60, 48, None, None])

	ax.xaxis.set_major_locator(ticker.MultipleLocator(25))
	ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05, offset = 0.0))

	fig.savefig('output/pbl-simulation')
	ppl.close(fig)

	with open('../4-linearity/output/samples.csv') as fid:
		lindata = [l.strip().split(',') for l in fid.readlines()]
	lindata = [{k: v for k, v in zip(lindata[0], l)} for l in lindata[1:]]
	for r in lindata:
		r['Sample'] = r['Sample'].replace('_', '-')
		for k in r:
			if k != 'Sample':
				if r[k] == '-':
					r[k] = None
				else:
					r[k] = float(r[k])

	with open('../4-linearity/output/results.csv') as fid:
		mixwaterresults = [l.strip().split(',') for l in fid.readlines()]
	mixwaterresults = [
		{
			k: float(v) if k != 'Sample' else v.replace('_', '-')
			for k, v in zip(mixwaterresults[0], l)
			if v and v != '-'
		}
		for l in mixwaterresults[1:]
	]
	mixwaterresults = {r['Sample']: r for r in mixwaterresults}

	for l in lindata:
		l['D17O_VSMOW_CO2_prediction'] = mixwaterresults[l['Sample']][
			'final_D17O_VSMOW_CO2_predicted'
		]

	lindata = {r['Sample']: r for r in lindata}

	with open('../7-ref-materials/output/samples.csv') as fid:
		rmdata = [l.strip().split(',') for l in fid.readlines()]
		rmdata = [{k: v for k, v in zip(rmdata[0], l)} for l in rmdata[1:]]
		for r in rmdata:
			r['Sample'] = r['Sample'].replace('_', '-')
			for k in r:
				if k != 'Sample':
					if r[k] == '-':
						r[k] = None
					else:
						r[k] = float(r[k])
		rmdata = {r['Sample']: r for r in rmdata}

	data = lindata | rmdata

	vsmowdata = [rmdata[s] for s in rmdata if s.startswith('VSMOW2')]
	N = sum([s['N'] for s in vsmowdata])
	data['VSMOW2'] = {
		'N': N,
		'd18O_VSMOW': sum([s['N'] / N * s['d18O_VSMOW'] for s in vsmowdata]),
		'SD_d18O_VSMOW': (
			(sum([s['SD_d18O_VSMOW'] ** 2 * (s['N'] - 1) for s in vsmowdata]))
			/ (N - len(vsmowdata))
		)
		** 0.5,
		'D17O_VSMOW': sum([s['N'] / N * s['D17O_VSMOW'] for s in vsmowdata]),
		'SD_D17O_VSMOW': (
			(sum([s['SD_D17O_VSMOW'] ** 2 * (s['N'] - 1) for s in vsmowdata]))
			/ (N - len(vsmowdata))
		)
		** 0.5,
		'D17O_VSMOW_CO2_prediction': sum(
			[s['N'] / N * s['D17O_VSMOW'] for s in vsmowdata]
		),
	}
	slapdata = [rmdata[s] for s in rmdata if s.startswith('SLAP2')]
	N = sum([s['N'] for s in slapdata])
	data['SLAP2'] = {
		'N': N,
		'd18O_VSMOW': sum([s['N'] / N * s['d18O_VSMOW'] for s in slapdata]),
		'SD_d18O_VSMOW': (
			(sum([s['SD_d18O_VSMOW'] ** 2 * (s['N'] - 1) for s in slapdata]))
			/ (N - len(slapdata))
		)
		** 0.5,
		'D17O_VSMOW': sum([s['N'] / N * s['D17O_VSMOW'] for s in slapdata]),
		'SD_D17O_VSMOW': (
			(sum([s['SD_D17O_VSMOW'] ** 2 * (s['N'] - 1) for s in slapdata]))
			/ (N - len(slapdata))
		)
		** 0.5,
		'D17O_VSMOW_CO2_prediction': sum(
			[s['N'] / N * s['D17O_VSMOW'] for s in slapdata]
		),
	}

	with open('output/samples.csv', 'w') as fid:
		fid.write(
			'Sample,Type,D17O_VSMOW_predicted,N,d18O_VSMOW,SD_d18O_VSMOW,95CL_d18O_VSMOW,D17O_VSMOW,SD_D17O_VSMOW,95CL_D17O_VSMOW,d13C_VPDB,SD_d13C_VPDB,95CL_d13C_VPDB,d18O_VPDB,SD_d18O_VPDB,95CL_d18O_VPDB'
		)
		for s in ['VSMOW2', 'SLAP2', 'HAWAI', 'OC4']:
			fid.write(f'\n{s}-CO2,W')
			fid.write(
				f",{data[s]['D17O_VSMOW_CO2_prediction']:.4f}"
				if 'D17O_VSMOW_CO2_prediction' in data[s]
				else ','
			)
			fid.write(f",{data[s]['N']:.0f}")
			fid.write(f",{data[s]['d18O_VSMOW']:.2f}*")
			fid.write(f",{data[s]['SD_d18O_VSMOW']:.2f}")
			fid.write(',')
			fid.write(f",{data[s]['D17O_VSMOW']:.4f}*")
			fid.write(f",{data[s]['SD_D17O_VSMOW']:.4f}")
			fid.write(',')
			fid.write(',')
			fid.write(',')
			fid.write(',')
			fid.write(',')
			fid.write(',')
			fid.write(',')

		for s in ['NEEM', 'MIX-NH', 'MIX-ONH', 'MIX-OH', 'GRESP']:
			fid.write(f'\n{s}-CO2,W')
			fid.write(
				f",{data[s]['D17O_VSMOW_CO2_prediction']:.4f}"
				if 'D17O_VSMOW_CO2_prediction' in data[s]
				else ','
			)
			fid.write(f",{data[s]['N']:.0f}")
			fid.write(f",{data[s]['d18O_VSMOW']:.2f}")
			fid.write(f",{data[s]['SD_d18O_VSMOW']:.2f}")
			fid.write(f",{data[s]['95CL_d18O_VSMOW']:.2f}")
			fid.write(f",{data[s]['D17O_VSMOW']:.4f}")
			fid.write(f",{data[s]['SD_D17O_VSMOW']:.4f}")
			fid.write(f",{data[s]['95CL_D17O_VSMOW']:.4f}")
			fid.write(',')
			fid.write(',')
			fid.write(',')
			fid.write(',')
			fid.write(',')
			fid.write(',')

		for s in ['NBS18']:
			fid.write(f'\n{s}-CO2,C')
			fid.write(f",,{data[s]['N']:.0f}")
			fid.write(f",{data[s]['d18O_VSMOW']:.2f}")
			fid.write(f",{data[s]['SD_d18O_VSMOW']:.2f}")
			fid.write(f",{data[s]['95CL_d18O_VSMOW']:.2f}")
			fid.write(f",{data[s]['D17O_VSMOW']:.4f}")
			fid.write(f",{data[s]['SD_D17O_VSMOW']:.4f}")
			fid.write(f",{data[s]['95CL_D17O_VSMOW']:.4f}")
			fid.write(f",{data[s]['d13C_VPDB']:.2f}")
			fid.write(f",{data[s]['SD_d13C_VPDB']:.2f}")
			fid.write(f",{data[s]['95CL_d13C_VPDB']:.2f}")
			fid.write(f",{data[s]['d18O_VPDB']:.2f}*")
			fid.write(f",{data[s]['SD_d18O_VPDB']:.2f}")
			fid.write(',')

		for s in ['NBS19', 'IAEA603']:
			fid.write(f'\n{s}-CO2,C')
			fid.write(f",,{data[s]['N']:.0f}")
			fid.write(f",{data[s]['d18O_VSMOW']:.2f}")
			fid.write(f",{data[s]['SD_d18O_VSMOW']:.2f}")
			fid.write(f",{data[s]['95CL_d18O_VSMOW']:.2f}")
			fid.write(f",{data[s]['D17O_VSMOW']:.4f}")
			fid.write(f",{data[s]['SD_D17O_VSMOW']:.4f}")
			fid.write(f",{data[s]['95CL_D17O_VSMOW']:.4f}")
			fid.write(f",{data[s]['d13C_VPDB']:.2f}*")
			fid.write(f",{data[s]['SD_d13C_VPDB']:.2f}")
			fid.write(f',')
			fid.write(f",{data[s]['d18O_VPDB']:.2f}*")
			fid.write(f",{data[s]['SD_d18O_VPDB']:.2f}")
			fid.write(',')

		for s in ['IAEA610', 'IAEA611', 'IAEA612']:
			fid.write(f'\n{s}-CO2,C')
			fid.write(f",,{data[s]['N']:.0f}")
			fid.write(f",{data[s]['d18O_VSMOW']:.2f}")
			fid.write(f",{data[s]['SD_d18O_VSMOW']:.2f}")
			fid.write(f",{data[s]['95CL_d18O_VSMOW']:.2f}")
			fid.write(f",{data[s]['D17O_VSMOW']:.4f}")
			fid.write(f",{data[s]['SD_D17O_VSMOW']:.4f}")
			fid.write(f",{data[s]['95CL_D17O_VSMOW']:.4f}")
			fid.write(f",{data[s]['d13C_VPDB']:.3f}*")
			fid.write(f",{data[s]['SD_d13C_VPDB']:.2f}")
			fid.write(f',')
			fid.write(f",{data[s]['d18O_VPDB']:.2f}")
			fid.write(f",{data[s]['SD_d18O_VPDB']:.2f}")
			fid.write(f",{data[s]['95CL_d18O_VPDB']:.2f}")
