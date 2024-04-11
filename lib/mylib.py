import isofunctions
import stdz
import co2irms
import numpy as np
from matplotlib import pyplot as ppl
from matplotlib import ticker
from matplotlib.patches import Ellipse
from matplotlib.figure import Figure as figureclass
from statistics import stdev
from loguru import logger
from scipy import stats

from isofunctions import (
	LAMBDA_17,
	R18_VSMOW,
	R17_VSMOW,
	CO2g_H2O_FRACTIONATION,
	alpha18_CO2g_H2O,
	theta17_CO2g_H2O,
	alpha17_CO2g_H2O,
)

import warnings

with warnings.catch_warnings():
	warnings.simplefilter('ignore')
	from allantools import adev
	# this is to suppress SyntaxWarnings
	# ("is" with a literal. Did you mean " = ="?)
	# when first importing allantools.

AUTHOR = 'M. Daëron (daeron@lsce.ipsl.fr)'

SAVEFIG_FORMATS = []
SAVEFIG_FORMATS.append('pdf')
# SAVEFIG_FORMATS.append('png')

_old_savefig = figureclass.savefig


def new_savefig(self, fname, *args, **kwargs):
	for f in SAVEFIG_FORMATS:
		_old_savefig(self, f'{fname}.{f}', *args, **({'format': f} | kwargs))


figureclass.savefig = new_savefig

ISO_COLORS = {
	'626': (0, 0, 0),
	'636': (0, 0.7, 0),
	'628': (0.9, 0, 0),
	'627': (1, 0, 1),
	'626b': (0.25, 0.75, 0.75),
}


def read_rawdata(
	list_of_csv_files,
	str_fields = ['UID', 'Run', 'Session', 'Sample'],
):
	from csv import DictReader

	data = []
	for f in list_of_csv_files:
		with open(f) as fid:
			data += [
				{k: l[k] if k in str_fields else float(l[k]) for k in l if l[k] != ''}
				for l in DictReader(fid)
			]

	return data
