import sys

sys.path.append('../lib')

from mylib import *

RAWDATA_FILENAME = 'input/carbonate-repeatability-rawdata.csv'

logger.remove(0)
logger.add(
	'output/carbonate-repeatability.log', mode = 'w', format = '[{line:3.0f}]    {message}'
)

ppl.style.use('../lib/mydefault.mplstyle')

if __name__ == '__main__':
	data = read_rawdata([RAWDATA_FILENAME])

	d627 = np.array([_['d627'] for _ in data])
	d628 = np.array([_['d628'] for _ in data])
	D17O = 1e6 * (np.log(1 + d627 / 1e3) - LAMBDA_17 * np.log(1 + d628 / 1e3))
	logger.info(
		f'Carbonate Δ17O repeatability = {D17O.std(ddof = 1):.1f} ppm (1SD, N = {len(data)})'
	)
