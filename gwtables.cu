struct fft {
	int max_exp;
	int fftlen;
		};

struct fft gwtable[] = {
	755,	32,	939,	40,	1113,	48,	1303,	56,	1499,	64,
	1857,	80,	2211,	96,	2585,	112,	2953,	128,	3663,	160,
	4359,	192,	5093,	224,	5833,	256,	7243,	320,	8639,	384,
	10085,448,	11537,512,	14301,640,	17047,768,	19881,896,
	22799,1024,	28295,1280,	33761,1536,	39411,1792,	45061,2048,
	55825,2560,	66519,3072,	77599,3584,	89047,4096,	110400,5120,
	131100,6144,152800,7168,175300,8192,217700,10240,258200,12288,
	301400,14336,346100,16384,430300,20480,511600,24576,596100,28672,
	683700,32768,848800,40960,1009000,49152,1177000,57344,1350000,65536,
	1678000,81920,1994000,98304,2324000,114688,2664000,131072,3310000,163840,
	3933000,196608,4593000,229376,5264000,262144,6545000,327680,7772000,393216,
	9071000,458752,10380000,524288,12890000,655360,15310000,786432,17890000,917504,
	20460000,1048576,25390000,1310720,30190000,1572864,35200000,1835008,40300000,2097152,
	50020000,2621440,59510000,3145728,69360000,3670016,79370000,4194304,0,0
	};

struct fft gwtablep[] = {
	755,	32,	1111,	48,	1485,	64,	2199,	96,	2947,	128,
	4345,	192,	5817,	256,	8607,	384,	11515,512,	17001,768,
	22701,1024,	33569,1536,	44951,2048,	66319,3072,	88747,4096,
	130600,6144,174000,8192,257700,12288,344700,16384,508600,24576,
	679400,32768,1006000,49152,1345000, 65536,1983000,98304,2652000,131072,
	3924000,196608,5242000,262144,7733000,393216,10320000,524288,15260000,786432,
	20360000,1048576,30070000,1572864,40110000,2097152,59360000,3145728,79100000,4194304,
	0,0};

