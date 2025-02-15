#!/usr/bin/env python3

import numpy as np


def make_data(args):
    """Make data dictionary for HD 189"""
    mass_MJ = 1.142
    radius_RJ = 1.138
    gravity_SI = 23.970 
    Rs_Rsun = 0.805
    inc = 85.71
    t0 = 2454037.612
    sma = 8.839304998   # semi major axis in stellar radiu
    orb_per = 2.21857545   #in days
    ecc = 0.0041
    w_peri = -24.1  # longiutude of periastron
    limbdark = "linear"
    
    u_limbdark = [0.35]
    
    num_transit = 1
    
    dates = [2458383.77055943, 2458383.77384704, 2458383.77707875,
       2458383.78030307, 2458383.78358918, 2458383.78681399,
       2458383.79004101, 2458383.79326712, 2458383.79655574,
       2458383.79984545, 2458383.80307906, 2458383.80629228,
       2458383.80958299, 2458383.8128124 , 2458383.81603942,
       2458383.81925973, 2458383.82248474, 2458383.82577195,
       2458383.82900097, 2458383.83223048, 2458383.8354501 ,
       2458383.83874811, 2458383.84196822, 2458383.84520053,
       2458383.84847654, 2458383.85170346, 2458383.85493727,
       2458383.85821578, 2458383.86144419, 2458383.86466921,
       2458383.86790322, 2458383.87118233, 2458383.87441074,
       2458383.87763435, 2458383.88092406, 2458383.88414957],
    #don't forget the coma at the end if there is only one transit !!!!!
   


    # Wmean = [2400.695909757236,2328.5343131275904,1972.9809993156186,
    #          1927.2107049022654,]
    # Wmean = [1634.5200937047302,1600.8109822367207],[1670.071564637037,1634.5459486709924,1600.8124596368639],
    Wmean = [2328.5343131275904],        
    orderstot = [33]
    orders = [33],
    # orderstot = [46,47,48]
    # orders = [47,48],[46,47,48],
        
    # Vfiles = ["Vcorr47_DRS2.txt",
    #           "Vcorr48_DRS2.txt",
    #           ],["Vcorr46_Jun19-1_DRS2.txt",
    #           "Vcorr47_Jun19-1_DRS2.txt",
    #           "Vcorr48_Jun19-1_DRS2.txt"
    #           ],
    Vfiles = ["V33_CO.txt"],             
        
    Ifiles = ["I33_CO.txt"],
    
    #  if Stdfiles are not needed, for example with the Brogi likelihood, 
    # uncomment the next line
    #Stdfiles = []
    Stdfiles = ["Std33_CO.txt"],
    
    lambdas = np.array([[   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [2291.84518119, 2362.55271775],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [1939.42197854, 1998.81548771],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [1758.50261646, 1812.39702422],
       [1718.50054581, 1771.64067835],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [1512.43747007, 1558.89713666],
       [1484.77586677, 1528.30354258],
       [1457.06015806, 1498.88570675],
       [1429.75333156, 1470.19096444],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [1306.967007  , 1343.21643463],
       [1285.02046052, 1320.56072659],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [1167.78440327, 1198.13940642],
       [1150.59417256, 1178.48372217],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ],
       [   0.        ,    0.        ]])

    return dict(
        mass_MJ=mass_MJ,
        radius_RJ=radius_RJ,
		gravity_SI = gravity_SI,
		Rs_Rsun = Rs_Rsun,
		inc = inc,
		t0 = t0,
		sma = sma,
		orb_per = orb_per,
		ecc = ecc,
		w_peri = w_peri,
        Wmean = Wmean,
		limbdark = limbdark,
		u_limbdark = u_limbdark,
		dates = dates,
		lambdas = lambdas,
        orders = orders,
        orderstot=orderstot,
        num_transit=num_transit,
		Vfiles = Vfiles,
		Ifiles = Ifiles,
		Stdfiles = Stdfiles
		    )
