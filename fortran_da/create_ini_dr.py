import numpy as np
import random
import sys

xmin = 0.5
xmax = 3.5
ymin = 0.5
ymax = 2.5
dr_num = int(sys.argv[1]) 
dr_lon = [random.uniform(xmin,xmax) for i in np.arange(dr_num)]
dr_lat = [random.uniform(ymin,ymax) for i in np.arange(dr_num)]

file = open("IC_dr_%0.3d.nml"%dr_num,"w")
file.write("!------------------------------------------------------------------------------!\n")
file.write("! Namelist file :\n")
file.write("! Initial condition.\n")
file.write("!------------------------------------------------------------------------------!\n")
file.write("\n")
file.write("&ICDRLIST\n")
for j in np.arange(dr_num):
    file.write("  IC_DR(%d) =    %f\n"%(j*2+1,dr_lon[j]))
    file.write("  IC_DR(%d) =    %f\n"%(j*2+2,dr_lat[j]))
file.write("&END")
file.write("\n")

file.write("!------------------------------------------------------------------------------!\n")
file.write("! Initialisation type.                                                         !\n")
file.write("!------------------------------------------------------------------------------!\n")
file.write("! type = 'read': use IC above (will generate a new seed);\n")
file.write("!        'rand': random state (will generate a new seed);\n")
file.write("!        'zero': zero IC (will generate a new seed);\n")
file.write("!        'seed': use the seed below (generate the same IC)\n")

file.write("&DRRAND\n")
file.write("  init_type= 'read'\n")
file.write("  size_of_random_noise =   0.0000000D+00\n")
file.write("  seed(1) =   2087639516\n")
file.write("  seed(2) =    590834459\n")
file.write("  seed(3) =    900301077\n")
file.write("  seed(4) =   -411484539\n")
file.write("  seed(5) =   2037166880\n")
file.write("  seed(6) =  -2127688463\n")
file.write("  seed(7) =  -1633681454\n")
file.write("  seed(8) =    398303250\n")
file.write("  seed(9) =   1362561686\n")
file.write("  seed(10) =   -644462581\n")
file.write("  seed(11) =   -991677581\n")
file.write("  seed(12) =  -1031490702\n")
file.write("&END\n")
file.close()
