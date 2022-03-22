import subprocess
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directory", help = "Data directory", required = True)
args = parser.parse_args()

data_paths = {
#    1: 'planets.csv',
#    2: 'rgb_tip_profile.data',
#    3: 'rgb_tip_profile.data',
#    4: 'successful_ejection.h5',
#    6: 'grid_1/',
#    7: 'drag_coefficients/drag.h5',
#    8: 'grid_2/',
#    9: './',
    10: './'
}

for i in data_paths.keys():
    print("Plotting figure %i" % i)
    subprocess.run("python %i.py -i %s" % ( i, os.path.join(args.directory, data_paths[i]) ), shell=True, check=True)
