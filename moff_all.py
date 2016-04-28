import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import glob
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import logging
import itertools
import argparse
import re
import os as os
import sys
import moff_mbr
import moff

parser = argparse.ArgumentParser(description='moFF match between run and apex module input parameter')

parser.add_argument('--map_file', dest='map_file', action='store',
                    help='specify a map file that contains input files   and raw file ', required=True)

parser.add_argument('--log_file_name', dest='log_label', action='store', default='moFF',
                    help='a label name to use for the log file', required=False)

parser.add_argument('--filt_width', dest='w_filt', action='store', default=2,
                    help='width value of the filter  k * mean(Dist_Malahobis)', required=False)

parser.add_argument('--out_filt', dest='out_flag', action='store', default=1,
                    help='filter outlier in each rt time allignment', required=False)

parser.add_argument('--weight_comb', dest='w_comb', action='store', default=0,
                    help='weights for model combination combination : 0 for no weight  1 weighted devised by trein err of the model.',
                    required=False)

# parser.add_argument('--input', dest='name', action='store',help='specify input list of MS2 peptides ', required=True)

parser.add_argument('--tol', dest='toll', action='store', type=float, help='specify the tollerance  parameter in ppm',
                    required=True)

parser.add_argument('--rt_w', dest='rt_window', action='store', type=float, default=3,
                    help='specify rt window for xic (minute). Default value is 3 min', required=False)

parser.add_argument('--rt_p', dest='rt_p_window', action='store', type=float, default=0.1,
                    help='specify the time windows for the peak ( minute). Default value is 0.1 ', required=False)

parser.add_argument('--rt_p_match', dest='rt_p_window_match', action='store', type=float, default=0.4,
                    help='specify the time windows for the matched peptide peak ( minute). Default value is 0.4 ',
                    required=False)

#parser.add_argument('--raw_repo', dest='raw', action='store', help='specify the raw file repository ', required=True)

parser.add_argument('--output_folder', dest='loc_out', action='store', default='', help='specify the folder output',
                    required=True)

args = parser.parse_args()

print 'Matching between run module (mbr)'

map_name =  pd.read_csv(args.map_file,sep="\t", header=None )

res_state = moff_mbr.run_mbr(args,map_name)
if res_state == -1:
    exit('An error is occurred during the writing of the mbr file')
folder = os.path.join( args.loc_out ,'mbr_output')

# os.chdir(folder)
print 'Apex module '
for file in glob.glob(folder + os.sep + "*.txt"):
    file_name = file
    tol = args.toll
    h_rt_w = args.rt_window
    s_w = args.rt_p_window
    s_w_match = args.rt_p_window_match
    #loc_raw = args.raw
    loc_output = args.loc_out
    moff.run_apex(file_name, tol, h_rt_w, s_w, s_w_match, map_name, loc_output)
print '-------------     -------------'  # moff.run_apex(args)
