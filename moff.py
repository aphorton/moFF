import numpy as np
import glob as glob
import pandas as pd
import os as os
import sys
import subprocess
import shlex
import logging
import argparse
import ConfigParser
import ast
from StringIO import StringIO
from sys import platform as _platform
import pymzml
import bisect
### input###
## - MS2 ID file
## - tol
## - half rt time window in minute
###### output
##  list of intensities..+


def pyMZML_xic_out( name, ppmPrecision, minRT, maxRT ,MZValue ):
    run = pymzml.run.Reader(name, MS1_Precision= ppmPrecision , MSn_Precision = ppmPrecision)
    timeDependentIntensities = []
    for spectrum in run: 
        if  spectrum['ms level'] == 1 and spectrum['scan start time'] > minRT and spectrum['scan start time'] < maxRT:
                        lower_index = bisect.bisect( spectrum.peaks, (float(MZValue  - ppmPrecision*MZValue  ),None  ))
                        upper_index = bisect.bisect( spectrum.peaks, (float(MZValue + ppmPrecision*MZValue )  ,None) )
                        maxI= 0.0
                        for sp  in spectrum.peaks[ lower_index : upper_index ] :
                            if sp[1] > maxI:
                             maxI= sp[1]
                        if maxI > 0:
                                timeDependentIntensities.append( [ spectrum['scan start time'], maxI ])	
        
    if len(timeDependentIntensities) > 5 :
        return (pd.DataFrame( timeDependentIntensities,  columns=['rt', 'intensity']) ,1  )
    else :
        return (pd.DataFrame( timeDependentIntensities, columns=['rt', 'intensity']) ,-1 )


def check_columns_name(col_list, col_must_have):
    for c_name in col_must_have:
        if not (c_name in col_list):
            # fail
            return 1
    # succes
    return 0


def run_apex(file_name, tol, h_rt_w, s_w, s_w_match, map_name, loc_output):
    # file_name=  args.name
    # tol= args.toll
    # h_rt_w = args.rt_window
    # s_w= args.rt_p_window
    # s_w_match= args.rt_p_window_match
    # loc_raw = args.raw
    # loc_output = args.loc_out
    # OS detect
    flag_windows = False
    if _platform == "linux" or _platform == "linux2" or _platform =='darwin' :
        flag_windows = False
    else:
        if _platform == "win32":
            flag_windows = True


        # flag_for matching
    mbr_flag = 0
    config = ConfigParser.RawConfigParser()
    ## get the  running path of moff
    moff_path= os.path.dirname( sys.argv[0])

    # it s always placed in same folder of moff.py

    config.read(os.path.join( moff_path, 'moff_setting.properties'))

    # case of moff_all more than one subfolderi
    name = os.path.basename(file_name).split('.')[0]
    if '_match' in name:
        ## in case of mbr , here i dont have evaluate the flag mbr
        start = name.find('_match')
        # extract the name of the file
        name = name[0:start]

    log = logging.getLogger('moFF apex module')
    log.setLevel(logging.INFO)
    if loc_output != '':
        if not (os.path.isdir(loc_output)):
            os.makedirs(loc_output)
            print "created output folder : ", loc_output

        ## outputname : name of the output
        ## it should be ok also in linux
        outputname = os.path.join(loc_output, name + "_moff_result.txt")
        fh = logging.FileHandler( os.path.join(loc_output,name + '__apex.log' ), mode='w')
    else:
        outputname = name + "_moff_result.txt"
        fh = logging.FileHandler( os.path.join(name + '__apex.log' ), mode='w')

    fh.setLevel(logging.INFO)
    log.addHandler(fh)
    
    loc =  str(map_name[map_name[1].str.contains(str(name))][0].values[0])
    flag_mzml= False
    if ('MZML' in loc.upper()):
        flag_mzml=True 
    if os.path.isfile(loc):
        log.info('-- raw file detected --')
    else:
        exit('ERROR:' + loc + ' wrong path or wrong file name  for the raw data')
    ## detect OS

    print 'moff Input file: %s  XIC_tol %s XIC_win %4.4f moff_rtWin_peak %4.4f ' % (file_name, tol, h_rt_w, s_w)
    print 'RAW file  :  %s' % (loc)
    log.info('moff Input file: %s  XIC_tol %s XIC_win %4.4f moff_rtWin_peak %4.4f ', file_name, tol, h_rt_w, s_w)
    log.info('Output_file in :  %s', outputname)
    log.info('RAW file and its location :  %s', loc)
    #print file_name
    ##read data from file
    data_ms2 = pd.read_csv(file_name, sep="\t", header=0)
    if check_columns_name(data_ms2.columns.tolist(), ast.literal_eval(config.get('moFF', 'col_must_have_x'))) == 1:
        exit('ERROR minimal field requested are missing or wrong')


    index_offset = data_ms2.columns.shape[0] - 1

    data_ms2["intensity"] = -1
    data_ms2["rt_peak"] = -1
    data_ms2["lwhm"] = -1
    data_ms2["rwhm"] = -1
    data_ms2["5p_noise"] = -1
    data_ms2["10p_noise"] = -1
    data_ms2["SNR"] = -1
    data_ms2["log_L_R"] = -1
    data_ms2["log_int"] = -1
    data_ms2["rt_peak"] = data_ms2["rt_peak"].astype('float64')
    data_ms2['intensity'] = data_ms2['intensity'].astype('float64')
    data_ms2['lwhm'] = data_ms2['lwhm'].astype('float64')
    data_ms2["rwhm"] = data_ms2['rwhm'].astype('float64')
    data_ms2["5p_noise"] = data_ms2['5p_noise'].astype('float64')
    data_ms2["10p_noise"] = data_ms2['10p_noise'].astype('float64')
    data_ms2["SNR"] = data_ms2['SNR'].astype('float64')
    data_ms2["log_L_R"] = data_ms2['log_L_R'].astype('float64')
    data_ms2["log_int"] = data_ms2['log_int'].astype('float64')


    ##set mbr_flag
    if 'matched' in data_ms2.columns:
        mbr_flag = 1
        log.info('Apex module has detected mbr peptides')
        log.info('moff_rtWin_peak for matched peptide:   %4.4f ', s_w_match)
    c = 0
    steps = data_ms2.shape[0] / 10
    #print 'Starting apex .........',
    print 'Starting  Apex [          ]',
    print '\b'*12,
    sys.stdout.flush()
    for index_ms2, row in data_ms2.iterrows():
        # log.info('peptide at line: %i',c)
        mz_opt = "-mz=" + str(row['mz'])
	##convert rt to sec to min
        time_w = row['rt'] / 60
        if mbr_flag == 0:
            log.info('peptide at line %i -->  MZ: %4.4f RT: %4.4f ', c, row['mz'], time_w)
            temp_w = s_w
        else:
            log.info('peptide at line %i -->  MZ: %4.4f RT: %4.4f matched(y/n): %i', c, row['mz'], time_w,row['matched'])
            if row['matched'] == 1:
                temp_w = s_w_match
            else:
                temp_w = s_w
        if row['rt'] == -1:
            log.warning('rt not found. Wrong matched peptide in the mbr step line: %i', c)
            c += 1
            continue
        try:
            if flag_mzml:
                # mzml raw file
                data_xic ,status = pyMZML_xic_out(  loc, float(tol / (10 ** 6)) ,   time_w - h_rt_w , time_w + h_rt_w , row['mz']  )
                if status==-1:
                    log.warning("WARNINGS: XIC not retrived line: %i", c)
                    log.warning('MZ: %4.4f RT: %4.4f Mass: %i', row['mz'], row['rt'], index_ms2)
                    c += 1
                    continue
            else:
                if flag_windows :
                    os.path.join('folder_name', 'file_name')
                    args_txic = shlex.split( os.path.join(moff_path ,"txic.exe") +" "+  mz_opt + " -tol=" + str(tol) + " -t " + str(time_w - h_rt_w) + " -t " + str(time_w + h_rt_w) + " " + loc, posix=False)
                else:
                    args_txic = shlex.split("./txic " + mz_opt + " -tol=" + str(tol) + " -t " + str(time_w - h_rt_w) + " -t " + str(time_w + h_rt_w) + " " + loc)
                    p = subprocess.Popen(args_txic, stdout=subprocess.PIPE)
                    output, err = p.communicate()
                    data_xic = pd.read_csv(StringIO(output.strip()), sep=' ', names=['rt', 'intensity'], header=0)
            
            
            if data_xic[(data_xic['rt'] > (time_w - temp_w)) & (data_xic['rt'] < (time_w + temp_w))].shape[0] >= 1:
                ind_v = data_xic.index
                pp = data_xic[data_xic["intensity"] ==  data_xic[(data_xic['rt'] > (time_w - temp_w)) & (data_xic['rt'] < (time_w + temp_w))]['intensity'].max()].index
                pos_p = ind_v[pp]
                if pos_p.values.shape[0] > 1:
                    log.warning(" RT gap for the time windows searched. Probably the ppm values is too small %i", c)
                    continue
                val_max = data_xic.ix[pos_p, 1].values
            # log.info(data_xic[(data_xic['rt']>   (time_w -1)   ) & ( data_xic['rt']<  ( time_w + 1   )    )]   )
            else:
                log.info("LW_BOUND window  %4.4f", time_w - temp_w)
                log.info("UP_BOUND window %4.4f", time_w + temp_w)
                log.info(data_xic[(data_xic['rt'] > (time_w - +0.60)) & (data_xic['rt'] < (time_w + 0.60))])
                log.info("WARNINGS: moff_rtWin_peak is not enough to detect the max peak line : %i", c)
                log.info('MZ: %4.4f RT: %4.4f Mass: %i', row['mz'], row['rt'], index_ms2)
                c += 1
                continue
            pnoise_5 = np.percentile(data_xic[(data_xic['rt'] > (time_w - (h_rt_w / 2))) & (data_xic['rt'] < (time_w + (h_rt_w / 2)))]['intensity'], 5)
            pnoise_10 = np.percentile(data_xic[(data_xic['rt'] > (time_w - (h_rt_w / 2))) & (data_xic['rt'] < (time_w + (h_rt_w / 2)))]['intensity'], 10)
        except (IndexError, ValueError, TypeError):
            log.warning(" size is not enough to detect the max peak line : %i", c)
            log.info('MZ: %4.4f RT: %4.4f index: %i', row['mz'], row['rt'], index_ms2)
            continue
            c += 1
        except pd.parser.CParserError:
            log.warning("WARNINGS: XIC not retrived line: %i", c)
            log.warning('MZ: %4.4f RT: %4.4f Mass: %i', row['mz'], row['rt'], index_ms2)
            c += 1
            continue
        else:
            # log.info("Intensisty at pos_p-1 %4.4f",data_xic.ix[(pos_p-1),1].values )
            log_time = [-1, -1]
            c_left = 0
            find_5 = False
            stop = False
            while c_left < (pos_p - 1) and stop != True:
                # print c_left

                if find_5 == False and (data_xic.ix[(pos_p - 1) - c_left, 1].values <= (0.5 * val_max)):
                    find_5 = True
                    # print "LWHM",c_left,data_xic.ix[(pos_p-1)-c_left,1]
                    # log_point[0] = np.log2(val_max)/np.log2(data_xic.ix[(pos_p-1)-c_left,1])
                    log_time[0] = data_xic.ix[(pos_p - 1) - c_left, 0].values * 60
                    stop = True
                c_left += 1
            find_5 = False
            stop = False
            r_left = 0
            while ((pos_p + 1) + r_left < len(data_xic)) and stop != True:
                if find_5 == False and data_xic.ix[(pos_p + 1) + r_left, 1].values <= (0.50 * val_max):
                    find_5 = True
                    # print "RWHM",r_left,data_xic.ix[(pos_p+1)+r_left,1]
                    # log_point[2] = np.log2(val_max) /np.log2(data_xic.ix[(pos_p+1)+r_left,1])
                    log_time[1] = data_xic.ix[(pos_p + 1) + r_left, 0].values * 60
                    stop = True
                r_left += 1

            data_ms2.ix[index_ms2, (index_offset + 1)] = val_max
            data_ms2.ix[index_ms2, (index_offset + 2)] = data_xic.ix[pos_p, 0].values * 60
            data_ms2.ix[index_ms2, (index_offset + 3)] = log_time[0]
            data_ms2.ix[index_ms2, (index_offset + 4)] = log_time[1]
            data_ms2.ix[index_ms2, (index_offset + 5)] = pnoise_5
            data_ms2.ix[index_ms2, (index_offset + 6)] = pnoise_10
            # conpute log_L_R SNR and log intensities
            if (pnoise_5 == 0 and pnoise_10 > 0):
                data_ms2.ix[index_ms2, (index_offset + 7)] = 20 * np.log10(data_xic.ix[pos_p, 1].values / pnoise_10)
            else:
                data_ms2.ix[index_ms2, (index_offset + 7)] = 20 * np.log10(data_xic.ix[pos_p, 1].values / pnoise_5)
            #### WARNING  time - log_time 0 / time -log_time 1
            data_ms2.ix[index_ms2, (index_offset + 8)] = np.log2(
                abs(data_ms2.ix[index_ms2, index_offset + 2] - log_time[0]) / abs(
                    data_ms2.ix[index_ms2, index_offset + 2] - log_time[1]))
            data_ms2.ix[index_ms2, (index_offset + 9)] = np.log2(val_max)
	    ## bar update
	    if c % steps == 0 :
		print '\b.',
        	sys.stdout.flush()  
            c += 1
    ## save  result i
    print '\b] Apex end'
    print 'Writing result in %s' % (outputname)
    data_ms2.to_csv(path_or_buf=outputname, sep="\t", header=True, index=False)
    fh.close()
    log.removeHandler(fh)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='moFF input parameter')

    parser.add_argument('--map_file', dest='map_file', action='store',
                    help='specify a map file that contains input files   and raw file ', required=True) 
   
    parser.add_argument('--tol', dest='toll', action='store', type=float,
                        help='specify the tollerance parameter in ppm', required=True)

    parser.add_argument('--rt_w', dest='rt_window', action='store', type=float, default=3,
                        help='specify rt window for xic (minute). Default value is 3 min', required=False)

    parser.add_argument('--rt_p', dest='rt_p_window', action='store', type=float, default=0.1,
                        help='specify the time windows for the peak ( minute). Default value is 0.1 ', required=False)

    parser.add_argument('--rt_p_match', dest='rt_p_window_match', action='store', type=float, default=0.4,
                        help='specify the time windows for the matched peptide peak ( minute). Default value is 0.4 ',
                        required=False)


    parser.add_argument('--output_folder', dest='loc_out', action='store', default='', help='specify the folder output',
                        required=False)

    args = parser.parse_args()
    map_name =  pd.read_csv(args.map_file,sep="\t", header=None )
    print 'Apex module '
    
    for file_name  in  map_name.ix[:,1].tolist():
    	
    	tol = args.toll
    	h_rt_w = args.rt_window
    	s_w = args.rt_p_window
    	s_w_match = args.rt_p_window_match
    	loc_output = args.loc_out
    	run_apex(file_name, tol, h_rt_w, s_w, s_w_match, map_name, loc_output)

