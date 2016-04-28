# moFF #
## A modest Feature Finder (but still robust) to extract apex MS1 intensity directly from Thermo raw file ##

[Introduction](#introduction)
[Requirement](#requirement)
[Sample data](#sample-data)
[Matching Between Runs](#matching-between-runs)
[Apex intensity](#apex-intensity)
[Entire workflow](#entire-workflow)


---

## Introduction  Pipeline Version##

moFF is a python tool that quantifies MS1 intensity peak starting from a list of MS2 idenfied peptide 
moFF works directly and only on Thermo Raw file thanks to a Go library that is able to read the data from Thersmos raw file without any kind of conversion in other formats.

moFF is composed by two stand alone modules :
- *moff_all.py* :  matching between run + Apex intensity;  this will run the entire workflow
- *moff.py* :  apex intensity; this will run aonly the ms2

To run  the entire workflow (mbr and apex ) you should  use  *moff_all.py*

[Top of page](#moff)

----

## Requirement ##

Required python libraries :
- Python 2.7
- pandas  > 0.17.
- numpy > 1.10.0
- argparse > 1.2.1 
- scikit-learn > 0.17

moFF uses *txic*/*txic.exe* to extract the XiC data from the raw files, so  *txic*  must be located in the same folder where you have all moFF scripts.

The txic program is compatibale with  the raw file of all the Orbitrap and triple quadrupole Thermo machines. 
For the moment it does not work with the Thermo Fusion machine.

The input files that contain the list of the MS2 identified peptides (you can use any search engines) must contains the information showed in *moFF_setting.property* for each peptide. The minimun specificic requirements of the input files are:
- tab delimited file
- the header of the input file should contain the following the fields  and columnns names :  
  - 'peptide' : sequence of the peptide
  - 'prot': protein ID 
  - 'rt': retention time of peptide   ( The retention time must be specified in second )
  - 'mz' : mass over charge
  - 'mass' : mass of the peptide
  - 'charge' : charge of the ionized peptide

see the sample input files in the folder *f1_folder* for more information .

[Top of page](#moff)

---


## Entire workflow: Matching between run + Apex ##

use `python moff_all.py -h`
```
        --map_file MAP_FILE  specify a map file that contains input files  and raw file     
        --log_file_name LOG_LABEL a label name to use for the log file  (not mandatory, moFF_mbr_.log is the default name)
        --filt_width W_FILT   width value of the filter k * mean(Dist_Malahobis)  Default val = 1.5
        --out_filt OUT_FLAG   filter outlier in each rt time allignment   Default val =1
        --weight_comb W_COMB  weights for model combination combination : 0 for no weight 1 weighted devised by trein err of the model. Default val =1
        --tol TOLL            specify the tollerance parameter in ppm
        --rt_w RT_WINDOW      specify rt window for xic (minute). Default value is  5  min
        --rt_p RT_P_WINDOW    specify the time windows for the peak ( minute). Default value is 0.1
        --rt_p_match RT_P_WINDOW_MATCH  specify the time windows for the matched peptide peak ( minute). Default value is 0.4
        --output_folder LOC_OUT         specify the folder output (mandatory)
```
`python moff_all.py  --map_file moff_map_file.txt  --output_folder test_moffgui/   --filt_width 1.5  --out_filt 1  --weight_comb 1   --tol 10  --rt_w 5 --rt_p 0.1 --rt_p_match 0.4  `

In  the  --output_folder is stored all the output :
```
	- in the mbr_output folder the output of the mbr module (also the log of the mbr is saved   moFF_mbr_.log )
	- all the moff result; (inputname_moff_result.txt)
	- a log file file for each apex task performed (inputname__apex.log .log )
```
[Top of page](#moff)


---

## Apex Intensity ##

use  `python moff.py -h`
````
  --map_file MAP_FILE  		      specify a map file that contains input files  and raw file
  --tol TOLL                          specify the tollerance parameter in ppm
  --rt_w RT_WINDOW                    specify rt window for xic (minute). Default value is 3 min
  --rt_p RT_P_WINDOW                  specify the time windows for the peak ( minute). Default value is 0.1
  --output_folder LOC_OUT             specifyO the folder output
```
`python moff.mbr  --map_file moff_map_file.txt  --output_folder test_moffgui/     --tol 10  --rt_w 5 --rt_p 0.1 --rt_p_match 0.4 ` 
 
It run the apex module on the input file specified in the map file where also the raw file location are specified.
In the output files moFF just adds the following fields to your origin input file:
- "intensity" intensity, taking the highest peak in the XIC
- "rt_peak" rt of the highest peak
- "lwhm" left width half maximun of the signal in seconds
- "rwhm" right width half maximun of the signal in seconds
- "SNR" signal-to-noise
- "log_L_R" log ratio of lwhm over rwhm (peak shape )
- "log_int" log 2 of the intesity 

It generates a __apex.log file (with same name of input file) that contains  detailesd information for each peak retrieved.
This module determines automaticaly if the input file contains matched peptides or not.

WARNING : the raw file names  MUST be the same of the input file otherwise the script give you an error !


[Top of page](#moff)

---


