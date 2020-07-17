# Greedy-Sensor-Selection-Algorithm
---
This repository contains Matlab R2013a code to reproduce results for a manuscript entitled __"Effect of Objective Function on Data-Driven Sparse Sensor Optimization"__ published in [arXiv](https://arxiv.org/abs/2007.05377).  
The sparse sensor selection problem is solved by the greedy method.  
To run the program, excute `P_greedy_demo`.  

## Directory  
---
- src: source code is stored  
- work: calculation results are stored (created automatically by running the program)  
- data: __NOAA Optimum Interpolation (OI) Sea Surface Temperature (SST) V2__ data is stored  
  - sst.wkmean.1990-present.nc  
  - lsmask.nc  
NOAA_OI_SST_V2 is provided by the NOAA/OAR/ESRL PSD, Boulder, Colorado, USA, from their Web site at https://www.esrl.noaa.gov/psd/.  
Due to GitHub file size limitations, a dataset is linked online: [NOAA Optimum Interpolation (OI) Sea Surface Temperature (SST) V2](https://www.esrl.noaa.gov/psd/data/gridded/data.noaa.oisst.v2.html)  

## Code  
---
### Main program  
- P_greedy_demo.m  

### Function  
#### Preprocessing  
- F_pre_read_NOAA_SST.m  
- F_pre_SVD_NOAA_SST.m  
- F_pre_truncatedSVD.m  

#### Sensor selection  
- F_sensor_random.m  
- F_sensor_convex.m  
  - F_sensor_convex_sub.m  
    - F_sensor_convex_approxnt_vec.m  
    - F_sensor_convex_approxnt.m  
    - F_sensor_convex_loc.m  
    - F_sensor_convex_locr.m  
- F_sensor_QR.m  
  - F_sensor_QR_pivot.m  
- F_sensor_DG.m  
  - F_sensor_DG_r.m  
  - F_sensor_DG_p.m  
- F_sensor_QD.m  
- F_sensor_AG.m  
  - F_sensor_AG_calc_trace.m  
- F_sensor_EG.m  
  - F_sensor_EG_calc_eigen.m  

#### Calculation
- F_calc_det.m  
- F_calc_trace.m  
- F_calc_eigen.m  
- F_calc_sensormatrix.  
- F_calc_error.m  
  - F_calc_reconst.m  
  - F_calc_reconst_error.m  

#### Data organization  
- F_data_ave1.m  
- F_data_ave2.m  
- F_data_arrange1.m  
- F_data_arrange2.m  
- F_data_arrange3.m  
- F_data_normalize.m  

#### Mapping
- F_map_original.m  
	- F_map_videowriter.m  
		- F_map_plot_sensors_forvideo.m  
- F_map_reconst.m  
	- F_map_plot_sensors.m  

## License  
---
[MIT-License](https://github.com/KumiNakai/Greedy-Sensor-Selection-Algorithm/blob/master/LICENSE)

## Author
---
Kumi Nakai  
[Experimental Aerodynamics Laboratory](http://www.aero.mech.tohoku.ac.jp/eng/)  
Department of Aerospace Engineering, Tohoku University  
Sendai, JAPAN  
E-mail: nakai@aero.mech.tohoku.ac.jp
