Update 09/17/2016


Vincent and Briggs cases were done by Young-Kwang Choi and Seungnam Seo
Choi and Seo provided all input files and post-processing files

folders
1) 01_eddy_breaking
   regular wave breaking case, eddy viscosity breaking

2) 02_shock_breaking
   regular wave breaking case, shock capturing

3) 03_nobreaking

postprocessing:

1) use Output_eta_save_v2.m to make .mat. Actually this process is not necessary but I think Young-Kwang try to use it for saving time.  .mat are saved in the same folder of results

2) use anal_brk_vis.m to compare with data along section 4. 

3) use plot_uvmean_simplefig.m to plot wave-averaged velocity

NOTE: 

1) The domain sizes used for breaking cases and non-breaking case are different since for breaking case, Delta_WK =2 is used for a larger wave height. The location of shoal is moved 5m forward to the shore. 

2) We recognized that the results using different time periods are slightly different. For example, for the viscosity breaking case, the analysis was conducted from 9.1s to 9.1+36.4s to get the best result. For the shock-capturing case, however, this period resulted in the over-prediction of wave height at the center region. So the time period chosen for shock-capturing case was from 24.9s to 24.9 + 36.4s.  

3) Choi and Seo did some comparisons between different numerical schemes for the Vincent and Braggs cases. Their results are very interesting.



 