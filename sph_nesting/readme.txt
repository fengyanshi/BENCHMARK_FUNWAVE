Updated /08/29/2016/  using Version 3.0

1) got to /src/, compile the code with an option (in Makefile) 
             FLAG_11 = -DSPHERICAL_IJ_STATION
   This allows the model output using (i,j) as stations specified in station.txt

2) go to /grid_A/, run the model

2) copy  your results: sta_0001, sta_0002, sta_0003 into /make_nest_file/
   use convert.f to get coupling.txt

3) goto /src/ again, use option FLAG_8 = -DCOUPLING to recompile the code (make clean first)

4) goto /grid_B/, run the model

5) goto /postprocessing/, use plot_ab.m to plot the figure
