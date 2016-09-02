
A) plunging case with shock capturing breaking
1) complile code with
            FLAG_3 = -DSAMPLES
            FLAG_4 = -DCARTESIAN
2) go to /pluging_shock/
   run the model
3) go to /postprocessing/ 
   use ht_setup_plunging_shock.m to plot result. 

B) plunging case with eddy viscosity breaking

   use ht_setup_plunging_visc.m to plot result.

   To compare shock capturing vs. viscosity, use comparison_visc_shock_pl.m


Other cases use the same procedure, but use different matlab files to plot.
