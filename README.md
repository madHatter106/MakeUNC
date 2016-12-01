Computes Monte-Carlo based uncertainty in ocean color remote sensing radiance to assess the impact of satellite-borne sensor noise and its propagation from at-sensor measurement, through atmospheric correction, to at-ocean-surface remote sensing radiance estimation and from there to derived ocean color products (e.g. as algal chlorophyll) that are vital to such concepts as climate change, fisheries yield, harmful algal blooms. A Monte-Carlo simulation boils down to reprocessing the same [scene](http://oceancolor.gsfc.nasa.gov/cgi/l3) about 1000 times to estimate the uncertainty resulting from of sensor noise and comparing the results to a perturbation-free baseline version of the scene.

Content:
The follwoing code can be called from the commandline 
  * RunL2genMC calls on [NASA's Ocean Biology Group's](http://oceancolor.gsfc.nasa.gov/cms/) remote sensing processing software l2gen, which is part of a very large code base written in C,  to run iterations of a Monte-Carlo process. I added a sensor-based noise model depending on a variety of inputs in l2gen.  One issue is that l2gen computes remote sensing data one pixel at a time. This makes for a very long wait when a simulation requires that a scene be reprocessed several hundred times. To mitigate the following, I added the option of running several processes in parallel by spawning concurrent l2gen processes.
The output is contained in netCDF4 files.
  * MakeUNC builds Rrs and optionally product uncertainty using data resulting from Monte-Carlo simulation runs  by statistically aggregating the results and adding them as uncertainty products to the baseline dataset, also contained in a netCDF4 file. 
