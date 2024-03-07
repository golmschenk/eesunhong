# File explanations

To run `eesunhong` fitting code, you will need a directory containing the data relevant to the event you are analyzing. The existing convention is to have a directory named after your event (e.g. MB110039). This directory will contain the following file structure.

## File structure

### Inputs
Files required for initial fitting run
1. Instrument parameter file - par{event_name}: A file containing a description of the instruments used to make observations. The naming will be par followed by the event name. For example, for an event called MB110039, the file will be named parMB110039.
2. Light curve files - lc{event_name}.{instrument_suffix}: The light curves files to the model to. Each light curve you want to fit requires that the instrument used to record the light curve is listed in the instrument parameter file. The instrument suffix (sfx) in that file must be extension used by the instrument in the instrument parameter file. For example, for the event MB110039 and the instrument suffix moa2r, the light curve file would be lcMB110039.moa2r.
3. Run input files - {run_prefix}_{run_number}.in: Files defining what fittings to run. Run input files with an incremented run_number will be created by earlier run files. However, the first run file for a sequence needs to be manually created.

### Outputs
Files resulting from runs
1. Run output files - {run_prefix}_{run_number}.out: The output of the fitting code for a specific run.
2. Model fit files - fit.lc_{run_prefix}_{run_number}: Files describing the model fit for a specific run.
3. Model fit residual files - resid.{run_prefix}_{run_number}: Files describing the model fit for a specific run.
If you are running an MCMC, you will also have:
4. mcmc_name.dat: [//]: # TODO 

## File content explanations
1. Instrument parameter file - par{event_name}
(TODO)

2. Light curve files - lc{event_name}.{instrument_suffix}
Light curve files columns use the convention HJD, Flux, Flux_err OR HJD, Magnitude, Magnitude_err. What decides one or the other convention is the row number corresponding to the instrument in the instrument parameter file, with the following rule:
If the column jclr of the instrument is (15 <= jclr <= 29) or (40 <= jclr <= 49), then the observable is the magnitude with a zero point equal to 21 (i.e., m = 21 - 2.5 * log(flux)).
If the column jclr of the instrument is (30 <= jclr <= 39) or (50 <= jclr <= 59), then the observable is the flux.
If the column jclr of the instrument is 9 <= jclr <= 14, it is a little special because the observable is the magnitude with a zero point equal to 0. In practice, you will not want to use this convention, unless you have a good reason to do so. Also, please keep in mind that if you use those line to define a new instrument (or replace it), then the code will expect more columns than 3.

3. Run input files - {run_prefix}_{run_number}.in
This file contains the description of a fitting run to be run. The first line of the file is a comment. Afterward, the file has two main parts: 
The parameter settings. This is the top part of the file with several parameters listed, one on each line.
The commands the fitting script should run. This is lower part of the file.

The parameter column values are:
The parameter index, to be referenced by the later commands.
The name of the parameter.
The initial value for the parameter.
(optional) A minimum value for the parameter.
(optional) A maximum value for the parameter.

The parameter meanings are:
1/t_E: inverse of Einstein radius crossing time = 2/t_hat.
t0: time of umin w.r.t. CM
umin: umin w.r.t. CM
sep: lens separation at time t_fix in units of Einstein radius
theta: angle between velocity vector and separation vector at time t_fix
eps1: mass fraction of mass 1
1/Tbin: inverse period of binary orbit
v_sep: velocity in separation direction.
Tstar: stellar radius in days
t_fix: A reference time to be used by the other fit parameters (should not be a fit parameter itself)
piEr: (related to parallax fitting)
piEtheta: (related to parallax fitting)

(TODO - write about the command section)


1. Run output files - {run_prefix}_{run_number}.out
(TODO)

2. Model fit files - fit.lc_{run_prefix}_{run_number}
This file has three sections:
the best model parameters for the corresponding run;
the source flux fs (A0) and blending fb (A2) for each dataset;
A0 and A2 provide a linear transformation from the raw light curve input data to the final magnification values for that dataset. Note that these values are different for each dataset, and are determined based on the fit model. That is, we don’t know the true magnification of a target in advance. Based on the current model being fit, we calculate the best linear transformation for the light curve to fit the model. This is then assumed to provide the magnification (for that given fit). This process aligns the data from different sources.
the light curve (see below).
The columns are given below.
Date (same convention as in input files)
Magnification (model) A(t)
Flux f(t)=A(t)fs + fb for dataset #1
... Idem for each of the N datasets
Source position xs (column N+4)
Source position ys (column N+5)

3. Model fit residual files - resid.{run_prefix}_{run_number}
This file has three sections:
the best model parameters for the corresponding run;
the source flux fs (A0) and blending fb (A2) for each dataset;
A0 and A2 provide a linear transformation from the raw light curve input data to the final magnification values for that dataset. Note that these values are different for each dataset, and are determined based on the fit model. That is, we don’t know the true magnification of a target in advance. Based on the current model being fit, we calculate the best linear transformation for the light curve to fit the model. This is then assumed to provide the magnification (for that given fit). This process aligns the data from different sources.
the light curve (see below).
The columns are given below.
Date (same convention as in input files)
Magnification (model) A(t)
Residual in magnification of each data point: A(t)−Adata(t)
Error in magnification σ/fs
Contribution of the data point to the χ2 in magnification
ID of the dataset
Name of the dataset


4. mcmc_name.dat
The name is the name of the input script file. The columns of this file are given below.
#1 chi2 #2 1/tE #3 t0 #4 umin #5 sep #6 theta #7 eps1 #8 1/Tbin #9 v_sep #10 Tstar #11 t_fix #12 A0 of observatory 1 (source flux) #13 A1 of observatory 1 (blending flux)
There are as many columns #12 and #13 as the total number of observatories.

