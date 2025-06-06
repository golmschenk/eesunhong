(Note: This content is taken directly from an old Slack post. It's possible context is missing, the contents are out of date, or the user who posted it was inexperienced with the code.)

# Explanation of the data files and file structure

To run David Bennett's fitting code, you will need a directory containing the data relevant to the event you are analyzing. The existing convention is to have a directory named after your event (e.g., `MO110039`). This directory will contain the following file structure.

## File structure

### Files required for initial fitting run

- **Instrument parameter file - `par{event_name}`**:  
  A file containing a description of the instruments used to make observations. The naming will be `par` followed by the event name. For example, for an event called `MO110039`, the file will be named `parMO110039`.

- **Light curve files - `lc{event_name}.{instrument_suffix}`**:  
  The light curve files to fit the model to. Each light curve you want to fit requires that the instrument used to record the light curve is listed in the instrument parameter file. The instrument suffix (`sfx`) in that file must be the extension used by the instrument in the instrument parameter file. For example, for the event `MO110039` and the instrument suffix `moa2`, the light curve file would be `lcMO110039_moa2`.

- **Run input files - `{run_prefix}{run_number}.in`**:  
  Files defining what fittings to run. Run input files with an incremented `run_number` will be created by earlier run files. However, the first run file for a sequence needs to be manually created.

### Files resulting from runs

- **Run output files - `{run_prefix}{run_number}.out`**:  
  The output of the fitting code for a specific run.

- **Model fit files - `fit.lc.{run_prefix}{run_number}`**:  
  Files describing the model fit for a specific run.

- **Model fit residual files - `resid.{run_prefix}{run_number}`**:  
  Files describing the model fit residual for a specific run.

- **(TODO) - `mcmc_name.dat` (TODO)**

## File content explanation

### Instrument parameter file - `par{event_name}`

### Light curve files - `lc{event_name}.{instrument_suffix}`

Light curve files columns use the convention `HJD`, `Flux`, `Flux_err`, OR `HJD`, `Magnitude`, `Magnitude_err`. What decides one or the other convention is the row number corresponding to the instrument in the instrument parameter file, with the following rule:

- If the column `jclr` of the instrument is `(15 <= jclr <= 29)` or `(40 <= jclr <= 49)`, then the observable is the magnitude with a zero point equal to 21 (i.e., `m = 21 - 2.5 * log(Flux)`).

- If the column `jclr` of the instrument is `(30 <= jclr <= 39)` or `(50 <= jclr <= 59)`, then the observable is the flux.

- If the column `jclr` of the instrument is `9 <= jclr <= 14`, it is a little special because the observable is the magnitude with a zero point equal to 0. In practice, you will not want to use this convention, unless you have a good reason to do so. Also, please keep in mind that if you use those line to define a new instrument (or replace it), then the code will expect more columns than 3.

### Run input files - `{run_prefix}_{run_number}.in`

This file contains the description of a fitting run to be run. The first line of the file is a comment. Afterward, the file has two main parts:

1. The parameter settings. This is the top part of the file with several parameters listed, one on each line.
2. The commands the fitting script should run. This is lower part of the file.

The parameter column values are:

- 1. The parameter index, to be referenced by the later commands.
- 2. The name of the parameter.
- 3. The initial value for the parameter.
- 4. (optional) A minimum value for the parameter.
- 5. (optional) A maximum value for the parameter.

The parameter names are:
- `t_E`: inverse of Einstein radius crossing time = 2/t_Earth.
- `u_0`: time of minimum w.r.t. CM
- `rho`: size of source in units of t_E
- `sep`: lens separation at time t_fix in units of Einstein radius
- `theta`: angle between velocity vector and separation vector at time t_fix
- `eps1`: mass fraction of mass 1
- `1/Tbin`: inverse period of binary orbit
- `v_sep`: velocity in separation direction.
- `Tstar`: stellar radius in days
- `t_fix`: A reference time to be used by the other fit parameters (should not be a fit parameter itself)
- `piEN`: (related to parallax fitting)
- `piEtheta`: (related to parallax fitting)

After specifying the model parameters, you may provide commands to the fitting scripts in the control section. Below is an example of the control section.
```
MB14360
run_
no limb
17 55 1.513 -30 51 2.09
0 run_1.dat
SET EPS 1.e-5
SET ERR 2.0
DSEEK 3000
EXIT
```
- `MB14360`: The first line on the control section should contain the name of the event as written in the parameter file name (`par{event_name}`).
- `run_`: `run_prefix` name. Should be the first part of the name of this file. 
- `no limb`: Tells eesunhong to assume the limb darkening constants provided in the instrument parameter file.
- `17 55 1.513 -30 51 2.09`: Coordinates of the event, precessed to the time of peak magnification.
- `0 run_1.dat`: This line is composed of two unrelated inputs. The numerical part may either be 0 or 1. If 1, the code will calculate and store the integration grid around the Einstein ring prior to doing the calculations. If 0, eesunhong will calculate this grid on the fly. So, setting it to 1 can save some time, but will not be helpful if orbital motion is included due to the changing geometry of the system. The second part is the name of the desired output file for an MCMC run. Only provide this name if doing an MCMC, or you will create an empty file.
- `SET EPS 1.e-5`: TODO
- `SET ERR 2.0`: TODO
- `DSEEK 3000`: DSEEK is one of the fitting commands. It performs $\chi^2$ minimization using the Metropolis-Hastings algorithm. Other fitting commands available are SCAN, which can be used to perform an initial condition grid search to find starting points for DSEEK, and OSEEK, which is used to launch an MCMC. The usage of these commands is described in more detail in [light curve analysis steps](light_curve_analysis_steps.md).
- `EXIT`: This command should always be called at the end of your control file. It prints the best fit lightcurve and residuals in the `fit.lc_{run_prefix}_{run_number}` and the `resid.{run_prefix}{run_number}` files described below.


### Run output files - `{run_prefix}_{run_number}.out`

(TODO)

### Model fit files - `fit.lc_{run_prefix}_{run_number}`

This file has three sections:

- the best model parameters for the corresponding run;
- the source flux fs (`A0`) and blending fb (`Az`) for each dataset;

    - `A0` and `Az` provide a linear transformation from the raw light curve input data to the final magnification values for that dataset. Note that these values are different for each dataset, and are determined based on the fit model. That is, we don't know the true magnification of a target in advance. Based on the current model being fit, we calculate the best linear transformation for the light curve to fit the model. This is then assumed to provide the magnification (for that given fit). This process aligns the data from different sources.

- the light curve (see below).

The columns are given below.

1. Date (same convention as in input files)
2. Magnification (model) A(t)
3. Residual in magnification of each data point: A(t)−Adata(t)
4. Error in magnification σ/fs
5. Contribution of the data point to the χ2 in magnification
6. ID of the dataset
7. Name of the dataset

### Model fit residual files - `resid.{run_prefix}{run_number}`

This file has three sections:

- the best model parameters for the corresponding run;
- the source flux fs (`A0`) and blending fb (`Az`) for each dataset;

    - `A0` and `Az` provide a linear transformation from the raw light curve input data to the final magnification values for that dataset. Note that these values are different for each dataset, and are determined based on the fit model. That is, we don't know the true magnification of a target in advance. Based on the current model being fit, we calculate the best linear transformation for the light curve to fit the model. This is then assumed to provide the magnification (for that given fit). This process aligns the data from different sources.

- the light curve (see below).

The columns are given below.

1. Date (same convention as in input files)
2. Magnification (model) A(t)
3. Residual in magnification of each data point: A(t)−Adata(t)
4. Error in magnification σ/fs
5. Contribution of the data point to the χ2 in magnification
6. ID of the dataset
7. Name of the dataset

### (TODO) - `mcmc_name.dat`

The name is the name of the input script file. The columns of this file are given below.

- `#1 tE #2 1/tE #3 t0 #4 umin #5 sep #6 theta #7 eps1 #8 1/Tbin #9 v_sep #10 Tstar #11 t_fix #12 A0 of observatory 1 (source flux) #13 A1 of observatory 1 (blending flux)`

There are as many columns `#12` and `#13` as the total number of observatories.
