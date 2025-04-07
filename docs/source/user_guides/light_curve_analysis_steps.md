(Note: This content is taken directly from an old Slack post. It's possible context is missing, the contents are out of date, or the user who posted it was inexperienced with the code.)

# Light Curve Analysis Steps

These are the main steps in general. You may have additional steps depending on your event. [Stela: I got this from Dave]

*(In italic I put things that are optional, depending on your data. If you got the data sets from Ian, you can skip step 2,3 and 4)*

1. Collect data sets

2. *Remove MOA data points with large error bars (typically > 2000), high seeing (typically > 5 pixels) and data more than a few years from the magnification. But, check the error bars and seeing for critical parts of the light curve first, so that you don't remove critical data.*

3. *Run fortran detrending code on MOA R-band data.*

4. *Run JD -> HJD conversion using MOA perl script on detrended MOA R-band data.*

5. Prepare parameter file that indicates data sets, error bar modification parameters, and finite source calculation parameters (particularly grid ratio).

6. Do scan runs to search for candidate initial conditions.

7. Run ~20 modeling runs based on best fit scan runs.

8. Possibly do higher density scan runs in the vicinity of good fits, possibly with a variety of t_* values.

9. Determine if you should include higher order effects like microlensing parallax or orbital motion, and check for possible degenerate solutions.

10. Do additional modeling runs if needed.

11. Once you have the best fit, renormalize the error bars and throw out outliers once you are sure you have a good fit. Be careful not to throw out good data that disagrees with your model. [See below for more details]

12. Rerun model with correct error bars to find the best fits, being sure to sample all degenerate models. [See below for more details]

13. Run code with MCMC call. This is normally a 2-step procedure to ensure that the posterior distribution is well sampled. [See below for more details]

14. Download OGLE-III photometry map, and find the red clump centroid. Use Nataf+16 to estimate the extinction. [See below for more details]

15. Determine the source magnitude and color and get the extinction correction magnitude and color to determine theta_*, mu_rel, and theta_E.

16. Run Galactic model code and then code to sum over light curve MCMC output to get distribution of detected parameters.

## I got the best-fit model! Now what? (Clarification regarding Steps 11 and 12 from 'Light Curve Analysis Steps')

Nice you got here. So you did a few scans, and later a few fittings using `DSEEK` and got some nice fit-models. It is time to move for the MCMC steps. However, before getting so excited, there is one more thing to do. [This file was done copying part of the conversion between Clément & Stela, and Stela & Greg]

### Refining your best fit-model:

1. Remove outliers data points (typically, data points where chi^2 > 16) - but check to see if any data points are removed from the vicinity of the planetary signal.

2. Renormalize your error bars by change your fudge factor. ( I recommend you calculate the chi^2/dof for each data set separately, and use the square root of this value `sqrt(chi^2/dof)` multiplied by the older fudge factor, to get a new fudge factor for each data set.

3. Re run another `DSEEK` with these new data sets and these new fudge factors in the parameter file.

Okay, now you have really your best fit model. And you are ready for starting your MCMC.

## MCMC in 2-steps (Clarification regarding Step 13 from 'Light Curve Analysis Steps')

### First MCMC Run
1. Create an MCMC without trying to optimize the exploration of the parameter space. *Notice that the coordinate of your event may be required or not depending on the version of Dave's code you are using.*

```
MB20135
p1d_
no limb
18 00 19.47 -28 08 56.0
0 mcmc_p1d_1.dat
SET EPS     1.e-5
SET ERR     2.0
DSEEK     200000
EXIT
```

(Example in file p1d_5.in)
Run it.

### Second MCMC Run
2. You can use the `mcmc_p1d_1.dat` (whatever name you called it) file to optimize the exploration, i.e., you use the covariance on the posterior to diagonalize the covariance matrix. The key line in the second file is:

```
DSEEK     500000 2.3 0.7 58000 60241
```

It means 500000 trials, 2.3 and 0.7 are parameters for the proposal, 58000 is the number of samples you keep to use to compute the covariance, 60241 is the total number of line in the mcmc output file from step 1. So the sub-steps are:

2.a Change your `DSEEK` line
2.b Before running it, copy/paste the file `mcmc_p1d_1.dat` with the name `mcmc_p1d_12.dat` because the run can append the MCMC samples at the end of the mcmc file you specify on the line:

```
0 mcmc_p1d_12.dat
```

(Example in files p1d_6.in . *I attached one without the coordinates, depending on your version of Dave's code, you should keep your coordinates there*)

#### Observations:
- `0 mcmc_p1d_12.dat` Dave said I could use 1 instead of 0 in my [Stela's] case to make it go faster. (TODO: figure out what this `1` does instead of `0`)

## How can I do my CMD plot? (Clarification regarding Step 14 from 'Light Curve Analysis Steps')

[This is a draft, we can write it better]

For the Color Magnitude Diagram you need:

1. An OGLE map with your target in the field
2. The color and brightness of your source star
3. The color and brightness of the red clump centroid

### OGLE Map:
You can find which map you want to download, by typing the coordinates events on this website.
[OGLE Field Finder](http://ogle.astrouw.edu.pl/radec2field.html)

The maps are here:
[Index of /ogle/ogle3/maps](http://www.astrouw.edu.pl/ogle/ogle3/maps/)

I will suggest to plot your target in the sky to make sure yours is not in the corner. (Mine was, so I combined two maps)

### Source Star:
For the source star, you should get the source flux of your best model, and convert it to magnitude. For MOA, you should have Ian's calibration values to converted it.
(I have a code that does it. I can share with you.)

### Red Clump Centroid:
You should calculate the color and magnitude of the red clump centroid. (I used a code from Dave which calculates it for stars within a distant radius)

### For the CMD plot
- Plot the magnitude I vs the color V-I .
- Do not include all the stars from the ogle map, but just ones within the radius used for the calculation of the red clump centroid.
