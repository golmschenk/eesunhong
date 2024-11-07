# How to run eesunhong

This document shows how to perform a basic run using `eesunhong`. For more details explaining what actually happened during this run, see the other user guides in the documentation.

First, you will need to gather some data files to run `eesunhong` on. There are several run configurations along with their expected results in the end to end tests. These can be found in the repository at `eesunhong/tests/end_to_end_tests`. To obtain a basic case, you can also just [click here](https://download-directory.github.io/?url=https%3A%2F%2Fgithub.com%2Fgolmschenk%2Feesunhong%2Ftree%2Fmain%2Ftests%2Fend_to_end_tests%2Fbinary_lens%2Fbasic_dseek%2Frun_directory). From this directory, with `eesunhong` already installed via `pip`, you can run:
```shell
eesunhong_main < run_1.in >& run_1.out
```

This will stream the input from the `run_1.in` file as input to `eesunhong`, and send the output stream to `run_1.out`. During this process, `eesunhong` will read the configuration file and any available light curve files, and will create several files, including a light curve fitting and residual.