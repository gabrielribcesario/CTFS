<h1 align="left">Continuous Time Fourier Series</h1>

<h1 align="left">Linux</h3>

Run the Bash script "fourier_routine" located in the bin folder. Optional arguments are:

- **-t** or **--period**: Period [s] (> 0.0) of the periodic function. Default is 1.0[s];
- **-d** or **--dc**: DC component of the periodic function. Default is 0.0;
- **-a** or **--amplitude**: Amplitude of the periodic function. Default is 1.0;
- **-p** or **--phase**: Phase shift [rad] of the periodic function. Default is 0.0[rad];
- **-f** or **--function**: Case insensitive name of the periodic function. Options are: SIN or SINE; SQR or SQUARE; SAW or SAWTOOTH; TRI or TRIANGLE. Default is SINE;
- **-e** or **--tolerance**: Tolerance (> 0.0) for the stopping criteria. Default is 0.1.

<h1 align="left">Windows</h3>

Run the Windows PowerShell script "fourier_routine.ps1" located in the bin folder. Optional arguments are the same as in the Linux version, only without the double dash (i.e. "--") for the long options.

<h1 align="left"></h3>

**NOTE 1**: The analysis.py script only animates the 50 first harmonics.

**NOTE 2**: As a safeguard, the Fortran algorithm truncates the series after 1000 coefficients.

![Sawtooth wave Fourier Series approximation gif](./misc/sawtooth.gif)
