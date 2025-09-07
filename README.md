<h1 align="left">Continuous-Time Fourier Series</h1>

Calculation of the Fourier coefficients of continuous-time periodic functions using Fortran.

<h1 align="left">Instructions</h1>

Run the Bash script "fourier_routine" located in the bin folder. Use a single dash for specifying a short option and a double dash for long option (e.g. "-o" for short, "--option" for long). Optional arguments are:

- **-h** or **--help**: Displays the optional arguments of the script;
- **-t** or **--period**: Period [s] (> 0.0) of the periodic function. Default is 1.0[s];
- **-d** or **--dc**: DC component of the periodic function. Default is 0.0;
- **-a** or **--amplitude**: Amplitude of the periodic function. Default is 1.0;
- **-p** or **--phase**: Phase shift [rad] of the periodic function. Default is 0.0[rad];
- **-f** or **--function**: Case insensitive name of the periodic function. Options are: SIN or SINE; SQR or SQUARE; SAW or SAWTOOTH; TRI or TRIANGLE. Default is SINE;
- **-e** or **--tolerance**: Tolerance (> 0.0) for the stopping criteria. Default is 0.1.

<h1 align="left"></h1></br>

**NOTE 1**: The analysis.py script only animates the 50 first harmonics.

**NOTE 2**: As a safeguard, the Fortran algorithm truncates the series at the 1000th harmonic.

![Sawtooth wave Fourier Series approximation gif](./misc/sawtooth.gif)
