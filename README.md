<h1 align="left">Continuous Time Fourier Series:</h1>

Run the bash script "fourier_routine" located in the bin folder. Optional arguments are:

- **-t** or **--period**: Period [s] of the periodic function. Default is 1.0[s];
- **-d** or **--dc**: DC component of the periodic function. Default is 0.0;
- **-a** or **--amplitude**: Amplitude of the periodic function. Default is 1.0;
- **-p** or **--phase**: Phase shift [rad/s] of the periodic function. Default is 0.0[rad/s];
- **-f** or **--function**: Case insensitive name of the periodic function. Options are: SIN or SINE; SQR or SQUARE; SAW or SAWTOOTH; TRI or TRIANGLE. Default is SINE;
- **-e** or **--tolerance**: Tolerance for the stopping criteria. Default is 0.1.

![](./sawtooth.gif)

**NOTE**: The analysis.py script only animates the 50 first harmonics.