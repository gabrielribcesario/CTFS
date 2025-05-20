#!/bin/pwsh

# Parse optional arguments
param(
    [switch]$HELP,

    [Alias("t")]
    [ValidateRange("Positive")]
    [double]$PERIOD=1.0,

    [Alias("d")]
    [double]$DC=0.0,

    [Alias("a")]
    [double]$AMPLITUDE=1.0,

    [Alias("p")]
    [double]$PHASE=0.0,

    [Alias("f")]
    [ValidateSet("SIN", "SINE", "SQR", "SQUARE", "SAW", "SAWTOOTH", "TRI", "TRIANGLE")]
    [string]$FUNCTION="SINE",

    [Alias("e")]
    [ValidateRange("Positive")]
    [double]$TOLERANCE=0.1
)

if ($HELP) {
    Write-Output "Options:"
    Write-Output "  -t, -period             Period [s] (> 0.0) of the periodic function. Default is 1.0[s]."
    Write-Output "  -d, -dc                 DC component of the periodic function. Default is 0.0."
    Write-Output "  -a, -amplitude          Amplitude of the periodic function. Default is 1.0."
    Write-Output "  -p, -phase              Phase shift [rad] of the periodic function. Default is 0.0[rad]."
    Write-Output "  -f, -function           Case insensitive name of the periodic function. Options are: SIN or SINE; SQR or SQUARE;"
    Write-Output "                          SAW or SAWTOOTH; TRI or TRIANGLE. Default is SINE."
    Write-Output "  -e, -tolerance          Tolerance (> 0.0) for the stopping criteria. Default is 0.1."
    exit 0
}

$FUNCTION=$FUNCTION.ToUpper()

# Store current location
Push-Location

# cd to script
Set-Location "$(Split-Path $MyInvocation.MyCommand.Path)"

# Compile executable
mkdir -p ../build
Set-Location ../build
if (-Not (Test-Path -Path "ctfs_win")) {
    Write-Output "Compiling executable..."
    gfortran -o ctfs_win ../src/functions.f95 ../src/fourier.f95 ../src/ctfs.f95 -O2
    if ($LASTEXITCODE -ne 0) {
        Write-Error "Internal error! Compilation failed."
        Pop-Location
        exit 1
    }
    Write-Output "Done."
}

# Run analysis
Write-Output "Running calculations..."
Write-Output "|   Period: $PERIOD [s]"
Write-Output "|   DC: $DC"
Write-Output "|   Amplitude: $AMPLITUDE"
Write-Output "|   Phase: $PHASE [rad/s]"
Write-Output "|   Function: $FUNCTION"
Write-Output "|   Tolerance: $TOLERANCE"
./ctfs_win "$PERIOD" "$DC" "$AMPLITUDE" "$PHASE" "$FUNCTION" "$TOLERANCE"
if ($LASTEXITCODE -ne 0) {
    Write-Error "Internal error! Failed to complete the calculations."
    Pop-Location
    exit 1
}
Write-Output "Done."

# Move output
mkdir -p "$FUNCTION/data" "$FUNCTION/figures"
rm -f "$FUNCTION/data/*.dat"
mv -t "$FUNCTION/data" *.dat
Set-Location ..

# Create and activate venv
if (-Not (Test-Path -Path "ctfs")) {
    ./bin/create_venv.ps1
    if ( $LASTEXITCODE -ne 0 ) {
        Pop-Location
        exit 1
    }
}
Write-Output "Activating venv..."
./ctfs/Scripts/Activate.ps1
if ($LASTEXITCODE -ne 0) {
    Write-Error "Internal error! Could not activate the virtual environment."
    Pop-Location
    exit 1
}
Write-Output "Done."

# Plot output of calculations
Set-Location bin
Write-Output "Running analysis..."
python3 -m analysis -t "$PERIOD" -d "$DC" -a "$AMPLITUDE" -p "$PHASE" -f "$FUNCTION" -e "$TOLERANCE"
if ($LASTEXITCODE -ne 0) {
    Write-Error "Internal error! Failed to plot the results."
    Pop-Location
    exit 1
}

# Restores the original location
Pop-Location

Write-Output "Done."

exit 0