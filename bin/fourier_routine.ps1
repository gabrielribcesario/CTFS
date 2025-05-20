#!/bin/pwsh

# Parse optional arguments
param(
    [switch]$HELP,

    [Alias("t")]
    [ValidateScript({
                    if ($_ -eq 0.0) {
                        throw "Cannot validate argument on parameter 'PERIOD': The argument '{0}' cannot be validated because its value is not greater than zero."
                    }
                    return $true
                    })]
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
    [ValidateScript({
                    if ($_ -eq 0.0) {
                        throw "Cannot validate argument on parameter 'TOLERANCE': The argument '{0}' cannot be validated because its value is not greater than zero."
                    }
                    return $true
                    })]
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
if ($FUNCTION -eq "SIN" -or $FUNCTION -eq "SINE") { $FUNCTION="SINE" }
elseif ($FUNCTION -eq "SQR" -or $FUNCTION -eq "SQUARE") { $FUNCTION="SQUARE" }
elseif ($FUNCTION -eq "SAW" -or $FUNCTION -eq "SAWTOOTH") { $FUNCTION="SAWTOOTH" }
elseif ($FUNCTION -eq "TRI" -or $FUNCTION -eq "TRIANGLE") { $FUNCTION="TRIANGLE" }

# Store current location
Push-Location

# cd to script
Set-Location -Path "$PSScriptRoot"

# Compile executable
New-Item -ItemType "Directory" -Force -Path ../build | Out-Null
Set-Location -Path ../build
if (-Not (Test-Path -Path "ctfs.exe")) {
    Write-Output "Compiling executable..."
    gfortran -o ctfs.exe ../src/functions.f95 ../src/fourier.f95 ../src/ctfs.f95 -O2
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
Write-Output "|   Phase: $PHASE [rad]"
Write-Output "|   Function: $FUNCTION"
Write-Output "|   Tolerance: $TOLERANCE"
./ctfs.exe "$PERIOD" "$DC" "$AMPLITUDE" "$PHASE" "$FUNCTION" "$TOLERANCE"
if ($LASTEXITCODE -ne 0) {
    Write-Error "Internal error! Failed to complete the calculations."
    Pop-Location
    exit 1
}
Write-Output "Done."

# Move output
New-Item -Force -ItemType "Directory" -Path "$FUNCTION/data" | Out-Null
New-Item -Force -ItemType "Directory" -Path "$FUNCTION/figures" | Out-Null
Remove-Item -Path "$FUNCTION/data/*.dat"
Move-Item -Path *.dat -Destination "$FUNCTION/data"
Set-Location -Path ..

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
Set-Location -Path bin
Write-Output "Running analysis..."
python -m analysis -t "$PERIOD" -d "$DC" -a "$AMPLITUDE" -p "$PHASE" -f "$FUNCTION" -e "$TOLERANCE"
if ($LASTEXITCODE -ne 0) {
    Write-Error "Internal error! Failed to plot the results."
    Pop-Location
    exit 1
}

# Restores the original location
Pop-Location

Write-Output "Done."

exit 0