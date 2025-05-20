#!/bin/pwsh

### Creates the Python virtual environment ###

Write-Output "Creating venv..."

# Store current location
Push-Location

# cd to script
Set-Location -Path "$PSScriptRoot/.."

# Create venv
python -m venv ctfs
if ($LASTEXITCODE -ne 0) {
    Write-Error "Internal error! Could not create virtual environment."
    Pop-Location
    exit 1
}
# Activate venv
./ctfs/Scripts/Activate.ps1

Write-Output "Installing required packages..."

# Install pip
python -m pip install pip --upgrade
if ($LASTEXITCODE -ne 0) {
    Write-Error "Internal error! Could not install pip."
    Pop-Location
    exit 1
}
# Install requirements.txt
pip install -r misc/requirements.txt
if ($LASTEXITCODE -ne 0) {
    Write-Error "Internal error! Could not install the required packages."
    Pop-Location
    exit 1
}

# Restores the original location
Pop-Location

Write-Output "Done."

exit 0