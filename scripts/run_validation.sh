#!/bin/bash
# Run crack1 + crack2 validation and produce per-crack CSVs.
# Usage (from project root): bash scripts/run_validation.sh [path/to/build]
set -e
BUILD=${1:-build}
cd "$BUILD"

echo "=== crack1 (diffraction.mac) ==="
./BBRSim diffraction.mac
cp diffraction_output.csv crack1_output.csv

echo "=== crack2 (diffraction_crack2.mac) ==="
./BBRSim diffraction_crack2.mac
cp diffraction_output.csv crack2_output.csv

echo ""
echo "Outputs: $BUILD/crack1_output.csv  $BUILD/crack2_output.csv"
echo "Plot:    conda run -n bbrsim python ../scripts/plot_validation.py crack1_output.csv crack2_output.csv"
