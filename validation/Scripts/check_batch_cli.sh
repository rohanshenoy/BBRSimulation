#!/usr/bin/env bash
# Exercise the real example executables with successful and failing batch macros.
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "usage: $0 <bbrsimTestWorld> <bbrsimLightPipe>" >&2
  exit 2
fi

scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT

printf '/control/verbose 0\n' > "$scratch/valid.mac"
printf '/bbr/thisCommandDoesNotExist\n' > "$scratch/invalid.mac"
printf '/bbr/thermal/setT definitely-not-a-number\n' > "$scratch/bad-parameter.mac"
printf '/control/execute %s\n/control/verbose 0\n' \
  "$scratch/invalid.mac" > "$scratch/nested-invalid.mac"
printf '/control/execute %s\n/control/verbose 0\n' \
  "$scratch/missing-child.mac" > "$scratch/nested-missing.mac"

for binary in "$@"; do
  if ! "$binary" "$scratch/valid.mac" > "$scratch/run.log" 2>&1; then
    echo "FAIL: $(basename "$binary") rejected a valid macro" >&2
    cat "$scratch/run.log" >&2
    exit 1
  fi
  for case in missing invalid bad-parameter nested-invalid nested-missing; do
    if "$binary" "$scratch/$case.mac" > "$scratch/run.log" 2>&1; then
      echo "FAIL: $(basename "$binary") accepted $case macro" >&2
      cat "$scratch/run.log" >&2
      exit 1
    fi
  done
done

echo "PASS: both example binaries report batch macro errors"
