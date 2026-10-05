"""
write_legacy_sidecars.py
Write the schema-1.0 sidecars (validation/README.md, Dataset sidecars) of the
two legacy 500 GHz datasets into data/waveguides/. The values come from the CSVs and from Blackbody-Simulations
configs/crack{1,2}_500GHz_reference.toml. The pose, the designs, the project
hash and the far-field steps were confirmed against InfParallelPlate.aedt on
2026-10-05. Only the PEC side walls remain inferred.

The script refuses to run if any CSV differs from the checksums recorded here,
so it cannot write a sidecar for changed data. Run from the repository root:
  conda run -n bbrsim python validation/Scripts/write_legacy_sidecars.py
"""
import os
import sys

from bbrsim import sidecar

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "waveguides")
EXPECTED_SHA256 = {
    "InfParallelPlate_crack1Rohan_500GHz_Ephi=0/far_field.csv": "e1eb1b28ff232d3b6ff217e67f0df208e796618833874187aa84216f4d5d67b9",
    "InfParallelPlate_crack1Rohan_500GHz_Ephi=0/waveguide.csv": "e75d467047bf5439a066732f1c6a02a0d10b355695130206940f6cec4e097bad",
    "InfParallelPlate_crack1Rohan_500GHz_Ephi=1/far_field.csv": "304faf1ec5a9c33e652c600c6c666e53e9957456d166a1504a165b638fc2a3b2",
    "InfParallelPlate_crack1Rohan_500GHz_Ephi=1/waveguide.csv": "a7bf2de1ad342975ad71ac68ea3a1f459ae44c2a2b89335783b26fb4c1f5f0e5",
    "InfParallelPlate_crack2_500GHz_Ephi=0/far_field.csv": "bee4fe0f50272a37f43a488717f9404c54f3f784526eb083bc2350b4de42537f",
    "InfParallelPlate_crack2_500GHz_Ephi=0/waveguide.csv": "98389c6d24fcdb4c12bdd214cccdfb72ef99d13adac04a2709b928a1d44be362",
    "InfParallelPlate_crack2_500GHz_Ephi=1/far_field.csv": "212b3da02dd9f4885e75e9239ecae568c7a46f38a812fcf5fd9f05ecd4cb51e5",
    "InfParallelPlate_crack2_500GHz_Ephi=1/waveguide.csv": "1e704751ba465cd7fafd8ead45af5c8087de760911b921fa8ebbfbd1a2fb54cd",
}
COMMON = {
    "producer": "legacy bbsim1freq.py", "hand_written": True,
    "written_by": "validation/Scripts/write_legacy_sidecars.py",
    "git_commit": None, "created_utc": None, "job_id": None, "manifest_sha256": None,
    "aedt_version": "2023 R2", "pyaedt_version": None, "project": "InfParallelPlate.aedt",
    "project_sha256": "a264705100ab80ff2a3b2315b0529a8505c2903b996c98cd0d65d7144097482f",
    "material": None, "model_units": "mm", "inferred": ["boundaries.walls"],
}
LEGACY = {   # id: (design, run_design, object, renamed design, gap [m], length [mm])
    "InfParallelPlate_crack1Rohan": ("crack1Rohan", "crack1Rohan_500GHz", "crack1", "parallel_plate_gap_50um", 5e-05, 1.0),
    "InfParallelPlate_crack2": ("crack2", "crack2_500GHz", "crack2", "parallel_plate_gap_100um", 1e-04, 1.5),
}

for dataset_id, (design, run_design, obj, renamed, gap, length) in LEGACY.items():
    stem = f"{dataset_id}_500GHz"
    got = {k: v["sha256"] for k, v in sidecar.file_entries(BASE, stem).items()}
    want = {k: v for k, v in EXPECTED_SHA256.items() if k.startswith(stem + "_")}
    if got != want:
        print(f"refusing: the CSVs of {stem} differ from the recorded checksums")
        sys.exit(2)
    gap_mm = gap * 1e3
    sc = sidecar.build_from_csvs(
        BASE, dataset_id, "500GHz",
        section={"shape": "rectangle", "y_e_half_m": 0.005, "z_e_half_m": gap / 2},
        extent_mm={"p": length, "l": 10.0, "g": gap_mm},
        provenance={**COMMON, "design": design, "run_design": run_design, "object": obj,
                    "renamed_to": {"design": renamed, "object": "gap"}},
        exit_origin_mm=[gap_mm / 2, 0.0, length], plane_wave_origin_mm=[gap_mm / 2, 0.0, 0.0],
        bounding_box_mm=[0.0, -5.0, 0.0, gap_mm, 5.0, length],
        pose_rule="legacy", resolution_mm=[0.0, 0.1, 0.001])
    path = sidecar.write(BASE, stem, sc)
    sidecar.check_full(BASE, dataset_id, stem, 500.0)
    print(f"wrote {os.path.relpath(path)}")
