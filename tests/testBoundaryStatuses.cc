#include "BBRTestSupport.hh"
#include "BBRAnalysis.hh"
#include "G4OpBoundaryProcess.hh"

#include <array>
#include <regex>

int main(int argc, char** argv) {
  return bbrtest::RunCase(argc, argv, {
    {"stock_statuses", [] {
      struct Expected { G4OpBoundaryProcessStatus status; const char* name; G4int type; };
      const std::array<Expected, 43> expected = {{
        {Undefined, "Undefined", 3}, {Transmission, "Transmission", 0},
        {FresnelRefraction, "FresnelRefraction", 0}, {FresnelReflection, "FresnelReflection", 1},
        {TotalInternalReflection, "TIR", 1}, {LambertianReflection, "LambertianReflection", 1},
        {LobeReflection, "LobeReflection", 1}, {SpikeReflection, "SpikeReflection", 1},
        {BackScattering, "BackScattering", 1}, {Absorption, "Absorption", 2},
        {Detection, "Detection", 2}, {NotAtBoundary, "NotAtBoundary", 3},
        {SameMaterial, "SameMaterial", 3}, {StepTooSmall, "StepTooSmall", 3},
        {NoRINDEX, "NoRINDEX", 2},
        {PolishedLumirrorAirReflection, "PolishedLumirrorAirReflection", 1},
        {PolishedLumirrorGlueReflection, "PolishedLumirrorGlueReflection", 1},
        {PolishedAirReflection, "PolishedAirReflection", 1},
        {PolishedTeflonAirReflection, "PolishedTeflonAirReflection", 1},
        {PolishedTiOAirReflection, "PolishedTiOAirReflection", 1},
        {PolishedTyvekAirReflection, "PolishedTyvekAirReflection", 1},
        {PolishedVM2000AirReflection, "PolishedVM2000AirReflection", 1},
        {PolishedVM2000GlueReflection, "PolishedVM2000GlueReflection", 1},
        {EtchedLumirrorAirReflection, "EtchedLumirrorAirReflection", 1},
        {EtchedLumirrorGlueReflection, "EtchedLumirrorGlueReflection", 1},
        {EtchedAirReflection, "EtchedAirReflection", 1},
        {EtchedTeflonAirReflection, "EtchedTeflonAirReflection", 1},
        {EtchedTiOAirReflection, "EtchedTiOAirReflection", 1},
        {EtchedTyvekAirReflection, "EtchedTyvekAirReflection", 1},
        {EtchedVM2000AirReflection, "EtchedVM2000AirReflection", 1},
        {EtchedVM2000GlueReflection, "EtchedVM2000GlueReflection", 1},
        {GroundLumirrorAirReflection, "GroundLumirrorAirReflection", 1},
        {GroundLumirrorGlueReflection, "GroundLumirrorGlueReflection", 1},
        {GroundAirReflection, "GroundAirReflection", 1},
        {GroundTeflonAirReflection, "GroundTeflonAirReflection", 1},
        {GroundTiOAirReflection, "GroundTiOAirReflection", 1},
        {GroundTyvekAirReflection, "GroundTyvekAirReflection", 1},
        {GroundVM2000AirReflection, "GroundVM2000AirReflection", 1},
        {GroundVM2000GlueReflection, "GroundVM2000GlueReflection", 1},
        {Dichroic, "Dichroic", 0},
        {CoatedDielectricReflection, "CoatedDielectricReflection", 1},
        {CoatedDielectricRefraction, "CoatedDielectricRefraction", 0},
        {CoatedDielectricFrustratedTransmission, "CoatedDielectricFrustratedTransmission", 0},
      }};
      for (const auto& row : expected) {
        const auto name = BBRAnalysis::BoundaryStatusName(row.status);
        CHECK(name == row.name);
        CHECK(BBRAnalysis::EventTypeForStatus(name) == row.type);
      }
      CHECK(BBRAnalysis::EventTypeForStatus("BBRDiffractionReflect") == 1);
      CHECK(BBRAnalysis::EventTypeForStatus("BBRReflect") == 1);
      CHECK(BBRAnalysis::EventTypeForStatus("BBRDiffractionTransmit") == 0);
    }},
    {"data_fingerprint", [] {
      // The metadata fingerprint hashes only <data root>/waveguides/*.dataset.json.
      bbrtest::TempDir tmp;
      const auto root = tmp.path() / "data";
      const auto wg = root / "waveguides";
      const auto fp = [&] { return BBRAnalysis::DataFingerprint(root.string()); };
      CHECK(BBRAnalysis::DataFingerprint((tmp.path() / "absent").string()) == "unavailable");
      bbrtest::WriteFile(root / "README", "no waveguides here\n");
      CHECK(fp() == "unavailable");
      bbrtest::WriteFile(wg / "c_500GHz_Ephi=0" / "waveguide.csv", "Freq\n500GHz\n");
      bbrtest::WriteFile(wg / "notes.txt", "not a sidecar\n");
      CHECK(fp() == "no-sidecars");
      bbrtest::WriteFile(wg / "c_500GHz.dataset.json", "{\"schema_version\": \"1.0\"}\n");
      const std::string first = fp();
      CHECK(std::regex_match(first, std::regex("^fnv1a64:[0-9a-f]{16}$")));
      CHECK(fp() == first);
      bbrtest::WriteFile(wg / "c_500GHz_Ephi=0" / "waveguide.csv", "Freq\n500GHz\n600GHz\n");
      CHECK(fp() == first);                                   // a CSV edit is not hashed
      bbrtest::WriteFile(wg / "refined_1.csv", "x\n");
      bbrtest::WriteFile(wg / "SHA256SUMS", "abc  x\n");
      bbrtest::WriteFile(wg / "c_500GHz_Ephi=1" / "extra.dataset.json", "{}\n");
      CHECK(fp() == first);                                   // nor a non-sidecar or nested file
      bbrtest::WriteFile(wg / "c_500GHz.dataset.json", "{\"schema_version\": \"1.1\"}\n");
      const std::string edited = fp();
      CHECK(std::regex_match(edited, std::regex("^fnv1a64:[0-9a-f]{16}$")));
      CHECK(edited != first);                                 // a sidecar edit is
      bbrtest::WriteFile(wg / "d_500GHz.dataset.json", "{}\n");
      CHECK(fp() != edited);                                  // and so is a new sidecar
    }},
  });
}
