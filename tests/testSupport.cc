// testSupport.cc — self-test of tests/BBRTestSupport.hh: the throwing handler,
// the scripted engine, TempDir, and that a failed check or an unexpected
// exception really fails a test (those two cases are registered WILL_FAIL).
#include "BBRTestSupport.hh"

#include "G4Exception.hh"
#include "G4StateManager.hh"
#include "Randomize.hh"

#include <filesystem>

int main(int argc, char** argv) {
  using namespace bbrtest;
  return RunCase(argc, argv, {
    {"handler_throws", [] {
      ExpectG4Exception("BBRT01", [] {
        G4Exception("testSupport", "BBRT01", FatalException, "fatal on purpose");
      }, "testSupport");
      ExpectG4Exception("BBRT01", [] {
        G4Exception("testSupport", "BBRT01", FatalException, "fatal on purpose");
      }, "testSupport", "on purpose");
    }},
    // Registered WILL_FAIL: the right code and origin with a message lacking messagePart.
    {"message_mismatch_fails", [] {
      ExpectG4Exception("BBRT01", [] {
        G4Exception("testSupport", "BBRT01", FatalException, "fatal on purpose");
      }, "testSupport", "not in the message");
    }},
    {"handler_records_warning", [] {
      G4Exception("testSupport", "BBRT02", JustWarning, "warning on purpose");
      G4Exception("testSupport", "BBRT02", JustWarning, "warning on purpose");
      CHECK(Handler().CountWarnings("BBRT02") == 2);
      CHECK(Handler().CountWarnings("BBRT01") == 0);
    }},
    {"scripted_engine", [] {
      ScriptedEngine eng{0.25, 0.5, 0.75};
      G4Random::setTheEngine(&eng);
      CHECK(G4UniformRand() == 0.25);
      CHECK(G4UniformRand() == 0.5);
      CHECK(eng.Draws() == 2);
      CHECK(eng.Remaining() == 1);
      CHECK(G4UniformRand() == 0.75);
      bool threw = false;
      try { G4UniformRand(); } catch (const std::runtime_error&) { threw = true; }
      CHECK(threw);
    }},
    {"tempdir", [] {
      std::filesystem::path kept;
      {
        TempDir d;
        kept = d.path();
        CHECK(std::filesystem::is_directory(kept));
        WriteFile(kept / "a" / "b" / "f.txt", "x\n");
        CHECK(std::filesystem::exists(kept / "a" / "b" / "f.txt"));
      }
      CHECK(!std::filesystem::exists(kept));
    }},
    // Registered WILL_FAIL: a failed CHECK must make the test exit non-zero.
    {"check_fails", [] { CHECK(1 == 2); }},
    // Registered WILL_FAIL: an unexpected fatal G4Exception must fail the test.
    {"unexpected_exception", [] {
      G4Exception("testSupport", "BBRT04", FatalException, "not expected by the case");
    }},
    // Run through ExpectAbort.cmake: with the stock (no) handler a fatal
    // G4Exception aborts the process, which is what the driver checks.
    {"abort_stock", [] {
      G4StateManager::GetStateManager()->SetExceptionHandler(nullptr);
      G4Exception("testSupport", "BBRT03", FatalException, "abort on purpose");
    }},
  });
}
