// BBRTestSupport.hh — shared helpers for the BBRsim C++ tests in tests/.
//
// Each test program holds several named cases; CTest runs one case per
// process (see tests/CMakeLists.txt), because the library's singletons
// (BBRCrackLibrary, BBRConfigManager) cannot be reset within a process.
//
//   int main(int argc, char** argv) {
//     return bbrtest::RunCase(argc, argv, {
//       {"some_case", [] { CHECK(1 + 1 == 2); }},
//     });
//   }
#ifndef BBRTestSupport_hh
#define BBRTestSupport_hh

#include "G4Electron.hh"
#include "G4ExceptionSeverity.hh"
#include "G4ThreeVector.hh"
#include "G4VExceptionHandler.hh"
#include "CLHEP/Random/RandomEngine.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>  // mkdtemp

// G4iosSetDestination arrived in Geant4 11.2. Before it the cout and cerr
// buffers took the destination directly, as G4UImanager::SetCoutDestination does.
#include "G4Version.hh"
#if G4VERSION_NUMBER < 1120
#  include "G4coutDestination.hh"
#  include "G4strstreambuf.hh"
inline void G4iosSetDestination(G4coutDestination* sink)
{
  G4coutbuf.SetDestination(sink);
  G4cerrbuf.SetDestination(sink);
}
#endif

namespace bbrtest {

// Before Geant4 11.4, G4OpticalPhysics::ConstructProcess always calls
// G4LossTableManager::Instance(), which creates e- after PreInit. That e- has no
// process manager, and G4OpticalPhysics aborts on it. An optical-only physics
// list therefore creates e- in PreInit, where it gets one.
inline void CreateElectronInPreInit()
{
#if G4VERSION_NUMBER < 1140
  G4Electron::Definition();
#endif
}

// A G4Exception turned into a C++ exception by ThrowingExceptionHandler.
struct G4ExceptionCaught : std::runtime_error {
  std::string code, origin;
  G4ExceptionSeverity severity;
  G4ExceptionCaught(const char* o, const char* c, G4ExceptionSeverity s, const char* d)
    : std::runtime_error(d ? d : ""), code(c ? c : ""), origin(o ? o : ""), severity(s) {}
};

// Throws on every G4Exception that is not a JustWarning, so a test can check
// the code in-process instead of letting the process abort. Warnings are
// printed and recorded. The G4VExceptionHandler base constructor registers the
// handler with this thread's G4StateManager; construct it (through Handler())
// before any run manager, which otherwise installs its own. Worker threads of
// an MT run manager get their own handler and still abort: test those fatals
// with tests/ExpectAbort.cmake.
class ThrowingExceptionHandler : public G4VExceptionHandler {
 public:
  std::vector<std::string> warnings;
  G4bool Notify(const char* origin, const char* code, G4ExceptionSeverity sev,
                const char* desc) override {
    if (sev == JustWarning) {
      warnings.emplace_back(code ? code : "");
      std::fprintf(stdout, "warning %s from %s\n", code ? code : "", origin ? origin : "");
      return false;
    }
    throw G4ExceptionCaught(origin, code, sev, desc);
  }
  int CountWarnings(const std::string& code) const {
    int n = 0;
    for (const auto& w : warnings) n += (w == code);
    return n;
  }
};

inline ThrowingExceptionHandler& Handler() {
  static ThrowingExceptionHandler* h = new ThrowingExceptionHandler;  // lives to exit
  return *h;
}

inline int& FailCount() { static int n = 0; return n; }

inline void Report(bool ok, const char* file, int line, const std::string& what) {
  if (ok) {
    std::fprintf(stdout, "ok   %s\n", what.c_str());
  } else {
    std::fprintf(stderr, "FAIL %s:%d: %s\n", file, line, what.c_str());
    ++FailCount();
  }
}

inline std::string Num(double v) {
  char b[64];
  std::snprintf(b, sizeof b, "%.12g", v);
  return b;
}

}  // namespace bbrtest

#define CHECK(cond) ::bbrtest::Report((cond), __FILE__, __LINE__, #cond)

// |got - exp| <= tol
#define CHECK_NEAR(got, exp, tol)                                              \
  do {                                                                         \
    const double g_ = (got), e_ = (exp), t_ = (tol);                           \
    ::bbrtest::Report(std::abs(g_ - e_) <= t_, __FILE__, __LINE__,             \
                      std::string(#got " = ") + ::bbrtest::Num(g_) +           \
                        ", expected " + ::bbrtest::Num(e_) + " +- " +          \
                        ::bbrtest::Num(t_));                                   \
  } while (0)

// |got - exp| <= rel * |exp|
#define CHECK_REL(got, exp, rel)                                               \
  do {                                                                         \
    const double g_ = (got), e_ = (exp), r_ = (rel);                           \
    ::bbrtest::Report(std::abs(g_ - e_) <= r_ * std::abs(e_), __FILE__,        \
                      __LINE__,                                                \
                      std::string(#got " = ") + ::bbrtest::Num(g_) +           \
                        ", expected " + ::bbrtest::Num(e_) + " (rel " +        \
                        ::bbrtest::Num(r_) + ")");                             \
  } while (0)

// |got - exp| <= tol for 3-vectors
#define CHECK_VEC(got, exp, tol)                                               \
  do {                                                                         \
    const G4ThreeVector g_ = (got), e_ = (exp);                                \
    std::ostringstream m_;                                                     \
    m_ << #got " = " << g_ << ", expected " << e_;                             \
    ::bbrtest::Report((g_ - e_).mag() <= (tol), __FILE__, __LINE__, m_.str()); \
  } while (0)

// Same axis up to sign (a polarization vector v and -v are one state).
#define CHECK_AXIS(got, exp)                                                   \
  do {                                                                         \
    const G4ThreeVector g_ = (got), e_ = (exp);                                \
    std::ostringstream m_;                                                     \
    m_ << #got " = " << g_ << ", expected +-" << e_;                           \
    ::bbrtest::Report(std::abs(std::abs(g_.unit().dot(e_.unit())) - 1.) < 1e-9, \
                      __FILE__, __LINE__, m_.str());                           \
  } while (0)

namespace bbrtest {

// Runs f and requires a G4Exception with the given code (and, if given, an
// origin containing originPart and a message containing messagePart, so a case
// pins the check it means to reach, not merely one with the same code that
// fires earlier). Any other outcome is a failure.
template <class F>
void ExpectG4Exception(const std::string& code, F&& f, const std::string& originPart = "",
                       const std::string& messagePart = "") {
  const std::string containing = messagePart.empty() ? "" : " containing \"" + messagePart + "\"";
  try {
    f();
    Report(false, __FILE__, __LINE__, "expected G4Exception " + code + containing + ", none raised");
  } catch (const G4ExceptionCaught& e) {
    const bool originOk = originPart.empty() || e.origin.find(originPart) != std::string::npos;
    const bool messageOk = messagePart.empty() || std::string(e.what()).find(messagePart) != std::string::npos;
    Report(e.code == code && originOk && messageOk, __FILE__, __LINE__,
           "expected " + code + (originPart.empty() ? "" : " from *" + originPart + "*") + containing +
             ", got " + e.code + " from " + e.origin + (messagePart.empty() ? "" : std::string(": ") + e.what()));
  } catch (const std::exception& e) {
    Report(false, __FILE__, __LINE__,
           "expected G4Exception " + code + ", got C++ exception: " + e.what());
  }
}

// A random engine that returns a scripted sequence, so a sampler's choice can
// be computed by hand. Install with G4Random::setTheEngine(&engine).
class ScriptedEngine : public CLHEP::HepRandomEngine {
 public:
  ScriptedEngine() = default;
  ScriptedEngine(std::initializer_list<double> v) : fQueue(v) {}
  void Push(double v) { fQueue.push_back(v); }
  std::size_t Remaining() const { return fQueue.size(); }
  std::size_t Draws() const { return fDraws; }
  double flat() override {
    if (fQueue.empty()) throw std::runtime_error("ScriptedEngine exhausted");
    const double v = fQueue.front();
    fQueue.pop_front();
    ++fDraws;
    return v;
  }
  void flatArray(const int n, double* v) override {
    for (int i = 0; i < n; ++i) v[i] = flat();
  }
  void setSeed(long, int) override {}
  void setSeeds(const long*, int) override {}
  void saveStatus(const char[]) const override {}
  void restoreStatus(const char[]) override {}
  void showStatus() const override {}
  std::string name() const override { return "ScriptedEngine"; }

 private:
  std::deque<double> fQueue;
  std::size_t fDraws = 0;
};

// A fresh directory under $TMPDIR (else /tmp), removed with its contents when
// the object goes out of scope. Never under the repository: data/ must stay
// exactly as tracked.
class TempDir {
 public:
  TempDir() {
    const char* base = std::getenv("TMPDIR");
    std::string tmpl = std::string(base && *base ? base : "/tmp") + "/bbrsim-test-XXXXXX";
    std::vector<char> buf(tmpl.begin(), tmpl.end());
    buf.push_back('\0');
    if (!mkdtemp(buf.data())) throw std::runtime_error("mkdtemp failed for " + tmpl);
    fPath = buf.data();
  }
  ~TempDir() {
    std::error_code ec;
    std::filesystem::remove_all(fPath, ec);
  }
  TempDir(const TempDir&) = delete;
  TempDir& operator=(const TempDir&) = delete;
  const std::filesystem::path& path() const { return fPath; }

 private:
  std::filesystem::path fPath;
};

// Writes text to p, creating parent directories.
inline void WriteFile(const std::filesystem::path& p, const std::string& text) {
  std::filesystem::create_directories(p.parent_path());
  std::ofstream f(p, std::ios::binary);
  if (!f) throw std::runtime_error("cannot write " + p.string());
  f << text;
}

using CaseMap = std::map<std::string, std::function<void()>>;

// Runs the case named by argv[1]. Returns 0 if every check passed, 1 if any
// failed or the case threw, 2 for a missing or unknown case name.
inline int RunCase(int argc, char** argv, const CaseMap& cases) {
  Handler();  // register the throwing handler before anything else runs
  if (argc < 2 || !cases.count(argv[1])) {
    std::fprintf(stderr, "usage: %s <case>; cases:", argc > 0 ? argv[0] : "test");
    for (const auto& kv : cases) std::fprintf(stderr, " %s", kv.first.c_str());
    std::fprintf(stderr, "\n");
    return 2;
  }
  try {
    cases.at(argv[1])();
  } catch (const G4ExceptionCaught& e) {
    Report(false, __FILE__, __LINE__,
           "unexpected G4Exception " + e.code + " from " + e.origin + ": " + e.what());
  } catch (const std::exception& e) {
    Report(false, __FILE__, __LINE__, std::string("unexpected C++ exception: ") + e.what());
  }
  std::fprintf(stdout, "%s: %d failure(s)\n", argv[1], FailCount());
  return FailCount() == 0 ? 0 : 1;
}

}  // namespace bbrtest

#endif
