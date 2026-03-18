# Plan: Remove Boost C++ Dependency

## Motivation

Boost is a large, heavy external dependency that requires platform-specific installation
(Homebrew on macOS, module load on HPC). Everything Boost provides in this project has a
standard library equivalent in C++17/C++20. Removing it simplifies builds, reduces
dependency surface, and makes the project more portable.

## Completed work

- **GSL removed** (2026-03-18): The only GSL use was `gsl_integration_qag` (21-point
  Gauss-Kronrod adaptive quadrature) in `src/libtao/base/utils.hh` (`redshift_to_age`).
  `redshift_to_age` was never called anywhere, so it was deleted entirely along with the
  `cosmology<T>` struct. GSL has been removed from `project.cmake`, `platform_config.cmake`,
  and `.cmake/3rd_party.cmake`. `boost::math::quadrature` is therefore also not needed.

## Current Boost Usage Inventory

Five compiled Boost components are declared in `project.cmake` and linked:

| Component | Where used | C++ stdlib replacement |
|-----------|-----------|------------------------|
| `program_options` | `src/libhpc/system/application.{cc,hh}` | No stdlib equivalent — use CLI11 (header-only) |
| `regex` | `src/libhpc/algorithm/multimatch.{cc,hh}` | `<regex>` (C++11) |
| `filesystem` | 7 files (see below) | `<filesystem>` (C++17) |
| `chrono` | `src/libhpc/system/timer.hh` | `<chrono>` (C++11) |
| `unit_test_framework` | Declared but **never used** | Remove from CMake |

Header-only Boost headers used (no linking required, but still a hard dependency):

| Header(s) | Where used | Replacement |
|-----------|-----------|-------------|
| `boost/algorithm/string.hpp`, `string/trim.hpp` | `kdapplication.cc`, `filter.cc`, `string.hh`, `group.cc`, `dataset.cc` | `<algorithm>` + `<cctype>` hand-rolled trim/lower |
| `boost/algorithm/string/join.hpp` | `multimatch.cc` | `std::ostringstream` loop |
| `boost/foreach.hpp` | `kdapplication.cc` | Range-based `for` (C++11) |
| `boost/tokenizer.hpp` | `kdapplication.cc` | `std::istringstream` or manual find/split |
| `boost/range/algorithm/fill.hpp` | `sage2kdtree.cc`, `batch.hh`, MPI files | `std::ranges::fill` (C++20) or `std::fill` |
| `boost/range/algorithm_ext/iota.hpp` | Multiple | `std::iota` (C++11) |
| `boost/range/algorithm.hpp`, `algorithm/copy.hpp` | MPI files, `kdtree.hh` | `std::ranges::copy`, `std::ranges::sort`, etc. |
| `boost/range/adaptors.hpp` | `csr.hh` | `std::views::` (C++20) or manual loop |
| `boost/range/adaptor/reversed.hpp` | `path_finder.cc` | `std::views::reverse` (C++20) or reverse iterators |
| `boost/range/numeric.hpp` | `partition.hh` | `std::accumulate` / `std::reduce` (C++17) |
| `boost/iterator/iterator_facade.hpp` | `range_map.hh`, `sort_by_key.hh`, `sort_permute_iter.hh`, `clip.hh`, `spline.hh`, `interp_iterator.hh`, `view.hh`, `matrix.hh` | Write iterators from scratch (C++17 style) — 8 iterator classes |
| `boost/iterator/iterator_adaptor.hpp` | `enumerator.hh` | Write from scratch |
| `boost/iterator/transform_iterator.hpp` | `factory.hh`, `clip.hh` | `std::views::transform` (C++20) or lambda wrapper |
| `boost/iterator/zip_iterator.hpp` | `clip.hh` (libtao), `clip.hh` (libhpc) | `std::views::zip` (C++23) or simple custom zip |
| `boost/tuple/tuple.hpp`, `tuple_io.hpp` | `standard.hh`, `derive.hh`, `stream.hh`, `functional.hh` | `std::tuple` / `std::tie` (C++11) |

---

## Replacement Strategy

### Step 1 — Drop `unit_test_framework` (trivial)

`unit_test_framework` is listed in `project.cmake` but no source file includes it.

- Remove `"unit_test_framework"` from `BOOST_COMPONENTS` in `project.cmake`
- No source changes needed

**Risk**: None. **Effort**: 5 minutes.

---

### Step 2 — Replace `boost::chrono` with `std::chrono` (easy)

`src/libhpc/system/timer.hh` is the only file. Drop-in replacement:

- `#include <boost/chrono.hpp>` → `#include <chrono>`
- `boost::chrono::` → `std::chrono::`
- `boost::chrono::process_real_cpu_clock` has no exact stdlib analogue — replace with
  `std::chrono::steady_clock` (wall time; good enough for the profiling use case here).
  If CPU time is ever needed, `std::clock()` from `<ctime>` can supplement.
- Remove `"chrono"` from `BOOST_COMPONENTS` in `project.cmake`

**Risk**: Low. **Effort**: ~30 minutes.

---

### Step 3 — Replace `boost::regex` with `std::regex` (easy)

Only two files: `src/libhpc/algorithm/multimatch.{cc,hh}`.

- `#include <boost/regex.hpp>` → `#include <regex>`
- `boost::regex` → `std::regex`
- `boost::regex_match()` → `std::regex_match()`
- `boost::regex_search()` → `std::regex_search()`
- `boost::smatch` (if present) → `std::smatch`
- Remove `"regex"` from `BOOST_COMPONENTS` in `project.cmake`

**Risk**: Low — API is nearly identical. **Effort**: ~30 minutes.

---

### Step 4 — Replace `boost::filesystem` with `std::filesystem` (easy)

Files affected:

- `src/apps/kdapplication.cc`
- `src/libtao/base/utils.hh`
- `src/libhpc/system/filesystem.hh`
- `src/libhpc/system/file_descriptor.hh`
- `src/libhpc/system/path_finder.hh`
- `src/libhpc/system/tmpfile.hh`
- `src/libhpc/mpi/init.cc`
- `src/libtao/modules/hdf5.hh`

Changes per file:

- `#include <boost/filesystem.hpp>` (and variants) → `#include <filesystem>`
- `namespace fs = boost::filesystem` → `namespace fs = std::filesystem`
- `using namespace boost::filesystem` → `using namespace std::filesystem` (or use `fs::`)
- `boost::filesystem::path` → `std::filesystem::path`
- `boost::filesystem::exists()` → `std::filesystem::exists()`
- `boost::filesystem::create_directories()` → `std::filesystem::create_directories()`
- `boost::filesystem::remove()` → `std::filesystem::remove()`

The API is nearly 1:1; `std::filesystem` was modelled on `boost::filesystem`.
Remove `"filesystem"` from `BOOST_COMPONENTS` in `project.cmake`.

**Risk**: Low. May need `-lstdc++fs` on older GCC — not an issue here (Apple Clang 17).
**Effort**: ~1 hour.

---

### Step 5 — Replace header-only range/algorithm Boost headers (medium)

These are all header-only and don't require a compiled component, but still need Boost
installed. Replace file by file:

#### ~~5a. `boost::math::quadrature`~~ — DONE (2026-03-18)

`redshift_to_age` and its `cosmology<T>` struct were deleted (never called). The
`boost/math/quadrature/gauss_kronrod.hpp` include is gone. No replacement needed.

#### 5c. String algorithms (`boost::algorithm::*`)

Files: `kdapplication.cc`, `filter.cc`, `string.hh`, `group.cc`, `dataset.cc`

- `boost::algorithm::to_lower_copy(s)` →
  ```cpp
  std::string lower = s;
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
  ```
- `boost::algorithm::trim(s)` → small inline helper (10 lines) in a utility header, or
  use `std::string`'s `find_first_not_of`/`find_last_not_of` idiom
- `boost::algorithm::join(container, sep)` → `std::ostringstream` loop
- `boost::algorithm::split(...)` → `std::istringstream` or manual `find`+`substr` loop

#### 5d. `boost::foreach` (one file: `kdapplication.cc`)

`BOOST_FOREACH(x, range)` → `for (auto& x : range)`. C++11 range-based for.

#### 5c. `boost::tokenizer` (one file: `kdapplication.cc`)

Replace with `std::istringstream` + `std::getline` using delimiter, or manual
`find`/`substr` loop.

#### 5d. Range algorithms (`boost::fill`, `boost::iota`, `boost::copy`, `boost::sort`, etc.)

All have direct `std::` equivalents:

| Boost | Replacement |
|-------|-------------|
| `boost::fill(r, v)` | `std::fill(r.begin(), r.end(), v)` or `std::ranges::fill(r, v)` |
| `boost::iota(r, v)` | `std::iota(r.begin(), r.end(), v)` |
| `boost::copy(src, dst)` | `std::copy(src.begin(), src.end(), dst)` or `std::ranges::copy` |
| `boost::sort(r)` | `std::sort(r.begin(), r.end())` or `std::ranges::sort(r)` |

#### 5e. Range adaptors (`boost::adaptors::reversed`, etc.)

- `boost::adaptors::reverse(r)` → `std::views::reverse(r)` (C++20), or use reverse
  iterators `rbegin()`/`rend()` with a range-based for.

#### 5f. `boost::tuple`

- `#include <boost/tuple/tuple.hpp>` → `#include <tuple>`
- `boost::tuple<A,B>` → `std::tuple<A,B>`
- `boost::get<N>(t)` → `std::get<N>(t)`
- `boost::make_tuple(...)` → `std::make_tuple(...)`
- `boost::tie(...)` → `std::tie(...)`
- Tuple I/O (`tuple_io.hpp`): write a small `operator<<` or use structured bindings.

**Risk**: Low — direct API matches. **Effort**: ~2–3 hours across all files.

---

### Step 6 — Replace `boost::iterator_facade` / `boost::iterator_adaptor` (hard)

This is the most labour-intensive step. Eight custom iterator classes use
`boost::iterator_facade` as their base:

- `src/libhpc/containers/range_map.hh`
- `src/libhpc/algorithm/sort_by_key.hh`
- `src/libhpc/algorithm/sort_permute_iter.hh`
- `src/libhpc/numerics/clip.hh`
- `src/libhpc/numerics/spline.hh`
- `src/libhpc/numerics/interp_iterator.hh`
- `src/libhpc/system/view.hh`
- `src/libhpc/system/matrix.hh`

Plus one using `boost::iterator_adaptor`:
- `src/libhpc/containers/enumerator.hh`

`boost::iterator_facade` provides boilerplate reduction for writing standard-compliant
iterators. The C++17/20 approach is to implement the iterator methods directly (they are
not that many). For each class, the migration is:

1. Remove `boost::iterator_facade<...>` base class
2. Add the five standard iterator typedef/using aliases:
   ```cpp
   using iterator_category = std::random_access_iterator_tag; // or forward/bidirectional
   using value_type        = ...;
   using difference_type   = std::ptrdiff_t;
   using pointer           = value_type*;
   using reference         = value_type&;
   ```
3. Implement the required operators directly (`++`, `--`, `+=`, `-`, `*`, `->`, `==`,
   `<`, `[]`). `iterator_facade` was generating these from a smaller set of primitives;
   now write them out explicitly. Most are one-liners.

For `enumerator.hh`, replace `boost::iterator_adaptor` with a small wrapper that holds
the wrapped iterator plus an index counter.

**Risk**: Medium — each iterator must be tested. Existing test suites in
`src/libhpc.project/tests/` cover some of these.
**Effort**: ~4–6 hours.

---

### Step 7 — Replace `boost::transform_iterator` and `boost::zip_iterator` (medium)

#### `boost::transform_iterator`

Used in `src/libtao/base/factory.hh` and `src/libhpc/numerics/clip.hh`.

With C++20: replace with `std::views::transform` at the call sites. If the iterator
object itself needs to be stored/passed, write a small templated `transform_iterator<It,
Fn>` struct (~30 lines) holding `It` and `Fn`, with `operator*` calling `fn(*it)`.

#### `boost::zip_iterator`

Used in `src/libtao/modules/clip.hh` and `src/libhpc/numerics/clip.hh`.

- C++23: `std::views::zip` — ideal, but C++23 is not yet required here.
- C++20 alternative: write a simple `zip_view` or use structured bindings with index
  loops. The usage in these files is limited; a manual index-based loop (`for (size_t
  i = 0; i < n; ++i)`) is the simplest and most readable replacement.

**Risk**: Low–medium. **Effort**: ~2 hours.

---

### Step 8 — Replace `boost::program_options` with CLI11 (medium)

`boost::program_options` is used in `src/libhpc/system/application.{cc,hh}` for all
command-line argument parsing in both `sage2kdtree` and `cli_lightcone`.

There is no stdlib equivalent. Recommended replacement: **CLI11**
(https://github.com/CLIUtils/CLI11) — single-header, MIT license, zero dependencies,
C++11 compatible, very similar API philosophy.

Migration:

1. Add `CLI11.hpp` (single header) to `src/libhpc/system/` or a `thirdparty/` directory.
   Alternatively vendor via CMake `FetchContent`.

2. In `application.hh`:
   - Remove `#include <boost/program_options.hpp>` and `namespace po = boost::program_options`
   - Replace `po::options_description` with `CLI::App`
   - Replace `po::variables_map` with direct typed variables bound via `app.add_option()`
   - Remove `po::store()` / `po::notify()` / `po::command_line_parser()` — CLI11 parses
     in a single `app.parse(argc, argv)` call

3. In `application.cc`:
   - Replace exception types: `po::required_option` / `po::error` →
     `CLI::RequiredError` / `CLI::ParseError`

4. Update all callers that add options to the application (both `sage2kdtree.cc` and
   `kdapplication.cc`). The `app.add_option("--name", var, "desc")` pattern in CLI11
   closely mirrors the Boost API; the transition is mostly mechanical.

5. Remove `"program_options"` from `BOOST_COMPONENTS` in `project.cmake`.

**Risk**: Medium — CLI11 and Boost program_options have subtle behavioural differences
(e.g. positional arguments, subcommands). Test both executables' `--help` output and
all argument combinations after migration.
**Effort**: ~3–4 hours.

---

### Step 9 — Remove Boost from CMake (final cleanup)

Once all components are replaced:

1. `project.cmake`: delete the entire `BOOST_COMPONENTS` block (lines 55–63) and the
   `set_3rd_party_required("Boost" ...)` call.
2. `.cmake/3rd_party.cmake`: remove the Boost `find_package` block (~lines 480–503).
3. `CMakeLists.txt`: the `CMP0167` policy workaround (lines 15–18) can be removed.
4. Build scripts (`setup_mac.sh`, `setup.sh`): remove any Boost path exports or module
   loads.
5. Run `cmake ..` fresh and confirm Boost is no longer searched for or linked.

---

## Recommended Execution Order

| Priority | Step | Effort | Risk |
|----------|------|--------|------|
| 1 | Step 1: drop unused `unit_test_framework` | 5 min | None |
| 2 | Step 2: `boost::chrono` → `std::chrono` | 30 min | Low |
| 3 | Step 3: `boost::regex` → `std::regex` | 30 min | Low |
| 4 | Step 4: `boost::filesystem` → `std::filesystem` | 1 hr | Low |
| 5 | Step 5: header-only string/range/tuple | 2–3 hr | Low |
| 6 | Step 7: transform/zip iterators | 2 hr | Medium |
| 7 | Step 6: iterator_facade / iterator_adaptor | 4–6 hr | Medium |
| 8 | Step 8: `boost::program_options` → CLI11 | 3–4 hr | Medium |
| 9 | Step 9: CMake cleanup | 30 min | Low |

**Total estimated effort: ~14–18 hours** of focused work across ~40 files.

Steps 1–5 are safe to do incrementally without a flag day. Each step can be built and
tested in isolation. Steps 6–8 are best done on a dedicated branch.

## Testing After Each Step

After each step, confirm:

```bash
./build_platform_aware.sh
cd tests/sage-model-tests
./run_test_sage_hdf5.sh
```

Output galaxy counts and lightcone content must match baseline values from
`plans/VALIDATE_SATELLITE.md`.

## Notes

- **C++20 helps significantly**: if the project moves to C++20 (`-DCMAKE_CXX_STANDARD=20`),
  `std::ranges::` and `std::views::` provide direct replacements for most
  `boost::range` and `boost::adaptors` usage, making Steps 5e and 7 trivial.
- **CLI11 is the only new dependency**: all other replacements use the stdlib. CLI11 is
  MIT-licensed and can be vendored as a single 8KB header — arguably not a "dependency"
  at all.
- **HPC portability**: `std::filesystem` requires GCC 8+ (linked with `-lstdc++fs` on
  GCC 8; automatic from GCC 9). HPC clusters with GCC 7 would need upgrading first.
  Apple Clang 17 and GCC 9+ have no issues.
