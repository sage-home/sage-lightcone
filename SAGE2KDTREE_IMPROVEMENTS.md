# Sage2KDTree Incremental Improvements Plan

## Overview

This document outlines small, incremental improvements to the `sage2kdtree` consolidated pipeline, aligned with the architectural principles defined in VISION.md. The goal is to enhance maintainability, reliability, and adherence to core design principles without disrupting the working implementation.

## Current State Analysis

### What Works Well
- **Consolidated Pipeline**: Successfully combines 4 phases (sageh5toh5 → sageimport → tree2snap → kdtree)
- **Working End-to-End**: `run_test_hdf5_one_step.sh` demonstrates full workflow
- **Clean Error Handling**: Uses `PipelineException` with phase tracking
- **Well-Structured**: Clear separation of concerns across 4 phases

### Architecture Alignment with VISION.md

| Principle | Current Status | Gap |
|-----------|---------------|-----|
| 1. SAGE HDF5 as single source of truth | ✅ Good | Cosmology metadata could be more comprehensive |
| 2. Runtime Modularity (C++) | ✅ Excellent | Platform-aware builds working |
| 3. Metadata-Driven Architecture | ⚠️ Partial | Field handling is mostly hardcoded, not fully introspected |
| 4. Single Source of Truth | ✅ Good | Data flows through consistent pipeline |
| 5. Unified Processing Model | ✅ Good | Single traversal approach |
| 6. Memory Efficiency | ⚠️ Needs attention | Some phases load full datasets into memory |
| 7. Format-Agnostic I/O (HDF5) | ✅ Excellent | Purely HDF5-based |
| 8. Type Safety & Validation | ⚠️ Partial | Limited input validation, raw H5 API in places |

---

## Proposed Incremental Improvements

### Priority 1: Validation & Error Handling (Principle 8)

**Goal**: Catch errors early with clear messages, validate inputs comprehensively.

#### 1.1 Input Validation Enhancement
**Effort**: Small | **Risk**: Low | **Impact**: High

**Changes**:
- Add validation for SAGE HDF5 structure at start of Phase 1:
  - Check for required groups (`Header/Simulation`, `Snap_*`)
  - Validate cosmology attributes exist (`BoxSize`, `Hubble_h`, `Omega_m`, `Omega_lambda`)
  - Check for required fields in `Snap_0` (Posx, Posy, Posz, SAGETreeIndex)
- Add file existence checks for parameter file and alist before processing
- Validate snapshot count matches redshift list length

**Files**: `sage2kdtree.cc`
**Implementation**: Add `validate_input()` method called before Phase 1

#### 1.2 Enhanced Error Messages
**Effort**: Small | **Risk**: Low | **Impact**: Medium

**Changes**:
- Wrap HDF5 operations in try-catch blocks with context
- Add file:line information to errors
- Report which galaxy/tree/snapshot caused issues when processing fails
- Add suggestions to error messages (e.g., "Check that SAGE parameter file has correct OutputFormat")

**Files**: `sage2kdtree.cc`

#### 1.3 Progress Reporting
**Effort**: Small | **Risk**: Low | **Impact**: Medium

**Changes**:
- Add progress indicators for long-running operations:
  - Galaxy count during Phase 1 loading
  - Tree traversal progress in Phase 2
  - Snapshot processing in Phase 3/4
- Use existing verbosity flag to control detail level
- Estimated time remaining for large datasets

**Files**: `sage2kdtree.cc`

---

### Priority 2: Metadata-Driven Architecture (Principle 3)

**Goal**: Reduce hardcoding, dynamically discover fields from input SAGE HDF5.

#### 2.1 Dynamic Field Discovery
**Effort**: Medium | **Risk**: Medium | **Impact**: High

**Changes**:
- In Phase 1, introspect SAGE HDF5 `Snap_0` to discover available fields
- Build field metadata list dynamically instead of hardcoding
- Preserve all discovered fields through pipeline (not just known ones)
- Maintain backward compatibility with current field expectations

**Files**: `sage2kdtree.cc` (Phase 1)
**Benefits**:
- Automatically handles new SAGE fields without code changes
- Aligns with "SAGE HDF5 as single source of truth" principle
- Simplifies maintenance when SAGE schema evolves

#### 2.2 Attribute Preservation
**Effort**: Small | **Risk**: Low | **Impact**: Medium

**Changes**:
- Copy all HDF5 attributes from input fields to output fields
- Preserve units, descriptions, labels from SAGE HDF5 if present
- Document expected vs optional attributes

**Files**: `sage2kdtree.cc` (Phases 3 & 4)

---

### Priority 3: Memory Efficiency (Principle 6)

**Goal**: Ensure memory usage remains bounded for large simulations.

#### 3.1 Memory Usage Profiling
**Effort**: Small | **Risk**: Low | **Impact**: Medium

**Changes**:
- Add optional memory usage reporting at phase boundaries
- Log peak memory consumption per phase
- Identify which phase dominates memory usage

**Files**: `sage2kdtree.cc`
**Implementation**: Use `/proc/self/status` on Linux, `task_info` on macOS

#### 3.2 Streaming Improvements (if needed)
**Effort**: Large | **Risk**: High | **Impact**: High

**Changes** (only if profiling shows issues):
- Process trees in batches rather than loading all into memory (Phase 2)
- Stream snapshot processing in chunks (Phase 3)
- Use HDF5 chunked I/O more extensively

**Files**: `sage2kdtree.cc` (Phases 2 & 3)
**Note**: Defer until profiling demonstrates necessity

---

### Priority 4: Testing & Validation

**Goal**: Ensure correctness, prevent regressions.

#### 4.1 Validation Against Baseline
**Effort**: Medium | **Risk**: Low | **Impact**: High

**Changes**:
- Create comparison script that validates `sage2kdtree` output against multi-step workflow:
  - Compare lightcone/snapshotNNN groups (bounds, splits, cell_counts)
  - Validate data field checksums match
  - Check subsize calculations are consistent
- Automate in CI or as part of test suite

**Files**: New script `tests/validate_sage2kdtree.sh` or Python script

#### 4.2 Unit Tests for Critical Functions
**Effort**: Medium | **Risk**: Low | **Impact**: High

**Changes**:
- Add unit tests for:
  - Depth-first ordering logic
  - BFS/DFS traversal metadata calculation
  - Subsize calculation for known tree structures
  - KDTree partitioning with known spatial distributions
- Use Google Test or Catch2 framework

**Files**: New directory `tests/unit/` with test cases

#### 4.3 Edge Case Testing
**Effort**: Small | **Risk**: Low | **Impact**: Medium

**Changes**:
- Test with edge cases:
  - Empty snapshots
  - Single galaxy
  - Snapshots with no mergers
  - Very deep trees
  - Very wide trees (many progenitors)

**Files**: `tests/sage-model-tests/`

---

### Priority 5: Code Quality & Documentation

**Goal**: Improve maintainability and developer onboarding.

#### 5.1 Inline Documentation
**Effort**: Small | **Risk**: Low | **Impact**: Medium

**Changes**:
- Add detailed comments explaining:
  - Why each phase exists (what transformation it performs)
  - Key algorithms (depth-first recursion, BFS/DFS traversal)
  - Data structure invariants (global_index assignment, subsize meaning)
- Document assumptions about SAGE HDF5 format

**Files**: `sage2kdtree.cc`

#### 5.2 Reduce Code Duplication with libhpc
**Effort**: Small | **Risk**: Low | **Impact**: Medium

**Changes**:
- Replace raw H5 API calls with libhpc wrappers where possible
- Create helper functions for common HDF5 operations:
  - Reading attributes with fallback
  - Creating extensible datasets
  - Copying attributes between groups

**Files**: `sage2kdtree.cc`, potentially new `src/libhpc/h5/utils.hh`

#### 5.3 Logging Consistency
**Effort**: Small | **Risk**: Low | **Impact**: Low

**Changes**:
- Use libhpc logging (`LOGILN`, `LOGBLOCKI`) consistently throughout
- Remove `std::cout` in favor of logging framework
- Add log levels: DEBUG, INFO, WARN, ERROR

**Files**: `sage2kdtree.cc`

---

## Implementation Strategy

### Recommended Sequence

1. **Start with Validation (1.1, 1.2)**: Immediate safety improvements, low risk
2. **Add Testing (4.1)**: Establish baseline for comparison before making changes
3. **Implement Metadata-Driven Fields (2.1, 2.2)**: High-impact architectural improvement
4. **Enhance Documentation (5.1)**: Helps team understand changes
5. **Memory Profiling (3.1)**: Data-driven decision for optimization
6. **Code Quality (5.2, 5.3)**: Cleanup and polish

### Phased Approach

**Phase A (Week 1-2)**: Safety & Validation
- 1.1 Input Validation Enhancement
- 1.2 Enhanced Error Messages
- 4.1 Validation Against Baseline
- 5.1 Inline Documentation

**Phase B (Week 3-4)**: Metadata-Driven Improvements
- 2.1 Dynamic Field Discovery
- 2.2 Attribute Preservation
- 4.3 Edge Case Testing

**Phase C (Week 5-6)**: Quality & Optimization
- 1.3 Progress Reporting
- 3.1 Memory Usage Profiling
- 5.2 Reduce Code Duplication
- 5.3 Logging Consistency
- 4.2 Unit Tests

---

## Success Metrics

1. **Correctness**: Validation script shows 100% match with baseline workflow
2. **Robustness**: Invalid inputs produce clear error messages (not crashes)
3. **Flexibility**: New SAGE fields automatically flow through pipeline without code changes
4. **Documentation**: New developers can understand each phase's purpose from comments
5. **Memory Efficiency**: Memory usage documented and remains bounded

---

## Non-Goals (Out of Scope)

- Changing the 4-phase pipeline structure (it works well)
- Optimizing for speed (correctness first)
- Supporting non-HDF5 formats (VISION.md specifies HDF5 only)
- Parallelizing within phases (current MPI structure adequate)
- GUI or interactive mode

---

## Risks & Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Breaking working pipeline | High | Validate against baseline after each change |
| Dynamic field discovery introduces bugs | Medium | Keep hardcoded fallback, extensive testing |
| Memory optimization too complex | Medium | Profile first, defer if not needed |
| Testing infrastructure overhead | Low | Start simple, expand incrementally |

---

## Open Questions

1. Should intermediate files be optional (command-line flag to skip writing)?
2. Is there a need for checkpoint/resume capability for large datasets?
3. Should we support partial pipeline execution (e.g., start from Phase 3)?
4. What level of MPI parallelism is needed for production datasets?

---

## Alignment with VISION.md

All improvements directly support core principles:

- **Principle 1 (SAGE HDF5 as truth)**: Enhanced by dynamic field discovery (2.1)
- **Principle 3 (Metadata-driven)**: Core focus of Priority 2
- **Principle 6 (Memory efficiency)**: Addressed by Priority 3
- **Principle 8 (Type safety & validation)**: Core focus of Priority 1
- **Quality Attributes (Reliability)**: Addressed by Priority 4 (testing)
- **Quality Attributes (Maintainability)**: Addressed by Priority 5 (documentation)

Every proposed change strengthens adherence to the architectural vision while maintaining the working pipeline.
