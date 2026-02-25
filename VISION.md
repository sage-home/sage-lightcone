# Sage-light-cone Architectural Vision

**Purpose**: Define the architectural principles and design philosophy for Sage-light-cone, an addition to the sage software ecosystem facilitating the creation of light-cone HDF5 datasets.

---

## Vision Statement

sage-light-cone is a **C++ addition**, to the sage software ecosystem.

This architecture enables researchers to reformat their sage hdf5 outputs and add kdtree indexing to facilitate the extraction of light-cone HDF5 datasets for further analysis.

---

## 8 Core Architectural Principles

These principles guide all design decisions and implementation choices in Sage-light-cone:

### 1. SAGE output hdf5 as the single point of truth

**Principle**: The process of reformating with kdtree indexing and eventual light-cone extraction maintains all the meta data available in the sage output hdf5 file.

**Requirements**:
- Core systems (memory management, I/O, tree processing) operate independently of physics.
- Inputs will be any conformant sage HDF5 output.
- The results can always be validated against the standard sage example.

**Benefits**: Enables independent development and simplifies testing, and reduces complexity in core systems.

**In Practice**: The reformatting to kd-tree indexed format is a simple step as is the extraction of light-cone datasets from the kd-tree indexed format.

### 2. Runtime Modularity

**Principle**: Fully implemented in C++.

**Requirements**:
- Builds on macOS, and HPC based on modules (nt.swin.edu.au).

**Benefits**: Provides scientific flexibility for different research questions, and simplifies testing of different sage conformant representations.

**In Practice**: Users can switch from one version of sage to another and apply the same workflow to extract light-cone datasets.

### 3. Metadata-Driven Architecture

**Principle**: The metadata contained in the sage hdf5 output is maintained throughout the transformations to kdtree format hdf5 file and ultimately to the light-cone hdf5 file.

**Requirements**:
- Galaxy properties (e.g., `StellarMass`, `ColdGas`) are defined in metadata the sage output hdf5 file.

**Benefits**: Reduces code duplication, eliminates manual synchronization between different representations, enables build-time optimization, and simplifies maintenance by creating a single source of truth.

**In Practice**: Can be tested on a variety of sage representations.

### 4. Single Source of Truth

**Principle**: Galaxy data has one authoritative representation with consistent access patterns.

**Requirements**:
- No dual property systems or synchronization code.
- All access to galaxy data goes through a unified approach where each property field is its own dataset.

**Benefits**: Eliminates synchronization bugs, simplifies debugging by having a single data path, reduces memory overhead, and improves performance through unified access patterns.

### 5. Unified Processing Model

**Principle**: Sage-light-cone has one consistent, well-understood method for processing spatially indexed galaxies grouped by snapshot.

**Requirements**:
- The sage output is assumed to have ordering of galaxies within snapshots to be based on merger trees.
- Consistent galaxy inheritance and property calculation methods.
- Robust orphan galaxy handling that prevents data loss.
- Clear separation between merger tree structure and physics calculations.
- Spatially indexing galaxies to optimise light-cone extraction

**Benefits**: Eliminates complexity from maintaining multiple processing modes, simplifies validation, reduces bug surface area, and makes the system easier to understand and modify.

### 6. Memory Efficiency and Safety

**Principle**: Memory usage is bounded, predictable, and safe.

**Requirements**:
- Memory usage is bounded and does not grow with the total number of forests processed.
- Memory management for galaxy arrays and properties is automatic where possible.
- Memory is allocated on a per-forest scope with guaranteed cleanup after processing.
- Tools for memory leak detection and prevention are built-in.

**Benefits**: Allows reliable processing of large simulations, reduces debugging overhead by preventing memory-related bugs, improves performance predictability, and enables processing of datasets larger than available RAM.

**In Practice**: Memory for a merger tree (halos, galaxies, module-specific data) is allocated at the start of processing and guaranteed to be freed upon completion.

### 7. Format-Agnostic I/O

**Principle**: Sage-light-cone input/output formats are universally based on hdf5.

**Requirements**:
- A common, abstract interface for all tree reading operations.
- A property-based output system that adapts to data available at runtime.
- Proper handling of cross-platform issues like endianness.
- Graceful fallback mechanisms for unsupported features in certain formats.

**Benefits**: Ensures scientific compatibility with different simulation codes and analysis tools, future-proofs against format changes, simplifies validation across formats, and eases integration with external tools.

### 8. Type Safety and Validation

**Principle**: Data access is type-safe with automatic validation.

**Requirements**:
- Type-safe property accessors (macros or functions) are generated from metadata.
- Automatic bounds checking and validation where appropriate.
- Fast failure with clear error messages upon invalid data access.

**Benefits**: Reduces runtime errors by catching problems at compile-time, improves debugging with clear messages, catches problems early, and increases confidence in scientific accuracy.

---

## Implementation Philosophy

### Metadata-Driven Development
- **Single Source of Truth**: The sage output hdf5 is the single source of truth.

### Type Safety First
- **Compile-Time Validation**: Errors are caught at compile-time rather than runtime.
- **Generated Access Patterns**: Type-safe property access is generated from metadata.
- **IDE Integration**: Full autocomplete, go-to-definition, and refactoring support in modern IDEs.

### Standard Tools
- **Industry Standards**: Proven tools (CMake, HDF5) rather than custom solutions.
- **Professional Workflow**: Modern development environment with IDE integration.
- **Debugging Support**: All standard debugging tools (GDB, Valgrind) work out of the box.

## System Architecture

### Component Structure


### Data Flow

## Quality Attributes

### Maintainability
- **Modularity**: Clear separation of concerns with well-defined interfaces.
- **Documentation**: Comprehensive documentation for developers and users.
- **Code Quality**: Professional coding standards with consistent style.  Adherence to formatting style as defined in ".clang-format"
- **Testing**: Automated testing covering all major functionality.

### Reliability
- **Error Handling**: Robust error detection and recovery mechanisms.
- **Validation**: Comprehensive input and state validation.
- **Memory Safety**: Automatic memory management preventing leaks and corruption.
- **Debugging**: Clear error messages and debugging capabilities.

### Usability
- **Debugging Support**: Clear diagnostics for troubleshooting problems.

---

## Design Philosophy

This architectural vision defines Sage-light-cone as a modern, maintainable scientific software framework. The core principles provide a clear foundation for all development decisions, ensuring the system remains focused on its primary goals: scientific accuracy, developer productivity, and long-term maintainability.

The key insight is that **scientific accuracy and architectural elegance are not mutually exclusive**. By applying proven software engineering principles to scientific computing, Sage-light-cone accelerates scientific discovery through improved flexibility, reliability, and maintainability.

Success is measured not just by technical metrics, but by the system's ability to enable new scientific insights through easier experimentation, more reliable results, and faster development of new physics models.
