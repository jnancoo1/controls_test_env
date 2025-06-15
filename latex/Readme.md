# C++ Controls Library

A modern C++ library for control systems engineering, providing tools for state space analysis, system transformations, and mathematical operations.

## 🚀 Features

- **State Space Analysis**
  - Controllability and observability analysis
  - Stability analysis for discrete and continuous-time systems
  - Stabilizability checks using PBH test

- **Canonical Forms**
  - Controllable canonical form
  - Observable canonical form
  - Phase variable form
  - Schur form
  - System diagonalization

- **Linear Algebra Operations**
  - Matrix operations
  - Vector operations
  - Linear system solvers
  - LU and QR decomposition

## 🛠️ Prerequisites

- Ubuntu Server 24.02 LTS
- G++ 13.3
- CMake 3.28.3
- Eigen (Linear algebra library)
- GTest (For unit tests)

## 📦 Installation

1. Clone the repository:
\`\`\`bash
git clone [your-repository-url]
cd controls_library
\`\`\`

2. Create a build directory:
\`\`\`bash
mkdir build
cd build
\`\`\`

3. Configure and build:
\`\`\`bash
cmake ..
make
\`\`\`

## 📚 Documentation

Detailed documentation is available in the generated \`refman.pdf\` file. This includes:
- Comprehensive API documentation
- Mathematical formulations
- Usage examples
- Implementation details

To generate fresh documentation:
\`\`\`bash
doxygen
cd latex
make
\`\`\`

## ⚡ Quick Example

\`\`\`cpp
#include "discrete_state_space.hpp"
#include "analysis.hpp"

int main() {
    // Create a state space system
    Discrete_StateSpace_System system;
    
    // Analyze system properties
    Analysis analyzer;
    bool isControllable = analyzer.is_controllable(system);
    bool isObservable = analyzer.is_observable(system);
    
    // Transform to canonical form
    Forms transformer;
    auto ccf = transformer.Cont_Cannonical_form(system);
    
    return 0;
}
\`\`\`

## 🚧 Work in Progress

This library is under active development. Upcoming features include:
- Additional system analysis tools
- More canonical forms
- Enhanced numerical stability
- Additional test coverage
- Performance optimizations

## ⚙️ Build Environment

- **OS**: Ubuntu Server 24.02 LTS
- **Compiler**: G++ 13.3
- **Build System**: CMake 3.28.3
- **Dependencies**:
  - Eigen (Matrix operations)
  - GTest (Testing framework)

## 📝 License

[Your chosen license]

## 👥 Contributing

This is a work in progress. Contributions, suggestions, and feedback are welcome.
