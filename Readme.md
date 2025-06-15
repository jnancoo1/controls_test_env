# C++ Controls Library

A C++ library for control systems engineering, providing tools for state space analysis, system transformations, and mathematical operations.

I was kinda doing this on kicks
## Features

- **State Space Analysis**
  - Controllability and observability analysis
  - Stability analysis for discrete and continuous-time systems
  - Stabilizability and detectabilities checks using PBH test
  - Lyapunov Solvers for grammians 

- **Canonical Forms**
  - Controllable canonical form
  - Observable canonical form
  - Phase variable form
  - Schur form
  - System diagonalization
  - Kalman Decom 

- **Linear Algebra Operations**
  - Matrix operations
  - Vector operations
  - Linear system solvers
  - LU and QR decomposition


- Ubuntu Server 24.02 LTS
- G++ 13.3
- CMake 3.28.3
- Eigen (Linear algebra library)
- GTest (For unit tests)

### Installing Dependencies

On Ubuntu 24.02 LTS, you can install Eigen and GTest using apt:

```bash
# Install Eigen
sudo apt-get update
sudo apt-get install libeigen3-dev

# Install GTest
sudo apt-get install libgtest-dev
```

## Installation

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


Detailed documentation is available in the generated \`refman.pdf\` file. This includes:
- Comprehensive API documentation
- Mathematical formulations
- Usage examples (comming soon)
- Implementation details

To generate fresh documentation:
```bash
doxygen
cd latex
make
```


```cpp
#include "discrete_state_space.hpp"
#include "analysis.hpp"

int main() {
    // Create a state space system
    Discrete_StateSpace_System system;
    System.A=(1,2;3,4)l
    System.B=(1;1);
    System.C=(1,0;1,0);
    // Analyze system properties
    Analysis analyzer;
    bool isControllable = analyzer.is_controllable(system);
    bool isObservable = analyzer.is_observable(system);
    
    // Transform to canonical form
    Forms transformer;
    auto ccf = transformer.Cont_Cannonical_form(system);
    
    return 0;
}
```


This library is under active development. Upcoming features include:
- Additional system analysis tools
- Controller Synthesis
- Frequency Domina tools
- Enhanced numerical stability
- Model Order Reduction
- Additional test coverage
- Performance optimizations


- OS: Ubuntu Server 24.02 LTS
- Compiler: G++ 13.3
- Build System**: CMake 3.28.3
- Dependencies:
  - Eigen (Matrix operations)
  - GTest (Testing framework)


## 👥 Contributing

This is a work in progress. Contributions, suggestions, and feedback are welcome.
