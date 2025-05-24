controls_test_env
A control systems testing environment developed on Ubuntu 22.04 LTS. This project leverages CMake as the build system, with dependencies on the Eigen linear algebra library and Google Test (gtest) for unit testing.

Overview
controls_test_env is a C++ project designed to provide a framework for control algorithm development and testing. It includes implementations, utilities, and unit tests to facilitate rapid prototyping and validation of control system components.

Features
Modular C++ codebase for control systems algorithms

Linear algebra support using Eigen

Unit tests powered by Google Test

CMake-based build system for easy compilation and integration

Requirements
Ubuntu 22.04 LTS (or compatible Linux distribution)

C++17 compatible compiler (e.g., g++ >= 9)

CMake >= 3.16

Eigen library (installed on the system or bundled)

Google Test framework (can be added as a submodule or installed separately)

Setup and Build Instructions
Clone the repository:

bash
Copy code
git clone https://github.com/jnancoo1/controls_test_env.git
cd controls_test_env
Create a build directory and navigate into it:

bash
Copy code
mkdir build && cd build
Configure the project with CMake:

bash
Copy code
cmake ..
Build the project:

bash
Copy code
make
Run the unit tests:

bash
Copy code
./tests/controls_test_env_tests
Project Structure
makefile
Copy code
controls_test_env/
├── CMakeLists.txt          # Build configuration
├── include/                # Header files
├── src/                    # Source files
├── tests/                  # Google Test unit tests
└── README.md               # This file
Dependencies
Eigen: A high-performance C++ linear algebra library

Google Test: C++ testing framework

You can install Eigen and Google Test on Ubuntu using:

bash
Copy code
sudo apt-get update
sudo apt-get install libeigen3-dev libgtest-dev
Note: For Google Test, after installing libgtest-dev, you may need to build the library manually:

bash
Copy code
cd /usr/src/gtest
sudo cmake .
sudo make
sudo cp *.a /usr/lib
