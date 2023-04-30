# Read Mapping Project

### Description

This project implements a fast and accurate read mapping algorithm for DNA sequences using C++.

The following algorithms are used:
- Suffix array construction using the SAIS algorithm.
- Smith-Waterman algorithm for string alignment. 

### Environment

I developed the project on an M1 Macbook Pro using Apple's clang g++ compiler (version 14.0.3). The project should be able to be compiled on any system with a C++ compiler that supports C++17.

### Compiling

To compile the project, run the following command in the project directory:

```bash
mkdir build
cd build
cmake ..
make
cd ..
```

Two binaries will be generated:

- `debug` - a debug build of the project, with debug symbols and no optimizations.
- `release` - a release build of the project, with -O3 optimizations. Please use this binary to process large inputs.

### Running

The program takes two arguments: the path to the reference genome, and the path to the reads file. Under the project directory, run the following command:

```bash
mkdir output
./build/release <path to reference genome> <path to reads file>
```

The output will be written to `output/output.txt`.

To reformat the output to the submission .zip format, run the following bash script:

```bash
./format_output.sh
```

### References and Acknowledgements

When developing the project, I acknowledge using the following sources in learning relevant algorithms and adapting codes from them.

- [SAIS algorithm](https://zork.net/~st/jottings/sais.html) for suffix array construction. 
- [SAIS algorithm C++ implementation](https://github.com/Tascate/Suffix-Array-Implementation).
- [C++17 parallelized std algorithms](https://github.com/mikekazakov/pstld) for Apple clang compiler.
- [Smith-Waterman Algorithm](https://www.wikiwand.com/en/Smith%E2%80%93Waterman_algorithm) for string alignment. 
