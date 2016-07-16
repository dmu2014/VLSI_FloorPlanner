The source code is written in C++ and consists of the following files - 
hfile.h - Header file 
classdefs.cpp - Contains all class definitions
floorplanner.cpp - main function

Compilation - 
Copy above 3 files in a working directory and invoke the g++ compiler.
g++ *.cpp -o ExecutableName
Example:
g++ *.cpp -o floorplanner

Usage - 
After compilation, run the executable created as follows. The benchmark files 
need to be in the working directory.
 
./floorplanner Benchmark_Name Stop_Temperature
Example: 
./floorplanner B10 0.1

The expected output .pl file is generated in the working directory.

Note: All results reported are with a stop temperature of 0.1. You can have a higher stop 
temperature for reduced execution time but that will generate a lower quality result.
Dependence of temperature on Quality of Solution is discussed in the report.

