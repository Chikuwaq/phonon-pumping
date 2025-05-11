# Phonon-pumping in magnetic|nonmagnetic bilayers
This program generates data file and gnuplot script that generate FMR power absorption spectra.

## Compilation
Build depends on 
* `cmake`
* some `make` tool
* C++ compiler and linker

On Mac, `cmake` can be installed by `Homebrew`. Xcode command line tools contains C++ build tools.
In the project root folder, run:
```sh
mkdir build
cd build
cmake -S .. -B .
make
```

On Windows, 
1. install `cmake`
2. download a zip of [MinGW binaries](https://github.com/niXman/mingw-builds-binaries/releases), extract them and move the folder `mingw64` under `C:\Program Files (x86)`. Add `C:\Program Files (x86)\mingw64\bin` to PATH.
In the project root folder, run:
```sh
mkdir build
cd build
cmake -S .. -B . -G "MinGW Makefiles"
mingw32-make
```
NOTE: The executable might require the dynamic libraries `mingw64\bin\libstdc++-6.dll` and `libgcc_s_seh-1.dll` to be present in the same directory. In that case, simply copy-paste the dlls.

## Usage
The material parameters are specified in a material parameter file. The layer thicknesses, magnetization angle, ac magnetic field strengths, and the ranges of static magnetic field and frequency are specified in a geometry file. You need to tell the program where these two files are.

1. Prepare material parameter file
2. Prepare geometry file
3. Run the program specifying the paths to those two files:
(Linux, Mac)
```sh
./PhononPumpingPlotter -m path/to/material_parameter_file -g path/to/geometry_file
```
(Windows)
```sh
.\PhononPumpingPlotter.exe -m path\to\material_parameter_file -g path\to\geometry_file
```
If the two input files are in the same directory as the executable, this simplifies to:
(Windows)
```sh
.\PhononPumpingPlotter.exe -m .\material_parameters.txt -g .\geometry.txt
```
4. Run the gnuplot script. A PDF file will be generated in the same directory.
5. If needed, you can modify the gnuplot script manually to optimize the visualization.
