# Phonon-pumping in magnetic|nonmagnetic bilayers
This program generates data file and [gnuplot](https://gnuplot.sourceforge.net) script for producing FMR power absorption spectra discussed in [our paper](https://doi.org/10.1103/PhysRevB.104.014403).
Specifically,
1. Loads material parameters and geometry from files. 'Geometry' includes plot settings such as magnetic field and frequency range.
2. Calculates the FMR power absorption spectrum according to Eq. (36).
3. Within the plot range, search for the frequencies of Fig. 3(c2) and (c3), which are expected to enhance Gilbert damping via magnetoelastic coupling.
4. Writes a data file of the FMR spectrum and a [gnuplot](https://gnuplot.sourceforge.net) script for producing Figs. 3(d)(e) and 4(b)(d).

## How to compile executable
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
2. download a zip containing the keyword `uvcrt` from [MinGW binaries](https://github.com/niXman/mingw-builds-binaries/releases), extract them and move the folder `mingw64` under `C:\Program Files (x86)`. Add `C:\Program Files (x86)\mingw64\bin` to PATH.
In the project root folder, run:
```sh
mkdir build
cd build
cmake -S .. -B . -G "MinGW Makefiles"
mingw32-make
```

`CMakeLists.txt` determines whether the build should be Debug or Release mode. Look for the variable `DEBUG_MODE`.
The executable will be saved under `build/bin/`.

NOTE: On Windows, the executable might require the dynamic libraries `mingw64\bin\libstdc++-6.dll` and `libgcc_s_seh-1.dll` to be present in the same directory. In that case, copy-paste the dlls to the directory.


## How to use
The material parameters are specified in a material parameter file. The layer thicknesses, magnetization angle, ac magnetic field strengths, and the ranges of static magnetic field and frequency are specified in a geometry file. You need to tell the program where these two files are.

1. Prepare material parameter file
2. Prepare geometry file
3. Run the program in Command Prompt/PowerShell (Windows) or Terminal (Linux/Mac), specifying the paths to those two files:

(Linux, Mac)
```sh
./PhononPumpingPlotter -m path/to/material_parameter_file -g path/to/geometry_file
```
(Windows)
```sh
.\PhononPumpingPlotter.exe -m path\to\material_parameter_file -g path\to\geometry_file
```

For example on Windows, run the following from the project root directory:
```sh
.\build\bin\release\PhononPumpingPlotter.exe -m .\examples\material_parameters_YIG-GGG.conf -g .\examples\geometry_Sato2021_fig3d.conf
```

4. Run the gnuplot script. A PDF file will be generated in the same directory.
5. If needed, you can modify the gnuplot script manually to optimize the visualization.
