# farm2d

A setup to run 2D RANS simulations of wind farm flows in OpenFOAM.

Since `farm2d` is relatively fast, you can run the simulations on a regular laptop. This simulation of Horns Rev 1 ran in 96 seconds on my laptop (using 4 cores).

![](imgs/hornsrev1_u_contour.png)

## Installation

Assuming that your git projects are located at `~/git` on your computer:
```
cd git
git clone https://github.com/mchba/farm2d.git
```

*Or download the repository as a ZIP-file or any other way you like*.

### OpenFOAM

The current setup has only been tested with v2206. You can see this [video](https://www.youtube.com/watch?v=CeEJS1eT9NE&t=477s) (start around 8:00) for a concise guide on how to download OpenFOAM.

To check that you have installed OpenFOAM correctly, type

```
simpleFoam -help
```

in the terminal, which, if successful, should print some information about the simpleFoam solver.

### Actuator Disk (AD) compilation

The built-in AD code of OpenFOAM uses a "monitor"-method to determine the freestream velocity, which does not make sense for wind farm studies. Therefore, `farm2d` instead uses an AD based on 1D momentum control.

To be able to use this AD, it must first be compiled and linked to your OpenFOAM installation.

1. Activate your OpenFOAM (if `simpleFoam -help` works, it is activated) and check if you have a OpenFOAM user directory:

   `cd $WM_PROJECT_USER_DIR`

    If you don't have this directory, create it with `mkdir -p $WM_PROJECT_USER_DIR/{run,applications,src}`. The user directory is used to add custom code, see more info [here](https://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2022/lectureNotes/01_initialPreparations.pdf).


2. `cp AD_calaf $WM_PROJECT_USER_DIR/src/ -r`.
3. `cd $WM_PROJECT_USER_DIR/src/AD_calaf`
4. `wmake`. This compiles the AD code and creates a library called `lib_ad_calaf.so` in the `$WM_PROJECT_USER_DIR/platforms` folder.

If you look at the examples of `farm2d`, you will find that they all have a reference to the `lib_ad_calaf.so`-file in their `system/controlDict`-file.

*This AD can also be used in 3D simulations by adjusting diskArea in constant/fvOptions and the turbine location in system/topoSetDict.*

### Python (optional)

It is recommended to use a virtual Python environment, for example through [Miniconda](https://docs.anaconda.com/miniconda/). Assuming you already have a conda installation, create a new environment as:

```
conda create -n farm2d python spyder numpy scipy xarray matplotlib pyvista
conda activate farm2d
```

Also add `farm2d` to your `~/.bashrc`-file:

```
# Add this somewhere in your bashrc-file:
export PYTHONPATH="${PYTHONPATH}:/home/<username>/git/farm2d/src"
```

## How to run

### Native OpenFOAM

You can run any of the examples with "pure" OpenFOAM. This does not require any Python installation.
```
cd examples
cp V80 my_first_test -r
cd my_first_test
./Allrun
```

This case runs in about 3 seconds on my laptop (using 4 cores). You can take a look at the results with for example `ParaFoam` or `Paraview`.

```
paraview test.foam &
```

![](imgs/paraview_u_contour.png)





### With Python

The main advantage of running through Python is that you can easily modify the example cases and automate various tasks.

```
cd examples
cp V80 my_second_test -r
cd my_second_test
python setup_case.py
```

A contour plot was automatically created with the python script:

![](imgs/farm2d_u_contour.png)


## License

This project is released under the MIT License.

You are free to use, modify, and distribute the code for academic purposes, commercial projects, personal experiments, or any other purpose. However, the code is provided "as-is," without any warranty of any kind. The developers of this project are not liable for any damages or issues arising from the use of this software. Users bear sole responsibility for the results obtained or consequences arising from its use.


