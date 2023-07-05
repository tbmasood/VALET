# VALET: Visual Analysis of Electronic Transitions

A python library for analysis of Natural Transition Orbitals (NTO) of molecules and molecular complexes. It supports computing electronic charges and charge transfer during electronic transitions at atomic as well as at subgroup level. The methods implemented in this library are described in the following paper:

Masood, T.B., Thygesen, S., Linares, M., Abrikosov, A.I., Natarajan, V. and Hotz, I. (2021), Visual Analysis of Electronic Densities and Transitions in Molecules. Computer Graphics Forum, 40: 287-298. [https://doi.org/10.1111/cgf.14307]

```
@article{MasoodCGF2021,
    author = {Masood, T. Bin and Thygesen, S.S. and Linares, M. and Abrikosov, A. I. and Natarajan, V. and Hotz, I.},
    title = {Visual Analysis of Electronic Densities and Transitions in Molecules},
    journal = {Computer Graphics Forum},
    volume = {40},
    number = {3},
    pages = {287-298},
    doi = {https://doi.org/10.1111/cgf.14307},
    url = {https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.14307},
    year = {2021}
}
```

# Installation

I recommend using a virtual environment to install this package and its dependencies. You can use environments like `miniconda`, `anaconda`, `virtualenv`, etc. Here, I demonstrate the procedure for `virtualenv`.

First check if you have `virtualenv` installed using the command `virtualenv --version`. It should show the currently installed version. Otherwise you have to install `virtualenv` according to your platform requirements or use `pip` to install `virtualenv` locally using the following command. I recommend the latter option for installing `virtualenv`.

`
$ pip install virtualenv
`

If `python3` and `virtualenv` are already installed on your system, then go ahead and create a virtual enviroment called `valet-venv` in the directory where you have downloaded this repo as follows:

`
$ virtualenv valet-venv
`

Activate the virtual environment using the command:

`
$ source ./valet-venv/bin/activate
`

You should see the prompt change accordingly. Then install `valet` package in the virtual environment along with its dependencies using the following command:

`
(valet-venv) $ python -m pip install -e .
`

The command above installs the package in the current directory which is `valet` in our case. Note that `valet` also needs `vtk` library if you are interested in generating the 3D output which can be loaded in VTK compatible software like Paraview. However, it is not installed by default. If you need the VTK file outputs, install `vtk` as follows, otherwise skip this step:

`
(valet-venv) $ python -m pip install -e .[vis]
`

You are then ready to use `valet` package in your python scripts or through command line, for example as:

`
(valet-venv) $ valet -d ./data/
`

Finally, you can deactivate and exit the virtual environment as follows:

`
(valet-venv) $ deactivate
`

We discuss the various ways to use this package in more detail below.

## Use as a python package
After installation, the whole package can be imported in your python scripts using the command `import valet`. However, in most cases only `from valet import transition_analysis_utils as tau` would be enough as this conatins all the utility fuctions required for charge computation and creating transition diagrams . See [`test_valet.py`](test_valet.py) for an example use case.

## Using through Jupyter notebook
The package can also be imported in a Jupyter lab and notebook. See [`plot_transition_diagram.ipynb`](plot_transition_diagram.ipynb) for an example use case. 

You can change the installation option to automatically install Jupyter lab in your virtual environment as follows:

`
(valet-venv) $ python -m pip install -e .[jupyter]
`
Then launch jupyter lab:
`
(valet-venv) $ jupyter-lab
`


## Command line execution
After installation, the library can also be executed from command line as follows:

`
python -m valet [-h] [-p] [-q] [-d] [-a] [-s] [-av] [-sv] [-t THREADS] input_dir
`

### Positional arguments:
|Argument | Description|
|---|---|
| input_dir | The directory containing the input cube files for transitions. The directory should a 'metadata.csv' file which contains 3 columns. The first specifies the name of the transition. The second and third columns specify hole and particle cube files respectively.|

### Optional arguments
|Option | Description|
|---|---|
|-h, --help | Show help message and exit. |
|-p, --different_atom_pos| If specified, Voronoi diagram will be computed for each cube file separately. By default it is assumed, the atomic positions remain the same for all the cubefiles and hence segmentation is computed only once for all the cubefiles. |
|-q, --use_quadratic_optimization | Use quadratic programming based optimization approach to compute the charge trasfer between subgroups. By default a hueristic approach is used for this purpose. |
|-d, --output_transition_diagram | Generate transition diagram and store all the diagrams in `<input_dir>/results/transition_diagrams/` directory. |
|-a, --output_atomic_charges | Save atomic charges for each transition as a csv file in the directory `<input_dir>/results/atomic_charges/`. |
|-s, --output_subgroup_charges |If this option is set, subgroup charges and amount of charge transfer for each transition would be output in a directory named `<input_dir>/results/subgroup_charges/`. |
|-av, --output_atoms_vtk | If this option is set, the atoms are saved as 3D model in VTK format. A file is generated for each transition which contains the data about the hole and particle charge, charge difference, subgroup, etc. These VTK files are saved in `<input_dir>/results/vtk/atoms/` and can be loaded in VTK compatible software like Paraview. |
|-sv, --output_segmentation_vtk | If specified, the computed Voronoi segmentation at atomic and subgroup level are saved as VTK files in the directory `<input_dir>/results/vtk/segmentation/`. These files can be loaded in VTK compatible software like Paraview. |
|-t THREADS, --threads THREADS | Specify the number of parallel threads to be used in computing the Voronoi diagram. By default four threads are used. |
