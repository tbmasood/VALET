# ElecTrans

A python library for analysis of natural transition orbitals of molecules and molecular complexes. It supports computing electronic charges and charge transfer during electronic transitions at atomic as well as at subgroup level.

# Installation
TODO

# Use as a python package
After installation, the whole package can be imported in your python scripts using the command `import electrans`. However, in most cases only `from electrans import transition_analysis_utils as tau` would be enough as this conatins all the utility fuctions required for charge computation and creating transition diagrams . See [`test_electrans.py`](test_electrans.py) for an example use case.

# Using through Jupyter notebook
The package can also be imported in a Jupyter lab and notebook. See [`plot_transition_diagram.ipynb`](plot_transition_diagram.ipynb) for an example use case. 

# Command line execution
After installation, the library can also be executed from command line as follows:

`
electrans [-h] [-p] [-q] [-d] [-a] [-s] [-av] [-sv] [-t THREADS] input_dir
`

### Positional arguments:
|Argument | Description|
|---|---|
| input_dir | The directory containing the input cube files for transitions. The directory should a 'metadata.csv' file which contains 3 columns. The first specifies the name of the transition. The second and third columns specify hole and particle cube files respectively.|

### Optional arguments
|Option | Description|
|---|---|
|-h, --help | Show help message and exit|
|-p, --different_atom_pos| If specified, Voronoi diagram will be computed for each cube file separately. By default it is assumed, the atomic positions remain the same for all the cubefiles.
|-q, --use_quadratic_optimization | Use quadratic programming based optimization approach to compute the charge trasfer between subgroups. By default a hueristic approach is used for this purpose.|
|-d, --output_transition_diagram | If specified, transition diagram will be generated forall the transitions in `<input_dir>/results/transition_diagrams/` directory.|
|-a, --output_atomic_charges | If specified, atomic charges for each transition will be output as a csv file in the directory `<input_dir>/results/atomic_charges/`.|
|-s, --output_subgroup_charges |If this option is set, subgroup charges and amount of charge transfer for each transition would be output in a directory named `<input_dir>/results/subgroup_charges/`.|
|-av, --output_atoms_vtk | If this option is set, the atoms are saved as 3D model in VTK format. A file is generated for each transition which contains the data about the hole and particle charge, charge difference, subgroup, etc. These VTK files are saved in `<input_dir>/results/vtk/atoms/` and can be loaded in VTK compatible software like Paraview.|
|-sv, --output_segmentation_vtk | If specified, the computed Voronoi segmentation at atomic and subgroup scales are saved as VTK compatible files in the directory `<input_dir>/results/vtk/segmentation/`. These files can be loaded in VTK compatible software like Paraview.|
|-t THREADS, --threads THREADS | Specify the number of parallel threads to be used in computing the Voronoi diagram. By default four threads are used.|