# Polymer Physics Demo
This folder contains demos I made for teaching polymer physics.

# Folders:
- `core`
    - `run_pe.py`: running polyethylene
    - `run_kremer_grest.py`: running coarse-grained kremer grest model
    - `run_fjc.py`: running freely jointed chain
    - `polymer.py`: provide helpful class
- `utils`
    - `align_trajectory.py`: align the trajectory to fix the first monomer, and maybe rotate the end to end vector to y-axis
- `tests`
    - `manual`
        - `FJC`: the `run.sh` is an example script for running FJC simulation
        - `PE`: the `run.sh` is an example script for running polyethylene simulation

# Thigs to set up
- install anaconda or miniconda from [here](https://conda.io/en/latest/miniconda.html)
- create a separate virtual environment `PP`:
```zsh
conda config --add channels conda-forge 
conda create --name pp python "numpy<1.24" scipy matplotlib pandas MDAnalysis seaborn
conda activate pp
conda install networkx -c anaconda
conda install -c conda-forge openmm cudatoolkit=YOUR_CUDA_VERSION 
```
for openmm installation, please see [documentation](http://docs.openmm.org/latest/userguide/application/01_getting_started.html#installing-openmm) for more details
- install [`git`](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
- navigate to a folder where you want to download this repository
- clone this repository `git clone https://github.com/yihengwuKP/polymer.git`
- set the `D_PP` variable to the path to this folder. `export D_PP=/path/to/PP`


# More guidance
1. `${D_PP}/core/run_pe.py -h` or `${D_PP}/core/run_fjc.py -h` will provide help on the command line options, same as the `${D_PP}/utils/align_trajectory.py -h`
