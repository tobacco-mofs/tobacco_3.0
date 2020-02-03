# tobacco_3.0

# Authors
- Ryther Anderson
- Yamil Colón
- Diego Gómez-Gualdrón

# Motivation
Our Topologically Based Crystal Constructor (ToBaCCo) was developed to rapidly produce molecular representations of porous crystals as crystallographic information (.cif) files, which can then be used for molecular simulation or for materials characterization. 

# Installation and Dependencies
We recommend building a Python 2.7 environment using Anaconda. You will need the following packages:
- numpy
- networkx
- scipy

To set up a Python 2.7 environment first install Anaconda (if you install Python 2.7 as your base environment there is no need to set up a new environment for ToBaCCo). After installing Anaconda, run:
```
conda create --name my_tobacco python=2.7
```
where my_tobacco is the name of the environment, which can be changed to whatever you like. Next, load your new environment:
```
conda activate my_tobacco
```
(this command changes for Windows, see https://docs.anaconda.com/anaconda/user-guide/tasks/switch-environment/)
Finally, install these packages just run the following with the environment loaded:
```
conda install numpy
conda install networkx
conda install scipy
```

# Usage
Execute the tobacco.py file to run ToBaCCo:
```
python tobacco.py
```
For more details on ToBaCCo inputs, outputs, and configuration see the ToBaCCo_manual.pdf included with the repository.

# License
GNU General Public License (can be viewed in the LICENSE file included in this repository)
