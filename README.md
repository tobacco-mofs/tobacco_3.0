# tobacco_3.0

# Authors
- Ryther Anderson
- Yamil Colón
- Diego Gómez-Gualdrón

# Motivation
Our Topologically Based Crystal Constructor (ToBaCCo) was developed to rapidly produce molecular representations of porous crystals as crystallographic information (.cif) files, which can then be used for molecular simulation or for materials characterization. 

# Installation and Dependencies
ToBaCCo has been tested for Python 2.7 and Python 3.7. We recommend building a Python 2.7/3.7 environment using Anaconda. You will need the following packages:
- numpy
- networkx
- scipy

To set up a Python X.7 environment first install Anaconda. After installing Anaconda, run:
```
conda create --name my_tobacco python=X.7
```
where my_tobacco is the name of the environment and X is 2 or 3, which can be changed to whatever you like. Next, load your new environment:
```
conda activate my_tobacco
```
(this command changes for Windows, see https://docs.anaconda.com/anaconda/user-guide/tasks/switch-environment/)
Finally, install the above packages, thus:
```
conda install numpy
conda install networkx
conda install scipy
python setup.py install 
```
Once these are installed you can clone or download the repository and start running ToBaCCo.

# Usage
Execute the tobacco command to run ToBaCCo:
```
tobacco -h
```
it will shows

```
usage: tobacco [-h] [-c CONFIG] [-n NODE] [-e EDGE] [-t TEMPLATE]

Topologically Based Crystal Constructor (ToBaCCo) was developed to rapidly
produce molecular representations of porous crystals as crystallographic
information (.cif) files, which can then be used for molecular simulation or
for materials characterization. To see the options, type "tobacco -h".

optional arguments:
  -h, --help            show this help message and exit
  -c CONFIG, --config CONFIG
                        configure file
  -n NODE, --node NODE  node directory
  -e EDGE, --edge EDGE  edge directory
  -t TEMPLATE, --template TEMPLATE
                        template directory
```

to run a real example , you may to to example directory, then use following command
```
mkdir cifs
tobacco 
```
or
```
mkdir cifs
tobacco -c config.yaml -n nodes -e edges -t templates
```

For more details on ToBaCCo inputs, outputs, and configuration see the ToBaCCo_manual.pdf included with the repository.

# License
GNU General Public License (can be viewed in the LICENSE file included in this repository)
