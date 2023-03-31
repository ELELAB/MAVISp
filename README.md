Cancer Systems Biology, Technical University of Denmark, 2800, Lyngby, Denmark  
Cancer Structural Biology, Danish Cancer Society Research Center, 2100, Copenhagen, Denmark  
Repository associated to the publication:  
> MAVISp: Multi-layered Assessment of VarIants by Structure for proteins  
> Matteo Arnaudi, Ludovica Beltrame, Kristine Degn, Mattia Utichi, Alberto Pettenella, Simone Scrima, Peter Wad Sackett, Matteo Lambrughi, Matteo Tiberti, Elena Papaleo  
> biorxiv, doi: https://doi.org/10.1101/2022.10.22.513328


# MAVISp web app

## Introduction

This is the web app of the MAVISp database for structural variants annotation.

If you use MAVISp, please cite our preprint:

> MAVISp: Multi-layered Assessment of VarIants by Structure for proteins  
> Matteo Arnaudi, Ludovica Beltrame, Kristine Degn, Mattia Utichi, Alberto Pettenella, Simone Scrima, Peter Wad Sackett, Matteo Lambrughi, Matteo Tiberti, Elena Papaleo  
> biorxiv, doi: https://doi.org/10.1101/2022.10.22.513328

Please see the [CHANGELOG.md](CHANGELOG.md) file on this repository for current 
and expected data releases as well as updates on the software.

Please see the [CURATORS.md](CURATORS.md) file in this repository for an up to
date list of our curators.

This web app uses by default a database based on a set of CSV files, available
on [our OSF repository](https://osf.io/ufpzm/). These were generated by using
the MAVISp Python package from raw data files. See below for details about
the Python package.

## Requirements

Running the MAVISp web app requires a working Python 3.7+ installation with the following
Python packages installed:

- streamlit
- streamlit-aggrid
- pandas
- osfclient

It has been tested on Linux and macOS, but in principle runs on any operating system
that supports Python.

The easiest way to run MAVISp is by means of a Docker container - see the next section.

If you would rather install and run it on your local hardware, please follow the instructions
below.

## Running the app as a Docker container

coming soon - please use one of the options below for now

## Installing requirements

These instructions apply to both Linux and macOS, using the terminal. You
will need to have a recent (>=3.7) Python distribution installed on your system
or [Anaconda](https://anaconda.org).

### Installing requirements using a virtualenv Python environment

1. if you have access to the `virtualenv` Python environment module, you can create a
virtual environment:

```
virtualenv -p python3.8 MAVISp_env
```

2. then, activate it:

```
source MAVISp_env/bin/activate
```

3. you can install the requirements in the environment using `pip`:

```
pip install pandas streamlit streamlit-aggrid pyyaml osfclient
```

### Installing requirements using a conda Python environment

1. if you have access to Anaconda or Miniconda (executable `conda`), you can use it
to create a virtual environment:

```
conda create -n MAVISp_env python==3.8 streamlit pandas pyyaml
```

2. then you can activate it:

```
conda activate MAVISp_env
```

3. you need to install two of the requirements, which is currently not on conda, using `pip`:

```
pip install streamlit-aggrid osfclient
```

## Downloading the required files

You will need to download the MAVISp web app from GitHub and the current MAVISp
database from our OSF repository:

1. If you haven't already, activate your Python environment (see previous steps)

2. create a local copy of the MAVISp repository in your system:

```
git clone https://github.com/ELELAB/MAVISp
```

3. Download the database files from [our OSF repository](https://osf.io/ufpzm/).
If you have installed the `osfclient` Python package (see requirements), just
enter the directory of the MAVISp repository you have created in the previous step
and download the files, as follows:

```
cd MAVISp
osf clone && mv ufpzm/osfstorage/database/ . && rm -r ufpzm
```

alternatively, you can download the whole OSF database from the web interface
and copy it in the `MAVISp` folder.

At the end of the process, you should have the OSF `database` folder and its
contents inside the `MAVISp` folder.

## Running the app

In the following instructions
  - `hostname` denotes the hostname of the server (i.e. the one you usually ssh to)
  - `user` denotes your username on the server

### Running the app locally

With your Python environment still active and from inside the `MAVISp` folder, run:

```
cd MAVISp
streamlit run app.py
```

a browser window displaying the MAVISp web app should open.

### Running the app remotely - if you have access to a host

the steps are the same as the last section, but a browser window will not open.
Instead, you will have to connect from your local browser to the
host printed out by streamlit in the terminal.

### Running the app remotely - with X forwarding

This option is useful for who wants to run MAVISp remotely, but doesn't have
direct access to the MAVISp web service (e.g. because the outbound port that
streamlit uses is blocked). It is however pretty slow and clunky, meaning it
is not very effective for in-depth data exploration or analysis. It also requires
to be able to connected with X forwarding to your server. For a faster 
alternative that doesn't use X-forwarding see instructions below.

1. connect to your server via ssh, with X forwarding:

```
ssh -XY user@hostname
```

2. Follow the previous instructions to install requirements and run the app locally.
A browser window should pop up.

### Running the app remotely - without X forwarding

You can use an ssh tunnel to connect to your streamlit instance.
This is slightly more complicated but leads to an overall much better
user experience.

1. ssh without X forwarding to the server:

```
ssh user@hostname
```

2. Follow the previous instructions to install the requirements for MAVISp

3. create a local copy of the repository:

```
git clone https://github.com/ELELAB/MAVISp
```

4. run the app in headless mode:

```
cd MAVISp
streamlit run --logger.level=info --server.headless=true
```

Please note down the port that Streamlit is using
(i.e. in this case it is `8501`)

5. **on your workstation** open a new terminal and open an SSH tunnel:

```
ssh -N -L 8080:hostname:8501 user@hostname
```

notice that you need to change the port number at the **right** of `hostname:`
with the one that Streamlit is providing service on 

6. **on your workstation** open a new browser window and visit the website `localhost:8080`. MAVISp should load.

# MAVISp Python package

the MAVISp Python package is not necessary to run the web app - it is however necessary to generate
the database file starting from raw input files. It includes a user-executable script (`mavisp`)
with this very purpose, which performs sanity checks on the input data, prints a report, and
performs the conversion. Please see the help text of the script itself for further details
(`mavisp -h`) or instructions below.

## Requirements

In order to install the package and all its requirements automatically, you will need to have a
working Python 3.7+ installation availalbe. We recommend installing the package in its own
virtual environment - please see previous instructions on how to perform this.

## Installation

Once your virtual environment is ready and active, you need to

1. create a local copy of the repository:

```
git clone https://github.com/ELELAB/MAVISp
```

2. install the Python package

```
cd MAVISp
pip install .
```

the `mavisp` executable will be available.

Notice that installing the Python package always installs all the requirements
for the web app, meaning it is ready to run.

## Running the `mavisp` script

A typical command line of the `mavisp` script looks like:

`mavisp -d input_data -o output_database`

where `input_data` is a folder containing the raw data in a specific format and
`output_database` is were the csv files database will be written. The script

1. performs all the parsing on the input file as well as basic sanity checks
2. prints a summary of all the available datasets and their status
3. prints a more detailed report of the status of each MAVISp module in each dataset
3. writes the output database. The output directories needs not to be present, unless
option `-f` is set, which forces writing or overwriting the database.

The MAVISp modules can have 4 possible states:

  - `OK`, colored in green, for when no warnings or errors were detected
  - `WARNING`, colored in yellow, for when prcessing the module generated some messages that
  can be useful for the user, but was still able to read and elaborate the data for the module correctly
  - `ERROR`, colored in red, for when processing the module resulted in a critical error
  i.e. a problem that couldn't be overcome
  - `NOT_AVAILABLE`, colored in the default terminal color, when the module was not present in the dataset

Additionally, if the mutation list or metadata file are not available, the whole entry is flagged
as being in a `CRITICAL` state - if this is the case, the corresponding modules will not be 
processed and the detected error will be displayed in the report. 

If any module generates at least one error or warning for a given dataset, the status of the corresponding
dataset in the summary table is set accordingly.

If any module generates an error, the script will not write the corresponding database file and exit.

It is possible to process only some of the available proteins in the dataset, or exclude some, by
using options `-p` or `-e` respectively, with a comma-separated list, e.g.

```
mavisp -d input_data -o output_database -p MAP1LC3B,BCL2
```

it is also possible to perform a check of the raw data without writing the database files by using
option `-n`.
