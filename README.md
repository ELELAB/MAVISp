Cancer Systems Biology, Technical University of Denmark, 2800, Lyngby, Denmark  
Cancer Structural Biology, Danish Cancer Institute, 2100, Copenhagen, Denmark  
Repository associated to the publication:  
> MAVISp: A Modular Structure-Based Framework for Genomic Variant Interpretation
> Matteo Arnaudi, Ludovica Beltrame, Kristine Degn, Mattia Utichi, Simone Scrima,
> Pablo Sanchez Izquierdo, Karolina Krzesinska, Francesca Maselli, Terezia Dorcakova,
> Jordan Safer, Alberte Heering Estad, Katrine Meldgard, Philipp Becker, Julie Bruun Brockhoff,
> Amalie Drud Nielsen, Valentina Sora, Alberto Pettenella, Jeremy Vinhas,
> Peter Wad Sackett, Claudia Cava, Anna Rohlin, Mef Nilbert, Sumaiya Iqbal, Matteo Lambrughi,
> Matteo Tiberti, Elena Papaleo. 
> bioRxiv https://doi.org/10.1101/2022.10.22.513328

# MAVISp web app

## Introduction

This is the web app of the MAVISp database for structural variants annotation.

If you use MAVISp, please cite our preprint:

> MAVISp: Multi-layered Assessment of VarIants by Structure for proteins  
> Matteo Arnaudi, Ludovica Beltrame, Kristine Degn, Mattia Utichi, Pablo Sánchez-Izquierdo, Simone Scrima, Francesca Maselli, Karolina Krzesińska, Terézia Dorčaková, Jordan Safer, Katrine Meldgård, Julie Bruun Brockhoff, Amalie Drud Nielsen, Alberto Pettenella, Jérémy Vinhas, Peter Wad Sackett, Claudia Cava, Sumaiya Iqbal,  View ORCID ProfileMatteo Lambrughi, Matteo Tiberti, Elena Papaleo 
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

Running the MAVISp web app requires a working Python 3.9+ installation with the following
Python packages:

- streamlit 1.28.2
- streamlit-aggrid 0.3.4.post3
- pandas 2.1.3
- matplotlib 3.7.4

In principle, it is compatible with all operating systems that support Python.

It has been last test on Linux (Ubuntu 18.04), and on macOS (13.5.2),
with Python 3.9.6 and the following package versions:

- streamlit 1.28.2
- streamlit-aggrid 0.3.4.post3
- pandas 2.1.3
- matplotlib 3.7.4

In order to download the full MAVISp dataset from OSF, you will also need the `wget` 
program (see below). We last tested the download with wget 1.21.3.

## Installing requirements

These instructions apply to both Linux and macOS, using the terminal. You
will need to have a recent (>=3.9) Python distribution installed on your system
or [Anaconda](https://anaconda.org).

### Installing requirements using a virtualenv Python environment

1. if you have access to the `virtualenv` Python environment module, you can create a
virtual environment:

```
virtualenv -p python3.9 MAVISp_env
```

2. then, activate it:

```
source MAVISp_env/bin/activate
```

3. you can install the requirements in the environment using `pip`:

```
pip install pandas==2.1.3 matplotlib==3.7.4 streamlit==1.28.2 streamlit-aggrid==0.3.4.post3
```

### Installing requirements using a conda Python environment

1. if you have access to Anaconda or Miniconda (executable `conda`), you can use it
to create a virtual environment:

```
conda create -n MAVISp_env python
```

2. then you can activate it:

```
conda activate MAVISp_env
```

3. you need to install the remaining requirements, using `pip`:

```
pip install pandas==2.1.3 matplotlib==3.7.4 streamlit==1.28.2 streamlit-aggrid==0.3.4.post3
```

Installation time is typically up to a few minutes.

## Running the app

In the following instructions
  - `hostname` denotes the hostname of the server (i.e. the one you usually ssh to)
  - `user` denotes your username on the server

### Running the app locally - full dataset

This requires downloading the full MAVISp dataset.

In order to run our web server locally with its full content, you will need to
download the full MAVISp dataset from OSF, as follow, as well as download our
web app from GitHub. If you'd rather test the web app on a small subset, please
follow the instructions in the "Running the app locally - test dataset" instead.

1. If you haven't already, activate your Python environment (see previous steps)

2. create a local copy of the MAVISp repository in your system:

```
git clone https://github.com/ELELAB/MAVISp
```

3. Download the database files from [our OSF repository](https://osf.io/ufpzm/).
This requires the `wget` program or similar. If it's not available, you can
manually download a zip file containing all the database files from [this link](https://files.de-1.osf.io/v1/resources/ufpzm/providers/osfstorage/65579865874c2e15e54e7d34/?zip=).

```
cd MAVISp
rm -rf ./database
mkdir database 
cd database
wget -O database.zip 'https://files.de-1.osf.io/v1/resources/ufpzm/providers/osfstorage/65579865874c2e15e54e7d34/?zip=' && unzip database.zip 
rm database.zip
cd ..
```

At the end of the process, you should have a `database` folder inside
the `MAVISp` folder including all the contents of the `database` folder on OSF.

4. With your Python environment still active and from inside the `MAVISp` repository
directory, run:

```
streamlit run Welcome.py
```

a browser window displaying the MAVISp web app should open.

### Running the app locally - test dataset

These instructions allow to run our web app on a minimal test dataset that
is included in the distribution.

1. If you haven't already, activate your Python environment (see previous steps)

2. create a local copy of the MAVISp repository in your system:

```
git clone https://github.com/ELELAB/MAVISp
```

3. Create a copy of the test dataset in the MAVISp directory

```
cd MAVISp
rm -rf ./database
cp -r test_data/mavisp_web_server database
```

4. With your Python environment still active and from inside the `MAVISp` repository
directory, run:

```
streamlit run Welcome.py
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
to be able to connected with X forwarding to your server, and to have a web
browser installed on your server. For a faster alternative that doesn't use
X-forwarding see instructions below.

1. connect to your server via ssh, with X forwarding:

```
ssh -XY user@hostname
```

2. Follow the previous instructions to install requirements, download required files
 and run the app locally. A browser window should pop up.

### Running the app remotely - without X forwarding

You can use an ssh tunnel to connect to your streamlit instance.
This is slightly more complicated but leads to an overall much better
user experience.

1. ssh without X forwarding to the server:

```
ssh user@hostname
```

2. Follow the previous instructions to install the requirements and download the
required files for MAVISp

3. run the app in headless mode:

```
cd MAVISp
streamlit run --logger.level=info --server.headless=true Welcome.py
```

Please note down the port that Streamlit is using
(i.e. in this case it is `8501`)

4. **on your workstation** open a new terminal and open an SSH tunnel:

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

The MAVISp Python package is designed to run on any operating system that supports
Python. It has been tested on Ubuntu Linux 18.04 and macOS (13.5.2) 

In order to install the package and all its requirements automatically, you will need to have a
working Python 3.9+ installation available. We recommend installing the package in its own
virtual environment - please see previous instructions on how to create a virtual environment.

The MAVISp Python package requires the following packages, and has been tested
with the following versions:

- pandas 2.1.3
- tabulate 0.9.0
- matplotlib 3.7.4
- numpy 1.26.2
- PyYAML 6.0.1 
- streamlit 1.28.2
- streamlit-aggrid 0.3.4.post3
- requests 2.31.0
- termcolor 2.3.0

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

Installation time is typically up to a few minutes.

Notice that installing the Python package always installs all the requirements
for the web app, meaning it is ready to run.

## Running the `mavisp` script on the example dataset

A typical command line of the `mavisp` script looks like:

`mavisp -d input_data -o output_database`

where `input_data` is a folder containing the raw data in a specific format and
`output_database` is were the csv files database will be written. The script

1. performs all the parsing on the input file as well as basic sanity checks
2. prints a summary of all the available datasets and their status
3. prints a more detailed report of the status of each MAVISp module in each dataset
3. writes the output database. The output directories needs not to be present, unless
option `-f` is set, which forces writing or overwriting the database.

To test the script on the MAVISp dataset included in the repository, from inside
the `MAVISp` folder you have cloned from GitHub (see previous instructions),
just run:

```
mavisp -d test_data/mavisp_python_package -o test_output
```

the `test_output` folder will now contain the output of the command. A reference
output for this command can be found in the `test_data/mavisp_web_server/` folder.
The execution should complete in a few seconds.

## Additional instructions

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
