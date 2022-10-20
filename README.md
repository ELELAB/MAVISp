# MAVISp

## Introduction

This is the web app of the MAVISp database for structural variants annotation.

If you use MAVISp, please cite our preprint:

PLACEHOLDER_REF_BIOARXIV

Please see the CHANGELOG.md file on this repository for current and expected
data releases as well as updates on the software.

Please see the CURATORS.md file for an up to date list of our curators.

## Requirements

Running the MAVISp web app requires a working Python 3.7+ installation with the following
Python packages installed:

- streamlit
- streamlit-aggrid
- pandas
- numpy
- pyyaml

It has been tested on Linux and macOS, but in principle runs on any operating system
that supports Python.

The easiest way to run MAVISp is by means of a Docker container - see the next section.

If you would rather install and run it on your local hardware, please follow the instructions
below. They apply to both Linux and macOS, using the terminal.

## Running the app as a Docker container

coming soon - please use one of the options below for now

## Installing requirements

### Installing requirements using a virtualenv Python environment

1. if you have access to the `virtualenv` Python environment module, you can create a
virtual environment:

```
virtualenv -p python3.8 MAVISp_env
```

2. then you can activate it:

```
source MAVISp_env/bin/activate
```

3. you can install the requirements in the environment using `pip`:

```
pip install pandas streamlit streamlit-aggrid
```

### Installing requirements using a conda Python environment

1. if you have access to Anaconda or Miniconda (executable `conda`), you can create a
virtual environment:

```
conda create -n MAVISp_env python==3.8 streamlit pandas
```

2. then you can activate it:

```
conda activate MAVISp_env
```

3. you can install one of the requirements, which is currently not on conda, using `pip`:

```
pip install streamlit-aggrid
```

## Running the app

In the following instructions
  - `hostname` denotes the hostname of the server (i.e. the one you usually ssh to)
  - `user` denotes your username on the server

### Running the app locally

Once you have installed your prerequisites, or created and activate your
virtual environment as described above, you will need a local copy of the
application and data:

```
git clone https://github.com/ELELAB/MAVISp
```

if the environment is active (see step 2. of the previous sections)
you can then run the app directly:

```
cd MAVISp
streamlit run app.py
```

a browser window redirecting you the web app should open.

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

3. Run the app locally in headless mode:

streamlit run --logger.level=info --server.headless=true app.py

streamlit will respond with a welcome message similar to:

```
  You can now view your Streamlit app in your browser.

  Network URL: http://0.0.0.0:8501
  External URL: http://0.0.0.0:8501
```

Please note down the port that Streamlit is using
(i.e. in this case it is `8501`)

4. **on your workstation** open a new terminal and open an SSH tunnel:

```
ssh -N -L 8080:hostname:8501 user@hostname
```

notice that you need to change the port number at the **right** of `hostname:`
with the one that Streamlit is providing service on 

5. **on your workstation** open a new browser window and visit the website `localhost:8080`. MAVISp should load.
