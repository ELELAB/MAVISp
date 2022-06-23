# MAVISp

the MAVISp database for structural variants annotation

## Running the app

In order to run the app on our server, please follow these instructions:

1) ssh with X forwarding to the server:

```
ssh -XY user@...
```

2) clone the repository somewhere on the server

```
cd /data/user/$USER
git clone https://github.com/ELELAB/MAVISp.git
cd MAVISp
```

3) activate our Python virtual environment

```
. /usr/local/envs/py37/bin/activate
```

4)

run the app:

```
streamlit run app.py
```

a Firefox window should pop up.
