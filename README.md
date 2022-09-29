# MAVISp

the MAVISp database for structural variants annotation

## Running the app - with X forwarding

This is the simplest way to run MAVISp. It is however pretty slow
and clunky, meaning it is not very effective for in-depth data exploration 
or analysis. For a slightly more complicated but much better performing
way without using X-forwarding see instructions below.

In order to run the app on our server, please follow these instructions:

1) ssh with X forwarding to the server:

```
ssh -XY user@hostname
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

## Running the app - without X forwarding

We can use an ssh tunnel to connect to our streamlit instance on our server.
This is slightly more complicated but leads to an overall much better
user experience.

1) ssh without X forwarding to the server:

```
ssh user@...
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
streamlit run app.py --server.headless=true
```

streamlit will respond with a welcome message similar to:

```
  You can now view your Streamlit app in your browser.

  Network URL: http://0.0.0.0:8501
  External URL: http://0.0.0.0:8501
```

Please note down the port that Streamlit is providing service on
(i.e. in this case it is `8501`)

5) **on your workstation** open a new terminal and open an SSH tunnel:

```
ssh -N -L 8080:hostname:8501 user@hostname
```

notice that you need to change:
- the hostname
- your username
- the port number a the **right** of `hostname:` with the one that Streamlit is providing service on

6) **on your workstation** open a newbrowser window and visit the website `localhost:8080`. MAVISp should load.
