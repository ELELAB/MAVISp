# MAVISp - Events page
# Copyright (C) 2025 Matteo Tiberti, Danish Cancer Society
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import streamlit as st

from streamlit_utils import *
import streamlit as st

st.set_page_config(layout="wide",
    page_title="Events",
    page_icon="🤝")

add_mavisp_logo("static/logo_small.png", image_width='50%')

add_affiliation_logo()

st.title("Events")

st.header("MAVISp BioCurator Training Workshop 2025")

st.subheader("Event details")

st.markdown("""
- When: 3rd to 5th of September 2025
- Where: fully online
- Participation fee: none (free of charge)
- Deadline for registering: ~June 27th 2025~ **Deadline extended!** Register by July 4th 2025
- Requirements:
  - Basic competence with Linux and its terminal
  - Basic knowledge of structural and computational biology
  - Reading: MAVISp preprint and MAVISp GitLab site
- We are able to accommodate up to 15 participants, on a first-come first-served basis
- Registration [here](https://docs.google.com/forms/d/e/1FAIpQLScQmPatyYt43JwA6wyQ2V4Pyh7nLVo0uWa9kAAk3kyZawvSlg/viewform?usp=dialog)

We are organizing the first **MAVISp BioCurator Training Workshop**!

During the event, you will learn the basics behind the MAVISp methodology and work towards becoming a curator for the MAVISp server and website, on your own proteins of interest. The program includes:
- An overview of the methodological framework behind MAVISp
- Hands on sessions, to learn how to use the tools supported by MAVISp to curate your own protein of choice
- Access to a computing server to perform the analyses

For a full programme, see below

The methodology you will learn can be used to generate data for novel proteins to be added to the MAVISp dataset, as well as standalone for your own research.

Read more about MAVISp on [our preprint](https://www.biorxiv.org/content/10.1101/2022.10.22.513328v6) and the [MAVISp gitbook](https://elelab.gitbook.io/mavisp/documentation/how-to-contribute-as-a-biocurator)

Please register by filling out [this form](https://docs.google.com/forms/d/e/1FAIpQLScQmPatyYt43JwA6wyQ2V4Pyh7nLVo0uWa9kAAk3kyZawvSlg/viewform?usp=dialog)
""")

st.subheader("Tentative programme")

st.markdown("""This is the tentative programme for the workshop. Dates and times displayed are expressed in Central European Time (CET, UTC+1), which is the time zone where Copenhagen (Denmark) is located.

#### Day 1 - September 3rd 2025

- **09:30 – 09:45** Welcome and overview
- **09:45 – 10:15** Introducing MAVISp: a community-driven framework for bio-curators and developers
- **10:15 – 10:30** Introduction to the workshop and practicalities
- **10:30 – 10:45** Coffee break
- **10:45 – 11:15** Structure trimming strategies and structure selection for variant analysis
- **11:15 – 11:45** *Practical:* Trimming your protein structure
- **11:45 – 12:45** Lunch break
- **12:45 – 13:15** Introduction to first Snakemake workflow: MutateX
- **13:15 – 13:45** *Practical:* Designing input for the MutateX Snakemake workflow
- **13:45 – 14:15** Long-range module with AlloSigMA2
- **14:15 – 14:45** *Practical:* submission of AlloSigMA2 jobs


#### Day 2 – September 4th 2025

- **09:30 – 09:45** Recap of day 1
- **09:45 – 10:45** Snakemake workflow for MAVISp automatization
- **10:45 – 11:15** *Practical:* Running MAVISp automatization (including coffee break)
- **11:15 – 11:45** AlloSigMA workflow overview
- **11:45 – 12:15** *Practical:* Running AlloSigMA workflow
- **12:15 – 13:15** Lunch break
- **13:15 – 13:45** How to request imports and validate aggregated CSV files
- **13:45 – 14:15** *Practical:* Requesting and running downstream analyses in MAVISp
- **14:15 – 14:45** First look at results obtained so far

#### Day 3 – September 5th 2025

- **09:30 – 09:45** Recap of day 2
- **09:45 – 11:15** Evaluating results and discussion
- **11:15 – 11:45** Introduction to GitBook reporting
- **11:45 – 12:45** 🍽 Lunch break
- **12:45 – 13:45** *Practical:* GitBook reporting (breakout rooms)
- **13:45 – 14:45** Overview of other MAVISp modules and ensemble mode
- **14:45 – 15:15** How to contribute as a curator or developer
- **15:15 – 15:45** Feedback and ideas for future development""")
