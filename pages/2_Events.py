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
    page_title="MAVISp - Events",
    page_icon="🤝")

add_mavisp_logo("static/logo_small.png", image_width='50%')

add_affiliation_logo()

st.title("Events")

st.header("Current events")

st.subheader("MAVISp BioCurator Training Workshop 2026")

st.markdown("""
### Essential information
- When: 7th to 9th of September 2026
- Where: Danish Cancer Institute, Strandboulevarden 49, 2100 Copenhagen, Denmark
- Participation fee: none (free of charge)
- Deadline for registering: June 28th 2026
- We are able to accommodate up to 20 participants
- Requirements:
  - Basic competence with Linux and its terminal
  - Basic knowledge of structural and computational biology
  - Reading: [MAVISp paper](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.70548) and [MAVISp GitBook site](https://elelab.gitbook.io/mavisp/)
- Registration [here](https://docs.google.com/forms/d/e/1FAIpQLSesi3Si0qLHlJuf1xJ3Sy4LoPUFy05uBpAzKK2FC1m3X_MzfA/viewform?usp=publish-editor)

### Full description
            
We are organizing the second MAVISp BioCurator Training Workshop!
            
During the event, you will learn the basics of the MAVISp methodology
and work towards becoming a curator for the MAVISp server and website.
            
The program includes:
            
- An overview of the methodological framework behind MAVISp
- Hands on sessions, to learn how to use the tools supported by MAVISp to curate a protein
- Access to a computing server to perform the analyses
- Exciting talks by invited speakers
            
The methodology you will learn can be used to generate data for novel proteins to be
added to the MAVISp dataset. 
            
Read more about MAVISp on [our publication](https://onlinelibrary.wiley.com/doi/full/10.1002/pro.70548) 
and the [MAVISp GitBook site](https://elelab.gitbook.io/mavisp/documentation/how-to-contribute-as-a-biocurator)

Please register by filling out [this form](https://docs.google.com/forms/d/e/1FAIpQLSesi3Si0qLHlJuf1xJ3Sy4LoPUFy05uBpAzKK2FC1m3X_MzfA/viewform?usp=preview)

A tentative programme for the three days follows

#### Day 1 – September 7th 2026 - DCI, 5.1.A.B
            
- **09:00 – 09:15** Welcome and overview
- **09:15 – 09:45** Introducing MAVISp: a community-driven framework for bio-curators and developers
- **09:45 – 10:00** Introduction to the workshop and practicalities
- **10:00 – 10:45** Structure trimming strategies and structure selection for variant analysis
- **10.45 - 11.00** Coffee  break
- **11:00 – 12:30** Practical: Trimming your protein structure and reporting results in small groups
- **12:30 – 13:30** Lunch break
- **13:30 – 14:15** Introduction to first Snakemake workflows: MutateX and ThermoMPNN
- **14:15 – 15:15** Practical: Designing input for the MutateX Snakemake workflow and ThermoMPNN workflow
- **15:15 - 15:30** Coffee break 
- **15:30 – 16:15** Long-range module in simple mode 
- **16:15 – 17:30** Practical: submission of jobs for long-range module (in small groups)

#### Day 2 – September 8th 2026 - DCI, 5.S.A

- **09:00 – 09:30** Recap of day 1
- **09:30 – 10:15** Snakemake workflow for MAVISp automatization
- **10:15 – 11:15** Practical: Running MAVISp automatization (in small groups)
- **11:15 - 11:30** Coffee break
- **11:30 - 12:00** How to request imports and validate aggregated CSV files 
- **12:00 - 12:30** Practical: Requesting and running downstream analyses in MAVISp (in small groups)
- **12:30 – 13:30** Lunch break
- **13:30 - 14:30** First look at results obtained so far and data analysis through the database (in small groups)
- **14:30 - 15:00** Evaluating results, discussion and Q&A 
- **15:00 - 15:15** Coffee break 
- **15:15 - 16:00** Introduction to GitBook reporting
- **16:00 - 17:00** Practical: GitBook reporting (in small groups)
- **17:00 - 19:00** Networking - Small bites and drinks 

#### Day 3 – September 9th 2026 - DCI, 4.1.A
            
- **09:00 – 09:30** Overview of other MAVISp modules and ensemble mode 
- **09:30 - 10:00** How to contribute as a curator or developer 
- **10:00 - 10:30** Closing training part and Feedback
- **10:30 - 10:45** Coffee break
- **10:45 - 11:00** Introduction of Final Session with Invited Talks
- **11:00 - 11:30** Invited Talk - Possible new tools for the framework - TBA
- **11:30 - 12:00** Invited Talk - Possible new tools for the framework - TBA
- **12:00 - 12:20** Invited Talk - Selected from the developer teams - TBA
- **12:20 - 13:20** Lunch break
- **13:20 - 14:00** Invited Talk - Selected from the developer teams - TBA
- **14:00 - 14:30** Invited Talk - Possible new sources of data and annotation for the framework -  TBA
- **14:30 - 15:00** Invited Talk - Possible new sources of data and annotation for the framework -  TBA
- **15:00 - 16.00** Closing remarks, Networking with coffee and cakes 

            
""")

st.divider()

st.header("Past events")

st.subheader("MAVISp BioCurator Training Workshop 2025")

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
