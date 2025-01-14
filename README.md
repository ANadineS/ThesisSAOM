# ThesisSAOM

Guide to the files in the code directory:
- General data exploration can be found in the data exploration markdown file.
- For both error scenarios, there is a separate main_analysis markdown file. Here, the simulation scenarios are prepared. Subsequently, the simulation studies are executed on the supercomputer. Then, in the same markdown files the output from the supercomputer is loaded, formatted and analzyed per scenario.
- The simulation run and execute simulation scripts are used to run the simulation study on the supercomputer.
- The format output and analyze output scripts are used to format and perform the most important analyses on the output of the supercomputer (they are called from the main analyses files).
- The helper functions file contains several useful functions used throughout the study, including the main analyses files and during the execution of the simulation study on the supercomputer.
- The graphs file contains functions to easily create density graphs for all parameters in each scenario.

The data files used to store the supercomputer output, as well as the formatted and analyzed versions of the output, are very large and cannot be uploaded to GitHub. They can be found here: https://drive.google.com/drive/folders/1_PBcyaZmHucm89LVYUTFN3_bnsWtyJZN?usp=sharing.
