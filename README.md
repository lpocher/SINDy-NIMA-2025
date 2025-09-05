# SINDy NIMA 2025
This GitHub folder represents much of the work done for the NIMA 2025 paper "Data-Driven Discovery of Centroid Dynamics". There are a variety of files uploaded herein. 

The various "umerlat_main....py" files are python files which are input decks to the Warp code (https://warp.lbl.gov) with source code found here https://bitbucket.org/berkeleylab/warp/src/master/. 
The Warp code uses a variety of support files: "work.gs", "work2.gs"
The script "umerlat_main.py" is the main script. It calls various other auxialiary files found in the "tune_util" and "currents" folder.
In addition, other *.py files are used by "umerlat_main.py": "util.py", "Umagnets_pkl.py", "beam.py", and "mphoto.py".
The various text files are Warp simulation data used in the post-processing and SINDy script written as Jupyter Notebook "Learning the Position.ipynb". 
The text file "package_SINDy_conda.txt" lists the conda environments packages used in for the aforementioned Jupyter Notebook.

The files uploaded here are uploaded as a convienence for fellow researchers. 
The data and processing files are uploaded as "working notebooks" and "as is".
We make not statements regarding how they will work on other machines.
