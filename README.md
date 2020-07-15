# clumpanalyser
Repo containing scripts to analyse fragments found in phantom dump files (https://github.com/danieljprice/phantom).
Data directory should contain all phantom dump files and output files from clumpfind routine .

To run the analysis:
python runanalysis.py

This will run the tracker.py routine which iterates through clumpfind files and traces individual clumps through their evolution history.

analyse.py will then analyse at plot clump properties and different stages of their evolution
