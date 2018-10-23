# hla-mm
This scripts read in output from HLAmatchmaker (Hmm) to quantify the number of mismatched epitopes in the graft-versus-host (GVH) and host-versus-graft (HVG) directions respectively. The dose of the mismatched epitopes are also taken into consideration. 
## Input files
* Class1 analysis result from Hmm, tab 4 and tab 5 exported as class1rec.csv and class1don.csv respectively in the same directory as the python script.
* Class2 analysis result from Hmm, tab 4 and tab 5 exported as class2rec.csv and class2don.csv respectively in the same directory as the python script.
## Running the script
```
python class1gvh_hvg_v2.py
python class2gvh_hvg_v2.py
```
## Output and further processing
"class1gvh_hvg" and "class2gvh_hvg" are generated by the scripts. The files can be imported into excel. Pivot using the "type of ep-mm". 
