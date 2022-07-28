# CutProba
CutProba-1.1.0
The tool for predicting probabilities of proteolytic events in proteins. It based on ML algorithms and data from CutDB database.
This is the second version of the tool. 

Running steps
1. Download directories "fitted_models", "scripts", "PWM" and the script "run.py"
2. Run the script "run.py". As input, you should provide:
  a) the protein ID from PDB database and the name of separate chain from this protein;
  b) the name of desired PWM;
  c) the name of desired protease;
  d) the name of desired ML algorithm (see models in directory "fitted_models".
  For help, use '-h' flag.
3. Output is the directory with name according to the protein ID + chain, which contains structural features file, file with structural and sequence scores calculated and plots with visualised scores.  
