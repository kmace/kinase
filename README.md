# kinase-analysis
Code base to analyze data generated from my thesis project

Directory structure is as follows:

/input:
	Directory to hold incoming resources
/input/raw:
	Raw data, directly from experiment
/input/meta:
	metadata associated with each experiment
/input/...:
	output from preprocessing steps

/src:
	Source code
/src/preprocessing:
	Code to transfer raw data, into forms for analysis
/src/utils:
	Code for doing general purpose transformations on data, including loading, ploting, filtering
/src/sci:
	Code for asking and answering specific scientific questions about the data. generally one off scripts to generate a single plot etc.
	these might be .R files, or Rmd files if there is specific formating desired. all output (such as pdf, or html) should go to output

/output:
	Directory for output from scripts, should follow the form that the filename follows the script that created it. should have 1:1 mapping
	to /src/sci
