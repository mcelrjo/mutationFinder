# Scripts for finding target-site mutations in herbicide-related genes.

abiMutationFinder_v.0.0.1.py is written to batch analyze abi files for mutations in ALS, ACCase, PPO, EPSPS, and psbA. Takes only a path to a directory of AB1 files.  

herbicideMutationFinder is written to look for target site mutations in the ALS, ACCase, psbA, alpha-tubulin, PPO, and EPSP synthase DNA/RNA coding sequences. I am writing different versions of the script to allow for different ways to enter the sequences (using command line arguments or entering sequences via copy and paste).  Check back for updates.  Please contact me if you are using and you find this useful -- just to help keep me motivated to work on this project.


### There are currently two versions.  See below for there different use.

#### abiMutationFinder 

##### abiMutationFinder_v.0.0.1.py is written to take a path to a diretory of .ab1 files which will be analyzed for possible mutations. Separate files will be written for each ab1 file. 

#### herbicideMutationFinder - Various Versions

##### herbicideMutationFinder_v.0.#.#.py is meant to run in command line and read in ABI files (.ab1).  Run in the command line as follows:

##### python herbicideMutationFinder_v.0.#.#.py YOURABIFILEHERE.ab1

The script is not written to accept flags at this time.  Give it a try and let me know what you think.  A batch version will be coming soon that can analyze all .ab1 files in a directory and output to a text file.

#### herbicideMutationFinder_v.1.#.#.py is meant to run in command line and take RAW INPUT.  Do not pass command line argument to this version.  Instead wait for the prompts and copy and paste your sequence into the command line (Remember CTRL-Shift-V to paste) and select your target site when prompted.

Run Version 1 as follows:

##### python herbicideMutationFinder_v.1.0.0.py 

This version will walk you through a series of questions to allow you to input your sequence (need to copy and paste the sequence) and select the target site of interest.  The raw input questions/prompts are as follows:

`Copy and paste your DNA search string here:`

Followed by

`For target site enter 'a' for alpha-tubulin, 
 'b' for acetolactate synthase, 
 'c' for acetylcoa carboxylase, 
 'd' for protoporphyrinogen oxidase, 
 'e' for epsp synthase, 
 'f' for psbA PSII D1 protein. 
 Enter a selection:`

 After a selection is made the output will immediately follow.


