This file describes in details how to use the program TmPy (Tm.py)

Overview
The program calculates various number of thermodynamics for a set of DNA duplexes using the nearest neighbor parameters (SantaLucia J Jr., (1998), PNAS, 95, 1460-65). It also corrects the monovalent and divalent salt effects (Ahsen et. al, (2001), Clinical Chemistry, 47(11), 1956-61). However, it does not correct for anything else, such as DMSO and glycerol. The DNA duplexes can have mismatches but cannot have consecutive mismatches.

Input
There are three ways DNA duplexes can be entered:
1. Through the program command line interface using -s1 and -s2 options.
     In this case, the duplex is named "YourSeq". The first strand (-s1) is required and must be in 5'->3' orientation. The second strand (-s2) can be in either 5'->3' or 3'->5' orientations. If it is in 5'->3' orientation, use -r or --revS2 to reverse the strand. If the second strand is not specified, the duplex is assumed to be a perfect match duplex. Note that -s2 cannot stand alone.
2. An input file.
     The file is a delimited file by comma or tab. It has a header and three columns. The first column is the name of the duplex and the second column is the first strand. Similar to the command line input, the second strand in the third column is optional for any given duplex. If it is not given, the duplex is assumed to be a perfect match duplex. If the second strand is given in 5'->3' orientation, use -r or --revS2 to reverse it.

     One example input file "example_input_file.txt" was given.

3. Combination of 1 and 2 above.

Output
The output has three verbose levels (use -v or -vv to increase the level) in terms of the amount of results obtained. The output is displayed to the screen and saved in a file. The output file name can be given by -o or --outfile, or inferred from the input. If inferred, and there is an input file, the output file is named as "TmPy.xxx.out.txt", where 'xxx' is the stem of the input file. If no input file is given, the output file is named as "TmPy.YourSeq.out.txt". The output is sorted by the name of the duplexes.  

An example output file "example_output_file.txt" corresponding to the example input file mentioned above is given. It was obtained using all the default parameters. One can compare his/her output to this file.

1. Default level
    Only the Tm under the given condition is returned.
2. Level 2 (specify by -v)
    Tm, perB, dG, dH, dS, and conditions (temperature, concentrations of primer, template, monovalent and divalent cations) are returned, where
           
           Tm     ---    temperature (C) at which half of the less strand is in duplex
           PerB  ---    percentage of the less strand in duplex
           dG     ---     free energy (kcal/mol)
           dH     ---     enthalpy (kcal/mol)
           dS      ---     entropy (e. u.)
3. Level 3 (specify by -vv)
    In addition to all the results at level 2, the following are given. One can use them to check against literature values:

           Tm_std   ---  Tm at standard condition (4 mM total strand concentration for non 
                                 self-complementary duplex and the two strand have the same concentrations.
                                1 mM for self-complementary duplex)
           dG_std, dS_std  --- for standard thermodynamics conditions (1 M for all strands and 1 M Na)

Installation
For a quick run, one can drop all the files in a working directory. For a long term, it is recommended to create a folder for the nearest neighbor parameter files 'nnSH.csv' and 'nnSS.csv', mark the files as read only and create an environment variable "NNDIR" for the folder they are in. It is also recommended to create a folder for the library files (thermo.py, error.py, util.py and utilSeq.py) and add its path to the environment variable "PYTHONPATH"

Disclaimer
This program is free to use but use at your own risk. No warranties or liabilities of any kind, explicit or implicit, are assumed.

Copyright
This program is free to use as a whole. Claim as your own, substantially copy into or otherwise used in your own tools, or similar acts, without public acknowledgment to the author is discouraged. Use in commercial activities, without express permission from the author, is prohibited.

Miscellaneous
Contact the author jerry.chen.ctnbio@gmail.com for bugs, improvements, comments, use permissions, or whatever you think is meaningful. 


