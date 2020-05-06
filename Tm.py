"""This tool calculates the nucleic acids thermodynamics using the 
nearest neighbor parameters. It calculates the melting temperatures
of DNA duplexes, and other properties based on user options.

The input can be a file ('-f'), or entered through the command line
argument ('-s1') and optionally '-s2'. One can use both an input
file and command line arguments as input. If neither is given a 
default duplex will be used as a demo.

The input file is a delimited file by comma or tab. It has a header
and three columns. The first column is the name of the duplex and 
the second column is the first strand. Similar to the command line
input, the second strand in the third column is optional for any
given duplex.

See the accompanying file 'README.txt' for more.

This script imports the following custom modules: 'util', 'utilSeq'
, 'thermo', 'error'. All the custom modules should come in one
distribution with this script.
 
"""

import sys, os
import argparse
from pathlib import Path 

import util, utilSeq
import thermo, error 

__version='R1.0.0.0'

def main():
    """DNA oligo duplex Tm calculator."""
    
    argParser=argparse.ArgumentParser(description=__doc__)
    msg="a sequence file in delimited ('[\t|,| ]') format"
    argParser.add_argument('-f', '--file', help=msg)
    argParser.add_argument('-o', '--outfile', help="output file")
    msg="the first strand in 5'->3' direction"
    argParser.add_argument('-s1', help=msg)
    msg="The second strand in either direction. See option -r for more"
    argParser.add_argument('-s2', help=msg)
    msg="A flag to reverse s2 if it is in 5'->3' direction"
    argParser.add_argument("-r", "--revS2", help=msg, action="store_true")
    argParser.add_argument('-v', help="verbosity level. e.g., -v, -vv."
                                             , action="count", default=0)
    argParser.add_argument('-t', '--temperature', type=float
                               , help="temperature in celsius degree")
    argParser.add_argument('-cp', help="primer concentration in nM (300)"
                                , type=float, default=300)
    msg="template concentration in nM (1.38e-15)"
    argParser.add_argument('-ct', help=msg
                                , type=float, default=1.38e-15)
    msg="monovalent salt concentration in mM (100)"
    argParser.add_argument("-n", "--na", help=msg, type=float, default=100)
    msg="divalent salt concentration in mM (0.0)"
    argParser.add_argument('-m', '--mg', help=msg, type=float, default=0.0)    
    argParser.add_argument('-V', '--version', action='version'
                                            , version=__version)
        
    args=argParser.parse_args()
    
    oligo={}

    if((not args.s1) and args.s2):
        print("\n*** Needs to specify s1 if s2 is specified.")
        print("\n*** See usage below\n\n")
        argParser.print_help()

        sys.exit(1)

    if((not args.file) and (not args.s1)):
        print("\n*** Used the following default DNA duplex as a demo.\n")
        args.s1="CGATCG"
        args.s2=None

    if args.file:
        print("Reading input file {} ...".format(args.file))
        
        oligo=util.readFileToDict(args.file)
        
    if args.s1 and args.s2:
        oligo['YourSeq']=args.s1+'\t'+args.s2
    elif args.s1:        
        oligo['YourSeq']=args.s1
        
    oligo=getOligoPair(oligo, args.revS2)

    # set up the condition
    cond={}
    
    cv=[args.temperature, args.cp, args.ct, args.na, args.mg]
    ck=['temper', 'cp', 'ct', 'na', 'mg']

    def _setkv(k, v):
        cond[k]=v
        
    [_setkv(ck[i], cv[i]) for i in range(5) if cv[i] != None]

    try:
        myThermo=thermo.Thermo(**cond)
    except (error.TemperatureRangeError, error.ConcentrationZeroError, error.ConcentrationOrderError) as e:
        print(str(e)+"\n")
        print("***Please see the usage below***\n\n")
        argParser.print_help()
        sys.exit(1)
    
    out=myThermo.thermoCal(oligo, args.v)
    
    # add header
    header1=["Name", "Duplex", "Tm(C)"]
    header2=["PerBound (%)", "dG(kcal/mol)", "dH(kcal/mol)", "dS(e.u.)"
                                            , "Temperature(C)"]
    header3=["Tm_std(C)", "dG_std(kcal/mol)", "dS_std(e.u.)"]
    headerCond=["Cprimer(nM)", "Ctemplate(nM)", "C_mono(mM)"
                                              , "C_divalent(mM)"]

    delimiter="\t"
    if args.v==0:
        out.insert(0, delimiter.join(header1))
        
    elif args.v==1:
        out.insert(0, delimiter.join(header1+header2+headerCond))
            
    elif args.v >=2:
        out.insert(0, delimiter.join(header1+header2+header3+headerCond))

    print("\n".join(out))

    # save the output
    if args.outfile:
        outFile=args.outfile
    else:
    
        outFile="TmPy."
    
        if args.file:
            outFile+=Path(args.file).stem
        else:
            outFile+="YourSeq"
        
        outFile+=".out.txt"
    
        outFile=Path(os.getcwd()+'/'+outFile)
                 
    util.saveListToFile(out, outFile)


def getOligoPair(oligo, r=True):
    """Explicitly match up the duplexes in an antiparallel fashion.
    
    All duplexes in the end should comprise of two explicit strands. 
    The two strands are concatenated as a string, separated by '/'.
    The first one in 5'-3' orientation and the second 3'-5'.

    If the second strand is not given, its reverse complement is created.     
        
    Parameters:
    oligo : dictionary --- dictionary holding the duplex info.  
    r     : int        --- if the second strand needs to be reversed.
                           (default 1)
    
    Returns:
    A dictionary. The keys and values are the duplex names and the explicit
                  duplexes, respectively.
    """
    
    oligoO={}
    for item in oligo:
        
        temp=oligo[item].split("\t")
        
        if len(temp)==1:

            value_t=temp[0]            
            oligoO[item]=value_t+'/'+utilSeq.seqComp(value_t)
            
        elif len(temp)==2:
            
            value_t=utilSeq.seqRev(temp[1]) if r==True else temp[1]
            
            oligoO[item]=temp[0]+'/'+value_t 


    return oligoO

    
if __name__ == '__main__':    
    """The program entry point"""
    
    main()
