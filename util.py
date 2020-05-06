"""This module contains a set of commonly used functions.

Functions:
getNshapes(*)     --- generate n shapes cycling through a short list.
getNcolors(*)     --- generate n distinct colors.
readFileToList(*) --- read a text file into a list.
saveListToFile(*) --- save a list to a text file.
printDict(*)      --- Print out a dictionary.
readFileToDict(*) --- read a delimited text file into a dictionary.
"""

import os, sys
import re
import colorsys
import itertools as it


def getNshapes(n, shape=None):
    """Generate n shapes cycling through a short list.
    
    Parameters:
    n : int      --- number of shapes to return.
    shape : list --- base shapes to cycle from.
                     if None, use a default set. 
    
    Returns:
    A list with the number of shapes.    
    """
    
    shapes=[]
 
    if shape==None or len(shape)==0:
        shape=['o', 'v', '^', "s", '+', '*', 'x', 'D', '|']
    
    itc=it.cycle(shape)
    [shapes.append(next(itc)) for _ in range(n)]

    return shapes


def getNcolors(n):
    """Generate n distinct colors.

    Parameters:
    n : int      --- number of colors to return.
    
    Returns:
    A list with the number of colors.
    """
    
    hsv=[(i/n, 1, 1) for i in range(n)]
    
    colors=[]
    for hsv_t in hsv:
        
        rgb_t=[int(c*255) for c in colorsys.hsv_to_rgb(*hsv_t)]

        print(rgb_t)

        rgb=["{:02x}".format(c) for c in rgb_t]
        
        colors.append("".join(rgb))
    
    return colors


def readFileToList(f, level=1):
    """A function reading a text file into a list.

    Parameters:
    f : str      --- the text file name.
    
    keyword arguments:
    level : int  --- keep empty or comments lines? (default 1).
                     If level=0, keep them.
    
    Returns:
    A list, each item is a line from the input file. 
    """ 

    try:
        fh=open(f, "r")
    
        ls=[]
        for line in fh:
            
            if level==1 and (re.match("^\s*$", line) or re.match("\s*#", line)):
                continue
            
            ls.append(line.strip())
    except IOError as e:
        print("\n*** error reading file {}***".format(f))
        print(e)
        
        sys.exit(1)
    
    fh.close()
        
    return ls


def saveListToFile(outL, f):
    """A function to save a list to a text file.
    
    Parameters:
    f : str      --- the text file name.
    outL : list  --- The list to be saved.
    """ 

    try:
        fh=open(f, "w")
    
        fh.writelines("\n".join(outL))
            
    except IOError as e:
        print("\n*** error saving file {}***".format(f))
        print(e)

        sys.exit(1)
        
    fh.close()
        

def printDict(dict, end="\n"):
    """Print out a dictionary.
    
    Parameters:
    dict : dictionary  --- dictionary to be saved.
    
    keyword arguments:
    end : str          --- end each key:value pair with (default "\\n")
    """
    
    for d in dict:        
        print("{} ===> {}".format(d, dict[d]), end=end)


def readFileToDict(f, keyC=1, valueC=0):
    """A function reading a delimited text file into a dictionary.
    
    Parameters:
    f : str     --- the file name of the delimited text file.
                    The delimiters can be '\\t', ',' or ' '
    
    keyword arguments:
    keyC : int   --- the column (1 based) used as the key of the
                     dictionary (default 1). Must >0.
    valueC : int --- the column (1 based) used as the value of
                     the dictionary (default 0).
                     When valueC is 0, take columns 2 up to the end.
                     When it is negative, take the whole line.
                     
    Returns:
    A dictionary 
    """  

    if keyC <1:        
        print("The key column must >1")
        sys.exit(1)

    try:
        fh=open(f, "r")
    except:
        print("\n*** error reading file {}***".format(f))
        sys.exit(1)
    else:
        ls=fh.readlines()

    fh.close()

    # get the headers
    line=ls.pop(0).strip()
    temp=re.split('[\t,]', line)
    
    nCol=len(temp)
    
    if nCol <keyC:
    
        msg=f"The key column {keyC} is too large."
        print(f"{msg}")
        
        sys.exit(1)
        
    if nCol <valueC:        

        print(f"The value column {valueC} is too large")
        
        sys.exit(1)
        
    lso={}
    for line in ls:
        
        line=line.strip()
        
        temp=re.split('[\t,]', line)
        
        key=temp[keyC-1]
        if valueC <0:
            lso[key]=line
            
        elif valueC==0:
            lso[key]="\t".join(temp[1:])
            
        else:            
            lso[key]=temp[valueC-1]
                        
    return lso
