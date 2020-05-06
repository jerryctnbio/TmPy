"""A set of common functions related to sequence analysis.

Functions:
seqRev(s)  --- reverse a DNA/RNA sequence.
seqComp(s) --- complement a DNA/RNA sequence.
seqRC(s)   --- reverse complement a DNA/RNA sequence.
isWC(b1, b2) --- check if the two bases form a canonical watson-crick pair.
matchUp(top, bottom) --- match the top strand to the bottom.
"""

import re


def isWC(b1, b2):
    """Check if the two bases form a canonical watson-crick pair.
    
    Parameters:
    b1 : str  --- a DNA/RNA base
    b2 : str  --- a DNA/RNA base
    
    Returns:
    A boolean. True when they form a pair; False when not.
    """
    
    pair=''.join(sorted((b1+b2).upper()))
    
    return True if pair=='CG' or pair=="AT" else False


def matchUp(top, bottom):
    """Match the top strand to the bottom.

    Generate a string for the match, where
    a match is denoted by '|' and mismatch by 'x'
    
    Parameters:
    top    : str --- top strand in 5'->3' orientation
    bottom : str --- bottom strand in 3'->5' orientation

    Returns:
    Three lines of text as shown in the following example.
    
    5' CGCAGT  3'
       |||x||
    3' GCGACA  5'
    """

    match=['|' if isWC(top[i], bottom[i]) else 'x'
                                    for i in range(len(top))]

    # turn into strings
    top=''.join(top)
    match=''.join(match)
    bottom=''.join(bottom)

    return f"5'{top}3'\n  {match}  \n3'{bottom}5'"
    

def seqRev(s):
    """Reverse a DNA/RNA sequence.
    
    Parameters:
    s : str     --- input sequence.
    
    Note:
    The sequences are case sensitive.
    
    Returns:
    A string reversed.
    """
    
    s_t="".join(reversed(s))
    
    return s_t


def seqComp(s):
    """Complement a DNA/RNA sequence.
    
    Parameters:
    s : str    --- input sequence.
    
    Note:
    The sequences are case sensitive.
    
    Returns:
    A string complementary to the input.
    """
    
    be="ACGTUacgtu"
    af="TGCAAtgcaa"
    
    s_t=s.translate(str.maketrans(be, af))
    
    return s_t


def seqRC(s):
    """Reverse complement a DNA/RNA sequence.
    
    Parameters:
    s : str    --- input sequence.
    
    Note:
    The sequences are case sensitive.
    
    Returns:
    A string reverse complementary to the input.
    """
    
    return seqComp(seqRev(s))
