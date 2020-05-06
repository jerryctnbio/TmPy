"""A set of custom exception handling classes.

This module imports a custom module 'utilSeq'. They should come
together in one distribution.

Classes:
Error(Exception) --- base of all the custom error class.
ConcentrationZeroError(Error) --- raised when concentration is zero.
ConcentrationOrderError(Error) --- raised when concentration of primer
                                   is lower than template.
TemperatureRangeError(Error) --- raised when temperature is out of range.
                                 It should be [-100, 200]
#FolderError(Error) --- raised when a file folder does not exist.
#FileNotExistError(Error) --- raised when a file does not exist.
NNnotExistError(Error) --- raised when an nearest neighbor is not
                           supported.
NNFileNotFoundError(Error) ---raised when the nearest neighbor parameters
                              files not found.
NotDNAError(Error) --- raised when a duplex contains letters other than
                       A, C, G and T.
"""

import utilSeq


class Error(Exception):
    """Base class for all the custom error class."""
    
    pass


class NotDNAError(Error):
    """Raised when a duplex contains letters other than A, C, G and T.
    
    Attributes:
    duplex  : str --- a duplex
    message : str --- explanation (Not DNA)    
    """ 
    
    def __init__(self, duplex):
        """constructor"""
        
        msg1="The duplex contains letters other than A, C, G and T!"
        self.message=f"\n***{msg1}***\n{duplex}"
        

    def __str__(self):
        """__str__"""
    
        return f"{__class__.__name__}:{self.message}"
        
        
class DuplexNotFlushError(Error):
    """Raised when a duplex is not flush.
    
    Attributes:
    duplex  : str --- a duplex
    message : str --- explanation (duplex is not flush)    
    """ 
    
    def __init__(self, duplex):
        """constructor"""
        
        msg1="The two strands of the duplex should have same length!"
        self.message=f"\n***{msg1}***\n{duplex}"
        

    def __str__(self):
        """__str__"""
    
        return f"{__class__.__name__}:{self.message}"
        
        
class ConcentrationZeroError(Error):
    """Raised when concentration is zero.
    
    Attributes:
    variable : str --- name of the variable
    message  : str --- explanation ('concentration should not be zero')    
    """ 
    
    def __init__(self, variable):
        """constructor
        
        Parameters:
        variable : str --- name of the variable
        """
    
        self.variable=variable
        self.message=f"\n***Concentration [{variable}] should not be zero!***"
        

    def __str__(self):
        """__str__"""
    
        return f"{__class__.__name__}:{self.message}"
        
        
class ConcentrationOrderError(Error):
    """Raised when concentration of primer is lower than template.
    
    Attributes:
    message : str --- explanation (primer concentration should be larger
                      than template)
    """ 
    
    def __init__(self):
        """constructor"""
    
        msg="Primer concentration should be larger than template!"
        self.message=f"\n***{msg}***"
        

    def __str__(self):
        """__str__"""
    
        return f"{__class__.__name__}:{self.message}"
        

class TemperatureRangeError(Error):
    """Raised when temperature is out of range.
    
    Attributes:
    message : str --- explanation (the temperature range is [-100, 200])
    """ 
    
    def __init__(self):
        """constructor"""
    
        msg="Temperature should be between -100 C to 200 C!"
        self.message=f"\n***{msg}***"


    def __repr__(self):
        """__repr__"""
    
        return f"{__class__.__name__}:{self.message}"

        
    def __str__(self):
        """__str__"""
    
        return f"{__class__.__name__}:{self.message}"


class NNnotExistError(Error):
    """ Raised when an nearest neighbor is not supported.
    
    Attributes:
    nn      : str --- nearest neighbor
    top     : str --- top strand in a duplex (5'->3')
    bottom  : str --- bottom strand in a duplex (3'->5')
    message : str --- explanation (nn is not supported)
    """ 

    def __init__(self, nn, top, bottom):
        """constructor"""
    
        msg=utilSeq.matchUp(top, bottom)
    
        msg1=f"The nearest neighbor '{nn}' is not supported!"
        msg2=f"Nearest neighbor with two mismatches not allowed."
        msg3='*********'
        self.message=f"\n{msg}\n\n{msg3}\n{msg1}\n{msg2}\n{msg3}"
        

    def __str__(self):
        """__str__"""
    
        return self.message


class NNFileNotFoundError(Error):
    """Raised when the nearest neighbor parameter files not found.
    
    Attribute:
    message : str --- explanation (the nearest neighbor parameter files not found)
    """ 
    
    def __init__(self):
        """constructor"""

        msg="The nearest neighbor parameter files 'nnSH.csv' and 'nnSS.csv'"
        msg+=" not found!\n"
        msg+="They should be in a folder pointed to by the environment"
        msg+=" variable 'NNDIR',\n"
        msg+="or in the current working directory.\n"
        
        self.message=f"\n***{msg}***"
                

    def __str__(self):
        """__str__"""
    
        return f"{__class__.__name__}:{self.message}"
