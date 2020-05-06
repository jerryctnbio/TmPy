"""This is a module for calculating nucleic acids thermodynamics
using nearest neighbor parameters.

The module contains only one class 'Thermo'. The class needs the
following custom modules: error, util and utilSeq. They should
come together with this module.

For more information about class 'Thermo', see the class for
details
""" 

import sys, os
from pathlib import Path
import re
import math

import error
import util, utilSeq


class Thermo(object):
    """A class to calculate thermodynamics using nearest neighbor
    parameters.
    
    Constants and default values
    _R_      --- The gas constant
    _cp_     --- default concentration of primer (300 nM)
    _ct_     --- default concentration of template (1.38e-15 nM)
    _na_     --- default concentration of monovalent salt (100 mM)
    _mg_     --- default concentration of divalent salt (0.0 mM)
    _temper_ --- default temperature (65 C)

    Attributes:
    _nndH    --- nearest neighbor parameters for enthalpy.
    _nndS    --- nearest neighbor parameters for entropy.
    
    Methods:
    thermoCal(*)  ---   calculates Tm, perB, dG, dH, dS, Tms, dGs
                        and dSs for all duplexes given in a dictionary.
                        perB is the percentage bound.
                        Tms, dGs and dSs are values under standard
                        conditions.                     
    thermoCal0(*) ---   calculates thermodynamics for one duplex.
    getMelting(*) ---   calculates the percentage bound for a duplex
                        at various temperatures given in a list.
   
    _get_dHdS(*)  ---   Calculates dH and dS for one duplex.
    _readNN(HorS) ---   read the nearest neighbor parameters for dH or dS.
    _perBcal(k, cp, ct)  --- calculate the percentage bound (perB).
    """

    _R_=1.987
    _cp_=300.0
    _ct_=1.38e-15
    _na_=100.0
    _mg_=0.0
    _temper_=65.0
    
    _nndH={}
    _nndS={}


    @staticmethod
    def _readNN(HorS):
        """Read the nearest neighbor parameters from a file.
        
        The file path is specified by an environment variable 'NNDIR'
        or at the current working directory.
        
        The parameters are made up of 4 6x6 blocks, ignoring blank and
        comment lines.
        
        The first row in each block is the top nearest neighbors
        in 3'->5'. The first columns in the following 5 rows are the
        bottom nearest neighbors in 5'->3' orientation. Together they
        form one complete nearest neighbor 5'-top-3'/3'-bottom-5'. 

        Parameters:
        HorS : str   --- flag if to read for dH or dS.
        
        Exceptions:
        NNFileNotFoundError  --- a custom error class, raised when
                                 neither of the parameter files
                                 "nnSH.csv" and "nnSS.csv" not found.
                                 This could happen when the environment
                                 variable 'NNDIR' not set or not pointing
                                 to the folder where the parameter files
                                 are or the parameter files are not at 
                                 the current working folder.
        
        Return:
        A dictionary. Its keys and values are the nearest neighbors and
                      the corresponding parameters, respectively 
        """
 
        if 'NNDIR' in os.environ:
            folder=os.environ['NNDIR']
        else:
            folder=os.getcwd()
            
        stem="nnSH.csv" if HorS=='dH' else "nnSS.csv"        
        name=Path(folder+'/'+stem)

        try:
            if not Path(name).exists():
                raise error.NNFileNotFoundError()
        except error.NNFileNotFoundError as e:
            print(e)
            sys.exit(1)
              
        ls=util.readFileToList(name)
        
        index=0
        NN={}
        while index < len(ls):
            
            # read six (6) rows
            # first row
            bottom=re.split("[\t,]", ls[index])
            
            # the following five (5) rows
            for i in range(5):
                
                index+=1
                top=re.split("[\t,]", ls[index])
                
                for j in range(1, 6):
                    
                    value=float(top[j]) if len(top[j]) >0 else 0.0
                    
                    NN[top[0]+"/"+bottom[j]]=value
                    
                    # get reverse complementary
                    nn_t=''.join(reversed(bottom[j]))+"/"+''.join(reversed(top[0]))
                    NN[nn_t]=value
                
            index+=1
        
        return NN


    @staticmethod
    def _perBcal(k, cp, ct):
        """Calculate the percentage bound.
        
        Paramaters:
        k  : float   --- equilibrium constant
        cp : float   --- primer concentration
        ct : float   --- template concentration
        
        Return:
        float  ---  percentage bound          
        """

        c=cp/ct
        b=-(c+1+1/k/ct)

        # since b is always negative, i.e., b<0
        x=math.sqrt(b*b-4.0*c)
        p1=(-b+x)/2.0
        p2=c/p1
        
        perB=p1 if p1>=0.0 and p1<=1.0 else p2
        
        return perB


    def _get_dHdS(self, pair):
        """Calculates dH and dS for one duplex.

        Parameters:
        pair : str  --- a duplex in "top/bottom" format. 
                        The top is in 5'->3' orientation.
                        The bottom is in 3'->5' orientation.
                        It should have the same length.

        Exceptions:
        NNnotExistError --- custom error class, raised when an nearest
                            neighbor parameter does not exist.
                            This happens when there are more
                            than one mismatches in one nearest neighbor.

        Return:
        a list with dH and dS
        """
             
        nndH=self._nndH
        nndS=self._nndS

        pair=pair.upper()

        try:
            if re.search("[^ACGT/]", pair):
                raise error.NotDNAError(pair)
        except error.NotDNAError as e:
            print(e)
            
            sys.exit(1)

        temp=pair.split("/")
        
        isSymm=1 if temp[0]==utilSeq.seqRC(temp[0]) else 0
        
        top=list(temp[0])
        bottom=list(temp[1])

        try:
            if len(top) !=len(bottom):        
                raise error.DuplexNotFlushError(pair)
        except error.DuplexNotFlushError as e:            
            print(e)

            sys.exit(1)

        # initiation  
        dH=0.2
        dS=-5.7
        
        # propagation
        for b in range(len(top)-1):
            
            nn=top[b]+top[b+1]+'/'+bottom[b]+bottom[b+1]
        
            try:    
                dH+=nndH[nn] 
                dS+=nndS[nn] 
            except KeyError:
                try:
                    raise error.NNnotExistError(nn, top, bottom)
                except error.NNnotExistError as e:
                    print(e)

                    sys.exit(1)
                    
        # symmetry correction
        dS+=-1.4 if isSymm==1 else 0.0

        # terminal AT correction
        if (top[0]+bottom[0] != 'GC') and (top[0]+bottom[0] != 'CG'):
            dH+=2.2
            dS+=6.9

        if (top[-1]+bottom[-1] != 'GC') and (top[-1]+bottom[-1] != 'CG'):
            dH+=2.2
            dS+=6.9

        return [dH, dS]
        

    def getMelting(self, pair, temp):
        """Calculates the percentage bound for a duplex at various 
        temperatures given in a list.
        
        Parameters:
        pair : str  --- a duplex in "top/bottom" format. 
                       The top is in 5'->3' orientation.
                       The bottom is in 3'->5' orientation.
                       It should have the same length.
        temp : list  --- temperatures in a list
        
        Return:
        A dictionary. Its key and value are temperature in celsius (c)
                      and percentage bound, respectively.
        """
        
        cp=self._cp
        ct=self._ct
        R=self._R_

        get_dHdS=self._get_dHdS
        perBcal=self._perBcal
        
        [dH, dS]=get_dHdS(pair)
        
        melt={}
        
        kelvin=[t+273.15 for t in temp]
        for k in kelvin:
            
            dG=dH-k*dS/1000
            keq=math.exp(-dG*1000/R/k)
        
            perB=perBcal(keq, cp, ct)
            
            melt[k-273.15]=perB
            
        return melt


    def thermoCal0(self, pair, verbose=0):
        """Thermodynamics calculation for one duplex.
        
        It calculates Tm, perB, dG, dH, dS, Tms, dGs and dSs for
        the given duplex.
        
        Parameters:
        pair : str  --- a duplex in "top/bottom" format. 
                        The top is in 5'->3' orientation.
                        The bottom is in 3'->5' orientation.
                        It should have the same length.
        verbose : int --- verbose level (default:0)
                          verbose=0: return only Tm.
                          verbose=1: verbose 0 + perB, dG, dH and dS.  
                          verbose=2: verbose 1 + Tms, dGs and dSs.  

        Return:
        *** everything is returned as string ***
        
        Tm   ---  melting temperature under the given condition.

        perB ---  percentage bound.
        dG, dH, dS --- thermodynamics calculated under the given condition         

        Tms --- Tm under conventional conditions.
                meaning total oligo concentraton is 0.1 mM for self
                complimentary duplexes; and 0.4 mM for non-self 
                complimentary duplexes and the two strands have the
                same concentrations. 
        dGs, dSs --- thermodynamics under standard condition.
        """
     
        R=self._R_

        temper=self._temper
        cp=self._cp
        ct=self._ct
        na=self._na
        mg=self._mg
        
        perBcal=self._perBcal
        get_dHdS=self._get_dHdS        
        
        [dH, dS]=get_dHdS(pair)

        # effective monovalent concentration
        # Ahsen et. al, (2001), Clinical Chemistry, 47(11), 1956-61
        na_eff=na+0.12*math.sqrt(mg*1000)

        # sale correction for dS        
        # SantaLucia J Jr., (1998), PNAS, 95, 1460-65 
        n=(len(pair)-1)/2-1
        dSeff=dS+0.368*n*math.log(na_eff)

        dG=dH-temper*dSeff/1000
        
        # at melting temperature, the equilibrium constant is
        ktm=1/(cp-ct/2)
        Tm=dH*1000/(dSeff-R*math.log(ktm))-273.15

        Tm_str="{:7.2f}".format(Tm)

        delimiter='\t'
        out=pair+delimiter+Tm_str

        if verbose>0:
            
            k=math.exp(-dG*1000/R/temper)
        
            perB=perBcal(k, cp, ct)
        
            perB_str="{:6.3e}%".format(perB*100.0)
            dG_str="{:8.3f}".format(dG)
            dH_str="{:8.3f}".format(dH)
            dSeff_str="{:8.3f}".format(dSeff)
            temper_str="{:6.2f}".format(temper-273.15)

            verb_str=delimiter.join([perB_str, dG_str, dH_str, dSeff_str
                                                             , temper_str])
        
            cp_str=f"{cp*1e9:7.4e}"
            ct_str=f"{ct*1e9:7.4e}"
            na_str=f"{na*1e3:7.4e}"
            mg_str=f"{mg*1e3:7.4e}"
            
            unit_str=delimiter.join([cp_str, ct_str, na_str, mg_str])           
        
        if verbose ==1:
            out+=delimiter+verb_str
            out+=delimiter+unit_str
            
        if verbose >=2:
            
            dGs=dH-310.15*dS/1000

            # under standard condition, the equilibirum constant is
            ks=1e4
            Tms=dH*1000/(dS-R*math.log(ks))-273.15

            Tms_str="{:7.2f}".format(Tms)
            dGs_str="{:8.3f}".format(dGs)
            dSs_str="{:8.3f}".format(dS)
    
            out+=delimiter+verb_str
            out+=delimiter+delimiter.join([Tms_str, dGs_str, dSs_str])
            out+=delimiter+unit_str

        return out
        
        
    def thermoCal(self, oligo, verbose=0):
        """Thermodynamics calculation for duplexes in a dictionary.
   
        Parameters:
        oligo : dictionary  ---  duplexes
        verbose : int       --- verbose level (default:0)
                                affects the amount of output.
                                see method 'thermoCal0' for details.
                                 
        Return:
        A list with thermodynamics calculated. The list is sorted by
        the keys of the duplexes.
        
        See method 'thermoCal0' for details.
        """

        delimiter="\t"
    
        # get an ordered list of duplexes
        oligo_order=sorted(oligo.keys())
        
        myThermo=[]
        for o in oligo_order:
            
            out_st=o+delimiter+self.thermoCal0(oligo[o], verbose)
            myThermo.append(out_st) 

        return myThermo
        

    def __init__(self, temper=_temper_, cp=_cp_, ct=_ct_, na=_na_, mg=_mg_):
        """Constructor.
        
        1. Check the conditions.
        2. Set the thermodynamics conditions including using defaults:
           Turn the concentration units into M, and the temperature
           from celcius to kelvin.
           read the nearest neighbor parameters.

        Parameters:
        temper : float  --- temperature (default: _temper_)
        cp     : float  --- primer concentration (default: _cp_)
        ct     : float  --- template concentration (default: _ct_)
        na     : float  --- monovalent salt concentration (default: _na_)
        mg     : float  --- divalent salt concentration (default: _mg_)        
        """

        if cp==0:
            raise error.ConcentrationZeroError('cp')

        if ct==0:
            raise error.ConcentrationZeroError('ct')
        
        if cp<ct:
            raise error.ConcentrationOrderError
        
        if temper>200 or temper<-100:
            raise error.TemperatureRangeError

        self._cp=float(cp)*1e-9
        self._ct=float(ct)*1e-9
        self._temper=float(temper)+273.15
        self._na=float(na)*1e-3
        self._mg=float(mg)*1e-3

        self._nndH=self._readNN("dH")
        self._nndS=self._readNN("dS")

        
    def __repr__(self):
        """A string representation of the class."""
        
        return "class:{}".format(__class__.__name__)
    