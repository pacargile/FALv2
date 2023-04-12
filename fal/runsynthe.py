import os,sys,glob,filecmp
from datetime import datetime
import subprocess
# import cStringIO
import shutil
from scipy import io as spIO
import numpy as np
import fal

class Synthe(object):
    """ 
    General class for running SYNTHE from python
    
    """
    def __init__(self, *args, **kwargs):
        super(Synthe, self).__init__()
        self.args = args
        self.kwargs = kwargs

        self.verbose = self.kwargs.get('verbose',True)

        # determine path to fortran exe files
        self.exedir = self.kwargs.get('exedir',fal.__abspath__+'/bin/')

        # define some general files
        self.molecules = self.kwargs.get('molecules','./data/molecules.dat')
        self.continuua = self.kwargs.get('continuua','./data/continuua.dat')
        self.he1tables = self.kwargs.get('he1tables','./data/he1tables.dat')

        # define atm file
        self.atmmod = self.kwargs.get('atmmod','./data/atmmod_sol.dat')
        
        # define spectrv input file
        self.spectrv_infile = self.kwargs.get('spectrv_infile','./data/spectrv.input')

        # string for rotate
        self.rotatevar = ("{NROT:5d}{NRADIUS:5d}\n{VROT:10.1f}\n")


    def xnfpelsyn(self,verbose_xnf=None):
        """
        Run XNFPELSYN code
        
        Reads In: 
            fort.2  (ascii)[molecules]
            fort.17 (bin)[continua]
            fort.18 (ascii)[he1lines]

        Writes Into:
            fort.10 (bin)
            
        """

        # write links to molecules, continua, and he1tables
        self._makesym(self.molecules,'fort.2')
        self._makesym(self.continuua,'fort.17')
        self._makesym(self.he1tables,'fort.18')

        if self.verbose:
            starttime_xnf = datetime.now()
            print("Running xnfpelsyn... [{0}]".format(starttime_xnf))
        self.xnfpelsynout = self._callpro("xnfpelsyn",
                                          inpipe=self.atmmod,
                                          verbose=verbose_xnf)
        if self.verbose:
            endtime_xnf = datetime.now()
            print("... Finished xnfpelsyn [{0}: {1}]".format(endtime_xnf,endtime_xnf-starttime_xnf))
        
        return
    
    def synthe(self,verbose_syn=True):
        """
        Run SYNTHE code

        Reads In: 
           fort.10
           fort.12
           fort.14
           fort.18
           fort.19
           fort.20
           fort.93
           
        Writes Into:
           fort.7
           fort.8
           fort.9
           fort.13
           fort.14
           fort.15
           fort.28
           fort.29
           
        Delete At Completion:
           fort.7
           fort.12
           fort.13
           fort.14
           fort.15
           fort.20
           fort.93

        """
        
        # write in information into input string
        self.synout = self._callpro("synthe",verbose=verbose_syn)
        # print("Finished Synthe")

        if self.verbose:
            starttime_syn = datetime.now()
            print("Running synthe... [{0}]".format(starttime_syn))
        self.synout = self._callpro("synthe",verbose=verbose_syn)
        if self.verbose:
            endtime_syn = datetime.now()
            print("... Finished synthe [{0}: {1}]".format(endtime_syn,endtime_syn-starttime_syn))
                
        return
    
    def spectrv(self,tau=False,verbose_sprv=True):
        """
        Run SPECTRV code

        Reads In: 
            fort.5 [model atm]
            fort.25 [spectrv.input]
            fort.9
            fort.10

        Writes Into:

        
        Delete At Completion
            fort.20
            fort.9
        """
        # make mod atm link to fort.5 and spectrv.input as fort.25
        self._makesym(self.atmomod,'fort.5')
        self._makesym(self.spectrv_infile,'fort.25')
        
        # write in information into input string
        if tau:
            if self.verbose:
                starttime_sprv = datetime.now()
                print("Running spectrv_tau... [{0}]".format(starttime_sprv))
            self.spectrvout = self._callpro("spectrv_tau",verbose=verbose_sprv)
            if self.verbose:
                endtime_sprv = datetime.now()
                print("... Finished spectrv_tau [{0}: {1}]".format(endtime_sprv,endtime_sprv-starttime_sprv))
        else:
            if self.verbose:
                starttime_sprv = datetime.now()
                print("Running spectrv... [{0}]".format(starttime_sprv))
            self.spectrvout = self._callpro("spectrv",verbose=verbose_sprv)
            if self.verbose:
                endtime_sprv = datetime.now()
                print("... Finished spectrv [{0}: {1}]".format(endtime_sprv,endtime_sprv-starttime_sprv))

        return
    
    def rotate(self,vrot=0.0,verbose_rot=True):
        """
        Run ROTATE code

        Reads In: 
             fort.1
        Writes Into:
             fort.19
        New Out: 
             ROTX (bin) X = # rotation velocities
        Delete At Completion:
        """
        # link fort.sol to fort.1
        self._makesym('fort.7','fort.1')            

        # write in information into input string
        if self.verbose:
            starttime_rot = datetime.now()
            print("Running rotate... [{0}]".format(starttime_rot))
        rotatestr = self.rotatevar.format(NROT=1,NRADIUS=0,VROT=vrot)
        self.rotateout = self._callpro("rotate",rotatestr,verbose=verbose_rot)
        if self.verbose:
            endtime_rot = datetime.now()
            print("... Finished rotate [{0}: {1}]".format(endtime_rot,endtime_rot-starttime_rot))

    def broaden(self,):
        pass

    def _callpro(self,function,inputstr=None,inpipe=None,verbose=None):
        """
        general function to call fortran code
        """
        # set stdout piping
        if verbose == None:
            fnull = self.FNULL
        elif verbose == True:
            fnull = sys.stdout
        elif verbose == False:
            fnull = open(os.devnull, 'w')
        else:
            fnull = self.FNULL

        # build the process
        if inpipe != None:
            pro = subprocess.Popen([self.exedir+function+".exe"],
                                   stdin=open(inpipe,'r'),stdout=fnull,encoding='ascii',
                                   universal_newlines=True)
        else:
            pro = subprocess.Popen([self.exedir+function+".exe"],
                                   stdin=subprocess.PIPE,stdout=fnull,encoding='ascii',
                                   universal_newlines=True)

        # if input to pro, then communicte
        if inputstr is not None:
            output = pro.communicate(input=inputstr)
        else:
            output = None

        # wait until pro is finished
        pro.wait()
        return output    
    
    def _fastcopy(self,src,dst,buffer_size=-1):#10485760):
        """
        Function that does a fast copy using buffers and binary fmt
        """
        # Optimize buffer for small files
        # buffer_size = min(buffer_size,os.path.getsize(src))
        # if(buffer_size == 0):
        #     buffer_size = 1024
        try:
            with open(src,'rb') as fsrc:
                with open(dst,'wb') as fdst:
                    shutil.copyfileobj(fsrc,fdst,buffer_size)
        except TypeError:
            shutil.copy(src,dst)

    def _makesym(self,src,outname):
        """
        Function that takes a file, copies to memory, and sets up a symlink
        """
        os.symlink(src,outname)
