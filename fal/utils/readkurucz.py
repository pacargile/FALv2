from ctypes import cdll, CDLL, POINTER, ARRAY, c_int, c_double, c_char_p, c_char, c_void_p, c_float, c_long
import numpy as np
import os
import fal

c_double_p = POINTER(c_double)
c_float_p = POINTER(c_float)
c_int_p = POINTER(c_int)
c_long_p = POINTER(c_long)
c_char_array = ARRAY(c_char,11)
c_char_array_p = POINTER(c_char_array)

class ReadKurucz(object):
    def __init__(self, *args, **kwargs):
        super(ReadKurucz, self).__init__()

        solib_path = kwargs.get('solib',fal.__abspath__+'/lib/readfort.so')

        try:
            self.rfort = cdll.LoadLibrary(solib_path)
        except IOError:
            print(f'ISSUE WITH READING IN SHARED OBJECT: {solib_path}')
            raise

        self.verbose = kwargs.get('verbose',False)

    def readfiles(self,
                  f12path='./fort.12',f14path='./fort.14',
                  f19path='./fort.19',f20path='./fort.20',
                  f93path='./fort.93'):

        
        # read fort.93 to grab needed line information
        if self.verbose:
            print('... reading in fort.93')
        self.f93in = self.readfort93(f93path)
        self.nlines12 = self.f93in['nlines']
        self.nlines19 = self.f93in['n19']

        # read fort.12
        if self.verbose:
            print('... reading in fort.12')
        self.f12in = self.readfort12(f12path)

        # read fort.14
        if self.verbose:
            print('... reading in fort.14')
        self.f14in = self.readfort14(f14path)

        # read fort.19
        if self.verbose:
            print('... reading in fort.19')
        self.f19in = self.readfort19(f19path)

        # read fort.20
        if self.verbose:
            print('... reading in fort.20')
        self.f20in = self.readfort20(f20path)
        

    def writefiles(self,
                   f12outpath='./fort_NEW.12',f14outpath='./fort_NEW.14',
                   f19outpath='./fort_NEW.19',f20outpath='./fort_NEW.20',
                   f93outpath='./fort_NEW.93'):

        # check if output file exisits, if so delete
        for ff in [f93outpath,f12outpath,f14outpath,f19outpath,f20outpath]:
            if os.path.isfile(ff):
                os.remove(ff)

        # write new fort.93
        if self.verbose:
            print('... writing fort.93')
        self.writefort93(f93outpath)
    
        # write new fort.12 and fort.14
        if self.verbose:
            print('... writing fort.12')
        self.writefort12(f12outpath)

        if self.verbose:
            print('... writing fort.14')
        self.writefort14(f14outpath)

        # write new fort.19 fort.20
        if self.verbose:
            print('... writing fort.19')
        self.writefort19(f19outpath)

        if self.verbose:
            print('... writing fort.20')
        self.writefort20(f20outpath)
        
    
    def readfort93(self,f93path):

        # turn path string into bytes
        s_i  = '{0}'.format(f93path)
        sb_i    = bytes(s_i,encoding='ascii')

        # initialize variables for ctypes to fill 
        NLINESi  = np.zeros_like(0.0,dtype=np.int32)
        LENGTHi  = np.zeros_like(0.0,dtype=np.int32)
        IFVACi   = np.zeros_like(0.0,dtype=np.int32)
        IFNLTEi  = np.zeros_like(0.0,dtype=np.int32)
        N19i     = np.zeros_like(0.0,dtype=np.int32)
        TURBVi   = np.zeros_like(0.0,dtype=np.float64)
        DECKJi   = np.zeros((7,99),dtype=np.float64)
        IFPREDi  = np.zeros_like(0.0,dtype=np.int32)
        WLBEGi   = np.zeros_like(0.0,dtype=np.float64)
        WLENDi   = np.zeros_like(0.0,dtype=np.float64)
        RESOLUi  = np.zeros_like(0.0,dtype=np.float64)
        RATIOi   = np.zeros_like(0.0,dtype=np.float64)
        RATIOLGi = np.zeros_like(0.0,dtype=np.float64)
        CUTOFFi  = np.zeros_like(0.0,dtype=np.float64)
        LINOUTi  = np.zeros_like(0.0,dtype=np.int32)

        self.rfort.readfile93(
            c_char_p(sb_i), 
            NLINESi.ctypes.data_as(c_int_p),
            LENGTHi.ctypes.data_as(c_int_p),
            IFVACi.ctypes.data_as(c_int_p),
            IFNLTEi.ctypes.data_as(c_int_p),
            N19i.ctypes.data_as(c_int_p),
            TURBVi.ctypes.data_as(c_double_p),
            DECKJi.ctypes.data_as(c_double_p),
            IFPREDi.ctypes.data_as(c_int_p),
            WLBEGi.ctypes.data_as(c_double_p),
            WLENDi.ctypes.data_as(c_double_p),
            RESOLUi.ctypes.data_as(c_double_p),
            RATIOi.ctypes.data_as(c_double_p),
            RATIOLGi.ctypes.data_as(c_double_p),
            CUTOFFi.ctypes.data_as(c_double_p),
            LINOUTi.ctypes.data_as(c_int_p),
        )
    
        output = ({
            'nlines':NLINESi, # Total number of lines
            'length':LENGTHi, # Length
            'ifvac':IFVACi, # If wl in vaccuum
            'ifnlte':IFNLTEi, # If include NLTE lines
            'n19':N19i, # Number of NLTE lines
            'turbv':TURBVi, # Vmic added to what is in the atm
            'deckj':DECKJi, # Additional information input in durng synbeg
            'ifpred':IFPREDi, # If predicted lines have been included
            'wlbeg':WLBEGi, # Begining wavelength
            'wlend':WLENDi, # End wavelength
            'resolution':RESOLUi, # resolution for opacity binning
            'ratio':RATIOi, # 1 + 1 / resolution 
            'ratiolg':RATIOLGi, # log(ratio)
            'cutoff':CUTOFFi, # Fraction of continuum opacity below which the wings of lines are computed
            'linout':LINOUTi, # Limit on lines printing out to STDOUT
        })
        
        return output

    def readfort12(self,f12path):

        # turn path string into bytes
        s_i  = '{0}'.format(f12path)
        sb_i    = bytes(s_i,encoding='ascii')

        # initialize variables for ctypes to fill 
        NBUFFi   = np.zeros(self.nlines12,dtype=np.int32)  
        CGFi   = np.zeros(self.nlines12,dtype=np.float32)
        NELIONi  = np.zeros(self.nlines12,dtype=np.int32)    
        ELOi     = np.zeros(self.nlines12,dtype=np.float32)
        GAMMARi  = np.zeros(self.nlines12,dtype=np.float32)   
        GAMMASi  = np.zeros(self.nlines12,dtype=np.float32)   
        GAMMAWi  = np.zeros(self.nlines12,dtype=np.float32)   

        self.rfort.readfile12(
            c_char_p(sb_i), 
            c_int(self.nlines12),
            NBUFFi.ctypes.data_as(c_int_p),    
            CGFi.ctypes.data_as(c_float_p),
            NELIONi.ctypes.data_as(c_int_p),    
            ELOi.ctypes.data_as(c_float_p),
            GAMMARi.ctypes.data_as(c_float_p),
            GAMMASi.ctypes.data_as(c_float_p),
            GAMMAWi.ctypes.data_as(c_float_p),
        )

        output = ({
            'nbuff':NBUFFi, # Number of lines in buffer
            'cgf':CGFi, # normalized GF in frequency
            'nelion':NELIONi, # Index for specific element
            'elo':ELOi, # Lower Energy level
            'gammar':GAMMARi, # Gamma R
            'gammas':GAMMASi, # Gamma S
            'gammaw':GAMMAWi, # Gamma W
        })

        return output

    def readfort14(self,f14path):

        # turn path string into bytes
        s_i  = '{0}'.format(f14path)
        sb_i    = bytes(s_i,encoding='ascii')

        # initialize variables for ctypes to fill 
        WLi       = np.zeros(self.nlines12,dtype=np.float64)
        Ei        = np.zeros(self.nlines12,dtype=np.float64)
        EPi       = np.zeros(self.nlines12,dtype=np.float64)

        LABELi    = np.zeros((self.nlines12,10),dtype='str')
        LABELPi   = np.zeros((self.nlines12,10),dtype='str')
        LABELx    = np.zeros((2,self.nlines12),dtype=np.float64)
        LABELPx   = np.zeros((2,self.nlines12),dtype=np.float64)

        # OTHER1i   = np.zeros((self.nlines12,10),dtype='str')
        # OTHER2i   = np.zeros((self.nlines12,10),dtype='str')

        ISHIFTi   = np.zeros(self.nlines12,dtype=np.int32)
        ISHIFTPi  = np.zeros(self.nlines12,dtype=np.int32)
        IXFIXFPi  = np.zeros(self.nlines12,dtype='str')
        LINESIZEi = np.zeros(self.nlines12,dtype=np.int32)
        AUTOi     = np.zeros(self.nlines12,dtype='str')

        OTHER1x   = np.zeros((2,self.nlines12),dtype=np.float64)
        OTHER2x   = np.zeros((2,self.nlines12),dtype=np.float64)

        WLVACi    = np.zeros(self.nlines12,dtype=np.float64)
        CENTERi   = np.zeros(self.nlines12,dtype=np.float64)
        CONCENi   = np.zeros(self.nlines12,dtype=np.float64)
        NELIONi   = np.zeros(self.nlines12,dtype=np.int32)
        GAMMARi   = np.zeros(self.nlines12,dtype=np.float32)
        GAMMASi   = np.zeros(self.nlines12,dtype=np.float32)
        GAMMAWi   = np.zeros(self.nlines12,dtype=np.float32)
        REFi      = np.zeros((self.nlines12,5),dtype='str')
        REFx      = np.zeros(self.nlines12,dtype=np.float32)
        NBLOi     = np.zeros(self.nlines12,dtype=np.int32)
        NBUPi     = np.zeros(self.nlines12,dtype=np.int32)
        ISO1i     = np.zeros(self.nlines12,dtype=np.int32)
        X1i       = np.zeros(self.nlines12,dtype=np.float32)
        ISO2i     = np.zeros(self.nlines12,dtype=np.int32)
        X2i       = np.zeros(self.nlines12,dtype=np.float32)
        GFLOGi    = np.zeros(self.nlines12,dtype=np.float32)
        XJi       = np.zeros(self.nlines12,dtype=np.float32)
        XJPi      = np.zeros(self.nlines12,dtype=np.float32)
        CODEi     = np.zeros(self.nlines12,dtype=np.float32)
        ELOi      = np.zeros(self.nlines12,dtype=np.float32)
        GFi       = np.zeros(self.nlines12,dtype=np.float32)
        GSi       = np.zeros(self.nlines12,dtype=np.float32)
        GRi       = np.zeros(self.nlines12,dtype=np.float32)
        GWi       = np.zeros(self.nlines12,dtype=np.float32)
        DWLi      = np.zeros(self.nlines12,dtype=np.float32)
        DGFLOGi   = np.zeros(self.nlines12,dtype=np.float32)
        DGAMMARi  = np.zeros(self.nlines12,dtype=np.float32)
        DGAMMASi  = np.zeros(self.nlines12,dtype=np.float32)
        DGAMMAWi  = np.zeros(self.nlines12,dtype=np.float32)
        DWLISOi   = np.zeros(self.nlines12,dtype=np.float32)
        ISOSHIFTi = np.zeros(self.nlines12,dtype=np.int32)
        EXTRA3i   = np.zeros(self.nlines12,dtype=np.float32)

        self.rfort.readfile14_20(
            c_char_p(sb_i), 
            c_int(self.nlines12),
            WLi.ctypes.data_as(c_double_p),    
            Ei.ctypes.data_as(c_double_p),    
            EPi.ctypes.data_as(c_double_p),    
            LABELi.ctypes.data_as(c_char_p),
            LABELx.ctypes.data_as(c_double_p),  
            LABELPi.ctypes.data_as(c_char_p),
            LABELPx.ctypes.data_as(c_double_p),  

            # OTHER1i.ctypes.data_as(c_int_p),
            # OTHER2i.ctypes.data_as(c_char_p),

            ISHIFTi.ctypes.data_as(c_int_p),
            ISHIFTPi.ctypes.data_as(c_int_p),
            IXFIXFPi.ctypes.data_as(c_char_p),
            LINESIZEi.ctypes.data_as(c_int_p),
            AUTOi.ctypes.data_as(c_char_p),

            OTHER1x.ctypes.data_as(c_double_p),    
            OTHER2x.ctypes.data_as(c_double_p),    
            WLVACi.ctypes.data_as(c_double_p),    
            CENTERi.ctypes.data_as(c_double_p),    
            CONCENi.ctypes.data_as(c_double_p),    
            NELIONi.ctypes.data_as(c_int_p),    
            GAMMARi.ctypes.data_as(c_float_p),    
            GAMMASi.ctypes.data_as(c_float_p),    
            GAMMAWi.ctypes.data_as(c_float_p),    
            REFi.ctypes.data_as(c_char_p),     
            REFx.ctypes.data_as(c_float_p),    
            NBLOi.ctypes.data_as(c_int_p),    
            NBUPi.ctypes.data_as(c_int_p),    
            ISO1i.ctypes.data_as(c_int_p),    
            X1i.ctypes.data_as(c_float_p),    
            ISO2i.ctypes.data_as(c_int_p),    
            X2i.ctypes.data_as(c_float_p),    
            GFLOGi.ctypes.data_as(c_float_p),    
            XJi.ctypes.data_as(c_float_p),    
            XJPi.ctypes.data_as(c_float_p),    
            CODEi.ctypes.data_as(c_float_p),    
            ELOi.ctypes.data_as(c_float_p),    
            GFi.ctypes.data_as(c_float_p),    
            GSi.ctypes.data_as(c_float_p),    
            GRi.ctypes.data_as(c_float_p),    
            GWi.ctypes.data_as(c_float_p),    
            DWLi.ctypes.data_as(c_float_p),    
            DGFLOGi.ctypes.data_as(c_float_p),    
            DGAMMARi.ctypes.data_as(c_float_p),    
            DGAMMASi.ctypes.data_as(c_float_p),    
            DGAMMAWi.ctypes.data_as(c_float_p),    
            DWLISOi.ctypes.data_as(c_float_p),    
            ISOSHIFTi.ctypes.data_as(c_int_p),    
            EXTRA3i.ctypes.data_as(c_float_p),
            )

        # # convert ctype char into python string
        # x = np.array([''.join(LABELi[i,:].tobytes('F').decode('ascii')) for i in range(self.nlines12)])
        # x = ''.join(x)
        # LABELi = np.array(list(map(''.join, zip(*[iter(x)]*10))))

        # x = np.array([''.join(LABELPi[i,:].tobytes('F').decode('ascii')) for i in range(self.nlines12)])
        # x = ''.join(x)
        # LABELPi = np.array(list(map(''.join, zip(*[iter(x)]*10))))

        # x = np.array([''.join(REFi[i,:].tobytes('F').decode('ascii')) for i in range(self.nlines12)])
        # x = ''.join(x)
        # REFi = np.array(list(map(''.join, zip(*[iter(x)]*5))))

        # IXFIXFPi = np.array([IXFIXFPi[ii].tobytes('F').decode('ascii') for ii in range(self.nlines12)])
        # AUTOi = np.array([AUTOi[ii].tobytes('F').decode('ascii') for ii in range(self.nlines12)])

        # for ii in range(10):
        #     print(ISHIFTi[ii], ISHIFTPi[ii],IXFIXFPi[ii],LINESIZEi[ii],AUTOi[ii])

        # x = np.array([''.join(OTHER1i[i,:].tobytes('F').decode('ascii')) for i in range(self.nlines12)])
        # x = ''.join(x)
        # OTHER1i = np.array(list(map(''.join, zip(*[iter(x)]*10))))

        # x = np.array([''.join(OTHER2i[i,:].tobytes('F').decode('ascii')) for i in range(self.nlines12)])
        # x = ''.join(x)
        # OTHER2i = np.array(list(map(''.join, zip(*[iter(x)]*10))))

        # for ii in range(10):
        #     print(OTHER2i[ii])

        output = ({
            'wl':WLi, # input wavelength
            'e':Ei, # first energy level
            'ep':EPi, # second energy level
            'label':LABELi, # first level label
            'labelx':LABELx, # unformatted first level label
            'labelp':LABELPi, # second level label
            'labelpx':LABELPx, # unformatted second level label
            # 'other1':OTHER1i, # Other1 info
            # 'other2':OTHER2i, # Other2 info
            'ishift':ISHIFTi,
            'ishiftp':ISHIFTPi,
            'ixfixfp':IXFIXFPi,
            'linesize':LINESIZEi,
            'auto':AUTOi,
            'other1x':OTHER1x, # Other1 info unformatted
            'other2x':OTHER2x, # Other2 info unformatted
            'wlvac':WLVACi, # wavelength in vacuum
            'center':CENTERi, # Intensity at center nu
            'concen':CONCENi, # Continuum at central nu
            'nelion':NELIONi, # Index for specific element
            'gammar':GAMMARi, # unlogged and shifted Gamma R
            'gammas':GAMMASi, # unlogged and shifted Gamma S
            'gammaw':GAMMAWi, # unlogged and shifted Gamma W
            'ref':REFi, # Reference for line
            'refx':REFx, # Unformatted reference for line
            'nblo':NBLOi, # Departure coeff for lower level
            'nbup':NBUPi, # Departure coeff for upper level
            'iso1':ISO1i, # Isotopic number for level 1
            'x1':X1i, # Log of fractional isotopic abundance for level 1
            'iso2':ISO2i, # Isotopic number for level 2
            'x2':X2i, # Log of fractional isotopic abundance for level 2
            'gflog':GFLOGi, # log(gf)
            'xj':XJi, # angular momentum of first level
            'xjp':XJPi, # angular momentum of second level
            'code':CODEi, # numeric code for line species
            'elo':ELOi, # lower energy level
            'gf':GFi, # log(gf) -> unlogged with shifts applied
            'gs':GSi, # log damping constant 
            'gr':GRi, # log damping constant 
            'gw':GWi, # log damping constant 
            'dwl':DWLi, # shift in wl
            'dgflog':DGFLOGi, # shift in log(gf)
            'dgammar':DGAMMARi, # shift in dampling constants 
            'dgammas':DGAMMASi, # shift in dampling constants
            'dgammaw':DGAMMAWi, # shift in dampling constants
            'dwliso':DWLISOi, # shift in isotopic wl shift
            'isoshift':ISOSHIFTi, # isotopic shift
            'extra3':EXTRA3i, # extra information
        })
        return output

    def readfort19(self,f19path):

        # turn path string into bytes
        s_i  = '{0}'.format(f19path)
        sb_i    = bytes(s_i,encoding='ascii')

        WLVACi   = np.zeros(self.nlines19,dtype=np.float64)
        ELOi     = np.zeros(self.nlines19,dtype=np.float32)
        GFi      = np.zeros(self.nlines19,dtype=np.float32)
        NBLOi    = np.zeros(self.nlines19,dtype=np.int32) 
        NBUPi    = np.zeros(self.nlines19,dtype=np.int32) 
        NELIONi  = np.zeros(self.nlines19,dtype=np.int32)    
        TYPEi    = np.zeros(self.nlines19,dtype=np.int32) 
        NCONi    = np.zeros(self.nlines19,dtype=np.int32) 
        NELIONXi = np.zeros(self.nlines19,dtype=np.int32)    
        GAMMARi  = np.zeros(self.nlines19,dtype=np.float32)   
        GAMMASi  = np.zeros(self.nlines19,dtype=np.float32)   
        GAMMAWi  = np.zeros(self.nlines19,dtype=np.float32)   
        NBUFFi   = np.zeros(self.nlines19,dtype=np.int32)  
        LIMi     = np.zeros(self.nlines19,dtype=np.int32)

        self.rfort.readfile19(
            c_char_p(sb_i), 
            c_int(self.nlines19),
            WLVACi.ctypes.data_as(  c_double_p),
            ELOi.ctypes.data_as(    c_float_p),
            GFi.ctypes.data_as(     c_float_p),
            NBLOi.ctypes.data_as(   c_int_p),    
            NBUPi.ctypes.data_as(   c_int_p),    
            NELIONi.ctypes.data_as( c_int_p),    
            TYPEi.ctypes.data_as(   c_int_p),    
            NCONi.ctypes.data_as(   c_int_p),    
            NELIONXi.ctypes.data_as(c_int_p),    
            GAMMARi.ctypes.data_as( c_float_p),
            GAMMASi.ctypes.data_as( c_float_p),
            GAMMAWi.ctypes.data_as( c_float_p),
            NBUFFi.ctypes.data_as(  c_int_p),    
            LIMi.ctypes.data_as(    c_int_p),        
        )

        output = ({
            'wlvac':WLVACi, # vacuum wl
            'elo':ELOi, # energy of lower level
            'gf':GFi, # log(gf) -> unlogged with shifts applied
            'nblo':NBLOi, # departure coeff for lower level
            'nbup':NBUPi, # departure coeff for upper level
            'nelion':NELIONi, # index for specific element
            'type':TYPEi, # type of line
            'ncon':NCONi, # number of continuum layers
            'nelionx':NELIONXi, # index for specific element
            'gammar':GAMMARi, # Gamma R unlogged with shifts applied
            'gammas':GAMMASi, # Gamma S unlogged with shifts applied
            'gammaw':GAMMAWi, # Gamma W unlogged with shifts applied
            'nbuff':NBUFFi, # index of line in buffer
            'lim':LIMi, # Limit for line size?
        })

        return output

    def readfort20(self,f20path):

        # turn path string into bytes
        s_i  = '{0}'.format(f20path)
        sb_i    = bytes(s_i,encoding='ascii')

        # initialize variables for ctypes to fill 
        WLi       = np.zeros(self.nlines19,dtype=np.float64)
        Ei        = np.zeros(self.nlines19,dtype=np.float64)
        EPi       = np.zeros(self.nlines19,dtype=np.float64)

        LABELi    = np.zeros((self.nlines19,11),dtype='str')
        LABELPi   = np.zeros((self.nlines19,11),dtype='str')

        LABELx    = np.zeros((2,self.nlines19),dtype=np.float64)
        LABELPx   = np.zeros((2,self.nlines19),dtype=np.float64)

        # OTHER1i   = np.zeros((self.nlines12,10),dtype='str')
        # OTHER2i   = np.zeros((self.nlines12,10),dtype='str')

        ISHIFTi   = np.zeros(self.nlines12,dtype=np.int32)
        ISHIFTPi  = np.zeros(self.nlines12,dtype=np.int32)
        IXFIXFPi  = np.zeros(self.nlines12,dtype='str')
        LINESIZEi = np.zeros(self.nlines12,dtype=np.int32)
        AUTOi     = np.zeros(self.nlines12,dtype='str')

        OTHER1x   = np.zeros((2,self.nlines19),dtype=np.float64)
        OTHER2x   = np.zeros((2,self.nlines19),dtype=np.float64)

        WLVACi    = np.zeros(self.nlines19,dtype=np.float64)
        CENTERi   = np.zeros(self.nlines19,dtype=np.float64)
        CONCENi   = np.zeros(self.nlines19,dtype=np.float64)
        NELIONi   = np.zeros(self.nlines19,dtype=np.int32)
        GAMMARi   = np.zeros(self.nlines19,dtype=np.float32)
        GAMMASi   = np.zeros(self.nlines19,dtype=np.float32)
        GAMMAWi   = np.zeros(self.nlines19,dtype=np.float32)
        REFi      = np.zeros((self.nlines19,5),dtype='str')
        REFx      = np.zeros(self.nlines19,dtype=np.float32)
        NBLOi     = np.zeros(self.nlines19,dtype=np.int32)
        NBUPi     = np.zeros(self.nlines19,dtype=np.int32)
        ISO1i     = np.zeros(self.nlines19,dtype=np.int32)
        X1i       = np.zeros(self.nlines19,dtype=np.float32)
        ISO2i     = np.zeros(self.nlines19,dtype=np.int32)
        X2i       = np.zeros(self.nlines19,dtype=np.float32)
        GFLOGi    = np.zeros(self.nlines19,dtype=np.float32)
        XJi       = np.zeros(self.nlines19,dtype=np.float32)
        XJPi      = np.zeros(self.nlines19,dtype=np.float32)
        CODEi     = np.zeros(self.nlines19,dtype=np.float32)
        ELOi      = np.zeros(self.nlines19,dtype=np.float32)
        GFi       = np.zeros(self.nlines19,dtype=np.float32)
        GSi       = np.zeros(self.nlines19,dtype=np.float32)
        GRi       = np.zeros(self.nlines19,dtype=np.float32)
        GWi       = np.zeros(self.nlines19,dtype=np.float32)
        DWLi      = np.zeros(self.nlines19,dtype=np.float32)
        DGFLOGi   = np.zeros(self.nlines19,dtype=np.float32)
        DGAMMARi  = np.zeros(self.nlines19,dtype=np.float32)
        DGAMMASi  = np.zeros(self.nlines19,dtype=np.float32)
        DGAMMAWi  = np.zeros(self.nlines19,dtype=np.float32)
        DWLISOi   = np.zeros(self.nlines19,dtype=np.float32)
        ISOSHIFTi = np.zeros(self.nlines19,dtype=np.int32)
        EXTRA3i   = np.zeros(self.nlines19,dtype=np.float32)

        self.rfort.readfile14_20(
            c_char_p(sb_i), 
            c_int(self.nlines19),
            WLi.ctypes.data_as(c_double_p),    
            Ei.ctypes.data_as(c_double_p),    
            EPi.ctypes.data_as(c_double_p),    
            LABELi.ctypes.data_as(c_char_p),
            LABELx.ctypes.data_as(c_double_p),  
            LABELPi.ctypes.data_as(c_char_p),
            LABELPx.ctypes.data_as(c_double_p),  
            # OTHER1i.ctypes.data_as(c_int_p),
            # OTHER2i.ctypes.data_as(c_char_p),

            ISHIFTi.ctypes.data_as(c_int_p),
            ISHIFTPi.ctypes.data_as(c_int_p),
            IXFIXFPi.ctypes.data_as(c_char_p),
            LINESIZEi.ctypes.data_as(c_int_p),
            AUTOi.ctypes.data_as(c_char_p),

            OTHER1x.ctypes.data_as(c_double_p),    
            OTHER2x.ctypes.data_as(c_double_p),    
            WLVACi.ctypes.data_as(c_double_p),    
            CENTERi.ctypes.data_as(c_double_p),    
            CONCENi.ctypes.data_as(c_double_p),    
            NELIONi.ctypes.data_as(c_int_p),    
            GAMMARi.ctypes.data_as(c_float_p),    
            GAMMASi.ctypes.data_as(c_float_p),    
            GAMMAWi.ctypes.data_as(c_float_p),    
            REFi.ctypes.data_as(c_char_p),     
            REFx.ctypes.data_as(c_float_p),    
            NBLOi.ctypes.data_as(c_int_p),    
            NBUPi.ctypes.data_as(c_int_p),    
            ISO1i.ctypes.data_as(c_int_p),    
            X1i.ctypes.data_as(c_float_p),    
            ISO2i.ctypes.data_as(c_int_p),    
            X2i.ctypes.data_as(c_float_p),    
            GFLOGi.ctypes.data_as(c_float_p),    
            XJi.ctypes.data_as(c_float_p),    
            XJPi.ctypes.data_as(c_float_p),    
            CODEi.ctypes.data_as(c_float_p),    
            ELOi.ctypes.data_as(c_float_p),    
            GFi.ctypes.data_as(c_float_p),    
            GSi.ctypes.data_as(c_float_p),    
            GRi.ctypes.data_as(c_float_p),    
            GWi.ctypes.data_as(c_float_p),    
            DWLi.ctypes.data_as(c_float_p),    
            DGFLOGi.ctypes.data_as(c_float_p),    
            DGAMMARi.ctypes.data_as(c_float_p),    
            DGAMMASi.ctypes.data_as(c_float_p),    
            DGAMMAWi.ctypes.data_as(c_float_p),    
            DWLISOi.ctypes.data_as(c_float_p),    
            ISOSHIFTi.ctypes.data_as(c_int_p),    
            EXTRA3i.ctypes.data_as(c_float_p),
            )

        # # convert ctype char into python string
        # x = np.array([''.join(LABELi[i,:].tobytes('F').decode('ascii')) for i in range(self.nlines19)])
        # x = ''.join(x)
        # LABELi = np.array(list(map(''.join, zip(*[iter(x)]*11))))

        # x = np.array([''.join(LABELPi[i,:].tobytes('F').decode('ascii')) for i in range(self.nlines19)])
        # x = ''.join(x)
        # LABELPi = np.array(list(map(''.join, zip(*[iter(x)]*11))))

        # x = np.array([''.join(REFi[i,:].tobytes('F').decode('ascii')) for i in range(self.nlines19)])
        # x = ''.join(x)
        # REFi = np.array(list(map(''.join, zip(*[iter(x)]*5))))

        # x = np.array([''.join(OTHER1i[i,:].tobytes('F').decode('ascii')) for i in range(self.nlines19)])
        # x = ''.join(x)
        # OTHER1i = np.array(list(map(''.join, zip(*[iter(x)]*10))))

        # x = np.array([''.join(OTHER2i[i,:].tobytes('F').decode('ascii')) for i in range(self.nlines19)])
        # x = ''.join(x)
        # OTHER2i = np.array(list(map(''.join, zip(*[iter(x)]*10))))

        output = ({
            'wl':WLi, # input wavelength
            'e':Ei, # first energy level
            'ep':EPi, # second energy level
            'label':LABELi, # first level label
            'labelx':LABELx, # unformatted first level label
            'labelp':LABELPi, # second level label
            'labelpx':LABELPx, # unformatted second level label
            # 'other1':OTHER1i, # Other1 info
            # 'other2':OTHER2i, # Other2 info
            'ishift':ISHIFTi,
            'ishiftp':ISHIFTPi,
            'ixfixfp':IXFIXFPi,
            'linesize':LINESIZEi,
            'auto':AUTOi,
            'other1x':OTHER1x, # Other1 info unformatted
            'other2x':OTHER2x, # Other2 info unformatted
            'wlvac':WLVACi, # wavelength in vacuum
            'center':CENTERi, # Intensity at center nu
            'concen':CONCENi, # Continuum at central nu
            'nelion':NELIONi, # Index for specific element
            'gammar':GAMMARi, # unlogged and shifted Gamma R
            'gammas':GAMMASi, # unlogged and shifted Gamma S
            'gammaw':GAMMAWi, # unlogged and shifted Gamma W
            'ref':REFi, # Reference for line
            'refx':REFx, # Unformatted reference for line
            'nblo':NBLOi, # Departure coeff for lower level
            'nbup':NBUPi, # Departure coeff for upper level
            'iso1':ISO1i, # Isotopic number for level 1
            'x1':X1i, # Log of fractional isotopic abundance for level 1
            'iso2':ISO2i, # Isotopic number for level 2
            'x2':X2i, # Log of fractional isotopic abundance for level 2
            'gflog':GFLOGi, # log(gf)
            'xj':XJi, # angular momentum of first level
            'xjp':XJPi, # angular momentum of second level
            'code':CODEi, # numeric code for line species
            'elo':ELOi, # lower energy level
            'gf':GFi, # log(gf) -> unlogged with shifts applied
            'gs':GSi, # log damping constant 
            'gr':GRi, # log damping constant 
            'gw':GWi, # log damping constant 
            'dwl':DWLi, # shift in wl
            'dgflog':DGFLOGi, # shift in log(gf)
            'dgammar':DGAMMARi, # shift in dampling constants 
            'dgammas':DGAMMASi, # shift in dampling constants
            'dgammaw':DGAMMAWi, # shift in dampling constants
            'dwliso':DWLISOi, # shift in isotopic wl shift
            'isoshift':ISOSHIFTi, # isotopic shift
            'extra3':EXTRA3i, # extra information
        })
        return output

    def writefort93(self,f93outpath):

        # turn path string into bytes
        s_o  = '{0}'.format(f93outpath)
        sb_o    = bytes(s_o,encoding='ascii')

        # extract variables from read in dict
        NLINESi  = self.f93in['nlines']
        LENGTHi  = self.f93in['length']
        IFVACi   = self.f93in['ifvac']
        IFNLTEi  = self.f93in['ifnlte']
        N19i     = self.f93in['n19']
        TURBVi   = self.f93in['turbv']
        DECKJi   = self.f93in['deckj']
        IFPREDi  = self.f93in['ifpred']
        WLBEGi   = self.f93in['wlbeg']
        WLENDi   = self.f93in['wlend']
        RESOLUi  = self.f93in['resolution']
        RATIOi   = self.f93in['ratio']
        RATIOLGi = self.f93in['ratiolg']
        CUTOFFi  = self.f93in['cutoff']
        LINOUTi  = self.f93in['linout']

        self.rfort.writefile93(
            c_char_p(sb_o), 
            c_void_p(NLINESi.ctypes.data),    
            c_void_p(LENGTHi.ctypes.data),
            c_void_p(IFVACi.ctypes.data),
            c_void_p(IFNLTEi.ctypes.data),
            c_void_p(N19i.ctypes.data),
            c_void_p(TURBVi.ctypes.data),
            c_void_p(DECKJi.ctypes.data),
            c_void_p(IFPREDi.ctypes.data),
            c_void_p(WLBEGi.ctypes.data),
            c_void_p(WLENDi.ctypes.data),
            c_void_p(RESOLUi.ctypes.data),
            c_void_p(RATIOi.ctypes.data),
            c_void_p(RATIOLGi.ctypes.data),
            c_void_p(CUTOFFi.ctypes.data),
            c_void_p(LINOUTi.ctypes.data),
        )

    def writefort12(self,f12outpath):

        # turn path string into bytes
        s_o  = '{0}'.format(f12outpath)
        sb_o    = bytes(s_o,encoding='ascii')

        # extract variables from read in dict
        NBUFFi    = self.f12in['nbuff']
        CGFi      = self.f12in['cgf']
        NELIONi   = self.f12in['nelion']
        ELOi      = self.f12in['elo']
        GAMMARi  = self.f12in['gammar']
        GAMMASi  = self.f12in['gammas']
        GAMMAWi  = self.f12in['gammaw']
                
        self.rfort.writefile12(
            c_char_p(sb_o), 
            c_int(self.nlines12),
            c_void_p(NBUFFi.ctypes.data),
            c_void_p(CGFi.ctypes.data),
            c_void_p(NELIONi.ctypes.data),
            c_void_p(ELOi.ctypes.data),
            c_void_p(GAMMARi.ctypes.data),
            c_void_p(GAMMASi.ctypes.data),
            c_void_p(GAMMAWi.ctypes.data),
        )
        
    def writefort14(self,f14outpath):

        # turn path string into bytes
        s_o  = '{0}'.format(f14outpath)
        sb_o    = bytes(s_o,encoding='ascii')
        
        WLi       = self.f14in['wl']
        Ei        = self.f14in['e']
        EPi       = self.f14in['ep']
        LABELx    = self.f14in['labelx']
        LABELPx   = self.f14in['labelpx']
        OTHER1x   = self.f14in['other1x']
        OTHER2x   = self.f14in['other2x']
        WLVACi    = self.f14in['wlvac']
        CENTERi   = self.f14in['center']
        CONCENi   = self.f14in['concen']
        NELIONi   = self.f14in['nelion']
        GAMMARi   = self.f14in['gammar']
        GAMMASi   = self.f14in['gammas']
        GAMMAWi   = self.f14in['gammaw']
        REFx      = self.f14in['refx']
        NBLOi     = self.f14in['nblo']
        NBUPi     = self.f14in['nbup']
        ISO1i     = self.f14in['iso1']
        X1i       = self.f14in['x1']
        ISO2i     = self.f14in['iso2']
        X2i       = self.f14in['x2']
        GFLOGi    = self.f14in['gflog']
        XJi       = self.f14in['xj']
        XJPi      = self.f14in['xjp']
        CODEi     = self.f14in['code']
        ELOi      = self.f14in['elo']
        GFi       = self.f14in['gf']
        GSi       = self.f14in['gs']
        GRi       = self.f14in['gr']
        GWi       = self.f14in['gw']
        DWLi      = self.f14in['dwl']
        DGFLOGi   = self.f14in['dgflog']
        DGAMMARi  = self.f14in['dgammar']
        DGAMMASi  = self.f14in['dgammas']
        DGAMMAWi  = self.f14in['dgammaw']
        DWLISOi   = self.f14in['dwliso']
        ISOSHIFTi = self.f14in['isoshift']
        EXTRA3i   = self.f14in['extra3']

        self.rfort.writefile14_20(
            c_char_p(sb_o), 
            c_int(self.nlines12),
            c_void_p(WLi.ctypes.data),    
            c_void_p(Ei.ctypes.data),    
            c_void_p(EPi.ctypes.data),    
            c_void_p(LABELx.ctypes.data),  
            c_void_p(LABELPx.ctypes.data),  
            c_void_p(OTHER1x.ctypes.data),    
            c_void_p(OTHER2x.ctypes.data),    
            c_void_p(WLVACi.ctypes.data),    
            c_void_p(CENTERi.ctypes.data),    
            c_void_p(CONCENi.ctypes.data),    
            c_void_p(NELIONi.ctypes.data),
            c_void_p(GAMMARi.ctypes.data),    
            c_void_p(GAMMASi.ctypes.data),    
            c_void_p(GAMMAWi.ctypes.data),    
            c_void_p(REFx.ctypes.data),    
            c_void_p(NBLOi.ctypes.data),
            c_void_p(NBUPi.ctypes.data),
            c_void_p(ISO1i.ctypes.data),
            c_void_p(X1i.ctypes.data),    
            c_void_p(ISO2i.ctypes.data),
            c_void_p(X2i.ctypes.data),    
            c_void_p(GFLOGi.ctypes.data),    
            c_void_p(XJi.ctypes.data),    
            c_void_p(XJPi.ctypes.data),    
            c_void_p(CODEi.ctypes.data),    
            c_void_p(ELOi.ctypes.data),    
            c_void_p(GFi.ctypes.data),    
            c_void_p(GSi.ctypes.data),    
            c_void_p(GRi.ctypes.data),    
            c_void_p(GWi.ctypes.data),    
            c_void_p(DWLi.ctypes.data),    
            c_void_p(DGFLOGi.ctypes.data),    
            c_void_p(DGAMMARi.ctypes.data),    
            c_void_p(DGAMMASi.ctypes.data),    
            c_void_p(DGAMMAWi.ctypes.data),    
            c_void_p(DWLISOi.ctypes.data),    
            c_void_p(ISOSHIFTi.ctypes.data),
            c_void_p(EXTRA3i.ctypes.data),
            )

    def writefort19(self,f19outpath):

        # turn path string into bytes
        s_o  = '{0}'.format(f19outpath)
        sb_o    = bytes(s_o,encoding='ascii')

        # extract variables from read in dict
        WLVACi   = self.f19in['wlvac']
        ELOi     = self.f19in['elo']
        GFi      = self.f19in['gf']
        NBLOi    = self.f19in['nblo']
        NBUPi    = self.f19in['nbup']
        NELIONi  = self.f19in['nelion']
        TYPEi    = self.f19in['type']
        NCONi    = self.f19in['ncon']
        NELIONXi = self.f19in['nelionx']
        GAMMARi  = self.f19in['gammar']
        GAMMASi  = self.f19in['gammas']
        GAMMAWi  = self.f19in['gammaw']
        NBUFFi   = self.f19in['nbuff']
        LIMi     = self.f19in['lim']

        self.rfort.writefile19(
            c_char_p(sb_o), 
            c_int(self.nlines19),
            c_void_p(WLVACi.ctypes.data),
            c_void_p(ELOi.ctypes.data),
            c_void_p(GFi.ctypes.data),
            c_void_p(NBLOi.ctypes.data),
            c_void_p(NBUPi.ctypes.data),
            c_void_p(NELIONi.ctypes.data),
            c_void_p(TYPEi.ctypes.data),
            c_void_p(NCONi.ctypes.data),
            c_void_p(NELIONXi.ctypes.data),
            c_void_p(GAMMARi.ctypes.data),
            c_void_p(GAMMASi.ctypes.data),
            c_void_p(GAMMAWi.ctypes.data),
            c_void_p(NBUFFi.ctypes.data),
            c_void_p(LIMi.ctypes.data),
        )
        
    def writefort20(self,f20outpath):

        # turn path string into bytes
        s_o  = '{0}'.format(f20outpath)
        sb_o    = bytes(s_o,encoding='ascii')
        
        WLi       = self.f20in['wl']
        Ei        = self.f20in['e']
        EPi       = self.f20in['ep']
        LABELx    = self.f20in['labelx']
        LABELPx   = self.f20in['labelpx']
        OTHER1x   = self.f20in['other1x']
        OTHER2x   = self.f20in['other2x']
        WLVACi    = self.f20in['wlvac']
        CENTERi   = self.f20in['center']
        CONCENi   = self.f20in['concen']
        NELIONi   = self.f20in['nelion']
        GAMMARi   = self.f20in['gammar']
        GAMMASi   = self.f20in['gammas']
        GAMMAWi   = self.f20in['gammaw']
        REFx      = self.f20in['refx']
        NBLOi     = self.f20in['nblo']
        NBUPi     = self.f20in['nbup']
        ISO1i     = self.f20in['iso1']
        X1i       = self.f20in['x1']
        ISO2i     = self.f20in['iso2']
        X2i       = self.f20in['x2']
        GFLOGi    = self.f20in['gflog']
        XJi       = self.f20in['xj']
        XJPi      = self.f20in['xjp']
        CODEi     = self.f20in['code']
        ELOi      = self.f20in['elo']
        GFi       = self.f20in['gf']
        GSi       = self.f20in['gs']
        GRi       = self.f20in['gr']
        GWi       = self.f20in['gw']
        DWLi      = self.f20in['dwl']
        DGFLOGi   = self.f20in['dgflog']
        DGAMMARi  = self.f20in['dgammar']
        DGAMMASi  = self.f20in['dgammas']
        DGAMMAWi  = self.f20in['dgammaw']
        DWLISOi   = self.f20in['dwliso']
        ISOSHIFTi = self.f20in['isoshift']
        EXTRA3i   = self.f20in['extra3']

        self.rfort.writefile14_20(
            c_char_p(sb_o), 
            c_int(self.nlines19),
            c_void_p(WLi.ctypes.data),    
            c_void_p(Ei.ctypes.data),    
            c_void_p(EPi.ctypes.data),    
            c_void_p(LABELx.ctypes.data),  
            c_void_p(LABELPx.ctypes.data),  
            c_void_p(OTHER1x.ctypes.data),    
            c_void_p(OTHER2x.ctypes.data),    
            c_void_p(WLVACi.ctypes.data),    
            c_void_p(CENTERi.ctypes.data),    
            c_void_p(CONCENi.ctypes.data),    
            c_void_p(NELIONi.ctypes.data),
            c_void_p(GAMMARi.ctypes.data),    
            c_void_p(GAMMASi.ctypes.data),    
            c_void_p(GAMMAWi.ctypes.data),    
            c_void_p(REFx.ctypes.data),    
            c_void_p(NBLOi.ctypes.data),
            c_void_p(NBUPi.ctypes.data),
            c_void_p(ISO1i.ctypes.data),
            c_void_p(X1i.ctypes.data),    
            c_void_p(ISO2i.ctypes.data),
            c_void_p(X2i.ctypes.data),    
            c_void_p(GFLOGi.ctypes.data),    
            c_void_p(XJi.ctypes.data),    
            c_void_p(XJPi.ctypes.data),    
            c_void_p(CODEi.ctypes.data),    
            c_void_p(ELOi.ctypes.data),    
            c_void_p(GFi.ctypes.data),    
            c_void_p(GSi.ctypes.data),    
            c_void_p(GRi.ctypes.data),    
            c_void_p(GWi.ctypes.data),    
            c_void_p(DWLi.ctypes.data),    
            c_void_p(DGFLOGi.ctypes.data),    
            c_void_p(DGAMMARi.ctypes.data),    
            c_void_p(DGAMMASi.ctypes.data),    
            c_void_p(DGAMMAWi.ctypes.data),    
            c_void_p(DWLISOi.ctypes.data),    
            c_void_p(ISOSHIFTi.ctypes.data),
            c_void_p(EXTRA3i.ctypes.data),
            )

    def readspecbin_init(self,specbinpath):
        # turn path string into bytes
        s_i     = '{0}'.format(specbinpath)
        sb_i    = bytes(s_i,encoding='ascii')
        
        TEFFi    = np.zeros_like(0.0,dtype=np.float64)
        GLOGi    = np.zeros_like(0.0,dtype=np.float64)
        WBEGINi  = np.zeros_like(0.0,dtype=np.float64)
        RESOLUi  = np.zeros_like(0.0,dtype=np.float64)
        NWLi     = np.zeros_like(0.0,dtype=np.int32)
        IFSURFi  = np.zeros_like(0.0,dtype=np.float64)
        NMUi     = np.zeros_like(0.0,dtype=np.float64)
        NEDGEi   = np.zeros_like(0.0,dtype=np.float64)
        WLEDGEi  = np.zeros_like(0.0,dtype=np.float64)
        NLINESi = np.zeros_like(0.0,dtype=np.int32)

        self.rfort.readspecbin_init(
            c_char_p(sb_i), 
            TEFFi.ctypes.data_as(c_double_p),
            GLOGi.ctypes.data_as(c_double_p),
            WBEGINi.ctypes.data_as(c_double_p),
            RESOLUi.ctypes.data_as(c_double_p),
            NWLi.ctypes.data_as(c_int_p),
            IFSURFi.ctypes.data_as(c_double_p),
            NMUi.ctypes.data_as(c_double_p),
            NEDGEi.ctypes.data_as(c_double_p),
            WLEDGEi.ctypes.data_as(c_double_p),
            NLINESi.ctypes.data_as(c_int_p),
        )

        output = ({
            'teff':TEFFi,
            'logg':GLOGi,
            'wbegin':WBEGINi,
            'resolu':RESOLUi,
            'nwl':NWLi,
            'ifsurf':IFSURFi,
            'nmu':NMUi,
            'nedge':NEDGEi,
            'wledge':WLEDGEi,
            'nlines':NLINESi,
        })
        
        return output

    def readspecbin(self,specbinpath):
        # first run the initializer to grab important info
        
        rsb_in = self.readspecbin_init(specbinpath)

        # turn path string into bytes
        s_i     = '{0}'.format(specbinpath)
        sb_i    = bytes(s_i,encoding='ascii')

        # initialize variables for ctypes to fill 
        NWLi      = np.int32(rsb_in['nwl'])
        NLINESi   = np.int32(rsb_in['nlines'])
        wavei     = np.zeros(rsb_in['nwl'],dtype=np.float64)    
        qmu1i     = np.zeros(rsb_in['nwl'],dtype=np.float64)
        qmu2i     = np.zeros(rsb_in['nwl'],dtype=np.float64)
        WLi       = np.zeros(rsb_in['nlines'],dtype=np.float64)   
        DWLi      = np.zeros(rsb_in['nlines'],dtype=np.float32)   
        CODEi     = np.zeros(rsb_in['nlines'],dtype=np.float32)
        GFLOGi    = np.zeros(rsb_in['nlines'],dtype=np.float32)   
        DGFLOGi   = np.zeros(rsb_in['nlines'],dtype=np.float32)   
        GRi       = np.zeros(rsb_in['nlines'],dtype=np.float32)   
        DGAMMARi  = np.zeros(rsb_in['nlines'],dtype=np.float32)   
        GSi       = np.zeros(rsb_in['nlines'],dtype=np.float32)   
        DGAMMASi  = np.zeros(rsb_in['nlines'],dtype=np.float32)   
        GWi       = np.zeros(rsb_in['nlines'],dtype=np.float32)   
        DGAMMAWi  = np.zeros(rsb_in['nlines'],dtype=np.float32)   
        RESIDi   = np.zeros(rsb_in['nlines'],dtype=np.float64)   


        self.rfort.readspecbin(
            c_char_p(sb_i), 
            c_int(NWLi),
            c_int(NLINESi),
            wavei.ctypes.data_as(c_double_p),    
            qmu1i.ctypes.data_as(c_double_p),    
            qmu2i.ctypes.data_as(c_double_p),    
            WLi.ctypes.data_as(c_double_p), 
            DWLi.ctypes.data_as(c_float_p), 
            GFLOGi.ctypes.data_as(c_float_p), 
            DGFLOGi.ctypes.data_as(c_float_p),
            CODEi.ctypes.data_as(c_float_p),
            GRi.ctypes.data_as(c_float_p),    
            DGAMMARi.ctypes.data_as(c_float_p),    
            GSi.ctypes.data_as(c_float_p),    
            DGAMMASi.ctypes.data_as(c_float_p),    
            GWi.ctypes.data_as(c_float_p),    
            DGAMMAWi.ctypes.data_as(c_float_p),    
            RESIDi.ctypes.data_as(c_double_p),    
        )

        output = ({
            'wave':wavei, 
            'qmu1':qmu1i, 
            'qmu2':qmu2i,
            'code':CODEi,
            'wl':WLi,
            'dwl':DWLi,
            'loggf':GFLOGi,
            'dloggf':DGFLOGi,
            'gammar':GRi,
            'dgammar':DGAMMARi, 
            'gammas':GSi,
            'dgammas':DGAMMASi, 
            'gammaw':GWi,
            'dgammaw':DGAMMAWi, 
            'resid':RESIDi,
        })

        return output

        
if __name__ == '__main__':
    RK = ReadKurucz()
    print('1')
    RK.readfiles()
    print('2')
    RK.writefiles()

