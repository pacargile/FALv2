import numpy as np
from jax import jit
from scipy import constants
speedoflight = constants.c / 1000.0
from datetime import datetime
from astropy.table import Table

from ..utils.runsynthe import Synthe
from ..utils.adjkurucz import AdjKurucz

class Like(object):
    """ 
    
    Calculate Likelihood Probability

    """
    
    def __init__(self, *args, **kwargs):
        super(Like, self).__init__()

        self.args = args
        self.kwargs = kwargs

        # pull out input data, this can either be arrays (if fitting one spectrum)
        # or list of arrays with each spectrum is a seperate star
        
        self.data_wave  = self.args[0]
        self.data_flux  = self.args[1]
        self.data_eflux = self.args[2]
        
        # if these are just a single array for one spectrum, wrap in list
        if not isinstance(self.data_wave,list):
            self.data_wave  = [self.data_wave]
            self.data_flux  = [self.data_flux]
            self.data_eflux = [self.data_eflux]
        
        # determine how many spec are being input
        
        self.nspec = len(self.data_wave)
        
        # get line indices for free parameters

        self.lineindex = self.kwargs.get('lineindex',[])
        
        # get index for free parameters associated with each line in lineindex
        
        self.linepars = self.kwargs.get('linepars',[])
        
        # pull out the transmission spectrum
        
        self.trans = self.kwargs.get('trans',None)
        
        # look for user defined info for each spectra in input list        
        self.inputspecinfo = self.kwargs.get('specinfo',[{}])

        self.verbose = self.kwargs.get('verbose',False)

        # check to see if user only input in one dictionary
        if isinstance(self.inputspecinfo,dict):
            self.inputspecinfo = [self.inputspecinfo for _ in range(self.nspec)]

        self.specinfo = {}
        self.specinfo['modatm']  = []
        self.specinfo['vmac']  = []
        self.specinfo['rotvel']  = []
        self.specinfo['R'] = []
        self.specinfo['isofrac'] = []
        self.specinfo['rvshiftbool'] = []
        self.specinfo['scalebool'] = []
        self.specinfo['weaklinemod'] = []
        self.specinfo['transmod'] = []
        
        for ii in range(self.nspec):
            inspecinfo = self.inputspecinfo[ii]
            if 'modatm' in inspecinfo.keys():
                self.specinfo['modatm'].append(inspecinfo['modatm'])
            else:
                self.specinfo['modatm'].append('./data/atmmod.dat')

            if 'vmac' in inspecinfo.keys():
                self.specinfo['vmac'].append(inspecinfo['vmac'])
            else:
                self.specinfo['vmac'].append(0.0)
        
            if 'rotvel' in inspecinfo.keys():
                self.specinfo['rotvel'].append(inspecinfo['rotvel'])
            else:
                self.specinfo['rotvel'].append(0.0)

            if 'R' in inspecinfo.keys():
                self.specinfo['R'].append(inspecinfo['R'])
            else:
                self.specinfo['R'].append(0.0)
            
            if 'isofrac' in inspecinfo.keys():
                self.specinfo['isofrac'].append(inspecinfo['isofrac'])
            else:
                self.specinfo['isofrac'].append(None)

            if 'rvshiftbool' in inspecinfo.keys():
                self.specinfo['rvshiftbool'].append(inspecinfo['rvshiftbool'])
            else:
                self.specinfo['rvshiftbool'].append(False)

            if 'scalebool' in inspecinfo.keys():
                self.specinfo['scalebool'].append(inspecinfo['scalebool'])
            else:
                self.specinfo['scalebool'].append(False)

            if 'weaklinemod' in inspecinfo.keys():
                weaklinepath = inspecinfo['weaklinemod']
            else:
                weaklinepath = './data/specWL_at12.fits'
            weaklinemod_i = Table.read(weaklinepath,format='fits')
            self.specinfo['weaklinemod'].append(weaklinemod_i)

        transmodpath = self.kwargs.get('transmod','./data/transmod.fits')            
        self.transmod = Table.read(transmodpath,format='fits')
                
        # init runsynthe for each spectrum in the input list
        
        self.RSarr = []

        self.numpix = 0

        for ii in range(self.nspec):
            # compute number of pixels
            self.numpix += len(self.data_wave[ii])
            
            inputdict = {}
            inputdict['verbose']   = False
            inputdict['exedir']    = './bin/'
            inputdict['molecules'] = './data/molecules.dat'
            inputdict['continuua'] = './data/continuua.dat'
            inputdict['he1tables'] = './data/he1tables.dat'
            inputdict['spectrv_infile'] = './data/spectrv.input'

            inputdict['isofrac']   = self.specinfo['isofrac'][ii]
            inputdict['atmmod'] = self.specinfo['modatm'][ii]
            inputdict['rotvel'] = self.specinfo['rotvel'][ii]
            inputdict['R'] = self.specinfo['R'][ii]
            inputdict['vmac'] = self.specinfo['vmac'][ii]
            inputdict['synspeed'] = 'fast'

            # init the class            
            RS_i = Synthe(**inputdict)
            
            RS_i.setfpaths(
                f12path='./ff/fort_tmp.12',
                f14path='./ff/fort_tmp.14',
                f19path='./ff/fort_tmp.19',
                f20path='./ff/fort_tmp.20',
                f93path='./ff/fort_tmp.93',
                )
            
            # jit it so that repeat calls will be much faster
            # JSrun = jit(RS_i.run)
            JSrun = RS_i.run
            self.RSarr.append(JSrun)

        # init the adjust kurucz class
        self.AK = AdjKurucz(
            f12path=self.kwargs.get('f12path','./ff/fort.12'),
            f14path=self.kwargs.get('f14path','./ff/fort.14'),
            f19path=self.kwargs.get('f19path','./ff/fort.19'),
            f20path=self.kwargs.get('f20path','./ff/fort.20'),
            f93path=self.kwargs.get('f93path','./ff/fort.93'),
            )

        if self.verbose:
            print(f'... Read in Number of lines: {len(self.AK.RK.f14in["wl"]+len(self.AK.RK.f20in["wl"]))}')
            print(f'... Aprox MIN Line WL {self.AK.RK.f14in["wl"].min()}')
            print(f'... Aprox MAX Line WL {self.AK.RK.f14in["wl"].max()}')

    def genmod(self,linepars={'dwl':[],'dloggf':[],'dgammaw':[],'dgammar':[],'dgammas':[]}):
        
        # self.AK.initfiles()
        
        # adjust the line parameters
        indictll = ({
            'lineind':self.lineindex,
            'dwl':linepars['dwl'],
            'dloggf':linepars['dloggf'],
            'dgammaw':linepars['dgammaw'],
            'dgammar':linepars['dgammar'],
            'dgammas':linepars['dgammas'],
            })
        # print(indictll)
        self.AK.adjll(lindict=indictll)
        # print('WL',self.AK.RK.f14in['wl'][indictll['lineind']])
        # write the tmp fortfiles
        self.AK.wfort(        
            f12path='./ff/fort_tmp.12',
            f14path='./ff/fort_tmp.14',
            f19path='./ff/fort_tmp.19',
            f20path='./ff/fort_tmp.20',
            f93path='./ff/fort_tmp.93',
            )
        
        # run synthe for each of the input spectra
        modarr = []
        for ii in range(self.nspec):
            starttime = datetime.now()
            mod_i = self.RSarr[ii]()            
                # f12path='./ff/fort_tmp.12',
                # f14path='./ff/fort_tmp.14',
                # f19path='./ff/fort_tmp.19',
                # f20path='./ff/fort_tmp.20',
                # f93path='./ff/fort_tmp.93',
                # )

            wave = mod_i['wave']
            qmu1 = mod_i['qmu1']
            qmu2 = mod_i['qmu2']
            flux = qmu1/qmu2

            wlmod = np.interp(wave,self.specinfo['weaklinemod'][ii]['wave'],self.specinfo['weaklinemod'][ii]['flux'])
            flux = flux * wlmod

            tmod = np.interp(wave,self.transmod['wave'],self.transmod['flux'])
            flux = flux * tmod

            modarr.append([wave,flux,qmu1,qmu2])
            if self.verbose:
                print(f'--->>> Total Model {ii} Eval Run Time: {datetime.now()-starttime}')

        return modarr
    
    def calcchisq(self,modarr):

        # compute chi-sq summing over each input spectrum
        chisq = 0.0
        for ii in range(self.nspec):
            # interpolate mod onto observed wave
            mod_i = np.interp(self.data_wave[ii],modarr[ii][0],modarr[ii][1])
            
            chisq += np.sum( ((mod_i-self.data_flux[ii])/self.data_eflux[ii])**2.0 )
        
        return chisq

    def run(self,pars_i):
        if isinstance(pars_i,dict):
            pars = []
            for ind in self.lineindex:
                pars.append(pars_i[f'dwl_{ind}'])
                pars.append(pars_i[f'dgf_{ind}'])
                pars.append(pars_i[f'dgw_{ind}'])
        else:
            pars = pars_i
        # last sets of pars are always RV shifts and then scale factors
        rvsnum = np.array(self.specinfo['rvshiftbool']).sum()
        scanum = np.array(self.specinfo['scalebool']).sum()

        if rvsnum > 0:
            rvpars = pars[-(rvsnum+scanum):-scanum]
        if scanum > 0:
            scpars = pars[-scanum:]
        
        # assemble the parameter dictionaries
        
        dwl     = np.zeros(len(self.lineindex),dtype=float)
        dloggf  = np.zeros(len(self.lineindex),dtype=float)
        dgammaw = np.zeros(len(self.lineindex),dtype=float)
        dgammar = np.zeros(len(self.lineindex),dtype=float)
        dgammas = np.zeros(len(self.lineindex),dtype=float)
        
        for ii,pind in enumerate(self.linepars):
            if pind[0] != -1:
                dwl[ii] = pars[pind[0]]
            if pind[1] != -1:
                dloggf[ii] = pars[pind[1]]
            if pind[2] != -1:
                dgammaw[ii] = pars[pind[2]]
            if pind[3] != -1:
                dgammar[ii] = pars[pind[3]]
            if pind[4] != -1:
                dgammas[ii] = pars[pind[4]]
        linepars = {'dwl':dwl,'dloggf':dloggf,'dgammaw':dgammaw,'dgammar':dgammar,'dgammas':dgammas}
        
        # generate the models
        modarr = self.genmod(linepars=linepars)
        
        # find which spectra need a RV shift and scaling
        for ii in range(self.nspec):
            jj = 0
            if self.specinfo['rvshiftbool'][ii]:
                modarr[ii][0] = modarr[ii][0] * (1.0 + (rvpars[jj]/speedoflight))
                jj += 1
            
            kk = 0
            if self.specinfo['scalebool'][ii]:
                modarr[ii][1] = modarr[ii][1] * scpars[kk]
                kk += 1
                
        # calculate the chi-sqaure
        chisq = self.calcchisq(modarr)
        
        # return the likelihood
        return (-0.5 * chisq, modarr)
    
    def compute_loss(self,pars):
        return (self.run(pars)[0]/-0.5)/self.numpix 
    
    def compute_like(self,pars):
        return self.run(pars)[0]
        