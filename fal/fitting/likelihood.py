import numpy as np
from jax import jit
from scipy import constants
speedoflight = constants.c / 1000.0
from datetime import datetime

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

        # check to see if user only input in one dictionary
        if isinstance(self.inputspecinfo,dict):
            self.inputspecinfo = [self.inputspecinfo for _ in range(self.nspec)]

        self.specinfo = {}
        self.specinfo['modatm']  = []
        self.specinfo['macvel']  = []
        self.specinfo['rotvel']  = []
        self.specinfo['R'] = []
        self.specinfo['C12/C13'] = []
        self.specinfo['rvshiftbool'] = []
        self.specinfo['scalebool'] = []
        
        for ii in range(self.nspec):
            inspecinfo = self.inputspecinfo[ii]
            if 'modatm' in inspecinfo.keys():
                self.specinfo['modatm'].append(inspecinfo['modatm'])
            else:
                self.specinfo['modatm'].append('./data/atmmod.dat')

            if 'macvel' in inspecinfo.keys():
                self.specinfo['macvel'].append(inspecinfo['macvel'])
            else:
                self.specinfo['macvel'].append(0.0)
        
            if 'rotvel' in inspecinfo.keys():
                self.specinfo['rotvel'].append(inspecinfo['rotvel'])
            else:
                self.specinfo['rotvel'].append(0.0)

            if 'R' in inspecinfo.keys():
                self.specinfo['R'].append(inspecinfo['R'])
            else:
                self.specinfo['R'].append(3E+5)
            
            if 'C12/C13' in inspecinfo.keys():
                self.specinfo['C12/C13'].append(inspecinfo['C12/C13'])
            else:
                self.specinfo['C12/C13'].append(None)

            if 'rvshiftbool' in inspecinfo.keys():
                self.specinfo['rvshiftbool'].append(inspecinfo['rvshiftbool'])
            else:
                self.specinfo['rvshiftbool'].append(False)

            if 'scalebool' in inspecinfo.keys():
                self.specinfo['scalebool'].append(inspecinfo['scalebool'])
            else:
                self.specinfo['scalebool'].append(False)

        
        # init runsynthe for each spectrum in the input list
        
        self.RSarr = []

        for ii in range(self.nspec):
            inputdict = {}
            inputdict['verbose']   = False
            inputdict['exedir']    = './bin/'
            inputdict['molecules'] = './data/molecules.dat'
            inputdict['continuua'] = './data/continuua.dat'
            inputdict['he1tables'] = './data/he1tables.dat'
            inputdict['spectrv_infile'] = './data/spectrv.input'

            inputdict['atmmod'] = self.specinfo['modatm'][ii]
            inputdict['rotvel'] = self.specinfo['rotvel'][ii]

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
            f12path='./ff/fort.12',
            f14path='./ff/fort.14',
            f19path='./ff/fort.19',
            f20path='./ff/fort.20',
            f93path='./ff/fort.93',
            )

        print(f'... READ IN LINES {len(self.AK.RK.f14in["wl"]+len(self.AK.RK.f20in["wl"]))}')
        print(f'... Aprox MIN WL {self.AK.RK.f14in["wl"].min()}')
        print(f'... Aprox MAX WL {self.AK.RK.f14in["wl"].max()}')

    def genmod(self,linepars={'dwl':[],'dloggf':[],'dgammaw':[],'dgammar':[],'dgammas':[]}):
        
        # adjust the line parameters
        
        indictll = ({
            'lineind':self.lineindex,
            'dwl':linepars['dwl'],
            'dloggf':linepars['dloggf'],
            'dgammaw':linepars['dgammaw'],
            'dgammar':linepars['dgammar'],
            'dgammas':linepars['dgammas'],
            })
        print(indictll)
        self.AK.adjll(lindict=indictll)
        print('WL',self.AK.RK.f14in['wl'][indictll['lineind']])
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
            mod_i = self.RSarr[ii]()            
                # # f12path='./ff/fort_tmp.12',
                # # f14path='./ff/fort_tmp.14',
                # # f19path='./ff/fort_tmp.19',
                # # f20path='./ff/fort_tmp.20',
                # # f93path='./ff/fort_tmp.93',
                # )
            wave = mod_i['wave']
            flux = mod_i['qmu1']/mod_i['qmu2']
            modarr.append([wave,flux])
        cond = (wave > 517.0) & (wave < 517.1)
        print(wave[cond])
        print(flux[cond])
    
        return modarr
    
    def calcchisq(self,modarr):

        # compute chi-sq summing over each input spectrum
        chisq = 0.0
        for ii in range(self.nspec):
            # interpolate mod onto observed wave
            mod_i = np.interp(self.data_wave[ii],modarr[ii][0],modarr[ii][1])
            
            chisq += np.sum( ((mod_i-self.data_flux[ii])/self.data_eflux[ii])**2.0 )
        
        return chisq

    def run(self,pars):
                        
        # last sets of pars are always RV shifts and then scale factors
        rvsnum = np.array(self.specinfo['rvshiftbool']).sum()
        scanum = np.array(self.specinfo['scalebool']).sum()
        
        rvpars = pars[-(rvsnum+scanum):-scanum]
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
        starttime = datetime.now()
        modarr = self.genmod(linepars=linepars)
        print(f'... model eval {datetime.now()-starttime}',flush=True)
        
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
        starttime = datetime.now()
        chisq = self.calcchisq(modarr)
        print(f'... chisq eval {datetime.now()-starttime}',flush=True)
        
        # return the likelihood
        return (-0.5 * chisq, modarr)