from . import readkurucz
import numpy as np
from collections import Counter

class AdjKurucz(object):
    def __init__(self, *args, **kwargs):
        super(AdjKurucz, self).__init__()

        self.args = args
        self.kwargs = kwargs

        self.verbose = kwargs.get('verbose',False)

        self.f12path=self.kwargs.get('f12path','./fort.12')
        self.f14path=self.kwargs.get('f14path','./fort.14')
        self.f19path=self.kwargs.get('f19path','./fort.19')
        self.f20path=self.kwargs.get('f20path','./fort.20')
        self.f93path=self.kwargs.get('f93path','./fort.93')

        self.RK = readkurucz.ReadKurucz(verbose=self.verbose)
        
        self.RK.readfiles(
            f12path=self.f12path,
            f14path=self.f14path,
            f19path=self.f19path,
            f20path=self.f20path,
            f93path=self.f93path)

    def wfort(self,*args,**kwargs):
        outf12path=kwargs.get('f12path','./fort_new.12')
        outf14path=kwargs.get('f14path','./fort_new.14')
        outf19path=kwargs.get('f19path','./fort_new.19')
        outf20path=kwargs.get('f20path','./fort_new.20')
        outf93path=kwargs.get('f93path','./fort_new.93')

        self.RK.writefiles(            
            f12outpath=outf12path,
            f14outpath=outf14path,
            f19outpath=outf19path,
            f20outpath=outf20path,
            f93outpath=outf93path)

    def adj93(self,newdict={}):
        """
        Function to adjust the global parameters in fort.93 
        and have that trickle down to the other files
        """

        if newdict == {}:
            return
                
        # adjust the wavelength range
        if 'wl' in newdict.keys():
            wlbeg = newdict['wl'][0]
            wlend = newdict['wl'][1]
        
            self.RK.f93in['wlbeg'] = np.array(wlbeg,dtype=np.float64)
            self.RK.f93in['wlend'] = np.array(wlend,dtype=np.float64)
        
            # redetermine the bounds for lines to be included
            ratiolg = self.RK.f93in['ratiolg']

            ixwlbeg = np.floor(np.log(wlbeg)/ratiolg)
            wbegin = np.exp(ixwlbeg*ratiolg)
            if wbegin < wlbeg:
                ixwlbeg = ixwlbeg + 1
            ixwlend = np.floor(np.log(wlend)/ratiolg)
            wlast = np.exp(ixwlend*ratiolg)
            if wlast >= wlend:
                ixwlend =  ixwlend - 1
            length = ixwlend-ixwlbeg+1

            # reset the length
            self.RK.f93in['length'] = np.array(length,dtype=np.int32)

            # trim the line lists based on new wavelength range
            if wlbeg > 500.0:
                delfactor = wlbeg / 500.0
            else:
                delfactor = 1.0

            DELIM_arr = [100.0,30.0,10.0,3.0,1.0,0.3,0.1]
            linsizep = 8.0-self.RK.f14in['linesize']
            lim = np.array(np.minimum(linsizep,7),dtype=int)
            delim = np.array([DELIM_arr[lim_i-1] for lim_i in lim])

            offset_beg = (wlbeg - delim*delfactor)
            offset_end = (wlend + delim*delfactor)
            
            condmol = self.RK.f14in['code'] > 100.0
            offset_beg[condmol] = (wlbeg - 0.1)
            offset_end[condmol] = (wlend + 1.0)

            condTiO = self.RK.f14in['code'] == 822.0
            offset_beg[condTiO] = (wlbeg - 1.145)
            offset_end[condTiO] = (wlend + 1.0)

            # select which stars are still included
            cond14 = (self.RK.f14in['wlvac'] >= offset_beg) & (self.RK.f14in['wlvac'] <= offset_end)

            # calculate NLINES
            self.RK.f93in['nlines'] = np.array(cond14.sum(),dtype=np.int32)

            # correct the variables nlines12/19
            self.RK.nlines12 = self.RK.f93in['nlines']
            self.RK.nlines19 = self.RK.f93in['n19']

            # parse the fort.12/14 lines
            for kk in self.RK.f12in.keys():
                self.RK.f12in[kk] = self.RK.f12in[kk][cond14]
            for kk in self.RK.f14in.keys():
                if self.RK.f14in[kk].ndim == 1:
                    self.RK.f14in[kk] = self.RK.f14in[kk][cond14]
                else:
                    if self.RK.f14in[kk].shape[0] == 2:
                        if kk == 'other1x':
                            print('2', kk)
                            x = self.RK.f14in[kk][:,np.flatnonzero(cond14)]
                            print(x[0][279],x[1][279])
                        self.RK.f14in[kk] = self.RK.f14in[kk][:,np.flatnonzero(cond14)]
                    elif self.RK.f14in[kk].shape[1] == 10:
                        self.RK.f14in[kk] = self.RK.f14in[kk][np.flatnonzero(cond14),:]
                    elif self.RK.f14in[kk].shape[1] == 5:
                        self.RK.f14in[kk] = self.RK.f14in[kk][np.flatnonzero(cond14),:]
                    else:
                        print(self.RK.f14in[kk].ndim,self.RK.f14in[kk].shape)
                        print('DID NOT WORK')
                        raise IOError
                        
            # recalculate the line buffer in fort.12
            nbuff_n = []
            for wl in self.RK.f14in['wl']:
                ixwl = np.log(wl)/ratiolg + 0.5
                nbuff_n.append(np.floor(ixwl - ixwlbeg + 1))
            self.RK.f12in['nbuff'] = np.array(nbuff_n,dtype=np.int32)

            # recalculate the line buffer for lines in fort.19
            nbuff_n = []
            for wl in self.RK.f20in['wl']:
                ixwl = np.log(wl)/ratiolg + 0.5
                nbuff_n.append(np.floor(ixwl - ixwlbeg + 1))
            self.RK.f19in['nbuff'] = np.array(nbuff_n,dtype=np.int32)


        # change the opacity binning resolution (lambda/delta(lambda))
        if 'res' in newdict.keys():
            res     = newdict['res']            
            ratio   = 1.0 + (1.0 / res)
            ratiolg = np.log(ratio)

            wlbeg = self.RK.f93in['wlbeg']
            wlend = self.RK.f93in['wlend']
        
            ixwlbeg = np.floor(np.log(wlbeg)/ratiolg)
            wbegin = np.exp(ixwlbeg*ratiolg)
            if wbegin < wlbeg:
                ixwlbeg = ixwlbeg + 1
            ixwlend = np.floor(np.log(wlend)/ratiolg)
            wlast = np.exp(ixwlend*ratiolg)
            if wlast >= wlend:
                ixwlend =  ixwlend - 1
            length = ixwlend-ixwlbeg+1

            self.RK.f93in['resolution'] = np.array(res,dtype=np.float64)
            self.RK.f93in['ratio']      = np.array(ratio,dtype=np.float64)
            self.RK.f93in['ratiolg']    = np.array(ratiolg,dtype=np.float64)
            self.RK.f93in['length']     = np.array(length,dtype=np.int32)

            # recalculate the line buffer in fort.12
            nbuff_n = []
            for wl in self.RK.f14in['wl']:
                ixwl = np.log(wl)/ratiolg + 0.5
                nbuff_n.append(np.floor(ixwl - ixwlbeg + 1))
            self.RK.f12in['nbuff'] = np.array(nbuff_n,dtype=np.int32)

            # recalculate the line buffer for lines in fort.19
            nbuff_n = []
            for wl in self.RK.f20in['wl']:
                ixwl = np.log(wl)/ratiolg + 0.5
                nbuff_n.append(np.floor(ixwl - ixwlbeg + 1))
            self.RK.f19in['nbuff'] = np.array(nbuff_n,dtype=np.int32)

        # adjust the microturbulence on top of what is in the atm files
        if 'vmic' in newdict.keys():
            self.RK.f93in['turbv'] = np.array(newdict['vmic'],dtype=np.float64)

        # adjust the cut off for including lines
        if 'cutoff' in newdict.keys():
            self.RK.f93in['cutoff'] = np.array(newdict['cutoff'],dtype=np.float64)

    def adjll(self,lindict={}):
        '''
        Function that takes in user def adjustments to line parameters.
        lindict can include either:
        
        linind -> int or array of int's with the index (indices) of line(s)
        pars   -> set of line parameter and code matches these to find line ind
        
        '''
        
        if lindict == {}:
            return
        
        if ('pars' in lindict.keys()) and not ('lineind' in lindict.keys()):
            # index within synthe line files (fort.12,14,etc.)
            lindict['lineind'] = np.array([])    
            lindict['rlte'] = []

            # user wants to match line using its parameters
            cond14 = np.ones(len(self.RK.f14in['wl']),dtype=bool)
            for kk in lindict['pars'].keys():
                if self.verbose:
                    print('Matching fort.14 on {}'.format(kk))
                cond14_i = np.in1d(self.RK.f14in[kk],lindict['pars'][kk])
                cond14 *= cond14_i
                if self.verbose:
                    print('Found {} matches'.format(cond14_i.sum()))

            if cond14.sum() == 1:
                lindict['rlte'] = lindict['rlte']+[False]
                lindict['lineind'] = np.append(lindict['lineind'],np.flatnonzero(cond14))
            elif cond14.sum() > 1:
                lindict['rlte'] = lindict['rlte']+[False for _ in range(cond14.sum())]
                lindict['lineind'] = np.append(lindict['lineind'],np.flatnonzero(cond14))
            else:
                pass

            cond20 = np.ones(len(self.RK.f20in['wl']),dtype=bool)
            for kk in lindict['pars'].keys():
                if self.verbose:
                    print('Matching fort.20 on {}'.format(kk))
                cond20_i = np.in1d(self.RK.f20in[kk],lindict['pars'][kk])
                cond20 *= cond20_i
                if self.verbose:
                    print('Found {} matches'.format(cond20_i.sum()))

            if cond20.sum() == 1:
                lindict['rlte'] = lindict['rlte']+[True]
                lindict['lineind'] = np.append(lindict['lineind'],np.flatnonzero(cond20))
            elif cond20.sum() > 1:
                lindict['rlte'] = lindict['rlte']+[True for _ in range(cond20.sum())]
                lindict['lineind'] = np.append(lindict['lineind'],np.flatnonzero(cond20))
            else:
                pass

            if cond14.sum() + cond20.sum() == 0:                
                if self.verbose:
                    print('! DID NOT FIND ANY LINES BASED ON USER PARS !')
                    print('! NOT ADJUSTING ANYTHING !')
                return
            if self.verbose:
                print('Found a total of {} matches'.format(cond14.sum() + cond20.sum()))
                
        if 'lineind' in lindict.keys():
            # user is using line indices to id lines
            if isinstance(lindict['lineind'],int):
                
                ind = lindict['lineind']

                fortfile = 12
                if 'rlte' in lindict.keys():
                    if lindict['rlte'][0]:
                        fortfile = 19
                        
                if 'dwl' in lindict.keys():
                    self.adjwl(ind,lindict['dwl'],fort=fortfile)
                if 'dloggf' in lindict.keys():
                    self.adjloggf(ind,lindict['dloggf'],fort=fortfile)
                if 'dgammaw' in lindict.keys():
                    self.adjgammaw(ind,lindict['dgammaw'],fort=fortfile)
                if 'dgammar' in lindict.keys():
                    self.adjgammar(ind,lindict['dgammar'],fort=fortfile)
                if 'dgammas' in lindict.keys():
                    self.adjgammas(ind,lindict['dgammas'],fort=fortfile)
            else:
                # MAY WANT TO VECTORIZE THIS STEP
                for ii,ind in enumerate(lindict['lineind']):
                    # make sure ind is int
                    ind = int(ind)
                    
                    # figure out if the user wants a global shift or 
                    # individual line shifts

                    deltapar = {}
                                        
                    for pp in ['dwl','dloggf','dgammaw','dgammar','dgammas']:
                        if pp in lindict.keys():
                            if isinstance(lindict[pp],float):
                                # global shift
                                deltapar[pp] = lindict[pp]
                            elif len(lindict[pp]) == len(lindict['lineind']):
                                deltapar[pp] = lindict[pp][ii]
                            else:
                                print('! NUMBER OF LINE SHIFTS MUST BE 1 (GLOBAL) OR EQUAL TO NUMBER OF LINES SELECTED !')
                                raise IOError
                        else:
                            deltapar[pp] = 0.0

                    fortfile = 12
                    if 'rlte' in lindict.keys():
                        if lindict['rlte'][ii]:
                            fortfile = 19
                    
                    if 'dwl' in lindict.keys():
                        self.adjwl(ind,deltapar['dwl'],fort=fortfile)
                    if 'dloggf' in lindict.keys():
                        self.adjloggf(ind,deltapar['dloggf'],fort=fortfile)
                    if 'dgammaw' in lindict.keys():
                        self.adjgammaw(ind,deltapar['dgammaw'],fort=fortfile)
                    if 'dgammar' in lindict.keys():
                        self.adjgammar(ind,deltapar['dgammar'],fort=fortfile)
                    if 'dgammas' in lindict.keys():
                        self.adjgammas(ind,deltapar['dgammas'],fort=fortfile)
                    
        else:
            print('! MUST INCLUDE EITEHR linind or pars as input !')
            print('! NOT ADJUSTING ANYTHING !')
            return
        
    def adjwl(self,linind,dwl,fort=None):
        # construct the shift in the line buffer
        wlbeg   = self.RK.f93in['wlbeg']
        ratiolg = self.RK.f93in['ratiolg']
        ixwlbeg = np.log(wlbeg)/ratiolg
        if np.exp(ixwlbeg*ratiolg) < wlbeg:
            ixwlbeg += 1
        
        if fort == 12:
            ixwl = np.log(self.RK.f14in['wlvac'][linind] + dwl)/ratiolg + 0.5
            nbuff = ixwl - ixwlbeg + 1
            
            self.RK.f12in['nbuff'][linind] = nbuff
            self.RK.f14in['wl'][linind]    = self.RK.f14in['wl'][linind] + dwl
            self.RK.f14in['wlvac'][linind] = self.RK.f14in['wlvac'][linind] + dwl

            self.RK.f14in['dwl'][linind]   = dwl


        if fort == 19:
            ixwl = np.log(self.RK.f19in['wlvac'][linind] + dwl)/ratiolg + 0.5
            nbuff = ixwl - ixwlbeg + 1

            self.RK.f19in['nbuff'][linind] = nbuff
            self.RK.f19in['wlvac'][linind] = self.RK.f19in['wlvac'][linind] + dwl
            self.RK.f20in['wl'][linind]    = self.RK.f20in['wl'][linind] + dwl
            self.RK.f20in['wlvac'][linind] = self.RK.f20in['wlvac'][linind] + dwl

            self.RK.f20in['dwl'][linind]   = dwl

            
    def adjloggf(self,linind,dloggf,fort=None):
        # construct a dlog(gf) -> 10^dlog(gf) as gf's are stored 
        # instead of log(gf)
        dgf = 10.0**dloggf

        if fort == 12:
            # shift log(gf) by dlog(gf) in terms of gf
            self.RK.f12in['cgf'][linind]    = self.RK.f12in['cgf'][linind] * dgf

            self.RK.f14in['gf'][linind]     = self.RK.f14in['gf'][linind] * dgf
            self.RK.f14in['gflog'][linind]  = self.RK.f14in['gflog'][linind] + dloggf
            self.RK.f14in['dgflog'][linind] = dloggf

        if fort == 19:
            self.RK.f19in['gf'][linind]     = self.RK.f19['gf'][linind] * dgf
            
            self.RK.f20in['gf'][linind]     = self.RK.f20in['gf'][linind] * dgf
            self.RK.f20in['gflog'][linind]  = self.RK.f20in['gflog'][linind] + dloggf
            self.RK.f20in['dgflog'][linind] = dloggf
    
    def adjgammaw(self,linind,dgammaw,fort=None):
        # construct a dgammw -> 10^dgammw as gammaw's are stored 
        # instead of log(gammaw)
        dgw = 10.0**dgammaw
        
        if fort == 12:
            self.RK.f12in['gammaw'][linind]  = self.RK.f12in['gammaw'][linind] * dgw

            self.RK.f14in['gammaw'][linind]  = self.RK.f14in['gammaw'][linind] * dgw
            self.RK.f14in['gw'][linind]      = self.RK.f14in['gw'][linind] + dgammaw
            self.RK.f14in['dgammaw'][linind] = dgammaw
            
        if fort == 19:
            self.RK.f19in['gammaw'][linind]  = self.RK.f19in['gammaw'][linind] * dgw

            self.RK.f20in['gammaw'][linind]  = self.RK.f20in['gammaw'][linind] * dgw
            self.RK.f20in['gw'][linind]      = self.RK.f20in['gw'][linind] + dgammaw
            self.RK.f20in['dgammaw'][linind] = dgammaw

    def adjgammas(self,linind,dgammas,fort=None):
        dgs = 10.0**dgammas
        
        if fort == 12:
            self.RK.f12in['gammas'][linind]  = self.RK.f12in['gammas'][linind] * dgs

            self.RK.f14in['gammas'][linind]  = self.RK.f14in['gammas'][linind] * dgs
            self.RK.f14in['gs'][linind]      = self.RK.f14in['gs'][linind] + dgammas
            self.RK.f14in['dgammas'][linind] = dgammas
            
        if fort == 19:
            self.RK.f19in['gammas'][linind]  = self.RK.f19in['gammas'][linind] * dgs

            self.RK.f20in['gammas'][linind]  = self.RK.f20in['gammas'][linind] * dgs
            self.RK.f20in['gs'][linind]      = self.RK.f20in['gs'][linind] + dgammas
            self.RK.f20in['dgammas'][linind] = dgammas

    def adjgammar(self,linind,dgammar,fort=None):
        dgr = 10.0**dgammar
        
        if fort == 12:
            self.RK.f12in['gammar'][linind]  = self.RK.f12in['gammar'][linind] * dgr

            self.RK.f14in['gammar'][linind]  = self.RK.f14in['gammar'][linind] * dgr
            self.RK.f14in['gr'][linind]      = self.RK.f14in['gr'][linind] + dgammar
            self.RK.f14in['dgammar'][linind] = dgammar
            
        if fort == 19:
            self.RK.f19in['gammar'][linind]  = self.RK.f19in['gammar'][linind] * dgr

            self.RK.f20in['gammar'][linind]  = self.RK.f20in['gammar'][linind] * dgr
            self.RK.f20in['gr'][linind]      = self.RK.f20in['gr'][linind] + dgammar
            self.RK.f20in['dgammar'][linind] = dgammar
    
if __name__ == '__main__':
    AK = AdjKurucz()

    linprt = 25

    sind = np.argsort(AK.RK.f14in['wl'])
    print(AK.RK.f14in['code'][sind][:linprt])
    print(AK.RK.f14in['wl'][sind][:linprt])
    print(AK.RK.f12in['nbuff'][sind][:linprt])
    print(AK.RK.f14in['code'][sind][-linprt:])
    print(AK.RK.f14in['wl'][sind][-linprt:])
    print(AK.RK.f12in['nbuff'][sind][-linprt:])
    for kk in AK.RK.f93in.keys():
        if kk != 'deckj':
            print(kk,AK.RK.f93in[kk])

    AK.adj93({'wl':[517.0,518.0],'res':50000.0})
    
    print('')
    print('')
    
    sind = np.argsort(AK.RK.f14in['wl'])
    print(AK.RK.f14in['code'][sind][:linprt])
    print(AK.RK.f14in['wl'][sind][:linprt])
    print(AK.RK.f12in['nbuff'][sind][:linprt])
    print(AK.RK.f14in['code'][sind][-linprt:])
    print(AK.RK.f14in['wl'][sind][-linprt:])
    print(AK.RK.f12in['nbuff'][sind][-linprt:])
    for kk in AK.RK.f93in.keys():
        if kk != 'deckj':
            print(kk,AK.RK.f93in[kk])
    
    AK.wfort()