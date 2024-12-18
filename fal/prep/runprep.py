import numpy as np
import os
import shutil
import glob
import h5py
import stat
from datetime import datetime
from astropy.table import Table

from ..utils import readkurucz, runsynthe, adjkurucz

from contextlib import contextmanager
@contextmanager
def cwd(path):
    oldpwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)

def is_unique(*lsts):
    arr = np.vstack(lsts)
    _, ind = np.unique(arr, axis=1, return_index=True)
    out = np.zeros(shape=arr.shape[1], dtype=bool)
    out[ind] = True
    return out

class RunPrep(object):
    def __init__(self, *args, **kwargs):
        super(RunPrep, self).__init__()
        self.args = args
        self.kwargs = kwargs

        # user input dir where master ll fort files are kept
        self.masterllfdir = self.kwargs.get('masterfort','./fortfiles/')

        # define line threshold for including as a fitted line, default 1% depth
        self.threshold = self.kwargs.get('threshold',0.01)

        # define path for files that are the same for all SYNTHE
        self.commonfilespath = self.kwargs.get('commonfiles',None)

        # # define the list of atm files to generate fitted lines
        # # the number of atm files here define how many different stars
        # # to consider.
        # self.atmflist = self.kwargs.get('atmlist',['./atm/atmmod_sol.dat'])

        # read in spec info for atm and spectra used to fit
        self.specinfo = self.kwargs.get('specinfo',{})

        self.atmflist = []
        for si_i in self.specinfo:
            self.atmflist.append(si_i['modatm'])

        # define path to SYNTHE exe binaries
        self.masterbinpath = self.kwargs.get('masterbin',None)

        # set the masterll file paths
        self.masterf12path=self.masterllfdir+'/fort.12'
        self.masterf14path=self.masterllfdir+'/fort.14'
        self.masterf19path=self.masterllfdir+'/fort.19'
        self.masterf20path=self.masterllfdir+'/fort.20'
        self.masterf93path=self.masterllfdir+'/fort.93'
        
        # set up a master info directory to store bulk data for masterLL
        if not os.path.exists('./masterinfo'):
            os.makedirs('./masterinfo/')

        # place copies of atm files into ./atm/
        if not os.path.exists('./masterinfo/atm'):
            os.makedirs('./masterinfo/atm')

        # copy atm files into masterLL info dir
        for ff in self.atmflist:
            filename = ff.split('/')[-1]
            dstfile = f'./masterinfo/atm/{filename}'
            srcfile = f'{ff}'
            if not os.path.exists(dstfile):
                shutil.copyfile(srcfile, dstfile)
        
        # copy commonfiles into masterinfo
        if not os.path.exists('./masterinfo/data'):
            os.makedirs('./masterinfo/data')

        # copy atm files into masterLL info dir    
        for ff in ['molecules.dat','continuua.dat','he1tables.dat','spectrv.input']:
            dstfile = f'masterinfo/data/{ff}'
            srcfile = f'{self.commonfilespath}/{ff}'
            if not os.path.exists(dstfile):
                shutil.copyfile(srcfile, dstfile)

        # define location of trans model
        self.transmodpath = self.kwargs.get('transmod',None)
        self.fulltransmod = {}
        with h5py.File(self.transmodpath,'r') as th5:
            self.fulltransmod['wave'] = th5['spec']['WAVE']
            self.fulltransmod['flux'] = th5['spec']['QMU1']/th5['spec']['QMU2']


    def readseg(self,segfile):
        # function that reads in list of segments, the ascii file format is:
        # segnum, starting wl, ending wl
        
        seginfo = np.loadtxt(segfile,
                             dtype={'names':('segnum','start_wl','end_wl'),
                                    'formats':(int,float,float)})
        segdict = {}
        try:
            numsegs = len(seginfo)
            for ii in range(numsegs):
                segdict[ii] = ({'segnum':seginfo['segnum'][ii],
                                                'startwl':seginfo['start_wl'][ii],
                                                'endwl':seginfo['end_wl'][ii]
                                                })
        except TypeError:
            segdict[0] = ({'segnum':seginfo['segnum'],
                                            'startwl':seginfo['start_wl'],
                                            'endwl':seginfo['end_wl']
                                            })

        return segdict

    def buildtree(self,segdict=None):
        # function to build directory tree for segment, copies in all necessary files
        # naming conv for subdirectories is seg_[segnum]_wl[startwl]_[endwl]
        
        # directories to be made:
        # - seg_[segnum]_wl[startwl]_[endwl]/
        #    - data/  
        #        contains common input files for synthe (e.g., molecules.dat)
        #    - samples/
        #        output from sampler
        #    - lineinfo/
        #        contains lineind and linepars files
        #    - mod/
        #        contains the specific model files for each spectrum
        #    - bin/
        #        contains the SYNTHE binaries, need to copy the binaries to this directory

        # loop through segs 
        for kk in segdict.keys():
            # build seg dir
            segdir = 'seg_{segnum}'.format(**segdict[kk])
            print(f'... Building Tree {segdir}',flush=True)

            if not os.path.exists(segdir):
                os.makedirs(segdir)

            # write wl info into text file in seg dir
            with open(f'{segdir}/wlinfo.txt','w') as wlf:
                wlf.write('segnum={segnum},startwl={startwl:.4f},endwl={endwl:.4f}'.format(**segdict[kk]))

            # build seg subdir
            for ff in ['data','samples','lineinfo','atm','bin','ff']:
                subsegdir = '{0}/{1}'.format(segdir,ff)
                if not os.path.exists(subsegdir):
                    os.makedirs(subsegdir)

            # copy common files
            for ff in ['molecules.dat','continuua.dat','he1tables.dat','spectrv.input']:
                dstfile = '{0}/data/{1}'.format(segdir,ff)
                srcfile = '{0}/{1}'.format(self.commonfilespath,ff)
                if not os.path.exists(dstfile):
                    shutil.copyfile(srcfile, dstfile)

            # copy binaries
            binlist = glob.glob('{}/*.exe'.format(self.masterbinpath))
            for bb in binlist:
                fname = bb.split('/')[-1]
                dstfile = '{0}/bin/{1}'.format(segdir,fname)
                if not os.path.exists(dstfile):
                    shutil.copyfile(bb, dstfile)
                    st = os.stat(dstfile)
                    os.chmod(dstfile,st.st_mode | stat.S_IEXEC)
                    
            # copy mod atm into subdir
            for aa in self.atmflist:
                fname = aa.split('/')[-1]
                dstfile = '{0}/atm/{1}'.format(segdir,fname)
                if not os.path.exists(dstfile):
                    shutil.copyfile(aa, dstfile)
            
            # cut transmission model to wavelength range and copy it to data/
            condwl = (self.fulltransmod['wave'] >= segdict[kk]['startwl']) & (self.fulltransmod['wave'] <= segdict[kk]['endwl'])
            transout = Table()
            transout['wave'] = self.fulltransmod['wave'][condwl]
            transout['flux'] = self.fulltransmod['flux'][condwl]
            transout.write('{0}/data/transmod.fits'.format(segdir),overwrite=True)


    def indexmasterll(self,):
        # function that takes the master LL fort.14 and fort.20, and builds an HDF5 files
        # with a master index

        # read the masterll files
        print('... Read in Master LL',flush=True)
        RK = readkurucz.ReadKurucz()
        RK.readfiles(
            f12path=self.masterf12path,
            f14path=self.masterf14path,
            f19path=self.masterf19path,
            f20path=self.masterf20path,
            f93path=self.masterf93path)
        
        # first stack f14 and f20 info
        mLL = {}
        for kk in RK.f14in.keys():
            mLL[kk] = np.append(RK.f14in[kk],RK.f20in[kk])

        # mLL['linsrc'] = np.array(
        #     [14 for _ in range(len(RK.f14in['wl']))] + 
        #     [20 for _ in range(len(RK.f20in['wl']))],dtype=int)

        # sort mLL based on wavelengths
        sort_ind = np.argsort(mLL['wl'])
        for kk in RK.f14in.keys():
            mLL[kk] = mLL[kk][sort_ind]
        
        # define index array based on length of mLL
        ind = range(len(mLL['wl']))
        
        print('... Create Master HDF5 file',flush=True)
        # create HDF5 file and fill it with info
        with h5py.File('./masterinfo/masterll.h5','w') as h5file:
            # first write a index table to simplify look-ups
            indexarr = np.vstack([ind,mLL['wl'],mLL['code']])
            h5file.create_dataset('index',data=indexarr, compression='gzip')

            # write individual line info into hdf5
            for kk in mLL.keys():
                dataset_name = f'{kk}'
                try:
                    h5file.create_dataset(dataset_name,data=mLL[kk], compression='gzip')
                except TypeError:
                    h5file.create_dataset(dataset_name,data=mLL[kk].astype('S1'), compression='gzip')

    
    def refactorll(self,segnum=0,startwl=0.0,endwl=np.inf):
        # function that takes masterll, parses it down to the specific wavelengths,
        # adds significant lines outside of wavelength range, and writes the new 
        # fort files out to segll subdirectory.
        
        # initialize AdjKurucz
        print('... Read Fort Files',flush=True)
        AK = adjkurucz.AdjKurucz(
            f12path=self.masterf12path,
            f14path=self.masterf14path,
            f19path=self.masterf19path,
            f20path=self.masterf20path,
            f93path=self.masterf93path)

        print(f'... Starting Number of Lines in fort.12: {AK.RK.f93in["nlines"]} ({len(AK.RK.f14in["wl"])})')
        print(f'... Starting Number of Lines in fort.19: {AK.RK.f93in["n19"]} ({len(AK.RK.f20in["wl"])})')
        # print('... First 10 lines')
        # print(AK.RK.f14in['wl'][:10])
        # print('... Last 10 lines')
        # print(AK.RK.f14in['wl'][-10:])


        print(f'... Adj the starting and ending wavelengths: {startwl} - {endwl}')
        AK.adj93(newdict={'wl':[startwl-0.05,endwl+0.05]})

        print(f'... Total Number of Lines in fort.12 after Ajd: {AK.RK.f93in["nlines"]} ({len(AK.RK.f14in["wl"])})')
        print(f'... Total Number of Lines in fort.19 after Ajd: {AK.RK.f93in["n19"]}  ({len(AK.RK.f20in["wl"])})')
        # print('... First 10 lines')
        # print(AK.RK.f14in['wl'][:10])
        # print('... Last 10 lines')
        # print(AK.RK.f14in['wl'][-10:])

        # write new fort files to segll directory
        print('... Write out new fort files to seg dir')
        AK.wfort(        
            f12path=f'seg_{segnum}/ff/fort.12',
            f14path=f'seg_{segnum}/ff/fort.14',
            f19path=f'seg_{segnum}/ff/fort.19',
            f20path=f'seg_{segnum}/ff/fort.20',
            f93path=f'seg_{segnum}/ff/fort.93',
            )
    
        # Now make expanded fort files for specfull synthesis
        AK = adjkurucz.AdjKurucz(
            f12path=self.masterf12path,
            f14path=self.masterf14path,
            f19path=self.masterf19path,
            f20path=self.masterf20path,
            f93path=self.masterf93path)
        AK.adj93(newdict={'wl':[startwl-5.0,endwl+5.0]})
        AK.wfort(        
            f12path=f'seg_{segnum}/ff/fort_specfull.12',
            f14path=f'seg_{segnum}/ff/fort_specfull.14',
            f19path=f'seg_{segnum}/ff/fort_specfull.19',
            f20path=f'seg_{segnum}/ff/fort_specfull.20',
            f93path=f'seg_{segnum}/ff/fort_specfull.93',
            )
        
    
    def selfitlines(self,segnum=0,startwl=0.0,endwl=np.inf):
        # function that selects which lines to be fit in segment (lineindex), assigns free parameter 
        # index for each line (linepars), and writes info files out for lineindex and linepars. 
        # lineindex also includes the translation from masterll to segll.
        
        # open masterindex HDF5 file
        mLL = h5py.File('./masterinfo/masterll.h5','r')
        # read in index array
        mindarr = mLL['index'][()]

        # parse mindarr to be generally around the wavelength range of fit region
        condwl_m = (mindarr[1,:] > startwl - 0.05) & (mindarr[1,:] < endwl + 0.05)
        mindarr = mindarr[:,condwl_m]

        print(f'... Number of Potentially matching lines in master LL: {len(mindarr[0,:])}')
                
        print('... Determining which lines need to be included in fit',flush=True)
        # temp change dir to seg_/ so that fortran is run there
        with cwd(f'seg_{segnum}/'):
            
            # # glob all atm in atm/ into list
            # atmlist_i = glob.glob('./atm/*.atm')

            for ii,atm_i in enumerate(self.atmflist):
                inputdict = {}
                inputdict['verbose']   = False
                inputdict['exedir']    = './bin/'
                inputdict['molecules'] = './data/molecules.dat'
                inputdict['continuua'] = './data/continuua.dat'
                inputdict['he1tables'] = './data/he1tables.dat'
                inputdict['spectrv_infile'] = './data/spectrv.input'

                inputdict['atmmod'] = self.specinfo[ii]['modatm']
                inputdict['rotvel'] = self.specinfo[ii]['rotvel']
                inputdict['R']      = self.specinfo[ii]['R']
                inputdict['vmac'] = self.specinfo[ii]['vmac']
                if 'isofrac' in self.specinfo[ii].keys():
                    inputdict['isofrac'] = self.specinfo[ii]['isofrac']
                else:
                    inputdict['isofrac'] = None
                
                inputdict['synspeed'] = 'slow'

                # first generate specfull spectrum for each atm input
                RS = runsynthe.Synthe(**inputdict)

                RS.setfpaths(
                    f12path='./ff/fort_specfull.12',
                    f14path='./ff/fort_specfull.14',
                    f19path='./ff/fort_specfull.19',
                    f20path='./ff/fort_specfull.20',
                    f93path='./ff/fort_specfull.93',
                    )

                print(f'---->>> specfull working on {atm_i}',flush=True)
                # set atm file path

                # run SYNTHE in seg directory
                synout_i = RS.run()

                # print(synout_i['qmu1'][:10])

                # write spectrum to seg_num/data/
                tmpspec = Table()
                tmpspec['wave'] = synout_i['wave']
                tmpspec['qmu1'] = synout_i['qmu1']
                tmpspec['qmu2'] = synout_i['qmu2']
                tmpspec['flux'] = synout_i['qmu1']/synout_i['qmu2']
                tmpspec.write(f'./data/specfull_{atm_i.split("/")[-1].replace(".atm",".fits")}',format='fits',overwrite=True)
                print(f'specfull min/max flux {min(tmpspec["flux"])}/{max(tmpspec["flux"])}',flush=True)

            

            # init list for fit lines
            code    = []
            wl      = []
            dwl     = []
            loggf   = []
            dloggf  = []
            gammar  = []
            gammas  = []
            gammaw  = []
            dgammar = []
            dgammas = []
            dgammaw = []
            resid   = []
            src     = [] # index for atm that flagged line

            # init list for all significant lines
            code_full    = []
            wl_full      = []
            dwl_full     = []
            loggf_full   = []
            dloggf_full  = []
            gammar_full  = []
            gammas_full  = []
            gammaw_full  = []
            dgammar_full = []
            dgammas_full = []
            dgammaw_full = []
            resid_full   = []
            src_full     = [] # index for atm that flagged line

            # set up initial SYNTHE runs for different atm files.
            # This is to determine which lines need to be included
            # as fit parameters.
            inputdict = {}
            inputdict['exedir']    = './bin/'
            inputdict['molecules'] = './data/molecules.dat'
            inputdict['continuua'] = './data/continuua.dat'
            inputdict['he1tables'] = './data/he1tables.dat'
            inputdict['spectrv_infile'] = './data/spectrv.input'
            inputdict['verbose'] = False
            inputdict['synspeed'] = 'fast'

            RS = runsynthe.Synthe(**inputdict)

            RS.setfpaths(
                f12path='./ff/fort.12',
                f14path='./ff/fort.14',
                f19path='./ff/fort.19',
                f20path='./ff/fort.20',
                f93path='./ff/fort.93',
                )

            # Do an inital synthesis for each atm saving resid info
            for ii,atm_i in enumerate(self.atmflist):
                starttime = datetime.now()
                print(f'---->>> Finding lines using {atm_i}',flush=True)
                # set atm file path
                RS.setatmpath(atmpath=atm_i)
                RS.setvrot(rotvel=self.specinfo[ii]['rotvel'])
                RS.setres(R=self.specinfo[ii]['R'])
                RS.setvmac(vmac=self.specinfo[ii]['vmac'])
                if 'isofrac' in self.specinfo[ii].keys():
                    isofrac = self.specinfo[ii]['isofrac']
                    RS.setisofrac(isofrac=isofrac)
                else:
                    RS.setisofrac(isofrac=None)
                    
                # run SYNTHE in seg directory
                synout_i = RS.run()

                code_fl    = synout_i['code']
                wl_fl      = synout_i['wl']
                dwl_fl     = synout_i['dwl']
                loggf_fl   = synout_i['loggf']
                dloggf_fl  = synout_i['dloggf']
                gammar_fl  = synout_i['gammar']
                gammas_fl  = synout_i['gammas']
                gammaw_fl  = synout_i['gammaw']
                dgammar_fl = synout_i['dgammar']
                dgammas_fl = synout_i['dgammas']
                dgammaw_fl = synout_i['dgammaw']
                resid_fl   = synout_i['resid']

                print(f'     ... Found {len(code_fl)} to consider.',flush=True)

                nonrepeatind = [True for _ in range(len(code_fl))]
                y = np.vstack([wl_full,loggf_full,gammar_full,gammas_full,gammaw_full]).T
                for jj,(x1,x2,x3,x4,x5) in enumerate(zip(wl_fl,loggf_fl,gammar_fl,gammas_fl,gammaw_fl)):
                    x = [x1,x2,x3,x4,x5]
                    if len((y == x).all(axis=1).nonzero()[0]) > 0:
                        nonrepeatind[jj] = False
                nonrepeatind = np.array(nonrepeatind,dtype=bool)
                
                print(f'     ... Adding {nonrepeatind.sum()}/{len(nonrepeatind)} to significant line list.',flush=True)

                # append the non-repeating to parent lists            
                code_full    = np.append(code_full,code_fl[nonrepeatind])
                wl_full      = np.append(wl_full,wl_fl[nonrepeatind])
                dwl_full     = np.append(dwl_full,dwl_fl[nonrepeatind])
                loggf_full   = np.append(loggf_full,loggf_fl[nonrepeatind])
                dloggf_full  = np.append(dloggf_full,dloggf_fl[nonrepeatind])
                gammar_full  = np.append(gammar_full,gammar_fl[nonrepeatind])
                gammas_full  = np.append(gammas_full,gammas_fl[nonrepeatind])
                gammaw_full  = np.append(gammaw_full,gammaw_fl[nonrepeatind])
                dgammar_full = np.append(dgammar_full,dgammar_fl[nonrepeatind])
                dgammas_full = np.append(dgammas_full,dgammas_fl[nonrepeatind])
                dgammaw_full = np.append(dgammaw_full,dgammaw_fl[nonrepeatind])
                resid_full   = np.append(resid_full,resid_fl[nonrepeatind])
                src_full     = np.append(src_full,ii*np.ones(len(code_fl[nonrepeatind]),dtype=int))

                # find lines where resid == 1.0
                lowlines = resid_fl > 0.9999
                print(f'     ... Found {lowlines.sum()} lines with RESID > 0.9999.',flush=True)
                
                # filter out lines less than threashold (resid -> continuum = 1.0)
                theshcond = resid_fl < (1.0 - self.threshold)

                print(f'     ... After threshold cut {theshcond.sum()}.',flush=True)
                
                # check to make sure there are lines to fit for this atm
                if theshcond.sum() == 0:
                    continue
                
                code_i    = code_fl[theshcond]
                wl_i      = wl_fl[theshcond]
                dwl_i     = dwl_fl[theshcond]
                loggf_i   = loggf_fl[theshcond]
                dloggf_i  = dloggf_fl[theshcond]
                gammar_i  = gammar_fl[theshcond]
                gammas_i  = gammas_fl[theshcond]
                gammaw_i  = gammaw_fl[theshcond]
                dgammar_i = dgammar_fl[theshcond]
                dgammas_i = dgammas_fl[theshcond]
                dgammaw_i = dgammaw_fl[theshcond]
                resid_i   = resid_fl[theshcond]
                
                # sort by wl
                sortcond = np.argsort(wl_i)
                
                code_i      = code_i[sortcond]
                wl_i      = wl_i[sortcond]
                dwl_i     = dwl_i[sortcond]
                loggf_i   = loggf_i[sortcond]
                dloggf_i  = dloggf_i[sortcond]
                gammar_i  = gammar_i[sortcond]
                gammas_i  = gammas_i[sortcond]
                gammaw_i  = gammaw_i[sortcond]
                dgammar_i = dgammar_i[sortcond]
                dgammas_i = dgammas_i[sortcond]
                dgammaw_i = dgammaw_i[sortcond]
                resid_i   = resid_i[sortcond]
                
                # filter out lines outside wavelength range
                wlcond = (wl_i > startwl) & (wl_i < endwl)

                print(f'     ... After wave cut {wlcond.sum()}.',flush=True)
                
                # check to make sure there are lines to fit for this atm
                if wlcond.sum() == 0:
                    continue
                
                code_i    = code_i[wlcond]
                wl_i      = wl_i[wlcond]
                dwl_i     = dwl_i[wlcond]
                loggf_i   = loggf_i[wlcond]
                dloggf_i  = dloggf_i[wlcond]
                gammar_i  = gammar_i[wlcond]
                gammas_i  = gammas_i[wlcond]
                gammaw_i  = gammaw_i[wlcond]
                dgammar_i = dgammar_i[wlcond]
                dgammas_i = dgammas_i[wlcond]
                dgammaw_i = dgammaw_i[wlcond]
                resid_i   = resid_i[wlcond]
                
                # # check to see if there are repeats
                # x = np.array([wl_i,loggf_i,gammar_i,gammas_i,gammaw_i])
                # y = np.array([wl,loggf,gammar,gammas,gammaw])
                # nonrepeatind = np.nonzero(np.all(~np.isin(x,y).T,axis=1))[0]

                nonrepeatind = [True for _ in range(len(code_i))]
                y = np.vstack([wl,loggf,gammar,gammas,gammaw]).T
                for jj,(x1,x2,x3,x4,x5) in enumerate(zip(wl_i,loggf_i,gammar_i,gammas_i,gammaw_i)):
                    x = [x1,x2,x3,x4,x5]
                    if len((y == x).all(axis=1).nonzero()[0]) > 0:
                        nonrepeatind[jj] = False
                nonrepeatind = np.array(nonrepeatind,dtype=bool)

                # append the non-repeating to parent lists            
                code    = np.append(code,code_i[nonrepeatind])
                wl      = np.append(wl,wl_i[nonrepeatind])
                dwl     = np.append(dwl,dwl_i[nonrepeatind])
                loggf   = np.append(loggf,loggf_i[nonrepeatind])
                dloggf  = np.append(dloggf,dloggf_i[nonrepeatind])
                gammar  = np.append(gammar,gammar_i[nonrepeatind])
                gammas  = np.append(gammas,gammas_i[nonrepeatind])
                gammaw  = np.append(gammaw,gammaw_i[nonrepeatind])
                dgammar = np.append(dgammar,dgammar_i[nonrepeatind])
                dgammas = np.append(dgammas,dgammas_i[nonrepeatind])
                dgammaw = np.append(dgammaw,dgammaw_i[nonrepeatind])
                resid   = np.append(resid,resid_i[nonrepeatind])
                src     = np.append(src,ii*np.ones(len(code_i[nonrepeatind]),dtype=int))

                print(f'     ... Registering {len(nonrepeatind)} lines -> Total: {len(wl)}',flush=True)
                print(f'     ... Finished {atm_i}: {datetime.now()-starttime}',flush=True)

            starttime = datetime.now()

            sortwl_full  = np.argsort(wl_full)
            code_full    = code_full[sortwl_full]
            wl_full      = wl_full[sortwl_full]
            dwl_full     = dwl_full[sortwl_full]
            loggf_full   = loggf_full[sortwl_full]
            dloggf_full  = dloggf_full[sortwl_full]
            gammar_full  = gammar_full[sortwl_full]
            gammas_full  = gammas_full[sortwl_full]
            gammaw_full  = gammaw_full[sortwl_full]
            dgammar_full = dgammar_full[sortwl_full]
            dgammas_full = dgammas_full[sortwl_full]
            dgammaw_full = dgammaw_full[sortwl_full]
            resid_full   = resid_full[sortwl_full]
            src_full     = src_full[sortwl_full]

            print(f'---->>> Total Number of Strong Lines {len(code_full)}')

            # initialize AdjKurucz
            AK = adjkurucz.AdjKurucz(
                f12path=f'./ff/fort.12',
                f14path=f'./ff/fort.14',
                f19path=f'./ff/fort.19',
                f20path=f'./ff/fort.20',
                f93path=f'./ff/fort.93',
                verbose=True,
                )
            # index the lines in just this segment
            # first stack f14 and f20 info
            sLL = {}
            for kk in AK.RK.f14in.keys():
                sLL[kk] = AK.RK.f14in[kk]
                # if kk in ['labelx','labelpx','other1x','other2x']:
                #     sLL[kk] = np.concatenate((AK.RK.f14in[kk]),axis=1)
                # else:
                #     sLL[kk] = np.concatenate((AK.RK.f14in[kk]),axis=0)
            sLL['index'] = np.array(range(len(sLL['wl'])))

            slindex = []
            for ii in range(len(code_full)):
                y = np.vstack([sLL['wl'],sLL['code'],sLL['gflog'],sLL['gr'],sLL['gs'],sLL['gw']]).T
                x = [wl_full[ii],code_full[ii],loggf_full[ii],gammar_full[ii],gammas_full[ii],gammaw_full[ii]]
                cond_sel = (y == x).all(axis=1)

                if cond_sel.sum() > 0:
                    for ss in sLL['index'][cond_sel]:
                        slindex.append(int(ss))
                                
                # if cond_sel.sum() == 1:
                #     slindex.append(int(sLL['index'][cond_sel]))
                # elif cond_sel.sum() > 1:
                #     for ss in sLL['index'][cond_sel]:
                #         slindex.append(ss)
                else:
                    print(f'Could not find match for this line:',flush=True)
                    print(wl_full[ii],code_full[ii],loggf_full[ii],gammaw_full[ii],flush=True)

                    # cond_fail = cond_sel
                    # print(sLL['wl'][cond_fail],sLL['code'][cond_fail],sLL['gflog'][cond_fail],sLL['gw'][cond_fail])
                    # raise IOError
                    pass

            slindex.sort()
            AK.filterll({'index':slindex})

            print(f'... Writing out rebuilt significant line fort files with {len(slindex)} lines',flush=True)
            AK.wfort(        
                f12path=f'./ff/fort_sl.12',
                f14path=f'./ff/fort_sl.14',
                f19path=f'./ff/fort_sl.19',
                f20path=f'./ff/fort_sl.20',
                f93path=f'./ff/fort_sl.93',
                )
            
            # read the segll files
            print('... Read SL fort files back in',flush=True)
            RK = readkurucz.ReadKurucz(verbose=True)
            RK.readfiles(
                f12path=f'./ff/fort_sl.12',
                f14path=f'./ff/fort_sl.14',
                f19path=f'./ff/fort_sl.19',
                f20path=f'./ff/fort_sl.20',
                f93path=f'./ff/fort_sl.93',
            )

            print(f'... Number of lines in fort files {RK.f93in["nlines"]+RK.f93in["n19"]}')

            # index the lines in just this segment
            # first stack f14 and f20 info
            sLL = {}
            for kk in RK.f14in.keys():
                if kk in ['labelx','labelpx','other1x','other2x']:
                    sLL[kk] = np.concatenate((RK.f14in[kk],RK.f20in[kk]),axis=1)
                else:
                    sLL[kk] = np.concatenate((RK.f14in[kk],RK.f20in[kk]),axis=0)

            print(f'... Number of lines in sLL dict {len(sLL["wl"])}')

            # sLL['linsrc'] = np.array(
            #     [14 for _ in range(len(RK.f14in['wl']))] + 
            #     [20 for _ in range(len(RK.f20in['wl']))],dtype=int)

            # sort sLL based on wavelengths
            # sort_ind = np.argsort(sLL['wl'])
            # for kk in RK.f14in.keys():
            #     sLL[kk] = sLL[kk][sort_ind]
                    
            # define index array based on length of sLL
            sLL['index'] = np.array(range(len(sLL['wl'])))

            # trim sLL to wavelength range
            condwl = (sLL['wl'] > startwl) & (sLL['wl'] < endwl)
            sLL_t = {}
            for kk in RK.f14in.keys():
                try:
                    if kk in ['labelx','labelpx','other1x','other2x']:
                        sLL_t[kk] = sLL[kk][...,condwl]
                    else:
                        sLL_t[kk] = sLL[kk][condwl,...]
                except:
                    print('PROBLEM',kk,len(sLL[kk]),sLL[kk].shape)
                    raise
            sLL_t['index'] = sLL['index'][condwl]
            
            print(f'... Number of lines in sLL dict after wl cut {len(sLL_t["wl"])}')
                    
            print('... Finding Master LL Matches',flush=True)
            # first find seg index and match with master index
            # Try using the index array first, maybe we'll get lucky and there is only one match
            starttime = datetime.now()
            lineindexarr = []
            for ii in sLL_t['index']:
                cond = (mindarr[1,:] == sLL['wl'][ii]) & (mindarr[2,:] == sLL['code'][ii])

                if cond.sum() == 1:
                    # found only one match, write line index
                    lineindexarr.append(int(mindarr[0,cond][0]))
                elif cond.sum() > 1:
                    # found more than one match, must sort out which one it is
                    # All columns must have matching values
                    potentiallines = mindarr[0,cond]
                    foundmat = False
                    for pp in potentiallines:
                        # mLL_i = mLL[f'{int(pp)}']
                        mat = True
                        for kk in sLL.keys():
                            if kk in (
                                ['index','linesrc','ref','auto','ixfixfp',
                                'labelp','label','labelx','labelpx',
                                'ishift','ishiftp','other1x','other2x']):
                                continue
                            # mat *= sLL[kk][ii] == mLL_i[kk]
                            try:
                                mat *= sLL[kk][ii] == mLL[kk][int(pp)]
                            except:
                                print(kk,int(pp))
                                raise
                        if mat:
                            lineindexarr.append(int(pp))
                            foundmat = True
                            break
                    if foundmat == False:
                        for kk in sLL.keys():
                            testpar = sLL[kk][ii]
                            matchedlist = np.array([mLL[kk][int(pp)] for pp in potentiallines])
                            condtest = matchedlist == testpar
                            print('A',kk,testpar,matchedlist,'->',matchedlist[condtest])

                        print(f'ISSUE WITH LINE {ii}, LOOKED AT POTENTIAL LINES AND DID NOT FIND MATCH',flush=True)
                        print(f'THIS SHOULD NOT HAPPEN',flush=True)
                        raise IOError
                        
                else:
                    print(sLL['wl'][ii])
                    print(mindarr[1,:].min())
                    print(mindarr[1,:].max())
                    print(f'ISSUE WITH LINE {ii}, COULD NOT FIND MATCH IN MASTERLL',flush=True)
                    print(f'THIS SHOULD NOT HAPPEN',flush=True)
                    raise IOError
            
            # close the masterLL HDF5 file
            mLL.close()

            print(f'... Master LL matched finished: {datetime.now()-starttime}',flush=True)
            
            # write in the master ll index for future use
            sLL_t['masterind'] = np.array(lineindexarr)

            starttime = datetime.now()                
            # write lineindex arrays to file
            with open(f'./lineinfo/lineindex.txt','w') as lif:
                lif.write('segind masterind wl code\n')
                for x,y,z,w in zip(sLL_t['index'],lineindexarr,sLL_t['wl'],sLL_t['code']):
                    lif.write(f'{x} {y} {z:.4f} {w:.2f}\n')
            print(f'... Fnished writing master LL match file: {datetime.now()-starttime}',flush=True)


            print(f'... setting up free parameter arrays',flush=True)
            # sort arrays by wl
            sortwl  = np.argsort(wl)
            code    = code[sortwl]
            wl      = wl[sortwl]
            dwl     = dwl[sortwl]
            loggf   = loggf[sortwl]
            dloggf  = dloggf[sortwl]
            gammar  = gammar[sortwl]
            gammas  = gammas[sortwl]
            gammaw  = gammaw[sortwl]
            dgammar = dgammar[sortwl]
            dgammas = dgammas[sortwl]
            dgammaw = dgammaw[sortwl]
            resid   = resid[sortwl]
            src     = src[sortwl]


            # int free par index for each line in fit line list
            wlind = []
            gfind = []
            gwind = []
            
            # init seg and master index
            segind_f    = []
            masterind_f = []

            # init final par columns
            code_f   = []
            wl_f     = []
            loggf_f  = []
            gammar_f = []
            gammas_f = []
            gammaw_f = []
            resid_f  = []
            src_f    = []


            parind_j = 0
            # wlind_j = 0
            # gfind_j = 1
            # gwind_j = 2
            
            for ii in range(len(code)):

                wlind_j,gfind_j,gwind_j = [parind_j,parind_j+1,parind_j+2]
                
                # figure out the indices
                wlind_i = wlind_j
                gfind_i = gfind_j
                if code[ii] > 99.0:
                    gwind_i = -1
                else:
                    gwind_i = gwind_j
                
                # figure out segind
                # cond_sel = (
                #     (sLL_t['wl'] == wl[ii]) & 
                #     (sLL_t['code'] == code[ii]) & 
                #     (sLL_t['gflog'] == loggf[ii]) &
                #     (sLL_t['gr'] == gammar[ii]) &
                #     (sLL_t['gs'] == gammas[ii]) &
                #     (sLL_t['gw'] == gammaw[ii]) 
                #     )

                y = np.vstack([sLL_t['wl'],sLL_t['code'],sLL_t['gflog'],sLL_t['gr'],sLL_t['gs'],sLL_t['gw']]).T
                x = [wl[ii],code[ii],loggf[ii],gammar[ii],gammas[ii],gammaw[ii]]
                cond_sel = (y == x).all(axis=1)

                # check to make sure if found the line
                if cond_sel.sum() == 1:
                    pass
                elif cond_sel.sum() > 1:
                    # just set the first line to be free, the other lines will likely 
                    # be added as HF/ISO lines
                    cond_sel = np.array(cond_sel).cumsum() == 1
                else:
                    print(f'Could not find  match to this line:',flush=True)
                    print(wl[ii],code[ii],loggf[ii],gammaw[ii],flush=True)

                    cond_fail = sLL_t['wl'] == wl[ii]
                    print(sLL_t['wl'][cond_fail],sLL_t['code'][cond_fail],sLL_t['gflog'][cond_fail],sLL_t['gw'][cond_fail])
                    raise IOError
                    
                segind_i    = sLL_t['index'][cond_sel][0]
                masterind_i = sLL_t['masterind'][cond_sel][0]
                
                # check to see if line already is in fit line list
                if segind_i in segind_f:
                    continue
                
                # book the initial line 
                segind_f.append(segind_i)
                masterind_f.append(masterind_i)
                wlind.append(wlind_i)
                gfind.append(gfind_i)
                gwind.append(gwind_i)
                code_f.append(code[ii])
                wl_f.append(wl[ii])
                loggf_f.append(loggf[ii])
                gammar_f.append(gammar[ii])
                gammas_f.append(gammas[ii])
                gammaw_f.append(gammaw[ii])
                resid_f.append(resid[ii])
                try:
                    src_f.append(src[ii])
                except:
                    print(src)
                    raise
                                
                # now look for connected line (molecular or HF/ISO)
                sLL_i = {}
                sLL_m = {}
                
                for kk in sLL_t.keys():
                    try:
                        if kk in ['labelx','labelpx','other1x','other2x']:
                            sLL_i[kk] = sLL_t[kk][...,cond_sel]
                        else:
                            sLL_i[kk] = sLL_t[kk][cond_sel,...]
                    except:
                        print('PROBLEM I',kk,len(cond_sel),sLL_t[kk].shape)
                        raise

                    # try:
                    #     sLL_i[kk] = sLL_t[kk][cond_sel]
                    # except:
                    #     print('I')
                    #     print(kk)
                    #     print(sLL_t[kk])
                    #     print(len(sLL_t[kk]))
                    #     print(len(cond_sel))
                    #     raise

                    try:
                        if kk in ['labelx','labelpx','other1x','other2x']:
                            sLL_m[kk] = sLL_t[kk][...,~cond_sel]
                        else:
                            sLL_m[kk] = sLL_t[kk][~cond_sel,...]
                    except:
                        print('PROBLEM M',kk,len(~cond_sel),sLL_t[kk].shape)
                        raise

                    # try:
                    #     sLL_m[kk] = sLL_t[kk][~cond_sel]
                    # except:
                    #     print('M')
                    #     print(kk)
                    #     print(sLL_t[kk])
                    #     raise

                # sLL_i = {sLL[kk][cond] for kk in sLL.keys()}
                # sLL_m = {sLL[kk][~cond] for kk in sLL.keys()}
                
                if code[ii] < 99.0:
                    # check for HF/ISO
                    HFISOcond = (
                        (sLL_m['code']  == sLL_i['code']) & 
                        (sLL_m['gflog'] == sLL_i['gflog']) & 
                        (sLL_m['e']  == sLL_i['e']) & 
                        (sLL_m['ep']  == sLL_i['ep'])
                    )

                    if HFISOcond.sum() > 0:
                        sLL_mm = {}
                        for kk in sLL_m.keys():
                            try:
                                if kk in ['labelx','labelpx','other1x','other2x']:
                                    sLL_mm[kk] = sLL_m[kk][...,HFISOcond]
                                else:
                                    sLL_mm[kk] = sLL_m[kk][HFISOcond,...]
                            except:
                                print('PROBLEM MM',kk,len(HFISOcond),sLL_m[kk].shape)
                                raise
                        # sLL_mm = {kk:sLL_m[kk][HFISOcond] for kk in sLL_m.keys()}
                        # sLL_mm = sLL_m[HFISOcond]
                        for jj in range(len(sLL_mm['index'])):
                            segind_f.append(sLL_mm['index'][jj])
                            masterind_f.append(sLL_mm['masterind'][jj])
                            wlind.append(wlind_i)
                            gfind.append(gfind_i)
                            gwind.append(gwind_i)
                            code_f.append(sLL_mm['code'][jj])
                            wl_f.append(sLL_mm['wl'][jj])
                            loggf_f.append(sLL_mm['gflog'][jj])
                            gammar_f.append(sLL_mm['gr'][jj])
                            gammas_f.append(sLL_mm['gs'][jj])
                            gammaw_f.append(sLL_mm['gw'][jj])
                            resid_f.append(0.0)
                            src_f.append(src[ii])
                else:
                    # check for other molecular lines in the same transition
                    MOLcond = (
                        (sLL_m['code']  == sLL_i['code']) & 
                        (sLL_m['iso1']  == sLL_i['iso1']) & 
                        (sLL_m['iso2']  == sLL_i['iso2']) & 
                        (sLL_m['x1']    == sLL_i['x1']) & 
                        (sLL_m['x2']    == sLL_i['x2'])  
                    )

                    if MOLcond.sum() > 0:
                        sLL_mm = {}
                        for kk in sLL_m.keys():
                            try:
                                if kk in ['labelx','labelpx','other1x','other2x']:
                                    sLL_mm[kk] = sLL_m[kk][...,MOLcond]
                                else:
                                    sLL_mm[kk] = sLL_m[kk][MOLcond,...]
                            except:
                                print('PROBLEM MM',kk,len(MOLcond),sLL_m[kk].shape)
                                raise

                        # sLL_mm = {kk:sLL_m[kk][MOLcond] for kk in sLL_m.keys()}
                        # sLL_mm = sLL_m[MOLcond]
                        for jj in range(len(sLL_mm['index'])):
                            segind_f.append(sLL_mm['index'][jj])
                            masterind_f.append(sLL_mm['masterind'][jj])
                            wlind.append(wlind_i)
                            gfind.append(gfind_i)
                            gwind.append(gwind_i)
                            code_f.append(sLL_mm['code'][jj])
                            wl_f.append(sLL_mm['wl'][jj])
                            loggf_f.append(sLL_mm['gflog'][jj])
                            gammar_f.append(sLL_mm['gr'][jj])
                            gammas_f.append(sLL_mm['gs'][jj])
                            gammaw_f.append(sLL_mm['gw'][jj])
                            resid_f.append(1.0)
                            src_f.append(src[ii])
                    
                parind_j = max([wlind_i,gfind_i,gwind_i])+1
                
                # wlind_j += 1 + wlind_j
                # gfind_j += 1 + gfind_j
                # if gwind_i != -1:
                #     gwind_j += 1
            
            print(f'... Writing line fit pars file',flush=True)
            # write fit pars file out
            with open('./lineinfo/linefitpars.txt','w') as lfp:
                lfp.write('segind masterind wlind gfind gwind code wl loggf gammar gammas gammaw resid src \n')
                for ii in range(len(segind_f)):
                    lfp.write(f'{segind_f[ii]} {masterind_f[ii]} {wlind[ii]} {gfind[ii]} {gwind[ii]} {code_f[ii]:.2f} {wl_f[ii]:.4f} {loggf_f[ii]:.3f} {gammar_f[ii]:.2f} {gammas_f[ii]:.2f} {gammaw_f[ii]:.2f} {resid_f[ii]:.4f} {src_f[ii]} \n')    
    
            print(f'... Making weak line spectrum',flush=True)
            # run SYNTHE with just the significant lines
            # this model is used to divide by the full specturm and 
            # the result is the weak line spectrum used in fitting
            RS = runsynthe.Synthe(
                exedir='./bin/',
                molecules='./data/molecules.dat',
                continuua='./data/continuua.dat',
                he1tables='./data/he1tables.dat',
                spectrv_infile='./data/spectrv.input',                
                verbose=False,
                synspeed='fast',
                )
            RS.setfpaths(
                f12path=f'./ff/fort_sl.12',
                f14path=f'./ff/fort_sl.14',
                f19path=f'./ff/fort_sl.19',
                f20path=f'./ff/fort_sl.20',
                f93path=f'./ff/fort_sl.93',
                )
            for ii,atm_i in enumerate(self.atmflist):
                starttime = datetime.now()
                print(f'... working on {atm_i}',flush=True)
                # set atm file path
                RS.setatmpath(atmpath=atm_i)
                RS.setvrot(rotvel=self.specinfo[ii]['rotvel'])
                RS.setres(R=self.specinfo[ii]['R'])
                RS.setvmac(vmac=self.specinfo[ii]['vmac'])
                if 'isofrac' in self.specinfo[ii].keys():
                    isofrac = self.specinfo[ii]['isofrac']
                    RS.setisofrac(isofrac=isofrac)
                else:
                    RS.setisofrac(isofrac=None)

                # run SYNTHE in seg directory
                synout_i = RS.run()

                # write spectrum to seg_num/data/
                tmpspec = Table()
                tmpspec['wave'] = synout_i['wave']
                tmpspec['qmu1'] = synout_i['qmu1']
                tmpspec['qmu2'] = synout_i['qmu2']
                tmpspec['flux'] = synout_i['qmu1']/synout_i['qmu2']
                tmpspec.write(f'./data/specSL_{atm_i.split("/")[-1].replace(".atm",".fits")}',format='fits',overwrite=True)

                # now read in specfull and specSL, divide them, and write out specWL
                specfull = Table.read(f'./data/specfull_{atm_i.split("/")[-1].replace(".atm",".fits")}',format='fits')
                specSL   = Table.read(f'./data/specSL_{atm_i.split("/")[-1].replace(".atm",".fits")}',format='fits')
                
                specfull_flux = np.interp(specSL['wave'],specfull['wave'],specfull['flux'])
                
                tmpspec = Table()
                tmpspec['wave'] = specSL['wave']
                tmpspec['flux'] = specfull_flux/specSL['flux']
                # set values > 1.0 to be == 1.0
                condfl = tmpspec['flux'] > 1.0
                tmpspec['flux'][condfl] = 1.0
                tmpspec.write(f'./data/specWL_{atm_i.split("/")[-1].replace(".atm",".fits")}',format='fits',overwrite=True)

    
    def defseg(self,):
        # function that defines segll given the masterll based on line density
        pass
    
    

    def __str__(self):
        return 'RunPrep@{:#x}: {} {}'.format(id(self), self.args, self.kwargs)