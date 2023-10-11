import numpy as np
import os
import shutil
import glob
import h5py

from ..utils import readkurucz, runsynthe, adjkurucz

class RunPrep(object):
    def __init__(self, *args, **kwargs):
        super(RunPrep, self).__init__()
        self.args = args
        self.kwargs = kwargs

        # user input dir where master ll fort files are kept
        self.masterllfdir = self.args[0]

        # define line threshold for including as a fitted line, default 1% depth
        self.threshold = self.kwargs.get('threshold',0.01)

        # define path for files that are the same for all SYNTHE
        self.commonfilespath = self.kwargs.get('commonfiles',None)

        # define the list of atm files to generate fitted lines
        # the number of atm files here define how many different stars
        # to consider.
        self.atmflist = self.kwargs.get('atmlist',['./data/atmmod_sol.dat'])

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

        # define path to SYNTHE exe binaries
        self.masterbinpath = self.kwargs.get('masterbin',None)

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
            segdir = 'seg_{segnum}'.format(**segdict)
            if not os.path.exists(segdir):
                os.makedirs(segdir)

            # write wl info into text file in seg dir
            with open(f'{segdir}/wlinfo.txt','w') as wlf:
                wlf.write('segnum={segnum},startwl={startwl:.4f},endwl={endwl:.4f}'.format(**segdict))

            # build seg subdir
            for ff in ['data','samples','lineinfo','mod','bin','ff']:
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

            # copy mod atm into subdir
            for aa in self.atmflist:
                fname = aa.split('/')[-1]
                dstfile = '{0}/mod/{1}'.format(segdir,fname)
                if not os.path.exists(dstfile):
                    shutil.copyfile(aa, dstfile)

    def indexmasterll(self,):
        # function that takes the master LL fort.14 and fort.20, and builds an HDF5 files
        # with a master index

        # read the masterll files
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
            mLL[kk] = np.hstack([RK.f14in[kk],RK.f20in[kk]])

        mLL['linsrc'] = np.array(
            [14 for _ in range(len(RK.f14in['wl']))] + 
            [20 for _ in range(len(RK.f20in['wl']))],dtype=int)

        # sort mLL based on wavelengths
        sort_ind = np.argsort(mLL['wl'])
        for kk in RK.f14in.keys():
            mLL[kk] = mLL[kk][sort_ind]
        
        # define index array based on length of mLL
        ind = range(len(mLL['wl']))
        
        # create HDF5 file and fill it with info
        with h5py.File('./masterinfo/masterll.h5','w') as h5file:
            # first write a index table to simplify look-ups
            indexarr = np.vstack([ind,mLL['wl'],mLL['code'],mLL['linsrc']])
            h5file.create_dataset('index',data=indexarr)

            # write individual line info into hdf5
            for ii in ind:
                for kk in mLL.keys():
                    dataset_name = f'{ii}/{kk}'
                    h5file.create_dataset(dataset_name,data=mLL[kk][ii])

    def readseg(self,segfile):
        # function that reads in list of segments, the ascii file format is:
        # segnum, starting wl, ending wl
        
        segfile = np.loadtxt(segfile,
                             dtype={'names':('segnum','start_wl','end_wl'),
                                    'formats':(int,float,float)})
        return segfile
    
    def refactorll(self,segnum=0,startwl=0.0,endwl=np.inf):
        # function that takes masterll, parses it down to the specific wavelengths,
        # adds significant lines outside of wavelength range, and writes the new 
        # fort files out to segll subdirectory.
        
        # initialize AdjKurucz
        AK = adjkurucz.AdjKurucz()
        AK.readfiles(
            f12path=self.masterf12path,
            f14path=self.masterf14path,
            f19path=self.masterf19path,
            f20path=self.masterf20path,
            f93path=self.masterf93path)

        AK.adj93(newdict={'wl':[startwl,endwl]})
    
        # write new fort files to segll directory
        AK.wfort(        
            f12path='seg_{segnum}/ff/fort.12'.format(segnum),
            f14path='seg_{segnum}/ff/fort.14'.format(segnum),
            f19path='seg_{segnum}/ff/fort.19'.format(segnum),
            f20path='seg_{segnum}/ff/fort.20'.format(segnum),
            f93path='seg_{segnum}/ff/fort.93'.format(segnum),
            )
    
    def selfitlines(self,segnum=0):
        # function that selects which lines to be fit in segment (lineindex), assigns free parameter 
        # index for each line (linepars), and writes info files out for lineindex and linepars. 
        # lineindex also includes the translation from masterll to segll.
        
        # open masterindex HDF5 file
        mLL = h5py.File('./masterinfo/masterll.h5','r')
        # read in index array
        mindarr = mLL['index'][()]

        # read the segll files
        RK = readkurucz.ReadKurucz()
        RK.readfiles(
            f12path='seg_{segnum}/ff/fort.12'.format(segnum),
            f14path='seg_{segnum}/ff/fort.14'.format(segnum),
            f19path='seg_{segnum}/ff/fort.19'.format(segnum),
            f20path='seg_{segnum}/ff/fort.20'.format(segnum),
            f93path='seg_{segnum}/ff/fort.93'.format(segnum),
        )

        # index the lines in just this segment
        # first stack f14 and f20 info
        sLL = {}
        for kk in RK.f14in.keys():
            sLL[kk] = np.hstack([RK.f14in[kk],RK.f20in[kk]])

        sLL['linsrc'] = np.array(
            [14 for _ in range(len(RK.f14in['wl']))] + 
            [20 for _ in range(len(RK.f20in['wl']))],dtype=int)

        # sort sLL based on wavelengths
        sort_ind = np.argsort(sLL['wl'])
        for kk in RK.f14in.keys():
            sLL[kk] = sLL[kk][sort_ind]
        
        # define index array based on length of sLL
        sLL['index'] = range(len(sLL['wl']))
        
        # first find seg index and match with master index
        # Try using the index array first, maybe we'll get lucky and there is only one match
        lineindexarr = []
        for ii in sLL['index']:
            cond = (sLL['wl'][ii] == mindarr[1,:]) & (sLL['code'][ii] == mindarr[2,:]) & (sLL['linsrc'][ii] == mindarr[3,:])

            if cond.sum() == 1:
                # found only one match, write line index
                lineindexarr.append(mindarr[0,cond])
            else:
                # found more than one match, must sort out which one it is
                # All columns must have matching values
                potentiallines = mindarr[0,cond]
                for pp in potentiallines:
                    mLL_i = mLL[pp]
                    mat = True
                    for kk in sLL.keys():
                        mat *= sLL[kk][ii] == mLL_i[kk]
                    if mat:
                        lineindexarr.append(pp)
                        break

        # write in the master ll index for future use
        sLL['masterind'] = lineindexarr
        
        # write lineindex arrays to file
        with open('seg_{segnum}/lineinfo/lineindex.txt','w') as lif:
            lif.write('segind masterind\n')
            for x,y in zip(sLL['index'],lineindexarr):
                lif.write(f'{x} {y}\n')
    
    
        # set up initial SYNTHE runs for different atm files.
        # This is to determine which lines need to be included
        # as fit parameters.
        RS = runsynthe.Synthe(
            exedir=self.masterbinpath,
            molecules='masterinfo/data/molecules.dat',
            continuua='masterinfo/data/continuua.dat',
            he1tables='masterinfo/data/he1tables.dat',
            spectrv_infile='masterino/data/spectrv.input',                
            verbose=False,
            )
        RS.setfpaths(
            f12path='seg_{segnum}/ff/fort.12'.format(segnum),
            f14path='seg_{segnum}/ff/fort.14'.format(segnum),
            f19path='seg_{segnum}/ff/fort.19'.format(segnum),
            f20path='seg_{segnum}/ff/fort.20'.format(segnum),
            f93path='seg_{segnum}/ff/fort.93'.format(segnum),
            )

        code    = np.array([],dtype=float)
        wl      = np.array([],dtype=float)
        dwl     = np.array([],dtype=float)
        loggf   = np.array([],dtype=float)
        dloggf  = np.array([],dtype=float)
        gammar  = np.array([],dtype=float)
        gammas  = np.array([],dtype=float)
        gammaw  = np.array([],dtype=float)
        dgammar = np.array([],dtype=float)
        dgammas = np.array([],dtype=float)
        dgammaw = np.array([],dtype=float)
        resid   = np.array([],dtype=float)
        src     = np.array([],dtype=int) # index for atm that flagged line
        
        
        # Do an inital synthesis for each atm saving resid info
        for ii,atm_i in enumerate(self.atmflist):
            # set atm file path
            RS.setatmpath(atmmod=atm_i)
            # run SYNTHE
            synout_i = RS.run()

            # filter out lines less than threashold (resid -> continuum = 1.0)
            theshcond = synout_i['resid'] < (1.0 - self.threshold)
            
            # check to make sure there are lines to fit for this atm
            if theshcond.sum() == 0:
                continue
            
            code_i    = synout_i['code'][theshcond]
            wl_i      = synout_i['wl'][theshcond]
            dwl_i     = synout_i['dwl'][theshcond]
            loggf_i   = synout_i['loggf'][theshcond]
            dloggf_i  = synout_i['dloggf'][theshcond]
            gammar_i  = synout_i['gammar'][theshcond]
            gammas_i  = synout_i['gammas'][theshcond]
            gammaw_i  = synout_i['gammaw'][theshcond]
            dgammar_i = synout_i['dgammar'][theshcond]
            dgammas_i = synout_i['dgammas'][theshcond]
            dgammaw_i = synout_i['dgammaw'][theshcond]
            resid_i   = synout_i['resid'][theshcond]
            
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
            
            # check to see if there are repeats
            x = np.array([wl_i,loggf_i,gammar_i,gammas_i,gammaw_i])
            y = np.array([wl,loggf,gammar,gammas,gammaw])
            nonrepeatind = np.nonzero(np.all(~np.isin(y,x).T,axis=1))[0]

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
            src     = np.append(src,[ii])

        # init seg and master index
        segind    = []
        masterind = []

        # int free par index for each line in fit line list
        wlind = []
        gfind = []
        gwind = []
        
        # init final par columns
        code_f   = []
        wl_f     = []
        loggf_f  = []
        gammar_f = []
        gammas_f = []
        gammaw_f = []
        resid_f  = []
        src_f    = []

        wlind_j = 0
        gfind_j = 0
        gwind_j = 0
        
        for ii in range(len(code)):
            
            # figure out the indices
            wlind_i = wlind_j
            gfind_i = gfind_j
            if code[ii] > 99.0:
                gwind_i = -1
            else:
                gwind_i = gwind_j
                        
            # figure out segind
            cond = (
                (sLL['wl'] == wl[ii]) & 
                (sLL['code'] == code[ii]) & 
                (sLL['gflog'] == loggf[ii]) &
                (sLL['gammar'] == gammar[ii]) &
                (sLL['gammas'] == gammas[ii]) &
                (sLL['gammaw'] == gammaw[ii]) 
                )

            # check to make sure if found the line
            assert cond.sum() == 1
            segind_i    = sLL['index'][cond]
            masterind_i = sLL['masterind'][cond]
            
            # book the initial line 
            segind.append(segind_i)
            masterind.append(masterind_i)
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
            src_f.append(src[ii])
            
            # now look for connected line (molecular or HF/ISO)
            sLL_i = sLL[cond]
            sLL_m = sLL[~cond]
            
            if code[ii] < 99.0:
                # check for HF/ISO
                HFISOcond = (
                    (sLL_m['code']  == sLL_i['code']) & 
                    (sLL_m['gflog'] == sLL_i['gflog']) & 
                    (sLL_m['e']  == sLL_i['e']) & 
                    (sLL_m['ep']  == sLL_i['ep'])
                )

                if HFISOcond.sum() > 0:
                    sLL_mm = sLL_m[HFISOcond]
                    for jj in range(len(sLL_mm)):
                        segind.append(sLL_mm['index'][jj])
                        masterind.append(sLL_mm['masterind'][jj])
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
                    sLL_mm = sLL_m[MOLcond]
                    for jj in range(len(sLL_mm)):
                        segind.append(sLL_mm['index'][jj])
                        masterind.append(sLL_mm['masterind'][jj])
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
                
            
            wlind_j += 1
            gfind_j += 1
            if gwind_i != -1:
                gwind_j += 1
        
        # write fit pars file out
        with open('seg_{segnum}/lineinfo/linefitpars.txt','w') as lfp:
            lfp.write('segind masterind wlind gfind gwind code wl loggf gammar gammas gammaw resid src \n')
            for ii in range(len(segind)):
                lfp.write(f'{segind[ii]} {masterind[ii]} {wlind[ii]} {gfind[ii]} {gwind[ii]} {code_f[ii]} {wl_f[ii]} {loggf_f[ii]} {gammar_f[ii]} {gammas_f[ii]} {gammaw_f[ii]} {resid_f[ii]} {src_f[ii]} \n')    
    
    def defseg(self,):
        # function that defines segll given the masterll based on line density
        pass
    
    

    def __str__(self):
        return 'RunPrep@{:#x}: {} {}'.format(id(self), self.args, self.kwargs)