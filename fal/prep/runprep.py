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

        # set the masterll file paths
        self.masterf12path=self.masterllfdir+'/fort.12'
        self.masterf14path=self.masterllfdir+'/fort.14'
        self.masterf19path=self.masterllfdir+'/fort.19'
        self.masterf20path=self.masterllfdir+'/fort.20'
        self.masterf93path=self.masterllfdir+'/fort.93'
        
        # define the list of atm files to generate fitted lines
        # the number of atm files here define how many different stars
        # to consider.
        self.atmflist = self.kwargs.get('atmlist',['./data/atmmod_sol.dat'])

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
        
        # define line threshold for including as a fitted line, default 1% depth
        self.threshold = self.kwargs.get('threshold',0.01)

        # define path for files that are the same for all SYNTHE
        self.commonfilespath = self.kwargs.get('commonfiles',None)

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
        # write lineindex arrays to file
        with open('seg_{segnum}/lineinfo/lineindex.txt','w') as lif:
            lif.write('segind masterind\n')
            for x,y in zip(sLL['index'],lineindexarr):
                lif.write(f'{x} {y}\n')
    
    
    
    # The following functions will eventually be used to determine seg info on the fly
    def genmastermod(self,):
        # function that runs SYNTHE and returns all of the significant lines with depths
        pass    
    def defseg(self,):
        # function that defines segll given the masterll based on line density
        pass
    
    

    def __str__(self):
        return 'RunPrep@{:#x}: {} {}'.format(id(self), self.args, self.kwargs)