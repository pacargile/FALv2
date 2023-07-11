import numpy as np
from ..utils import readkurucz

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

        # read the masterll files
        self.RK = readkurucz.ReadKurucz()
        self.RK.readfiles(
            f12path=self.masterf12path,
            f14path=self.masterf14path,
            f19path=self.masterf19path,
            f20path=self.masterf20path,
            f93path=self.masterf93path)
        
        # define the list of atm files to generate fitted lines
        # the number of atm files here define how many different stars
        # to consider.
        self.atmflist = self.kwargs.get('atmlist',['./data/atmmod_sol.dat'])

        # define line threshold for including as a fitted line, default 1% depth
        self.threshold = self.kwargs.get('threshold',0.01)


    def __str__(self):
        return 'RunPrep@{:#x}: {} {}'.format(id(self), self.args, self.kwargs)