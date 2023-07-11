

class Prior(object):
    def __init__(self, *args, **kwargs):
        super(Prior, self).__init__()

        self.args = args
        self.kwargs = kwargs

    def run(self,pars):
        return 0.0
