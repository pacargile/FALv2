import numpy as np

def setpars(fortll,lineindex=None):
    
    # map through lineindex
    # - identify if the line is atom or molec
    # - if atom:
    #   - see if line has hyperfine
    #   - see if line has any iotopic splits
    # - if molec
    #   - look for any other lines with from same transition
    
    fitlineindex = []
    fitlinepind  = []

    for int_ind in lineindex:
        # check to see if int_ind is in fitlineindex
        # meaning it has already be assigned
        if int_ind in fitlineindex:
            continue
    
    
    return (fitlineindex,fitlinepind)
