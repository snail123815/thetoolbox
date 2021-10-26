############################################
# by Chao DU
# Institute of Biology Leiden
# Leiden University
# c.du@biology.leidenuniv.nl
# durand[dot]dc[at]hot[no space]mail.com
############################################
# Designed in mind that all types of data can be calculated
# resulting the same hash across platform.
# Should be safe with nested dict but no guarantee
################################################################

import json
from collections import OrderedDict
from hashlib import md5

# Optional imports
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression as PLS


def calHash(*args) -> str:
    """Produce 6 digit string with unlimited number of arguments passed in
    Designed in mind that all types of data can be calculated
    resulting the same hash across platform.
    Should be safe with nested dict but no guarantee
    """
    def orderDict(di):
        try:
            d = di.copy()
        except:
            d = di
        if isinstance(d, dict):
            od = OrderedDict(sorted(d.items()))
            for k, v in od.items():
                v = orderDict(v)
                od[str(k)] = v
            d = od
        elif isinstance(d, (list, tuple, set)):
            d = [orderDict(el) for el in d]
        else:
            d = str(d).strip()
        return d

    def hashDict(di):
        od = orderDict(di)
        ha = md5(
            json.dumps(
                od,
                sort_keys=True,
                ensure_ascii=True,
                default=str
            ).encode()
        ).digest()
        return ha

    hastr = ''.encode()
    for arg in args:
        if isinstance(arg, str):
            hastr += arg.encode()
        elif isinstance(arg, set):
            hastr += str(sorted(list())).encode()
        elif isinstance(arg, dict):
            hastr += hashDict(arg)
        elif isinstance(arg, pd.core.frame.DataFrame) or isinstance(arg, pd.core.series.Series):
            hastr += md5(arg.to_json().encode()).digest()
        elif isinstance(arg, PCA):
            hastr += arg.components_.tobytes()
        elif isinstance(arg, PLS):
            hastr += arg.x_loadings_.tobytes()
        else:
            hastr += str(arg).encode()
    return md5(hastr).hexdigest()[:6]


# TEST
if __name__=="__main__":
    import numpy as np
    a = pd.DataFrame(dict(a=[1,2],b=[3,4]), index=[1,2])
    b = dict(x='iv', y='83', z=dict(i=0, u=5, k=dict(kj='kj', hs='qf')))
    c = 'avb'
    d = 1
    e = PLS()
    e.fit(np.linspace(1,100,50).reshape(-1,2).T, [[1],[1.5]])
    f = PCA()
    f.fit(np.linspace(4,5,100).reshape(-1,5))
    ha = calHash(a,b,c,d,e,f) 
    correctHa = 'b41ecd'
    if ha == correctHa:
        print("Test OK")
    else:
        print(f"Test fail, previous ha is {correctHa}, current ha is {ha}")
        print(f'a now {calHash(a)} - previous 177438')
        print(f'b now {calHash(b)} - previous 9536de')
        print(f'c now {calHash(c)} - previous 539183')
        print(f'd now {calHash(d)} - previous c4ca42')
        print(f'e now {calHash(e)} - previous 8889c7')
        print(f'f now {calHash(f)} - previous 5473de')
