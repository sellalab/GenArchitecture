from collections import Sequence
from random import getrandbits
from datetime import datetime
import os

# easy-to-use writer for various statistics
# statistics are not defined in advance;
# usage: w.write('statName',*vals)

class statWriter(object):
    
    # baseDir:      base statistics directory. Results would be written to:
    #               baseDir/id_DIR/statistic/runFile.txt;
    #               where id_DIR is shared between runs with the same identifiers;
    #               statistic is a directory for each saved statistic 'statName'
    #               and runFile is unique to each run regardless of identifiers.
    # identifiers:  dict by keywoards. all values must be hashable.
    
    def __init__(self, baseDir, **identifiers):
        
        # map the identifiers specified to a unique directory (sorting required for deterministic hashing order)
        hashID = hash(tuple(sorted(identifiers.items())))
        self._dir = baseDir + '/' + str(abs(hashID)) + '/'
        
        # if the directory doesn't exist, create it and save the identifiers
        if not os.path.exists(self._dir):
            os.makedirs(self._dir)
            with open(self._dir+'identifiers.txt', "w") as f:
                f.write('creation time : ' + str(datetime.now()) + '\n')
                f.writelines([str(k) + '\t: ' + str(v) + '\n' for k,v in identifiers.items()])
        
        self._writers = {}
        self._runID = str(getrandbits(64))        
    
    def write(self, key, *vals):
        
        # if there is no writer for this key, create it
        if key not in self._writers:
            dirKey = self._dir + str(key) + '/'
            # if there is no directory for this key, create it
            if not os.path.exists(dirKey):
                os.makedirs(dirKey)
            # create writer
            self._writers[key] = open(dirKey + key + self._runID + '.csv','w')
        
        # convert values to string and write to file. 
        self._writers[key].write(self._recuToStr(vals) + '\n')
    
    # conversion to string: sequences recursively converted to tuples.
    # for example: (1,2,[3,4,[5,6]]) -> '1,2,...,6'
    def _recuToStr(self, inp):
        if isinstance(inp, str):
            return inp
        res = ''
        for x in inp:
            if isinstance(x,Sequence):
                res+= self._recuToStr(x) + ','
            else:
                res+= str(x) + ','
        return res[:-1]

    def closeWriters(self):
        for _,w in self._writers.iteritems():
            w.close()
            
