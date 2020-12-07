from pykmapcore import KmapDB
import numpy as np
db = KmapDB('build/rdbtest')

NDATASETS = 50
print(db.get('[[[[M',NDATASETS))
