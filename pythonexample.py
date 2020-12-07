from pykmapcore import KmapDB

NDATASETS = 50

db = KmapDB('build/rdbtest')
print(db.get('[[[[M',NDATASETS))
