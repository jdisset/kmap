# KMAP
## Installation instructions
Kmap can be built using the included CMakeLists.txt. Kmap depends on rocksDB, which will be automatically pulled from the official repository and compiled. RocksDB depends on libBz2, which should be available from most system package managers.

To install kmap, clone this repository, cd into it and

```
mkdir build
cd build
cmake ..
make -j8
```

## Usage
