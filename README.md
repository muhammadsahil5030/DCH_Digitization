# DCH_Digitization
To run the digitization of the Drift chamber first we need:
### 1. EDM4hep
### 2. FCCDetector
### 3. k4RecTracker

## EDM4hep:
Source the key4hep setup
```
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
```

#### Clone EDM4hep repository and then build and install
```
git clone https://github.com/key4hep/EDM4hep
cd EDM4hep
mkdir build && cd build

# Configure with install prefix
cmake .. -DCMAKE_INSTALL_PREFIX=/afs/cern.ch/user/m/msaiel/FCC_study/EDM4hep/install

# Build and install
make -j4
make install
```
### Export the Environment variables
```
export EDM4HEP_DIR=/afs/cern.ch/user/m/msaiel/FCC_study/EDM4hep/install
export PATH=$EDM4HEP_DIR/bin:$PATH
export LD_LIBRARY_PATH=$EDM4HEP_DIR/lib:$LD_LIBRARY_PATH
export CMAKE_PREFIX_PATH=$EDM4HEP_DIR:$CMAKE_PREFIX_PATH
```
## FCCDetector
```
git clone https://github.com/HEP-FCC/FCCDetectors
cd FCCDetectors
mkdir build install
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
cd ../../
export FCCDETECTORS=$PWD/FCCDetectors/;PATH=$PWD/FCCDetectors/install/bin/:$PATH;CMAKE_PREFIX_PATH=$PWD/FCCDetectors/install/:$CMAKE_PREFIX_PATH;LD_LIBRARY_PATH=$PWD/FCCDetectors/install/lib:$LD_LIBRARY_PATH;export PYTHONPATH=$PWD/FCCDetectors/install/python:$PYTHONPATH;LD_LIBRARY_PATH=$PWD/FCCDetectors/install/lib64:$LD_LIBRARY_PATH
```
## k4RecTracker
```
git clone git@github.com:key4hep/k4RecTracker.git
cd k4RecTracker
mkdir build install
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
cd ../../
export K4RECTRACKER=$PWD/k4RecTracker/install/share/k4RecTracker; PATH=$PWD/k4RecTracker/install/bin/:$PATH; CMAKE_PREFIX_PATH=$PWD/k4RecTracker/install/:$CMAKE_PREFIX_PATH; LD_LIBRARY_PATH=$PWD/k4RecTracker/install/lib:$PWD/k4RecTracker/install/lib64:$LD_LIBRARY_PATH; export PYTHONPATH=$PWD/k4RecTracker/install/python:$PYTHONPATH
```
The following should be called in the folder hosting both k4RecTracker and FCCDetectors each time you start a new session (even if you do not need to rebuild):
```
export FCCDETECTORS=$PWD/FCCDetectors/;PATH=$PWD/FCCDetectors/install/bin/:$PATH;CMAKE_PREFIX_PATH=$PWD/FCCDetectors/install/:$CMAKE_PREFIX_PATH;LD_LIBRARY_PATH=$PWD/FCCDetectors/install/lib:$LD_LIBRARY_PATH;export PYTHONPATH=$PWD/FCCDetectors/install/python:$PYTHONPATH;LD_LIBRARY_PATH=$PWD/FCCDetectors/install/lib64:$LD_LIBRARY_PATH
export K4RECTRACKER=$PWD/k4RecTracker/install/share/k4RecTracker; PATH=$PWD/k4RecTracker/install/bin/:$PATH; CMAKE_PREFIX_PATH=$PWD/k4RecTracker/install/:$CMAKE_PREFIX_PATH; LD_LIBRARY_PATH=$PWD/k4RecTracker/install/lib:$PWD/k4RecTracker/install/lib64:$LD_LIBRARY_PATH; export PYTHONPATH=$PWD/k4RecTracker/install/python:$PYTHONPATH
```
