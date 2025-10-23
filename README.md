# DCH_Digitization
To run the digitization of the Drift chamber, first we need:
### 1. EDM4hep
### 2. FCCDetector
### 3. k4RecTracker

## Dependencies:
1. ROOT
2. PODIO
3. EDM4HEP
4. Gaudi
5. k4FWCore
6. DD4HEP

First we need to source the key4hep setup
```
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
```
## EDM4hep
#### Clone the EDM4hep repository and then build and install
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
#### Export the Environment variables
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
#### Run the digitization
Finally, we are ready to run the digitization.
To run the simple digitization, we can follow the steps given in the k4RecTracker repository, but to run the actual digitization, we need to follow some steps
1. Running the simulated steering
   ```
   ddsim --steeringFile sim_steering.py --outputFile 'dch_proton_10GeV.root' -N 10 --runType batch --random.seed 42
   ```
2. Download the file for cluster counting
   ```
   https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/IDEA/DataAlgFORGEANT.root
   ```
3. Run the digitizer for position smearing and cluster counting calculation
   ```
   k4run runDCHdigi.py
   ```
4. Check the distribution of distances from the hit position to the wire
   ```
   python3 check_DCHdigi_output.py
   ```
Now, if we want to make some changes in the DCHdigi_v01.cpp file, then we need to:
1. Rebuild the k4RecTracker
   ```
   build -j4
   ```
2. Export
   ```
   export LD_LIBRARY_PATH=/afs/cern.ch/user/m/msaiel/Digitization/k4RecTracker/build/DCHdigi:$LD_LIBRARY_PATH
   ```
3. Run the Python code again
   ```
   k4run DCHdigi_v01.py
   ```
4. And then follow your analysis codes.
