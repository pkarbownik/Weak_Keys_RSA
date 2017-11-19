# CUDA_weak_RSA_keys_breaker
  <h3>Install Nvidia drivers and CUDA toolkit 6.5 on Ubuntu 16.04 for GT216M [NVS 5100M] :
Driver:</h3>

  1. Check name of nvidia product that you have: lspci | grep -i nvidia
  2. Select appropiate driver based on above result of command http://www.nvidia.com/Download/index.aspx
  3. Download driver (version 340.102 for GT216M [5100M] - http://www.nvidia.com/download/driverResults.aspx/114719/en-us).
  4. Go to console mode using Alt+Ctrl+F1 and type: sudo service lightdm stop to diable ubuntu GUI
  5. In directory where downloaded driver is located type: chmod +x NVIDIA[TAB] (for GT216M: chmod +x NVIDIA-Linux-x86_64-340.102.run)
  6. Install by typing: ./NVIDIA[TAB] (or GT216M: ./NVIDIA-Linux-x86_64-340.102.run)
  7. Type: sudo service lightdm start to enable GUI and Go to GUI mode using Alt+Ctrl+F7

Information about Nvidia products and drivers: https://docs.parrotsec.org/doku.php/nvidia-drivers
  
  
  <h3>CUDA:</h3>
  
  Check which version of CUDA you need: https://en.wikipedia.org/wiki/CUDA in section named: "GPUs supported".
  For example for GT216M [NVS 5100M] CUDA version is 6.5, because compute capability is 1.2. 
  In below steps CUDA 6.5 will be installed:
  
  1. Download CUDA 6.5 for linux 64-bit from this site: https://developer.nvidia.com/cuda-toolkit-65
  2. In directory where downloaded CUDA installer type: chmod +x cuda_6.5.14_linux_64.run or 
  3. Install by typing: ./cuda_6.5.14_linux_64.run -toolkit -samples -silent -override
  4. Add libraries to PATH: export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64
  5. Add cuda binaries to PATH: export PATH=$PATH:/usr/local/cuda/bin
  6. Check nvcc version to verify is everything is right: nvcc -V
  
  Your changes in PATH and LD_LIBRARY_PATH enviroment are memorized only in current session. You can change enviroment for any session    by adding:
  PATH="$PATH:/usr/local/cuda/bin"
  LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64
  in ~/.profile file.
  
  <h3>Make and run example:</h3>
  
  Type to change samples owner: sudo chown user_name /usr/local/cuda/samples -R (user_name must be changed)
  
  1. cd /usr/local/cuda/samples/1_Utilities/deviceQuery
  2. make
  3. ./deviceQuery
  
  When everything is OK you see:</br></br>
  
  ./deviceQuery Starting...</br></br>

 CUDA Device Query (Runtime API) version (CUDART static linking)</br></br>

Detected 1 CUDA Capable device(s)</br></br>

Device 0: "NVS 5100M"</br>
  CUDA Driver Version / Runtime Version          6.5 / 6.5</br>
  CUDA Capability Major/Minor version number:    1.2</br>
  Total amount of global memory:                 1024 MBytes (1073414144 bytes)</br>
  ( 6) Multiprocessors, (  8) CUDA Cores/MP:     48 CUDA Cores</br>
  GPU Clock rate:                                1210 MHz (1.21 GHz)</br>
  Memory Clock rate:                             790 Mhz</br>
  Memory Bus Width:                              128-bit</br>
  Maximum Texture Dimension Size (x,y,z)         1D=(8192), 2D=(65536, 32768), 3D=(2048, 2048, 2048)</br>
  Maximum Layered 1D Texture Size, (num) layers  1D=(8192), 512 layers</br>
  Maximum Layered 2D Texture Size, (num) layers  2D=(8192, 8192), 512 layers</br>
  Total amount of constant memory:               65536 bytes</br>
  Total amount of shared memory per block:       16384 bytes</br>
  Total number of registers available per block: 16384</br>
  Warp size:                                     32</br>
  Maximum number of threads per multiprocessor:  1024</br>
  Maximum number of threads per block:           512</br>
  Max dimension size of a thread block (x,y,z): (512, 512, 64)</br>
  Max dimension size of a grid size    (x,y,z): (65535, 65535, 1)</br>
  Maximum memory pitch:                          2147483647 bytes</br>
  Texture alignment:                             256 bytes</br>
  Concurrent copy and kernel execution:          Yes with 1 copy engine(s)</br>
  Run time limit on kernels:                     Yes</br>
  Integrated GPU sharing Host Memory:            No</br>
  Support host page-locked memory mapping:       Yes</br>
  Alignment requirement for Surfaces:            Yes</br>
  Device has ECC support:                        Disabled</br>
  Device supports Unified Addressing (UVA):      No</br>
  Device PCI Bus ID / PCI location ID:           1 / 0</br>
  Compute Mode:</br>
     < Default (multiple host threads can use ::cudaSetDevice() with device simultaneously) ></br></br>

deviceQuery, CUDA Driver = CUDART, CUDA Driver Version = 6.5, CUDA Runtime Version = 6.5, NumDevs = 1, Device0 = NVS 5100M</br>
Result = PASS</br>

 <h3>Problems which may occur:</h3>
  
1.  /usr/local/cuda-5.5//bin/..//include/host_config.h:82:2: error: #error -- unsupported GNU version! gcc 4.9 and up are not supported! - 
Solution: </br>
sudo apt-get install gcc-4.8 g++-4.8</br>
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 50</br>
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.8 50</br></br>

or </br></br>

delete or comment this line in "host_config.h" file.</br>

Source: https://github.com/cbuchner1/CudaMiner/issues/155 

2. Missing libraries:
Missing recommended library: libGLU.so</br>
Missing recommended library: libXi.so</br>
Missing recommended library: libXmu.so</br></br>

Just type: apt-get install libglu1-mesa libxi-dev libxmu-dev libglu1-mesa-dev</br>
Source: https://stackoverflow.com/questions/22360771/missing-recommended-library-libglu-so</br>

3. Xorg system error:
 sudo echo -e 'Section "Device"\n\tIdentifier "My GPU"\n\tDriver "nvidia"\nEndSection' > /etc/X11/xorg.conf.d/20-nvidia.conf 

4. After login system return login screen again. [LOGIN LOOP]

This problem occurs when you update linux. You should check nvidia driver using nvidia-smi. When nvidia's drivers are lost you need install drivers and CUDA again. 
  
<h3>To unistall nvidia drivers and cuda - useful commands:</h3> 
sudo apt-get purge nvidia*</br>
sudo apt-get purge nvidia-cuda*</br>
sudo apt-get remove nvidia-cuda-toolkit</br>
sudo apt-get remove --auto-remove nvidia-cuda-toolkit</br>
sudo apt-get purge --auto-remove nvidia-cuda-toolkit</br></br>



<h3>System details:</h3></br></br>
$uname -a:</br></br>
Linux pkarbownik-EliteBook 4.4.0-79-lowlatency #100-Ubuntu SMP PREEMPT Wed May 17 20:56:57 UTC 2017 x86_64 x86_64 x86_64 GNU/Linux </br></br>
$gcc --version:</br></br>
gcc (Ubuntu 5.4.1-2ubuntu1~16.04) 5.4.1 20160904
Copyright (C) 2015 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
