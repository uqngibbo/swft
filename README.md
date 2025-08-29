# swft -- A Q1D supersonic flowpath analysis code

swft is a quasi-onedimensional compressible flow code for analysing supersonic flowpaths. It is written in the D programming language and relies on the Gasdynamic Toolkit (GDTk).

## Build Instructions
Install dependancies:

    sudo apt install build-essential libreadline-dev libncurses5-dev python3-cffi

Install the LLVM D compiler:

    wget https://github.com/ldc-developers/ldc/releases/download/v1.40.1/ldc2-1.40.1-linux-x86_64.tar.xz
    tar -xvf ldc2-1.41.1-linux-x86_64.tar.xz
    echo 'export PATH=${PATH}:${HOME}/ldc2-1.40.1-linux-x86_64/bin' >> .bashrc
    echo 'export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/ldc2-1.40.1-linux-x86_64/lib' >> .bashrc

Clone the GDTk repository:

     git clone https://github.com/gdtk-uq/gdtk.git gdtk
     echo 'export DGD_REPO=$HOME/gdtk' >> .bashrc

Build the code:

    git clone git@github.com:uqngibbo/swft.git
    cd swft/source
    make install

## Author:
Nick Gibbons (n.gibbons@uq.edu.au)

## License:
This Code Should Not Be Used By Anyone Under Any Circumstances
