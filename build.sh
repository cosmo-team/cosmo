    set -xe
    mkdir 3rd_party_src
    mkdir -p 3rd_party_inst/boost
    INSTDIR=`pwd`/3rd_party_inst    
    
    cd 3rd_party_src
    git clone https://github.com/refresh-bio/KMC
    cd KMC
    git checkout f090276855a3f7c0b14e9f3abc8c99d3213247b3
    cd ..
    git clone https://github.com/cosmo-team/sdsl-lite.git
    cd sdsl-lite
    git checkout 9fa981958a9d2ddade12d083548f2b09939514fb
    cd ..
    git clone https://github.com/stxxl/stxxl
    cd stxxl
    git checkout 5b9663e6b769748f3b3d3a9a779b4b89e24d7a27
    cd ..
    git clone https://github.com/eile/tclap
    cd tclap
    git checkout f41dcb5ce3d063c9fe95623193bba693338f3edb
    cd ..
    git clone https://github.com/emil-e/rapidcheck
    cd rapidcheck
    git checkout b2032e6e7029a6e5183ddc9f0c1f6edf5d619a4a
    cd ..
    git clone https://github.com/mmuggli/sdreader
    wget http://sourceforge.net/projects/boost/files/boost/1.54.0/boost_1_54_0.tar.bz2
    tar -xjf boost_1_54_0.tar.bz2

    # Build the  dependencies
    cd boost_1_54_0
    ./bootstrap.sh --prefix=../../3rd_party_inst/boost
    ./b2 install
    cd ..
    
    cd sdsl-lite/
    /usr/bin/time sh install.sh ${INSTDIR}
    cd ..

    cd stxxl
    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${INSTDIR} -DBUILD_STATIC_LIBS=ON
    make
    make install
    cd ../..

    cd KMC
    make
    cd ..

    cd tclap/
    autoreconf -fvi
    ./configure --prefix=${INSTDIR}
    make
    make install
    cd ..

    cd rapidcheck
    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${INSTDIR}
    make
    make install
    cd ../..

    cd sdreader
    make
    cd ..
    
    # Build VARI
    cd ..
    make

