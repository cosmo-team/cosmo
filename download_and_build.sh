    set -xe

    # Fetch software and setup directories
    git clone https://github.com/cosmo-team/cosmo/
    cd cosmo/
    git checkout VARI-merge
    source build.sh    
