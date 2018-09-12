**INSTALLATION**

The following commands will install the core dependencies for Debian-based systems (e.g., Ubuntu):

    sudo apt-get install -y -qq build-essential zlib1g-dev unzip libncurses5 libncurses5-dev libdb-dev git cpanminus libexpat1 libexpat1-dev

For RHEL-based systems (e.g., CentOS/Fedora):

    sudo yum groupinstall -y "Development Tools"
    sudo yum install -y perl-App-cpanminus ncurses ncurses-devel libdb-devel expat expat-devel zlib-devel java-1.7.0-openjdk

The next two commands install BioPerl, and these can be skipped if BioPerl is installed:
    
    echo "n" | cpanm -n  Data::Stag DB_File Bio::Root::Version Bio::SearchIO::blastxml

Finally, download the [latest release](https://github.com/sestaton/tephra/releases/latest) and run the following commands from the root directory:

    cpanm --installdeps .
    perl Makefile.PL
    make test
    make install

Please note, the above instructions will install Tephra for a single user. If you would like to configure Tephra to be installed for all users on a cluster, you will need to set the TEPHRA_DIR environment variable. For example,

    export TEPHRA_DIR=/usr/local/tephra
    perl Makefile.PL
    make test
    make install

will configure the software for all users. Please note that if Tephra is configured in a custom location this way it will be necessary to set this variable prior to using Tephra so the configuration can be found. In this case, just export the variable the same way. For a regular user, this can be done with a single line as below (note that this is the same command used to install/configure Tephra):

    export TEPHRA_DIR=/usr/local/tephra

Now you can type any command to use the usage, for example:

    tephra findltrs -h

For developers, please run the tests with:

    export TEPHRA_ENV='development' && make test

Please report any test failures or installation issues with the [issue tracker](https://github.com/sestaton/tephra/issues).