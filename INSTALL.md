**INSTALLATION OF CORE DEPENDENCIES**

The following commands will install the core dependencies for Debian-based systems (e.g., Ubuntu):

    sudo apt-get install -y -qq build-essential zlib1g-dev libgd-dev unzip libncurses5 libncurses5-dev libdb-dev git cpanminus libexpat1 libexpat1-dev

For RHEL-based systems (e.g., CentOS/Fedora):

    sudo yum groupinstall -y "Development Tools"
    sudo yum install -y perl-App-cpanminus ncurses ncurses-devel gd-devel libdb-devel expat expat-devel zlib-devel java-1.7.0-openjdk

**INSTALLATION OF PERL DEPENDENCIES**

Tephra requires Perl to be installed and though it is installed by default on Unix-based systems I highly recommend that users  configure you own Perl and leave the system Perl alone. 

For regular users, I recommend using [Perlbrew](https://perlbrew.pl/) or [Plenv](https://github.com/tokuhirom/plenv) (see also this [blog post](https://weblog.bulknews.net/plenv-alternative-for-perlbrew-7b5bf00a419e)). FWIW, I use Perlbrew on all Linux machines and Plenv works better on Mac in my experience with the frequent updates and permission changes. 

For System Adminstrators wanting to install Tephra for all users, please see the note below.

---

Assuming you have a working Perl and App::cpanminus, the next two commands install BioPerl (this can be skipped if BioPerl is installed):
    
    echo "n" | cpanm -n  Data::Stag DB_File Bio::Root::Version Bio::SearchIO::blastxml Bio::SearchIO::hmmer

Finally, download the [latest release](https://github.com/sestaton/tephra/releases/latest) and run the following commands from the root directory:

    cpanm --installdeps .
    perl Makefile.PL
    make test
    make install

**FOR SYSTEM ADMINISTRATORS**

Please note, the above instructions will install Tephra for a single user. If you would like to configure Tephra to be installed for all users on a cluster, you will need to set the TEPHRA_DIR environment variable. For example,

    export TEPHRA_DIR=/usr/local/tephra
    perl Makefile.PL
    make test
    make install

will configure the software for all users. It is important to note that if Tephra is configured in a custom location this way it will be necessary to set this ENV variable prior to a user running Tephra; that way the configuration can be found. In this case, just export the variable the same way. For a regular user, this can be done with a single line as below (note that this is the same command used to install/configure Tephra):

    export TEPHRA_DIR=/usr/local/tephra

Therefore, it only requires that one extra variable (and the permissions) to install and run Tephra from a custom location. Now you can type any command to see the usage, for example:

    tephra findltrs -h

**FOR DEVELOPMENT**

For developers, please run the tests with:

    export TEPHRA_ENV='development' && make test

Be warned that this will automatically fetch the latest *Arabidopsis thaliana* genome from TAIR and run much more extensive tests, as well as the full pipeline, which will take a number of hours. See the test files in the [/t](https://github.com/sestaton/tephra/tree/master/t) directory for information about the test suite, and don't hesitate to ask if me if you have any questions/concerns.

**FOR ISSUES AND ASSISTANCE**

Please report any test failures or installation issues with the [issue tracker](https://github.com/sestaton/tephra/issues).