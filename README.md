# Tephra
Transposable Element Paleontology by Historical Reconstruction from Annotations

### What is Tephra?

Tephra is a command line application to annotate [transposable elements](http://en.wikipedia.org/wiki/Transposable_element) from a genome assembly. The goal is to provide a high quality set of de novo annotations for all transposon types, describe the structure and evolution of those sequences, and do it without a reference set of transposon sequences (therefore being unbiased as possible).

**DEPENDENCIES**

Part of the utility of Tephra is to provide family-level TE classifications and infer patterns of molecular evoltion. To be efficient as possible, these tasks require a few external programs. Specifically, you will need to download [MUSCLE](http://http://drive5.com/muscle/) and [Vmatch](http://vmatch.de), both of which require a license so I cannot distribute them (but they are free). If you are only interested in TE identification, you can skip the need to download these programs. However, you do need to install [HMMER](http://hmmer.org/) version 3 and [EMBOSS](http://emboss.sourceforge.net/) at this time.

**INSTALLATION**

After cloning the repo or downloading the latest release, the following commands will build and install the package:

    perl Makefile.PL
    make test
    make install

The test step is not strictly necessary, but it is highly encouraged to ensure that everything is configured correctly.

**SUPPORT AND DOCUMENTATION**

You can get usage information at the command line with the following command:

    perldoc tephra

The `tephra` program will also print a diagnostic help message when executed with no arguments, and display the available subcommands.

You can also look for information at:

    Tephra wiki
            https://github.com/sestaton/tephra/wiki

    Tephra issue tracker
            https://github.com/sestaton/tephra/issues

**USAGE**

Tephra is a command-line program only for now. The command `tephra` itself controls all the action of the subcommands, which perform specific tasks. Typing the command `tephra` will show the available commands. Here is an example,

    $ tephra
    tephra <command> 
    
    Available commands:
    
           commands: list the application's commands
               help: display a command's help screen

       classifyltrs: Classify LTR retrotransposons into superfamilies and families.
       classifytirs: Classify TIR transposons into superfamilies.
      findhelitrons: Find Helitons in a genome assembly.
           findltrs: Find LTR retrotransposons in a genome assembly.
        findnonltrs: Find non-LTR retrotransposons in a genome assembly.
           findtirs: Find TIR transposons in a genome assembly.
          findtrims: Find TRIM retrotransposons in a genome assembly.
          illrecomb: Characterize the distribution of illigetimate recombination in a genome.
             ltrage: Calculate the age distribution of LTR retrotransposons.
            maskref: Mask a reference genome with transposons.
            sololtr: Find solo-LTRs in a genome assembly.

Typing a subcommand will show the usage of that command, for example:

    $ tephra findnonltrs
    ERROR: Required arguments not given.
    
    USAGE: tephra findnonltrs [-h] [-m]
        -m --man      :   Get the manual entry for a command.
        -h --help     :   Print the command usage.
    
    Required:
        -g|genome     :   The genome sequences in FASTA format to search for non-LTR-RTs. 

    Options:
        -o|outdir     :   The location to place the results.
        -p|pdir       :   Location of the HMM models (Default: configured automatically).


**CITATION**

A manuscript is in preparation, which includes a description of the all the methods and their uses, a comparison to other programs, and results from model systems. These will be provided in some form ahead of publication, as soon as they are available.

For now, please cite the github URL of this repo if you use Tephra. Thank you. 

**LICENSE AND COPYRIGHT**

Copyright (C) 2015 S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package.
If not, it can be found here: http://www.opensource.org/licenses/mit-license.php

