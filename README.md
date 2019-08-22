# Tephra
A tool for discovering transposable elements and describing patterns of genome evolution

[![Build Status](https://travis-ci.org/sestaton/tephra.svg?branch=master)](https://travis-ci.org/sestaton/tephra) [![Coverage Status](https://coveralls.io/repos/github/sestaton/tephra/badge.svg?branch=master)](https://coveralls.io/github/sestaton/tephra?branch=master) [![GitHub version](https://badge.fury.io/gh/sestaton%2Ftephra.svg)](https://badge.fury.io/gh/sestaton%2Ftephra)

### What is Tephra?

Tephra is a command line application to annotate [transposable elements](http://en.wikipedia.org/wiki/Transposable_element) from a genome assembly. The goal is to provide a high quality set of *de novo* annotations for all transposon types, describe the structure and evolution of those sequences, and do it without a reference set of transposon sequences (therefore being unbiased as possible).

**RECOMMENDED USAGE**

With [Singularity](https://sylabs.io/docs/), you will want to pull the latest container version and then run the container:

    singularity pull library://sestaton/default/tephra
    LC_ALL=C singularity run --bind $PWD/db:/db tephra_latest.sif

That assumes you have your data files for analysis (genome, config file, etc.) in a directory called `db` in the current working directory. After that, change to the `/db` directory and run your analysis.

---

With [Docker](https://www.docker.com/), you can create a container to run Tephra with the following command:

    docker run -it --name tephra-con -v $(pwd)/db:/db:Z sestaton/tephra

That will create a container called `tephra-con` and start an interactive shell. The above assumes you have a directory called `db` in the working directory that contains your database files and the Tephra configuration. To run the full analysis, change to the mounted directory with `cd /db` in your container and run the following command:

    tephra all -c tephra_config.yml

I recommend using `nohup` and then logging out, which will allow you to leave the container running in the background (I will expand this part of the install by v0.13.0).

If you cannot use Singularity or Docker, please see the [INSTALL](https://github.com/sestaton/tephra/blob/master/INSTALL.md) file included with this distribution to install Tephra on various operating systems.

**BASIC USAGE**

Tephra is a command-line program only for now. The command `tephra` itself controls all the action of the subcommands, which perform specific tasks. Typing the command `tephra` will show the available commands. Here is an example,

    $ tephra 

    Tephra version 0.12.5
    
    Copyright (C) 2015-2019 S. Evan Staton
    LICENSE -- MIT

    Citation: Staton, SE. 2019. https://github.com/sestaton/tephra

    Name:
         Tephra - A tool for discovering transposable elements and describing
         patterns of genome evolution
    
    Description:
         This is an application to find transposable elements based on structural and sequence similarity features,
         group those elements into recognized (superfamilies) and novel (families) taxonomic groups,
         and infer patterns of evolution.
    
    -------------------------------------------------------------------------------------------
    USAGE: tephra <command> [options]
    
    Available commands:
         
                age: Calculate the age distribution of LTR or TIR transposons.
                all: Run all subcommands and generate annotations for all transposon types.
       classifyltrs: Classify LTR retrotransposons into superfamilies and families.
       classifytirs: Classify TIR transposons into superfamilies.
      findfragments: Search a masked genome with a repeat database to find fragmented elements.
      findhelitrons: Find Helitons in a genome assembly.
           findltrs: Find LTR retrotransposons in a genome assembly.
        findnonltrs: Find non-LTR retrotransposons in a genome assembly.  
           findtirs: Find TIR transposons in a genome assembly.
          findtrims: Find TRIM retrotransposons in a genome assembly.
          illrecomb: Characterize the distribution of illegitimate recombination in a genome.
            maskref: Mask a reference genome with transposons.
         reannotate: Transfer annotations from a reference set of repeats to Tephra annotations.
            sololtr: Find solo-LTRs in a genome assembly.
    
    Most common usage:
    
        tephra all -c tephra_config.yml
    
     That will produce a FASTA and GFF3 of all intact and fragmented transposons in the genome,
     and generate a table of annotation results.
    
    To get the configuration file, run:
    
        wget https://raw.githubusercontent.com/sestaton/tephra/master/config/tephra_config.yml
    
    To see information about a subcommand, run:
    
        tephra <command> --help
    
    To get more detailed information, run:
    
        tephra <command> --man


Typing a subcommand will show the usage of that command, for example:

    $ tephra findnonltrs

    [ERROR]: Required arguments not given.
    
    Name:
         tephra findnonltrs - Find non-LTRs retrotransposons in a genome assembly.
    
    Description:
         Find non-LTR retrotransposons in a reference genome, classify them into known superfamilies, 
         and generate a GFF file showing their locations and properties.
    
    USAGE: tephra findnonltrs [-h] [-m]
        -m --man      :   Get the manual entry for a command.
        -h --help     :   Print the command usage.
    
    Required:
        -g|genome     :   The genome sequences in FASTA format to search for non-LTR-RTs. 
        -o|gff        :   The GFF3 outfile to place the non-LTRs found in <genome>.
    
    Options:
        -r|reference  :   The non-masked reference genome for base correction.
        -d|outdir     :   The location to place the results.
        -p|pdir       :   Location of the HMM models (Default: configured automatically).
        -t|threads    :   The number of threads to use for BLAST searches (Default: 1).
        -v|verbose    :   Display progress for each chromosome (Default: no).


**SUPPORT AND DOCUMENTATION**

You can get usage information at the command line with the following command:

    perldoc tephra

The `tephra` program will also print a diagnostic help message when executed with no arguments, and display the available subcommands.

You can also look for information at:

    Tephra wiki
            https://github.com/sestaton/tephra/wiki

    Tephra issue tracker
            https://github.com/sestaton/tephra/issues

 
**CITATION**

A manuscript is in preparation, which includes a description of the all the methods and their uses, a comparison to other programs, and results from model systems. These will be provided in some form ahead of publication, as soon as they are available.

For now, please cite the github URL of this repo if you use Tephra. Thank you. 

**CURRENT STATUS**

Please check the [wiki](https://github.com/sestaton/tephra/wiki) for progress updates.

**CONTRIBUTING**

I welcome any comments, bug reports, feature requests, or contributions to the development of the project. Please submit a new issue (preferred) or send me an email and I would be happy to talk about Tephra or transposons.

**LICENSE AND COPYRIGHT**

Part of this project uses code from [MGEScan-nonLTR](http://darwin.informatics.indiana.edu/cgi-bin/evolution/nonltr/nonltr.pl), which is released under the GPL license. With permission of the authors, this code is packaged with Tephra. Below is the copyright for MGEScan-nonLTR:

    Copyright (C) 2015. See the LICENSE file for license rights and limitations (GPL v3).

    This program is part of MGEScan.

    MGEScan is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

The license for Tephra is below:

Copyright (C) 2015-2019 S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package.
If not, it can be found here: http://www.opensource.org/licenses/mit-license.php

