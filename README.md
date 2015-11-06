# Tephra
Transposable Element Paleontology by Historical Reconstruction from Annotations

### What is Tephra?

Tephra is a command line application to annotate [transposable elements](http://en.wikipedia.org/wiki/Transposable_element) from a genome assembly. The goal is to provide a high quality set of de novo annotations for all transposon types, describe the structucture and evolution of those sequences, and do it without a reference set of transposon sequences (therefore being unbiased as possible)..

**DEPENDENCIES**

One of the final stages of Tephra is to infer patterns of molecular evoltion, and this requires a few external programs. Specifically, you will need to download [ClustalW](http://clustal.org/download/2.1/) and [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html#download). The versions tested are Clustal v2.1 and PAML v4.8. Unfortunately, I cannot guarantee everything will work, especially with different versions of PAML, so please aim to have those versions in your PATH. If this is a problem, please let me know.

**SUPPORT AND DOCUMENTATION**

You can get usage information at the command line with the following command:

    perldoc tephra

The `tephra` program will also print a diagnostic help message when executed with no arguments.

You can also look for information at:

    Tephra wiki
            https://github.com/sestaton/tephra/wiki

    Tephra issue tracker
            https://github.com/sestaton/tephra/issues


**LICENSE AND COPYRIGHT**

Copyright (C) 2015 S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package.
If not, it can be found here: http://www.opensource.org/licenses/mit-license.php

