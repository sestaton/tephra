# Tephra
Transposable Element Paleontology by Historical Reconstruction from Annotations

### What is Tephra?

Tephra is a command line application to annotate [transposable elements](http://en.wikipedia.org/wiki/Transposable_element) from a genome assembly. The goal is to provide a high quality set of de novo annotations for all transposon types, describe the structure and evolution of those sequences, and do it without a reference set of transposon sequences (therefore being unbiased as possible).

**DEPENDENCIES**

Part of the utility of Tephra is to provide family-level TE classifications and infer patterns of molecular evoltion. To be efficient as possible, these tasks require a few external programs. Specifically, you will need to download [MUSCLE](http://http://drive5.com/muscle/) and [Vmatch](http://vmatch.de), both of which require a license so I cannot distribute them (but they are free). If you are only interested in TE identification, you can skip the need to download these programs. At this time, you also need to install [HMMER](http://hmmer.org/) version 3 and [EMBOSS](http://emboss.sourceforge.net/), which are freely available, as well.

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

