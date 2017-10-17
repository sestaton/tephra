# TEPHRA TODO

This file is for logging feature requests and bugs during development. Hopefully, having one list will make it easier to keep track of proposed changes. It would be nice to rank the lists in order to prioritize tasks. It should be noted this list is for development purposes and it may go away once a stable release is made.

## Command `tephra classifyltrs`
 - [x] Classify 'best' LTR-RT elements into superfamilies based on domain content and organization
 - [x] Classify elements into families based on cluster organization; Generate FASTA for each family
 - [x] Create singletons file of ungrouped LTR-RT sequences
 - [x] Create FASTA files of exemplars for each LTR-RT family
 - [x] clean up logs and intermediate files from vmatch (dbluster-*)
 - [x] combine domain organization from both strands (if the same)?
 - [x] add family classifications to GFF Name attribute
 - [x] incorporate legacy annotations from input GFF/reference to family classification
 - [x] merge overlapping hits in chain of protein matches, and contatenate the rest for each element
 - [ ] mark unclassified elements with no protein domains as LARDs
 - [x] combine exemplars for efficiently comparing to a reference set
 - [x] identify fragmented elements with refined full-length elements (handled in v0.08.0+ in 'getfragments'
       command)
 - [x] include measure of similarity within/between families
 
## Command `tephra findtirs`
 - [x] Find all non-overlapping TIR elements passing thresholds
 - [x] Generate combined GFF3 of high-quality TIRs
 - [x] Check for index (if given)
 - [ ] Add optional test for the presence of coding domains similar to 'LTRRefine' class. This should reduce the
       number of DTX elements. Add this to the configuration file for the 'all' command the same as for LTRs.


## Command `tephra sololtr`
 - [x] Create HMM of LTRs for each LTR-RT
 - [x] Search masked ref with LTR HMM
 - [x] Create GFF with SO terms of solo-LTRs
 - [x] parallelize hmmsearch to speed things up. likely this is faster than multiple cpus for one model at time
 - [x] check if input directory exists
 - [x] make sure to set path to correct version of hmmer
 - [x] add family name to GFF output (the family name is now in the Parent tag)
 - [x] add option to pick on the top 20 families to speed up execution
 - [ ] consider preprocessing all LTR files so we don't block on one superfamily waiting for threads to finish
 - [x] if the soloLTR sequence file is empty, delete all other files and warn no soloLTRs were found

## Command `tephra classifytirs`
 - [x] Classify 'best' TIR elements into superfamilies based on domain content, TSD, and/or motif
 - [x] Group TIR elements into families based on TIR similarity and/or cluster-based method used for LTR-RT classification 
 - [x] in tests, skip if empty output (none found). This is not a good test honestly, need a new reference
 - [x] write fasta of each superfamily, and combined library
 - [x] identify	fragmented elements with refined full-length elements

## Command `tephra findltrs`
 - [x] Find all non-overlapping LTR-RTs under strict and relaxed conditions
 - [x] Filter elements by quality score, retaining the best elements
 - [x] Generate combined GFF3 of high-quality LTR-RTs
 - [x] Check for index (if given)
 - [x] change header format to be ">id_source_range"
 - [x] adjust filtering command to not increment if element has been deleted (inflated filtering stats)
 - [x] reporting of superfamilies after ltr search?.. better to do that at classification stage
 - [x] add options for LTR size parameters
 - [ ] add LTR_Finder
 - [x] add config file to handle the multitude of LTR-RT constraints
 - [x] clean up ltrharvest and ltrdigest intermediate files
 - [x] Add optional test for the presence of coding domains to 'LTRRefine' class. This should reduce the
       number of RLX elements.

 - Domain matches 
   - [ ] adjust duplicate domain filtering to consider strand and range of matches
   - [x] fix reporting of overlapping domain matches by LTRdigest? (issue reported: https://github.com/genometools/genometools/issues/706)
   - [x] add e-value threshold option and domain filtering method 

## Command `tephra findhelitrons`
 - [x] Find helitrons in reference sequences with HelitronScanner
 - [x] Generate GFF3 of full-length helitrons
 - [ ] Annotate coding domains in helitrons and include domains in GFF 
 - [x] Adjust header for full length elements to match output of other commands
 - [ ] Remove strand from FASTA header for consistency with other commands

## Command `tephra findtrims`
 - [x] Find all non-overlapping TRIMs under strict and relaxed conditions
 - [x] Filter elements by quality score, retaining the best elements
 - [x] Generate combined GFF3 of high-quality TRIMs
 - [x] Create a feature type called 'TRIM_retrotransposon' to distinguish these elements from other LTR-RTs
 - [x] create developer tests to operate on a larger data set to positively identify elements rather than just
       operation of the command

## Command `tephra findnonltrs`
 - [ ] break chromosomes to reduce memory usage in hmmsearch (only applies to HMMERv3)
 - [x] check HMMER2 var and program version
 - [x] remove backticks and shell exec of hmmer
 - [x] remove nasty regex parsing in favor or bioperl reading of report
 - [x] use list form of system to not fork
 - [ ] run domain searches in parallel
 - [ ] use multiple CPUs (make option) for domain searches
 - [x] write GFF of results
 - [x] add verbose option so as to not print progress when there are 5k scaffolds
 - [x] write combined file of all elements
 - [x] take a multifasta as input and create directories for input/output to methods
 - [x] use complete elements to find truncated nonLTRs after masking (do this with complete file at the end
       on masked genome to get fragments for all types)
 - [x] use domain/blast based method for classifying elements into families
 - [ ] investigate issues related to why most elements reported on negative strand and contain
       many gaps

## Command `tephra ltrage`
 - [x] Calculate age for each LTR-RT
 - [x] Take substitution rate as an option
 - [x] check if input directory exists
 - [ ] write age to GFF file
 - [x] Clean up results if requested

## Command `tephra maskref`
 - [x] Generate masked reference from custom repeat library 
 - [x] Add outfile option instead of creating filename
 - [x] Make some kind of statistical report about masking percentage. It would be helpful to format
       the output like RepeatMasker to give a global view of what was masked.
 - [x] Clean up the intermediate folders for each chromosome when masking the genome
 - [x] Create overlapping windows for masking subsets to solve the issue of reduced representation when
       generating smaller chunks

## Command `tephra illrecomb`

 - [x] Add correct sequence IDs to report
 - [ ] Investigate the apparent disagreement between the query/subject string and homology strings
 - [ ] Summarize the stats in a more intuitive way so it is clear what the gap summaries mean

## Command `tephra tirage`

 - [x] Update menu for all available options.
 - [x] Add 3-letter code to age file IDs
 - [x] Clean up results if requested
 - [x] Add method to select the top families instead of --all (requires generating families first)

## Command `tephra all`

 - [ ] Allow the user to pass a genome and repeat database, along with a species name instead of configuration file.
 - [x] Generate summary statistics for TE types (domain content, length distribution, diversity, etc.) See
       (sesbio/transposon_annotation/count_families.pl) for starters.
 - [ ] Generate HTML output for all command. Will need to store JSON data for graphs and tables.
 - [x] Add tirage options to configuration file.
 - [ ] Remove FASTA/GFF3 files of unclassified elements once the classification process is complete. 

*** 

## Meta
 - [x] logging results/progress (need to log progress and errors to the correct location)
 - [x] add debug options for seeing commands (done for LTR search)
 - [ ] documentation of algorithms, in addition to usage
 - [ ] reduce LTRs/TRIMs....perhaps when combining all GFFs
 - [ ] save tnp matching ltr-rts and search for cacta tes..or just add as putative classII
 - [ ] add kmer mapping command (see tallymer2gff.pl)
 - [x] create config role for setting paths
 - [x] change config module to be an Install namespace
 - [x] add subcommand to merge all GFFs (can be done with `gt gff3 sort`, though we want to be careful
       not to bake in too many subcommands for things that are easily done at the command line already,
       as this will make the package harder to use and maintain).
 - [ ] handle compressed input/output
 - [X] add fasta-handling classes from Transposome, which are faster than BioPerl (Won't do: Added kseq.h methods from HTSlib)
 - [ ] add verbose option for quickly debugging the installation of dependencies
 - [x] add command to get TIR ages
 - [x] investigate why tests fail with Perl version 5.12 or lower (Bio::DB::HTS needs 5.14.2, so that's why)
 - [x] add subcommand to run/log all methods as a pipeline
 - [ ] document the configuration file format and usage
 - [x] add 'findfragments' subcommand to be run after final masking prior to complete GFF generation
 - [x] add classification method for TRIMs
 - [ ] move 'classify[ltr|tir]' commands to 'find[ltr|tir]' commands to simplify the process similar to the methods
       for the commands for helitrons and tirs
 - [ ] modify header to include element number in family. The element number should be listed numerically according 
       to chromosome position? (Wicker et al., 2007)
