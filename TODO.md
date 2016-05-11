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
 - [ ] incorporate legacy annotations from input GFF/reference to family classification
 - [x] merge overlapping hits in chain of protein matches, and contatenate the rest for each element

## Command `tephra findtirs`
 - [x] Find all non-overlapping TIR elements passing thresholds
 - [x] Generate combined GFF3 of high-quality TIRs
 - [x] Check for index (if given)

## Command `tephra sololtr`
 - [x] Create HMM of LTRs for each LTR-RT
 - [x] Search masked ref with LTR HMM
 - [x] Create GFF with SO terms of solo-LTRs
 - [ ] parallelize hmmsearch to speed things up. likely this is faster than multiple cpus for one model at time
 - [x] check if input directory exists
 - [x] make sure to set path to correct version of hmmer
 - [ ] add family name to GFF output

## Command `tephra classifytirs`
 - [x] Classify 'best' TIR elements into superfamilies based on domain content, TSD, and/or motif
 - [ ] Group TIR elements into families based on TIR similarity and/or cluster-based method used for LTR-RT classification 
 - [x] in tests, skip if empty output (none found). This is not a good test honestly, need a new reference
 - [ ] write fasta of each superfamily, and combined library

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

 - Domain matches 
   - [ ] adjust duplicate domain filtering to consider strand and range of matches
   - [x] fix reporting of overlapping domain matches by LTRdigest? (issue reported: https://github.com/genometools/genometools/issues/706)
   - [x] add e-value threshold option and domain filtering method 

## Command `tephra findhelitrons`
 - [x] Find helitrons in reference sequences with HelitronScanner
 - [x] Generate GFF3 of full-length helitrons
 - [ ] Annotate coding domains in helitrons and include domains in GFF 
 - [ ] Adjust header for full length elements to match output of other commands

## Command `tephra findtrims`
 - [x] Find all non-overlapping TRIMs under strict and relaxed conditions
 - [x] Filter elements by quality score, retaining the best elements
 - [x] Generate combined GFF3 of high-quality TRIMs
 - [ ] Create a feature type called 'TRIM_retrotransposon' to distinguish these elements from other LTR-RTs

## Command `tephra findnonltrs`
 - [ ] break chromosomes to reduce memory usage in hmmsearch (only applies to HMMERv3)
 - [x] check HMMER2 var and program version
 - [x] remove backticks and shell exec of hmmer
 - [x] remove nasty regex parsing in favor or bioperl reading of report
 - [x] use list form of system to not fork
 - [ ] run domain searches in parallel
 - [ ] use multiple CPUs (make option) for domain searches
 - [x] write GFF of results
 - [ ] add verbose option so as to not print progress when there are 5k scaffolds
 - [x] write combined file of all elements
 - [x] take a multifasta as input and create directories for input/output to methods
 - [ ] use complete elements to find truncated nonLTRs after masking

## Command `tephra ltrage`
 - [x] Calculate age for each LTR-RT
 - [x] Take substitution rate as an option
 - [x] check if input directory exists
 - [ ] write age to GFF file

## Command `tephra maskref`
 - [x] Generate masked reference from custom repeat library 
 - [x] Add outfile option instead of creating filename

*** 

## Meta
 - [ ] logging results/progress
 - [ ] add debug options for seeing commands (done for LTR search)
 - [ ] documentation of algorithms, in addition to usage
 - [ ] reduce LTRs/TRIMs....perhaps when combining all GFFs
 - [ ] save tnp matching ltr-rts and search for cacta tes..or just add as putative classII
 - [ ] add kmer mapping command (see tallymer2gff.pl)
 - [x] create config role for setting paths
 - [x] change config module to be an Install namespace
 - [ ] add subcommand to merge all GFFs (can be done with `gt gff3 sort`, though we want to be careful
       not to bake in too many subcommands for things that are easily done at the command line already,
       as this will make the package harder to use and maintain).
 - [ ] handle compressed input/output
 - [ ] add fasta-handling classes from Transposome, which are faster than BioPerl
