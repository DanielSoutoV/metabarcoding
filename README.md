# metabarcoding
Analyzing metabarcoding data using the packages Metacoder and MicrobiotaProcess

This project uses COI metabarcoding data from two nights of light-trapping at Barro Colorado Island in two seasons (Wet/Dry).

The project aims at describing biodiversity classified through BOLD's mBrave platform based on BINs, and then compare it to biodiversity classified through the ForestGeo Arthrpod Initiative (tradiditional taxonomy and barcoding).

Raw metabarcoding data is processed through mBrave (see parameters.doc) and downloaded as .csv file which includes organism (in ranked classification) and the number of reads for that organism on each light-trap sample for each of the two sampling nights/seasons.

To follow this pipeline, go through the scripts in order 01 - 02 - 03 etc. 

NOTE 1: SCRIPT 01 IS FOR INFORMATION PURPOSES ONLY. You should not re-run these scripts as you will overwrite the existing files (which have already been merged and cleaned).

NOTE 2: Most critical merging/cleaning steps, however, for the metabarcoding data a very meticulous process was later done in excel where 'arbitrarily' I had to decide which species was that BIN associated to, for BINs with more than one species name. This happened much more frequently than I would hope, but its a consequence of BOLD. As a general rule, I used the first species name for that BIN, wich normally corresponds to the name of first database which the barcodes were compared against (in this case, DS-BCIARTH; look at params.doc). 

NOTE 3: For the MicrobiotaProcess, some further cleaning will be required, as of now, script 03 still contains some preliminary figures (cladograms, specifically) because the data needs to be further refined to include regular expressions to define the taxonomic rank (e.g. p__;c__;o__) in order for us to be able to properly filter the data to focal groups and/or focus only on a specific rank to visualize the differences between seasons. To bypass this issue, I filteered the datasets using Metacoder and created new 'files' for each focal order.


