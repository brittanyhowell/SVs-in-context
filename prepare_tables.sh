#!/bin/bash
## Usage: Prepare reference tables and sv tables for plotting in specified windows
## Author: Brittany Howell (bh10@sanger.ac.uk)
## Date: 4th February 2019


# Script home:
scriptDIR=/Users/bh10/Documents/scripts/brie_scripts/sv-view-multi-track/
# Where the annotation tables are:
refDIR=${scriptDIR}/reference/ucsc/
dummy_4=${scriptDIR}/reference/dummy_line.txt
dummy_6=${scriptDIR}/reference/dummy_six.txt

# This should be a bed3+1 format (chrom, start, stop, svtype)
tab_SVs=/Users/bh10/Documents/data/sv-analyse/QC/sorted_all.txt
filter_SVs=${refDIR}/SVs.tmp


tab_clinvar_long=${refDIR}/clinvar_long
tab_clinvar_short=${refDIR}/clinvar_short
tab_dgv_SV=${refDIR}/dgv_SV
tab_gencode_v28=${refDIR}/gencode_v28_comprehensive
tab_refSeq_func_elements=${refDIR}/refSeq_funcElems
tab_seg_dups=${refDIR}/seg_dups
tab_simple_repeats=${refDIR}/simple_repeats

filter_clinvar_long=${refDIR}/clinvar_long_filter.tmp
filter_clinvar_short=${refDIR}/clinvar_short_filter.tmp
filter_dgv_SV=${refDIR}/dgv_SV_filter.tmp
filter_gencode_v28=${refDIR}/gencode_v28_comprehensive_filter.tmp
filter_refSeq_func_elements=${refDIR}/refSeq_funcElems_filter.tmp
filter_seg_dups=${refDIR}/seg_dups_filter.tmp
filter_simple_repeats=${refDIR}/simple_repeats_filter.tmp




#### Set the required regions
    ## These thresholds will be used to set the x limits on the final plot.
    chrom=chr2 #$1
    start=77850000 #$2
    end=78000000  #$3


#### Prepare tables

    ## refSeq_funcElems
        cat ${tab_refSeq_func_elements} | cut -f1-4,11,15 | awk -v chr="${chrom}"  -F"\t" '$1 == chr { print $0 }' | awk -v start=${start} -v end=${end} '{ if (($2 > start) && ($3 < end)) { print $0 } }' > ${filter_refSeq_func_elements}.1
        # Add dummy line
            cat ${filter_refSeq_func_elements}.1 ${dummy_6} > ${filter_refSeq_func_elements}
            rm ${filter_refSeq_func_elements}.1

    ## seg_dups
        cat ${tab_seg_dups} | awk '{print $1 "\t" $2 "\t" $3 "\tSegmental Duplication" }' | awk -v chr="${chrom}"  -F"\t" '$1 == chr { print $0 }' | awk -v start=${start} -v end=${end} '{ if (($2 > start) && ($3 < end)) { print $0 } }' > ${filter_seg_dups}.1
        # Add dummy line
        cat ${filter_seg_dups}.1 ${dummy_4} > ${filter_seg_dups}
        rm ${filter_seg_dups}.1

    ## SVs 
    cat ${tab_SVs}  | awk '{print $2 "\t" $3 "\t" $11 "\t" $9}'| awk -v chr="${chrom}"  -F"\t" '$1 == chr { print $0 }' | awk -v start=${start} -v end=${end} '{ if (($2 > start) && ($3 < end)) { print $0 } }' >  ${filter_SVs}.1 
    cat ${filter_SVs}.1 ${dummy_4} > ${filter_SVs}
    rm ${filter_SVs}.1

    ## Simple Repeats
    ## Extract chr, start, end, period_copynumber_GCcontent
    cat ${tab_simple_repeats} |  awk '{print $2 "\t" $3 "\t" $4 "\t" $6"_"$7"_"$13+$14}' |   awk -v chr="${chrom}"  -F"\t" '$1 == chr { print $0 }' | awk -v start=${start} -v end=${end} '{ if (($2 > start) && ($3 < end)) { print $0 } }' >  ${filter_simple_repeats}.1
    # Add dummy line
    cat ${filter_simple_repeats}.1 ${dummy_4} > ${filter_simple_repeats}
    rm ${filter_simple_repeats}.1

    ## Database of genomic variants
    ## Extract chr, start, end, period_copynumber_GCcontent
    cat ${tab_dgv_SV} |  awk '{print $2 "\t" $3 "\t" $4 "\t" $5"_"$11}' |   awk -v chr="${chrom}"  -F"\t" '$1 == chr { print $0 }' | awk -v start=${start} -v end=${end} '{ if (($2 > start) && ($3 < end)) { print $0 } }' >  ${filter_dgv_SV}.1
    # Add dummy line
    cat ${filter_dgv_SV}.1 ${dummy_4} > ${filter_dgv_SV}
    rm ${filter_dgv_SV}.1


## Run the plotting script:

Rscript ${scriptDIR}/load_and_view.R  ${filter_SVs} ${filter_refSeq_func_elements} ${filter_seg_dups} ${filter_simple_repeats} ${start} ${end} ${scriptDIR} ${filter_dgv_SV}

# sv.file <- read.table(file = args[1], sep = "\t", stringsAsFactors = F, header = T)
# ucsc.file <- read.table(file =args[2], sep = "\t", stringsAsFactors = F)
# segdup.file <- read.table(file = args[3], sep = "\t", stringsAsFactors = F)
# sr.file <- read.table(file = args[4], sep = "\t", stringsAsFactors = F, header = T)

# dgv <- read.table(file = args[8], sep = "\t", stringsAsFactors = F, header = F)

# start.coord <- args[5]
# end.coord <- args[6] 

# plotDIR is args[7]

rm ${filter_clinvar_long}
rm ${filter_clinvar_short}
rm ${filter_dgv_SV}
rm ${filter_gencode_v28}
rm ${filter_refSeq_func_elements}
rm ${filter_seg_dups}
rm ${filter_simple_repeats}