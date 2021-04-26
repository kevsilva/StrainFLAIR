#!/bin/bash
#*****************************************************************************
#  StrainFLAIR (STRAIN-level proFiLing using vArIation gRaph) is a tool for 
#  strain identification and quantification that uses variation graph 
#  representation of genes sequences. The input is a collection of complete 
#  genomes, draft genomes or metagenome-assembled genomes from which genes 
#  will be predicted. StrainFLAIR is sub-divided into two main parts: first, 
#  an indexing step that stores clusters of reference genes into variation 
#  graphs, and then, a query step using mapping of metagenomic reads to 
#  infere strain-level abundances in the queried sample.
#
#  Authors: Kevin Da Silva - Pierre Peterlongo
#
#  AGPL-3.0 License
#*****************************************************************************


#-----------------------------------------------------------------------------
# Bash colors 
#-----------------------------------------------------------------------------
red=`tput setaf 1`
green=`tput setaf 2`
yellow=`tput setaf 3`
cyan=`tput setaf 6`
bold=`tput bold`
reset=`tput sgr0`

#-----------------------------------------------------------------------------
# Various functions
#-----------------------------------------------------------------------------
die() {
    printf '%s\n' "$1" >&2
    exit 1
}


#-----------------------------------------------------------------------------
# Check Dependencies
#-----------------------------------------------------------------------------
check_bin(){
    which $1 > /dev/null
    if [ $? -ne 0 ]
    then
        echo "$red $0 requires $1$reset"
        exit 1
    fi
}
#
check_bin "prodigal"
check_bin "cd-hit-est"
check_bin "minimap2"
check_bin "seqwish"
check_bin "vg"

#-----------------------------------------------------------------------------
# Variables common to index and query
#-----------------------------------------------------------------------------
Tstart="$(date +%s)"
version="0.0.2"


#-----------------------------------------------------------------------------
# Index or Query
#-----------------------------------------------------------------------------
if [ $# -lt 1 ]
then
    echo "${red}Usage: $0 [query/index]$reset"
    exit 1
fi

main=$1
if [[ $main !=  "index" && $main != "query" ]]
then
    echo "${red}Usage: $0 [query/index]$reset"
    exit 1
fi



#-----------------------------------------------------------------------------
# INDEX
#-----------------------------------------------------------------------------
if [[ $main ==  "index" ]]
then
    input_data=""               # sequences to be indexed. Either a fasta file containing one genome per line, or a file of file indicating for each line the absolute path to a fasta file.
    len_extension=75            # Len of the sequences on the left and right part of each predicted gene, added to the indexation graph.
    directory_output=""         # Name of the directory in which all files are output.
    
    cdhit_c=0.95                # CD-HIT Sequence identity threshold for clustering
    cdhit_aS=0.90               # CD-HIT Alignment coverage for the shorter sequence
    cdhit_g=1                   # CD-HIT Clustering algorithm (sequence clustered to the first vs the best cluster that meet the threshold)
    cdhit_d=0                   # CD-HIT Length of description in .clstr file
    cdhit_M=0                   # CD-HIT Memory limit
    cdhit_T=0                   # CD-HIT Number of threads
    cdhit_G=0                   # CD-HIT Local or global sequence identity


    function help_index {
        echo " ******************"
        echo " *** HELP index ***"
        echo " ******************"
        echo "$0 index: indexation of the genes of a set of bacterial genomes."
        echo "Version "$version
        echo "Usage: $0 index -i file_of_files -o directory_output_name [OPTIONS]"
        echo -e "\nMANDATORY"
        echo -e "\t -i <file name of a file of file(s) or of a fasta file>"
        echo -e "\t    In case of a fasta file: each fasta input line is considered as a genome"
        echo -e "\t    In case of a .txt file: each line contains a fasta file, and each of these fasta is considered as a genome. In this case a genome can span several line, for instance for perfectly assembled genomes"
        echo -e "\t -o <directory_output_name>. This directory must not exist. It is created by the program. All results are stored in this directory"

        echo -e "\nOPTIONS"
	echo -e "\ngenes prediction:"
        echo -e "\t -l value <int value>. Set the length of the sequences on the left and right part of each predicted gene, added to the indexation graph. [default: 75]"
	echo -e "\ngenes clustering:"
        echo -e "\t -c value <float value>. Sequence identity threshold [default: 0.95]"
	echo -e "\t -aS value <float value>. Alignment coverage for the shorter sequence [default: 0.90]"
	echo -e "\t -g 0 or 1. [default: 1]"
	echo -e "\t    If set to 0, a sequence is clustered to the first cluster that meet the threshold."
	echo -e "\t    If set to 1, a sequence is clustered to the most similar cluster that meet the threshold."
	echo -e "\t -d value <int value>. Length of description in .clstr file [default: 0]"
	echo -e "\t -M value <int value>. Memory limit (in MB) ; 0 for unlimited. [default: 0]"
	echo -e "\t -T value <int value>. Number of threads ; with 0, all CPUs will be used. [default: 0]"
	echo -e "\t -G 0 or 1. [default: 0]"
	echo -e "\t    If set to 0, use local sequence identity."
	echo -e "\t    If set to 1, use global sequence identity."
	echo -e "\nglobal:"
        echo -e "\t -h Prints this message and exit\n"
    
        echo "Any further question: read the readme file or contact the development team"
    }

    shift # pass "index"
    while :; do
        case $1 in
        -i) 
            if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
                input_data=$2
                shift
            else
                die 'ERROR: "'$1'" option requires a non-empty option argument.'
            fi
            ;;
        -o) 
            if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
                directory_output=$2
                shift
            else
                die 'ERROR: "'$1'" option requires a non-empty option argument.'
            fi
            ;;

        -l) 
            if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
                len_extension=$2
                shift
            else
                die 'ERROR: "'$1'" option requires a non-empty option argument.'
            fi
            ;;
        -h|-\?|--help)
            help_index
            exit 
            ;;

        -?*)
            printf 'WARN: Unknown option (exit): %s\n' "$1" >&2
            exit 1
            ;;

        :)
            echo "Option $1 requires an argument." >&2
            exit 1
            ;;
            --)              # End of all options.
                shift
                break
                ;;
            -?*)
                printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
                ;;
            *)               # Default case: No more options, so break out of the loop.
                break
        esac
        shift
    done
    
    # --------------
    # CHECK OPTIONS
    # --------------
    if [ -z "$input_data" ]; then
        echo "$red Error: You must provide an input file (-i)"
        help_index
        echo $reset
        exit 1
    fi
    
    if [ -z "${directory_output}" ]; then
        echo "$red Error: You must provide an output directory name (-o)"
        help_index
        echo $reset
        exit 1
    fi
    
    if [ -f ${directory_output} ] || [ -d ${directory_output} ]; then
        echo "$red Error: ${directory_output} already exists"
        help_index
        echo $reset
        exit 1
    fi
    mkdir ${directory_output}
    
    # --------------
    # RECAP OPTIONS
    # --------------
    echo -e "$yellow"
    echo "INDEXATION OPTIONS:"
    echo "-Index file ${input_data}"
    echo "-Output results in directory ${directory_output}"
    echo "-Extension len ${len_extension}"
    echo "-cd-hit-est options:"
    echo "  c=${cdhit_c}"
    echo "  aS=${cdhit_aS}"
    echo "  g=${cdhit_g}"
    echo "  d=${cdhit_d}"
    echo "  M=${cdhit_M}"
    echo "  T=${cdhit_T}"
    echo "  G=${cdhit_G}"
    echo -e "$reset"
    
    
    # --------------
    # INDEX CREATION
    # --------------
    
    # Gene Prediction
    echo "${yellow}GENE PREDICTION$reset"

    cmd="mkdir ${directory_output}/predicted_genes"
    echo $green$cmd$cyan
    $cmd
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with the gene prediction$reset"
        exit 1
    fi

    cmd="genes_prediction -s ${input_data} -o ${directory_output}/predicted_genes -l ${len_extension}"
    echo $green$cmd$cyan
    T="$(date +%s)"
    $cmd
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with the gene prediction$reset"
        exit 1
    fi
    T="$(($(date +%s)-T))"
    echo "$yellow Gene Prediction time in seconds: ${T}$reset"
    
    
    
    # Gene Clustering
    echo "${yellow}GENE CLUSTERING$reset"
    
    cmd="mkdir ${directory_output}/clusters"
    echo $green$cmd$cyan
    $cmd
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with the gene clustering$reset"
        exit 1
    fi
    
    cmd="cd-hit-est -i ${directory_output}/predicted_genes/all_genes.fasta -o ${directory_output}/clusters/all_genes_clusters -c ${cdhit_c} -aS ${cdhit_aS} -g ${cdhit_g} -d ${cdhit_d} -M ${cdhit_M} -T ${cdhit_T} -G ${cdhit_G}"
    echo $green$cmd$cyan
    T="$(date +%s)"
    $cmd
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with the gene clustering$reset"
        exit 1
    fi
    T="$(($(date +%s)-T))"
    echo "$yellow Gene Clustering time in seconds: ${T}$reset"
    
    
    
    # Graphs Construction
    echo "${yellow}GRAPHS CONSTRUCTION$reset"
    
    cmd="mkdir ${directory_output}/graphs"
    echo $green$cmd$cyan
    $cmd
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with the graphs construction$reset"
        exit 1
    fi
    
    cmd="graphs_construction -s ${directory_output}/predicted_genes/all_genes_extended.fasta -c ${directory_output}/clusters/all_genes_clusters.clstr -o ${directory_output}/graphs"
    echo $green$cmd$cyan
    T="$(date +%s)"
    $cmd
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with the graphs construction$reset"
        exit 1
    fi
    T="$(($(date +%s)-T))"
    echo "$yellow Graphs Construction time in seconds: ${T}$reset"
    
    
    
    
    # Graphs Concatenation
    echo "${yellow}GRAPHS CONCATENATION$reset"
    cmd="concat_graphs -i ${directory_output}/graphs -s 1000"
    echo $green$cmd$cyan
    T="$(date +%s)"
    $cmd
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with the Graphs Concatenation$reset"
        exit 1
    fi
    T="$(($(date +%s)-T))"
    echo "$yellow Graphs Concatenation time in seconds: ${T}$reset"
    
    # From vg format to gfa format. 
    echo "${yellow}VG TO GFA FORMAT$reset"
    cmd="vg view ${directory_output}/graphs/all_graphs.vg" #> ${directory_output}/graphs/all_graphs.gfa
    echo "$green$cmd > ${directory_output}/graphs/all_graphs.gfa $cyan"
    T="$(date +%s)"
    $cmd > ${directory_output}/graphs/all_graphs.gfa
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with the Graphs format modification$reset"
        exit 1
    fi
    T="$(($(date +%s)-T))"
    echo "$yellow Graphs format modification time in seconds: ${T}$reset"
    
    # Graph Index.
    T="$(date +%s)"
    cmd="vg prune ${directory_output}/graphs/all_graphs.vg" #| vg index -g ${directory_output}/graphs/all_graphs.gcsa -
    echo "$green$cmd | vg index -g ${directory_output}/graphs/all_graphs.gcsa - $cyan"
    $cmd | vg index -g ${directory_output}/graphs/all_graphs.gcsa -
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with the Graphs Indexation (GCSA file)$reset"
        exit 1
    fi
    cmd="vg index -x ${directory_output}/graphs/all_graphs.xg ${directory_output}/graphs/all_graphs.vg"
    $cmd
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with the Graphs Indexation (XG file)$reset"
        exit 1
    fi
    cmd="vg snarls ${directory_output}/graphs/all_graphs.vg" # > ${directory_output}/graphs/all_graphs.snarls
    echo "$green$cmd > ${directory_output}/graphs/all_graphs.snarls $cyan"
    $cmd > ${directory_output}/graphs/all_graphs.snarls
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with the Graphs Indexation (SNARLS file)$reset"
        exit 1
    fi
fi # END INDEX









#-----------------------------------------------------------------------------
# QUERY
#-----------------------------------------------------------------------------
if [[ $main ==  "query" ]]
then
    graph=""                    # variation graph directory and name (no extensions)
    reads1=""                   # sequencing reads (pair 1)
    reads2=""                   # sequencing reads (pair 2)
    vg_t=16                     # number of threads for mapping
    clusters_data=""            # pickle containing the dictionary for the clusters
    directory_output=""         # Name of the directory in which all files are output.
    filename_output=""          # Name of the output files.
    threshold=0.5               # threshold on the proportion of detected genes
	
    function help_query {
        echo " ******************"
        echo " *** HELP query ***"
        echo " ******************"
        echo "$0 query: query of reads on a variation graph."
        echo "Version "$version
        echo "Usage: $0 query -g graph -f1 reads1 -f2 reads2 -t threads -p dict_clusters -d output_directory_name -o output_files_name [OPTIONS]"
        echo -e "\nMANDATORY"
        echo -e "\t -g <file name of a graph (no format)>"
	echo -e "\t -f1 <single-end reads or pair1 of paired-end reads (fastq or fastq.gz)>"
        echo -e "\t -t <number of threads to use for mapping>"
	echo -e "\t -p <pickle file containing the dictionary of clusters and their genes>"
        echo -e "\t -d <output_directory_name>. Name of the directory in which all files are output."
        echo -e "\t -o <output_files_name>. Specific name for the output files."

        echo -e "\nOPTIONS"
        echo -e "\t -f2 <pair2 of paired-end reads (fastq or fastq.gz)>"
        echo -e "\t -s value <float value between [0-1]>. Set the threshold on proportion of detected specific genes. [default=0.5]"
        echo -e "\t -h Prints this message and exit\n"
    
        echo "Any further question: read the readme file or contact the development team"
    }

    shift # pass "query"
    while :; do
        case $1 in
        -g) 
            if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
                graph=$2
                shift
            else
                die 'ERROR: "'$1'" option requires a non-empty option argument.'
            fi
            ;;
        -f1) 
            if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
                reads1=$2
                shift
            else
                die 'ERROR: "'$1'" option requires a non-empty option argument.'
            fi
            ;;
        -f2) 
            if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
                reads2=$2
                shift
            else
                die 'ERROR: "'$1'" option requires a non-empty option argument.'
            fi
            ;;
        -t) 
            if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
                vg_t=$2
                shift
            else
                die 'ERROR: "'$1'" option requires a non-empty option argument.'
            fi
            ;;
        -p) 
            if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
                clusters_data=$2
                shift
            else
                die 'ERROR: "'$1'" option requires a non-empty option argument.'
            fi
            ;;
	-d) 
            if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
                directory_output=$2
                shift
            else
                die 'ERROR: "'$1'" option requires a non-empty option argument.'
            fi
            ;;
	-o) 
            if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
                filename_output=$2
                shift
            else
                die 'ERROR: "'$1'" option requires a non-empty option argument.'
            fi
            ;;
       -s) 
            if [ "$2" ] && [ ${2:0:1} != "-" ] ; then # checks that there exists a second value and its is not the start of the next option
                threshold=$2
                shift
            else
                die 'ERROR: "'$1'" option requires a non-empty option argument.'
            fi
            ;;
        -h|-\?|--help)
            help_query
            exit 
            ;;

        -?*)
            printf 'WARN: Unknown option (exit): %s\n' "$1" >&2
            exit 1
            ;;

        :)
            echo "Option $1 requires an argument." >&2
            exit 1
            ;;
            --)              # End of all options.
                shift
                break
                ;;
            -?*)
                printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
                ;;
            *)               # Default case: No more options, so break out of the loop.
                break
        esac
        shift
    done
    
    # --------------
    # CHECK OPTIONS
    # --------------
    if [ -z "${graph}" ]; then
        echo "$red Error: You must provide a graph (no format) (-g)"
        help_query
        echo $reset
        exit 1
    fi
    
    if [ -z "${reads1}" ]; then
        echo "$red Error: You must provide sequencing reads (fastq or fastq.gz format) (-f1)"
        help_query
        echo $reset
        exit 1
    fi

    if [ -z "${clusters_data}" ]; then
        echo "$red Error: You must provide a clusters dictionary (pickle file) (-p)"
        help_query
        echo $reset
        exit 1
    fi

    if [ -z "${directory_output}" ]; then
        echo "$red Error: You must provide output directory name (-o)"
        help_query
        echo $reset
        exit 1
    fi
    
    if [ ! -d ${directory_output} ]; then
        mkdir ${directory_output}
    fi
    if [ ! -d ${directory_output}/mapping ]; then
        mkdir ${directory_output}/mapping
    fi
    if [ ! -d ${directory_output}/results ]; then
        mkdir ${directory_output}/results
    fi

    if [ -f ${directory_output}/mapping/mapping_${filename_output}.gamp ] || [ -f ${directory_output}/mapping/mapping_${filename_output}.gamp ] || [ -f ${directory_output}/results/genelevel_${filename_output}.csv ] || [ -f ${directory_output}/results/strainsprofile_${filename_output}.csv ] ; then
        echo "$red Error: ${filename_output} is already used!"
        help_query
        echo $reset
        exit 1
    fi
    
    # --------------
    # RECAP OPTIONS
    # --------------
    echo -e "$yellow"
    echo "QUERY OPTIONS:"
    echo "-Graph file ${graph}"
    echo "-Reads file ${reads1} ${reads2}"
    echo "-Number of threads ${vg_t}"
    echo "-Clusters dictionary ${clusters_data}"
    echo "-Output results in directory ${directory_output}"
    echo "-Threshold on proportion of detected genes ${threshold}"
    echo -e "$reset"
    
    # --------------
    # MAPPING
    # --------------
    
    echo "${yellow}MAPPING$reset"
    if [ -z "${reads2}" ]
    then
    	cmd="vg mpmap -x ${graph}.xg -g ${graph}.gcsa -s ${graph}.snarls -f ${reads1} -t ${vg_t} -M 10 -m -L 0" # > ${directory_output}/mapping/mapping_${filename_output}.gamp"
    else
    	cmd="vg mpmap -x ${graph}.xg -g ${graph}.gcsa -s ${graph}.snarls -f ${reads1} -f ${reads2} -t ${vg_t} -M 10 -m -L 0" # > ${directory_output}/mapping/mapping_${filename_output}.gamp"
    fi
    echo "$green$cmd > ${directory_output}/mapping/mapping_${filename_output}.gamp $cyan"
    T="$(date +%s)"
    $cmd > ${directory_output}/mapping/mapping_${filename_output}.gamp
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with the mapping$reset"
        exit 1
    fi
    T="$(($(date +%s)-T))"
    echo "$yellow Mapping computation time in seconds: ${T}$reset"

    echo "${yellow}GAMP TO JSON FORMAT$reset"
    cmd="vg view -j -K ${directory_output}/mapping/mapping_${filename_output}.gamp"
    echo "$green$cmd > ${directory_output}/mapping/mapping_${filename_output}.json $cyan"
    T="$(date +%s)"
    $cmd > ${directory_output}/mapping/mapping_${filename_output}.json
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with the JSON conversion$reset"
        exit 1
    fi
    T="$(($(date +%s)-T))"
    echo "$yellow JSON conversion time in seconds: ${T}$reset"

    
    # --------------
    # GENE-LEVEL
    # --------------
    
    echo "${yellow}GENE-LEVEL ABUNDANCES$reset"
    cmd="json2csv -g ${graph}.gfa -m ${directory_output}/mapping/mapping_${filename_output}.json -p ${clusters_data} -o ${directory_output}/results/genelevel_${filename_output}"
    echo $green$cmd$cyan
    T="$(date +%s)"
    $cmd
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with the gene-level table generation$reset"
        exit 1
    fi
    T="$(($(date +%s)-T))"
    echo "$yellow Gene-level abundances computation time in seconds: ${T}$reset"

    # --------------
    # STRAIN-LEVEL
    # --------------
    
    echo "${yellow}STRAIN-LEVEL ABUNDANCES$reset"
    cmd="compute_strains_abundance -i ${directory_output}/results/genelevel_${filename_output}.csv -o ${directory_output}/results/strainsprofile_${filename_output} -t ${threshold}"
    echo $green$cmd$cyan
    T="$(date +%s)"
    $cmd
    if [ $? -ne 0 ]
    then
        echo "$red there was a problem with the strain-level table generation$reset"
        exit 1
    fi
    T="$(($(date +%s)-T))"
    echo "$yellow Strain-level abundances computation time in seconds: ${T}$reset"
    
fi # END QUERY
