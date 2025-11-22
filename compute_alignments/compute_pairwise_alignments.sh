#! /bin/bash

usage() {
cat << EOF
Pipeline for computing pairwise alignments for all species pairs in a list of given species

Usage: ${0##*/} -s SPECIES -t TARGETS [optional: -d DATA_DIR -f FORCE -@ NTHREADS -n DRY_RUN]

required (one or the other):
	-s SPECIES			Comma-separated list or file (first column) with species (genome builds) for which to compute pairwise alignments
	-t TARGETS			Comma-separated list or file with target files (e.g. a list of chain files) to be produced

optional:
	-c CREATE_PWALN_COLLECTION	Flag for creating a pairwise alingment collection file. Used as input for IPP
	-r REFERENCE			Reference species. Necessary for naming the collection file
	-q QUERY			Query species. Necessary for naming the collection file
	-d DATA_DIR			Output directory. The script will create a predefined folder structure in that directory
	   				(e.g. separate folders for alignment files, fasta files, etc.)
	-f FORCE			Force execution and overwriting existing files. Possible values: none (default) | space-separated rule names | all
	-@ NTHREADS			Number of parallel threads
	-n DRY_RUN			Flag for a dry run, i.e. only print the commands instead of executing
	-! SNAKEMAKE_ARGS  		Optional arguments for Snakemake in quotes, i.e. -! "--debug-dag --unlock -p"
EOF
}

# parse arguments (more info at https://ndench.github.io/bash/parsing-bash-flags)
create_pwaln_collection=False
ref=None
qry=None
data_dir="./ipp-data"
force="none"
nthreads="--cores 1"
script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # directory of this script
snakefile="$script_dir/Snakefile"
snakemake_args="-p --rerun-incomplete"
while getopts ":s:t:cr:q:d:f:@:!:n" OPTION; do
 	case $OPTION in
		s) species=$OPTARG ;;
		t) targets=$OPTARG ;;
		c) create_pwaln_collection=True ;;
		r) ref=$OPTARG ;;
		q) qry=$OPTARG ;;
		d) mkdir -p "$OPTARG" 2>/dev/null; data_dir=$(realpath "$OPTARG" 2>/dev/null || echo "$(cd "$(dirname "$OPTARG")" && pwd)/$(basename "$OPTARG")") ;;
		f) force=$OPTARG ;;
		@) nthreads="--cores $OPTARG" ;;
		n) nthreads="$nthreads -n" ;;
		!) snakemake_args="$snakemake_args $OPTARG" ;;
		?) usage; exit 1 ;;
 	esac
done

# throw error if not all mandatory arguments are passed
if [ -z ${species+x} ] && [ -z ${targets+x} ]; then usage && exit 1; fi
[[ $create_pwaln_collection == True ]] && [[ -z $ref || -z $qry ]] && usage && exit 1

# define which rules to force run
if [[ $force == "all" ]]; then
	force_flag="--forceall"
elif [[ $force == "none" ]]; then
	force_flag=""
else
	force_flag="--forcerun $force"
fi

# read species / targets file if one was passed
if [ -n "$species" ]; then
	[[ -e $species ]] && species_list=$(awk -vORS=, '/^[^#]/ {print $1}' $species | sed 's/,$/\n/') || species_list=$species
	species_or_targets_arg="species_list=${species_list}"
elif [ -n "$targets" ]; then
	[[ -e $targets ]] && target_list=$(awk -vORS=, '/^[^#]/ {print $1}' $targets | sed 's/,$/\n/') || target_list=$targets
	species_or_targets_arg="target_list=${target_list}"
fi

# define command
cmd="""
snakemake \
	--snakefile $snakefile \
	--config \
	$species_or_targets_arg \
	create_pwaln_collection=$create_pwaln_collection \
	ref=$ref \
	qry=$qry \
	data_dir=$data_dir \
	$force_flag \
	$snakemake_args \
	$nthreads
"""

if [[ $nthreads == *-n ]]; then echo $cmd; fi # print command during dry run
eval $cmd # run snakefile
