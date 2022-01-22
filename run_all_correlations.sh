#!/bin/bash

set -o errexit # exit on error

# make sure that modules are loaded and directories are set
if [[ -z "$ROOT_DIR" ]]; then
	echo "ROOT not loaded; try setupgj first"
	exit 1
fi

if [[ -z "$HDF5_DIR" ]]; then
	echo "HDF5 not loaded; try setupgj first"
	exit 1
fi

if [[ -z "$HDF5_USE_FILE_LOCKING" ]]; then
	echo "\$HDF5_USE_FILE_LOCKING not loaded; try setupgj first"
	exit 1
fi

if [[ -z "$OMP_NUM_THREADS" ]]; then
	echo "\$OMP_NUM_THREADS not loaded; try setupgj first"
	exit 1
fi

if [[ -z "$YAMLCPP" ]]; then
	echo "\$YAMLCPP not defined"
	exit 1
fi

# read input parameters
runconfig=$1
itername=$2
maxevents=$3

if [[ -z $runconfig ]]; then
	echo "Usage: ./run_all_correlations.sh [runconfig] [itername] [maxevents (optional)]"
	exit 1
fi

if [[ -z $itername ]]; then
	echo "Usage: ./run_all_correlations.sh [runconfig] [itername] [maxevents (optional)]"
	exit 1
fi

# by default, infinite max events, which is just this very large number
if [[ -z $maxevents ]]; then
	maxevents=999999999999999
fi

workingdir=iterations/$itername
if [[ -d $workingdir ]]; then
    echo "${workingdir} already exists"
    exit 1
fi
mkdir $workingdir

# this yaml parsing solution comes from https://stackoverflow.com/a/47791935
yaml() {
    python3 -c "import yaml; print(yaml.safe_load(open('$1'))$2)"
}

main() {
	date
	echo "runconfig: ${runconfig}, itername: ${itername}, maxevents: ${maxevents}"

	# make working directory, copy in configs, and move there
	systemconfig=$(yaml $runconfig "['systemconfig']")
	globalconfig=$(yaml $systemconfig "['globalconfig']")

	gjdir=$(pwd)
	echo "Moving configs"

	mkdir $workingdir/config
	mkdir $workingdir/root
	cp ${runconfig} ${workingdir}/config/
	cp ${systemconfig} ${workingdir}/config/
	cp ${globalconfig} ${workingdir}/config/
	# also make a copy with generic names
	cp ${runconfig} ${workingdir}/config/runconfig.yaml
	cp ${systemconfig} ${workingdir}/config/systemconfig.yaml
	cp ${globalconfig} ${workingdir}/config/globalconfig.yaml

	cp example-gj-correlations.ipynb ${workingdir}/gj-correlations-${itername}.ipynb
	responsematrix=$(yaml $runconfig "['filelists']['correlations']['responsematrix']")
	cp $responsematrix ${workingdir}/${responsematrix}

	cd $workingdir

	# run same event
	echo "Running same event"
	$gjdir/cpp/build/same_event $runconfig $maxevents

	# define mixed event labels to run
	mixlabels=()
	mixlabels+=("18q_int7_1")
	mixlabels+=("18q_int7_2")
	mixlabels+=("18q_int7_3")

	mixlabels_skimcent5090=()
	mixlabels_skimcent5090+=("skimcent5090_18q_int7_1")
	mixlabels_skimcent5090+=("skimcent5090_18q_int7_2")
	mixlabels_skimcent5090+=("skimcent5090_18q_int7_3")

	# get the name of the mixed event correlation from the config
	mixcorrelation=$(yaml $runconfig "['filelists']['correlations']['mixedevent']")
	mixcorrbase=${mixcorrelation%.*}

	for mixlabel in "${mixlabels[@]}"
	do
		echo
		date
		echo "Running mixed event: ${mixlabel}"
		$gjdir/parallel_run_mix.sh 4 $runconfig $mixlabel $maxevents $gjdir
		# for now, keep the individual files from each mixlabel juuuust in case
		mixcorrlabelbase=${mixcorrbase}_${mixlabel}
		hadd -f ${mixcorrlabelbase}.root ${mixcorrlabelbase}.rootmix_*
		rm ${mixcorrlabelbase}.rootmix_*
	done

	# merge all mixed event correlations
	echo
	hadd -f ${mixcorrelation} ${mixcorrbase}*.root

	# do the same for the skimcent5090 mixed event correlations
	# since they're small, we can just merge them all together
	for mixlabel in "${mixlabels_skimcent5090[@]}"
	do
		echo
		date
		echo "Running mixed event: ${mixlabel}"
		$gjdir/parallel_run_mix.sh 4 $runconfig $mixlabel $maxevents $gjdir
	done
	
	echo
	hadd -f ${mixcorrbase}_skimcent5090.root ${mixcorrbase}_skimcent5090*rootmix*
	rm ${mixcorrbase}_skimcent5090*rootmix*

	date
}

main | tee "${workingdir}/output.log"
exit 0
