#!/usr/bin/bash

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

echo date

# make working directory, copy in configs, and move there
yaml() {
	python -c “import yaml; print(yaml.safe_load(open(‘$1’))$2)”
}

systemconfig=$(yaml $runconfig “[‘systemconfig’]”)
globalconfig=$(yaml $systemconfig “[‘globalconfig’]”)

workingdir=iterations/$itername
if [[ -d $workingdir ]]; then
	echo "${workingdir} already exists"
	exit 1
fi

echo "Making ${workingdir} and moving configs"

mkdir $workingdir
cp ${runconfig} ${workingdir}/runconfig.yaml
cp ${systemconfig} ${workingdir}/systemconfig.yaml
cp ${globalconfig} ${workingdir}/globalconfig.yaml
cd $workingdir

# run same event
date
echo "Running same event"
./../cpp/build/same_event $runconfig maxevents

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
	date
	echo "Running mixed event: ${mixlabel}"
	./parallel_run_mix.sh 4 $runconfig $mixlabel
	# for now, keep the individual files from each mixlabel juuuust in case
	mixcorrlabelbase=${mixcorrbase}_${mixlabel}
	hadd -f ${mixcorrlabelbase}.root ${mixcorrlabelbase}.rootmix_*
	rm ${mixcorrlabelbase}.rootmix_*
done

# merge all mixed event correlations
hadd -f ${mixcorrelation} ${mixcorrbase}*.root

# do the same for the skimcent5090 mixed event correlations
# since they're small, we can just merge them all together

for mixlabel in "${mixlabels_skimcent5090[@]}"
do
	date
	echo "Running mixed event: ${mixlabel}"
	./parallel_run_mix 4 $runconfig $mixlabel
done
hadd -f ${mixcorrbase}_skimcent5090.root ${mixcorrbase}_skimcent5090*rootmix*
rm ${mixcorrbase}_skimcent5090*rootmix*

echo date
exit 0
