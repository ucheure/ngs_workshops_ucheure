#!/usr/bin/env bash
set -o nounset 
set -o errexit
set -o pipefail

MY_PERSONAL_DIRECTORY=${1}

mkdir ${HOME}/share/projects/${MY_PERSONAL_DIRECTORY}

cd ${HOME}/share/projects/${MY_PERSONAL_DIRECTORY}

mkdir ngs_project

cd ngs_project

mkdir fastq alignments variant_calls variant_annotations reports scripts temp
