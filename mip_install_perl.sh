#!/usr/bin/env/bash
set -e

usage="$(basename "$0") [-e <environment_name> -p <conda_path> -x <existing_conda_env> -d <install dev tools>]

where:
    -h  Show this help text
    -e  Conda environment name [string]
    -p  Conda path [string]
    -x  Install in already existing conda environment [flag]
    -d  Install developer tools"

unset OPTARG
unset OPTIND

CONDA_PATH="$HOME/miniconda3"
ENV_NAME='mip_rd-dna'
EXISTING_ENV='false'
DEV='false'

while getopts ':he:p:xd' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    e) ENV_NAME=$OPTARG
       ;;
    p) CONDA_PATH=$OPTARG
       ;;
    x) EXISTING_ENV='true'
       ;;
    d) DEV='true'
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND-1))

## Enable trap for signal(s) DEBUG
trap 'previous_command="$BASH_COMMAND"' DEBUG

error() {

        local program="$1"

        local return_code="$2"

        ## Display error message and exit
        echo "${program}: ${return_code}: Unknown Error - ExitCode=$return_code" 1>&2
        exit 1
}

## Enable trap for signal(s) ERR
trap '$(error "$previous_command" "$?")' ERR

## Source conda
source "$CONDA_PATH"/etc/profile.d/conda.sh

## Create or install conda env
if [ "$EXISTING_ENV" = true ] && [ -d ${CONDA_PATH}/envs/${ENV_NAME} ]
then
  conda install --name "$ENV_NAME" --yes  -c conda-forge libgcc-ng gxx_linux-64=7.3.0 python=3.7.7
elif [ -d ${CONDA_PATH}/envs/${ENV_NAME} ]
then
    echo "Environment already exists. Please supply flag -x to command"
    exit 1
else
  conda create --name "$ENV_NAME" --yes  -c conda-forge libgcc-ng gxx_linux-64=7.3.0 python=3.7.7
fi

conda install --name "$ENV_NAME" --yes -c bioconda -c conda-forge perl=5.26.2=h516909a_1006 perl-app-cpanminus perl-log-log4perl perl-moosex-app perl-file-copy-recursive perl-timedate perl-set-intervaltree perl-params-validate

conda activate "$ENV_NAME"

cpanm --installdeps .

if [ "$DEV" = true ]; then 
    
    ## Install dev modules from cpan
    cpanm Perl::Tidy@20200110 Perl::Critic@1.138 Data::Printer@0.40

    ## Install yamllint 
    conda install --name "$ENV_NAME" --yes -c conda-forge yamllint=1.20.0
fi

conda deactivate

