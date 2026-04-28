#!/bin/bash

# Locate project root and work from the test data directory
_F="${BASH_SOURCE[0]:-$0}"
_SCRIPTS_DIR="$(cd "$(dirname "$_F")" && pwd)"
MY_ROOT="$(cd "${_SCRIPTS_DIR}/.." && pwd)"
cd "${_SCRIPTS_DIR}" || exit 1

echo "Welcome the Semi Analytic Galaxy Evolution model!"
echo "Always keep up to date by visiting our main Github page."
echo "https://github.com/sage-home"
echo ""

echo "Making the 'output/millennium' directory to hold the model output."

# By using `-p`, `mkdir` won't fail if the directory already exists.
#mkdir -p output/millennium
echo "Done."
echo ""

echo "Creating the 'input/millennium/trees' directory."
mkdir -p ../RAW/trees/minimillennium/trees
echo "Done."
echo ""

pushd ../RAW/trees/minimillennium/trees

if [ ! -f trees_063.7 ]; then

    # To download the trees, we use either `wget` or `curl`. By default, we want to use `wget`.
    # However, if it isn't present, we will use `curl` with a few parameter flags.
    echo "First checking if either 'wget' or 'curl' are present in order to download trees."

    clear_alias=0
    command -v wget

    if [[ $? != 0 ]]; then
        echo "'wget' is not available. Checking if 'curl' can be used."
        command -v curl

        if [[ $? != 0 ]]; then
            echo "Neither 'wget' nor 'curl' are available to download the Mini-Millennium trees."
            echo "Please install one of these to download the trees."
            exit 1
        fi

        echo "Using 'curl' to download trees."

        # `curl` is available. Alias it to `wget` with some options.
        alias wget="curl -L -O -C -"

        # We will need to clear this alias up later.
        clear_alias=1
    else
        echo "'wget' is present. Using it!"
    fi

    # Now that we have confirmed we have either `wget` or `curl`, proceed to downloading the trees.
    echo "Fetching Mini-Millennium trees."
    wget "https://www.dropbox.com/s/l5ukpo7ar3rgxo4/mini-millennium-treefiles.tar?dl=0" -O "mini-millennium-treefiles.tar"
    if [[ $? != 0 ]]; then
        echo "Could not download tree files from the Manodeep Sinha's Dropbox."
        echo "Please open an issue on the 'sage-model' repository and we will assist ASAP :)"
        exit 1
    fi
    echo "Done."
    echo ""

    # If we used `curl`, remove the `wget` alias.
    if [[ $clear_alias == 1 ]]; then
        unalias wget
    fi

    tar -xvf mini-millennium-treefiles.tar
    if [[ $? != 0 ]]; then
        echo "Could not untar the Mini-Millennium tree files."
        echo "Please open an issue on the 'sage-model' repository and we will assist ASAP :)"
        exit 1
    fi
    echo "Done."
    echo ""

    rm -f mini-millennium-treefiles.tar
    echo "Mini-Millennium trees successfully gathered and placed into '${PWD}'"

else
    echo "Mini-Millennium trees already present in 'input/millennium/trees'."
fi
