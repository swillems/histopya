#!/bin/sh

PYTHON_COMMAND=$1

echo "################################################################################"
echo "Installating HistoPyA..."
echo "################################################################################"


DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $DIR
cd ..
mkdir data
mkdir lib/databases
echo "Setting directory to" $(pwd)


echo "################################################################################"
if [ -z $PYTHON_COMMAND ] ; then
  PYTHON_COMMAND=python3.6
fi

if hash $PYTHON_COMMAND 2>/dev/null; then
  echo "Found python3.6 at" $(command -v PYTHON_COMMAND)
else
  echo "Python3.6 could not be found, aborting HistoPyA installation!"
  exit
fi

echo "Creating python venv..."
$PYTHON_COMMAND -m venv venv

echo "Setting python venv path to HistoPyA..."
pwd > venv/lib/python3.6/site-packages/paths.pth

echo "Upgrading venv pip..."
venv/bin/python -m pip install --upgrade pip

echo "Installing python dependancies..."
venv/bin/python -m pip install -r install/pip_requirements.txt


echo "################################################################################"
echo "Unpacking databases"
tar xzvf install/databases.tar.gz -C lib/databases

echo "Unpacking test data"
tar xzvf install/test_data.tar.gz -C data


echo "################################################################################"
echo "Installation complete!"
echo "################################################################################"
