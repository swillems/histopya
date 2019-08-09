#!/bin/sh

PYTHON_COMMAND=$1

echo "################################################################################"
echo "Installating HistoPyA..."
echo "################################################################################"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $DIR
cd ..
mkdir data
mkdir projects
echo "Setting directory to" $(pwd)

if [ -z $PYTHON_COMMAND ] ; then
  PYTHON_COMMAND=python3.6
fi

if ! hash $PYTHON_COMMAND 2>/dev/null; then
  echo "Python3.6 could not be found, aborting HistoPyA installation!"
  exit
fi

echo "Creating python venv..."
$PYTHON_COMMAND -m venv venv
pwd > venv/lib/python3.6/site-packages/paths.pth

echo "Installing python dependancies..."
venv/bin/python -m pip install --upgrade pip
venv/bin/python -m pip install -r install/pip_requirements.txt

# TODO set download paths to proteomeXchange instead of dropbox
echo "Unpacking databases"
wget --output-document install/databases.tar.gz "https://filesender.belnet.be/download.php?token=bee888e9-6cc3-46b2-9610-7dcb28db3059&files_ids=48998"
tar xzvf install/databases.tar.gz -C lib
rm install/databases.tar.gz

# TODO set download paths to proteomeXchange instead of dropbox
echo "Unpacking test data"
wget --output-document install/test_lfq.tar.gz "https://filesender.belnet.be/download.php?token=bee888e9-6cc3-46b2-9610-7dcb28db3059&files_ids=48999"
tar xzvf install/test_lfq.tar.gz -C data
rm install/test_lfq.tar.gz

echo "################################################################################"
echo "Installation complete!"
echo "################################################################################"
