#!/bin/sh

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $DIR
#export PYTHONSTARTUP=.interactive_startup
venv/bin/python src/run_gui.py "$@"
