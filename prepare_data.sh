#! /bin/bash

# fail automatically if any of the scripts fail
set -e
cd prepare
bash download.sh
python3 gen_norm.py
python3 gen_uniform.py --many