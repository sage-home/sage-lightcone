#!/bin/bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r sage-model/requirements.txt
pip install mpi4py
