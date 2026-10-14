PYTHON ?= python3
PIPELINE := scripts/run_gyre/pipeline/run_gyre.py
MODEL ?= models/my_star.mesa
OUTPUT ?= results/run_gyre/my_star

.PHONY: help check dry-run run notebook

help:
	@printf '%s\n' \
	  'make check                         validate the Python entry point' \
	  'make dry-run MODEL=path/to/model  generate inputs without running GYRE' \
	  'make run MODEL=path/to/model      run GYRE with the documented defaults' \
	  'make notebook                     open the guided Jupyter notebook'

check:
	$(PYTHON) -m py_compile scripts/run_gyre/pipeline/*.py
	$(PYTHON) $(PIPELINE) --help >/dev/null

dry-run:
	$(PYTHON) $(PIPELINE) $(MODEL) --output-dir $(OUTPUT) --dry-run --overwrite

run:
	$(PYTHON) $(PIPELINE) $(MODEL) --output-dir $(OUTPUT) --overwrite

notebook:
	jupyter lab scripts/run_gyre/notebook/run_gyre.ipynb
