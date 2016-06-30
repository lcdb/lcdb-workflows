SHELL=/bin/bash

env:
	python test/run_test.py --build-env --no-test workflows/mapping/Snakefile

tests:
	source activate lcdb-workflows-$(USER)-env && python -m unittest discover -s 'lcdb' -p "*_test.py"
	source activate lcdb-workflows-$(USER)-env && python lcdb/interface.py

workflowtest: $(workflow)
	test/run_test.py . --build-env --clean --workflow=$(workflow)

osx-dev:
	clear
	fswatch -o -0 -r lcdb | xargs -0 -n 1 bash -c "clear && make tests"
