SHELL=/bin/bash

env:
	python test/run_test.py --build-env --no-test workflows/mapping/Snakefile

tests:
	source activate lcdb-workflows-$(USER)-env && python -m unittest discover -s 'lcdb' -p "*_test.py"
	source activate lcdb-workflows-$(USER)-env && python lcdb/interface.py

report-test:
	source activate lcdb-workflows-$(USER)-env && python -m unittest lcdb/test/reporting_test.py
	source activate lcdb-workflows-$(USER)-env && python lcdb/reporting.py 

workflowtest: $(workflow)
	test/run_test.py . --build-env --clean --workflow=$(workflow)

dev:
	# pyinotify can be installed with: conda install -c conda-forge pyinotify
	clear
	python -m pyinotify -e IN_MODIFY -r lcdb/ -c "clear && make report-test"
