
tests:
	python -m unittest discover -s 'lcdb' -p "*_test.py"

workflowtest: $(workflow)
	test/run_test.py . --build-env --clean --workflow=$(workflow)
