
tests:
	python -m unittest discover -s 'lcdb' -p "*_test.py"

workflowtest: $(workflow)
	test/run_test.py . --build-env --clean --workflow=$(workflow)

osx-dev:
	clear
	fswatch -o -0 -r lcdb | xargs -0 -n 1 bash -c "clear && make tests"
