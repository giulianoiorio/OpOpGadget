PYTHON_DIR=/Users/Giuliano/anaconda/envs

install:
	$(PYTHON_DIR)/py33/bin/./python setup.py install

33:
	$(PYTHON_DIR)/py33/bin/./python setup.py install

clean-pyc:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -f -r {} +

clean-build:
	rm -f -r build/
	rm -f -r dist/
	rm -f -r *.egg-info

clean:
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -f -r {} +
	rm -f -r build/
	rm -f -r dist/
	rm -f -r *.egg-info

.PHONY: clean-pyc clean-build clean 33

help:
	@echo "    clean-pyc"
	@echo "        Remove python artifacts."
	@echo "    clean-build"
	@echo "        Remove build artifacts."
	@echo "    clean"
	@echo "        Remove build and python artifacts."