.PHONY: install test

install:
	python -m pip install -e .

test:
	PYTHONPATH=src pytest tests

