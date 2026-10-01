.PHONY: test lint typecheck build

test:
	python -m unittest discover -s tests -v
lint:
	ruff check src tests
typecheck:
	mypy src
build:
	python -m build
