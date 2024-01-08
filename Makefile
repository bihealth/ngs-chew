.PHONY: default
default: help

.PHONY: help
help:
	@echo "Usage: make <target>"
	@echo ""
	@echo "Targets:"
	@echo "  format        Format code with black and isort"
	@echo "  black         Format code with black"
	@echo "  black-check   Check code formatting with black"
	@echo "  isort         Format imports with isort"
	@echo "  isort-check   Check imports with isort"
	@echo "  flake8        Run flake8"
	@echo "  mypy          Run mypy"
	@echo "  lint          Run all linting"
	@echo "  test          Run tests (excluding slow)"
	@echo "  test-all      Run all tests"

.PHONY: format
format: black isort

.PHONY: black
black:
	black -l 100 .

.PHONY: black-check
black-check:
	black -l 100 --check .

.PHONY: isort
isort:
	isort --force-sort-within-sections --profile=black .

.PHONY: isort-check
isort-check:
	isort --force-sort-within-sections --profile=black --check .

.PHONY: flake8
flake8:
	flake8

.PHONY: mypy
mypy: export MYPYPATH=stubs
mypy:
	mypy chew tests

.PHONY: lint
lint: flake8 isort-check black-check mypy

.PHONY: test
test:
	pytest . -m "not slow"

.PHONY: test-all
test-all:
	pytest .