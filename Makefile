.DEFAULT_GOAL := help
.PHONY: deps lint format release test help

REPO_ROOT:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

deps: ## install deps (library & development)
	python3 -m pip install --upgrade pip
	python3 -m pip install -r requirements/dev.txt

lint: ## run linters, formatters for current python versions
	python3 -m flake8 redpic tests setup.py
	python3 -m pylint redpic tests setup.py

format: ## autoformat code with black and isort
	python3 -m isort redpic tests setup.py
	python3 -m black redpic tests setup.py

release: ## release package on pypi
	python3 -m setup sdist bdist_wheel
	python3 -m twine upload dist/*

test: ## run test
	python3 -m pytest tests/functional/src

help: ## Show help message
	@grep -E '^[a-zA-Z0-9 -]+:.*#'  Makefile | sort | while read -r l; do printf "\033[1;32m$$(echo $$l | cut -f 1 -d':')\033[00m:$$(echo $$l | cut -f 2- -d'#')\n"; done