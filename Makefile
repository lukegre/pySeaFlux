.PHONY: help

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-12s\033[0m %s\n", $$1, $$2}'

clean:  ## Clean up the temporary build and test directories
	@echo Cleaning up build and test dirs
	rm -rf *.egg* .eggs .coverage .pytest_cache build dist site

lint: clean  ## Perform checks on your syntax (also cleans up builds)
	@echo Linting with black and flake
	black --check .
	flake8 .

test:  ## Run tests located in the tests folder and return a coverage report
	@echo Testing
	pytest --cov=seaflux ./tests/

build:  ## Build the files required to submit to Pypi (will first test)
	@echo Testing and then building
	python setup.py sdist bdist_wheel

dev_env: clean  ## set up your python for developing
	conda create -y --name seaflux-dev python==3.7
	conda install -y -c conda-forge --name seaflux-dev --file requirements-dev.txt

docs_deploy:  ## Comples and submits the documentation to GitHub pages
	mkdocs gh-deploy
