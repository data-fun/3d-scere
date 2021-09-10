default: help


run: ## Run app on port 8050
	@test "${CONDA_DEFAULT_ENV}" = "dashboard" && echo "Conda env ${CONDA_DEFAULT_ENV} found" || { echo "Conda env not activated"; exit 1; }
	python app.py
.PHONY: run


run-gunicorn:  ## Run app with gunicorn on port 8000
	@test "${CONDA_DEFAULT_ENV}" = "dashboard" && echo "Conda env ${CONDA_DEFAULT_ENV} found" || { echo "Conda env not activated"; exit 1; }
	gunicorn --config gunicorn.py app:server
.PHONY: run-gunicorn


help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: help