.PHONY: checks fmt lint type tests

# Synchronize environment and build extension
sync-build: uv-sync build-ext

# Synchronize Python environment and dependencies
uv-sync:
	uv sync

# Build the C++ extension module
build-ext:
	python setup.py build_ext --inplace

# Run pre-commit checks
checks:
	uvx pre-commit run --all-files

# Format code
fmt:
	uv run ruff format src tests

# Lint code
lint:
	uv run ruff check --fix src tests

# Run tests
tests:
	if [ "$(MAKECMDGOALS)" != "tests" ]; then \
		uv run pytest $(MAKECMDGOALS) -s -vv; \
	else \
		uv run pytest tests/ -s -vv; \
	fi
