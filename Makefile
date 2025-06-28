.PHONY: checks fmt lint type tests

# Run pre-commit checks
checks:
	uvx pre-commit run --all-files

fmt:
	uv run ruff format src tests

lint:
	uv run ruff check --fix src tests

type:
	uv run mypy src --install-types --non-interactive --show-traceback

tests:
	if [ "$(MAKECMDGOALS)" != "tests" ]; then \
		uv run pytest --cov=src --cov-report=term-missing $(MAKECMDGOALS) -s -vv; \
	else \
		uv run pytest --cov=src --cov-report=term-missing tests/ -s -vv; \
	fi