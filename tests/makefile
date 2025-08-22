# Look for pytest under different names
PYTEST := $(shell which pytest)
PYTEST := $(or $(PYTEST), pytest-3)

test:
	PYTHONDONTWRITEBYTECODE=1 $(PYTEST) -v

clean:
	- rm -rf *.bin temp.yaml

