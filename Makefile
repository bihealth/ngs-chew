.PHONY: default black

default: black

black:
	black -l 100 .
