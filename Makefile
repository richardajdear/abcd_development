# ABCD pipeline. One variable selects the run; every step picks it up.
#
#   export ABCD_CONFIG=ct_70_genetic
#   make all
#
# or per-invocation, without exporting:
#
#   make all ABCD_CONFIG=ct_70_genetic
#
# No ABCD_ROOT needed: 7.0 is vendored in the repo (gitignored) and 5.1 is
# found under ~/Git/ABCD.  See src/abcd/paths.py.

export ABCD_CONFIG

# Work whether or not `pip install -e .` has been run. Without this every
# target fails with ModuleNotFoundError, and `help` silently prints a blank
# run dir because $(shell ...) discards stderr.
#
# The export covers recipes; the $(shell ...) below needs PYTHONPATH passed
# explicitly as well, because GNU Make 3.81 (the macOS default) does not
# propagate exported variables into the shell function.
export PYTHONPATH := src:$(PYTHONPATH)
SRCPATH = PYTHONPATH=src:$(PYTHONPATH)

PY ?= python
RSCRIPT ?= Rscript
CORES ?= 4
RUN = $(shell $(SRCPATH) $(PY) -m abcd.run_dir --allow-missing 2>/dev/null)

.PHONY: all config assemble fit phenotype gcta test clean-run help

help:
	@echo "ABCD_CONFIG=$(ABCD_CONFIG)"
	@echo "run dir     =$(RUN)"
	@echo
	@echo "targets: assemble -> fit -> phenotype -> gcta   (all = every step)"
	@echo "configs: $(notdir $(basename $(wildcard configs/*.yaml)))"

# Fail early with a readable message rather than midway through a step.
config:
	@$(PY) -m abcd.run_dir --allow-missing >/dev/null

all: gcta

assemble: config
	$(PY) -m abcd.assemble

fit: assemble
	$(RSCRIPT) R/fit_lmm.R --cores $(CORES)

phenotype: fit
	$(PY) -m abcd.phenotype

gcta: phenotype
	$(PY) -m abcd.gcta_export

test:
	$(PY) -m pytest tests -q

# Delete only the active run, never all of out/.
clean-run: config
	@test -n "$(RUN)" && rm -rf "$(RUN)" && echo "removed $(RUN)"
