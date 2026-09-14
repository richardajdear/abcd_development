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

# The interpreter with abcd's dependencies.  A bare `python` on this machine is a
# pyenv shim without pandas, so prefer the project env when it exists; override
# with PY=... on the command line.
PY ?= $(shell test -x $(HOME)/mambaforge/envs/abcd/bin/python && echo $(HOME)/mambaforge/envs/abcd/bin/python || echo python)
# R/fit_lmm.R shells out to $PY to resolve the run directory.  The export must
# come AFTER the assignment: in GNU make 3.81 a bare `export PY` defines PY as
# empty and the `?=` above is then skipped.
export PY
RSCRIPT ?= Rscript
CORES ?= 4
RUN = $(shell $(SRCPATH) $(PY) -m abcd.run_dir --allow-missing 2>/dev/null)

.PHONY: all config assemble fit phenotype gcta test clean-run help audit-artifacts

help:
	@echo "ABCD_CONFIG=$(ABCD_CONFIG)"
	@echo "run dir     =$(RUN)"
	@echo
	@echo "targets: assemble -> fit -> phenotype -> gcta   (all = every step)"
	@echo "         test, audit-artifacts, clean-run"
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

# Reconcile repo deliverables against the Claude Science artifact tray.
#
# One-directional by design: every repo deliverable must be current in the tray,
# but the tray may hold plenty that is not in the repo (exploratory figures,
# checkpoints, one-off analyses) and that is fine.
#
# This target writes the manifest only.  The comparison needs host.artifacts(),
# a kernel global that no subprocess can import, so an agent finishes the audit
# by reading .artifact_audit.json in its python kernel:
#
#   import json, os
#   rows = json.load(open(".artifact_audit.json"))
#   art  = {a["filename"]: a for a in host.artifacts(limit=500)["artifacts"]}
#   for r in rows:
#       a = art.get(r["artifact_name"])
#       if a is None:                            print("UNSAVED", r["path"])
#       elif a["size_bytes"] != r["size_bytes"]: print("STALE  ", r["path"], a["id"])
#
# then re-saves the offenders with version_of={name: artifact_id} for the stale
# ones.  Sizes are the cheap check; .artifact_audit.json also carries sha256 for
# when a file changes without changing length.
audit-artifacts:
	$(PY) tools/audit_artifacts.py

# Delete only the active run, never all of out/.
clean-run: config
	@test -n "$(RUN)" && rm -rf "$(RUN)" && echo "removed $(RUN)"
