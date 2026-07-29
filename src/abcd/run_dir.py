"""Print the run directory for the active config, and nothing else.

Exists so non-Python steps (``R/fit_lmm.R``, ``hpc/*.sh``) can resolve the same
run as the Python steps without reimplementing the config hash.  Duplicating
that hash in R or bash would drift the moment a config field is added -- and the
failure would be silent, pointing a fit at a stale run directory.  One
implementation, shelled out to::

    RUN=$(python -m abcd.run_dir)

Writes only the path to stdout (errors go to stderr, exit 1), so it is safe to
capture in a shell variable.

Usage:
    python -m abcd.run_dir                    # uses $ABCD_CONFIG
    python -m abcd.run_dir ct_70_genetic      # explicit
    python -m abcd.run_dir --allow-missing    # path even if not yet assembled
"""

from __future__ import annotations

import argparse
import sys

from .config import ConfigError, active_run_dir


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("config", nargs="?", default=None,
                    help="config name or path; omit to use $ABCD_CONFIG")
    ap.add_argument("--allow-missing", action="store_true",
                    help="print the path even if the run has not been assembled")
    a = ap.parse_args(argv)
    try:
        print(active_run_dir(a.config, must_exist=not a.allow_missing))
    except ConfigError as exc:
        print(f"abcd.run_dir: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
