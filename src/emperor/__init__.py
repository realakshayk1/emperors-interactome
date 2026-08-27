"""The Emperor's Interactome — conformal FDR audit of an AI protein interactome."""
__version__ = "0.1.0"

# Several modules print Greek letters (sigma, delta, Delta) and em dashes in their progress
# output. On Windows the default console encoding is cp1252, which raises UnicodeEncodeError
# and aborts the step -- `make reproduce` died at .validate and `make audit-self` at
# .sensitivity. Force UTF-8 on our own streams so the pipeline is portable; harmless where
# UTF-8 is already the default.
import sys as _sys

for _stream in (_sys.stdout, _sys.stderr):
    try:
        _stream.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):  # non-reconfigurable stream (piped/captured)
        pass
