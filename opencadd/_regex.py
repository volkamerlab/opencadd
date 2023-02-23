"""
General regular expressions used in the package.
"""

import re

ANY_NUM = re.compile(r"^-?\d+(?:\.\d+)?$")
"""Match any number; positive or negative; decimal or integer.

This regular expression matches numbers that may or may not have a minus sign at the beginning, 
and may or may not have a decimal point with one or more digits after it.
Breakdown:
* '^': Matches the start of the string.
* '-?': Matches an optional minus sign.
* '\d+': Matches one or more digits (0-9).
* '(?:\.\d+)?': Matches an optional decimal point followed by one or more digits. 
  The '?:' syntax creates a non-capturing group.
* '$': Matches the end of the string.
"""