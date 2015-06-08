#!/usr/bin/env python
import sys
try:
  template, testing = sys.argv[1:]
except:
  sys.exit("Expecting 2 args: a template file and an output")
with open(template, 'rU') as temp_in:
  pats = [i.strip() for i in temp_in if i.strip()]
if not pats:
  sys.exit(0)
with open(testing, 'rU') as test_in:
  for line in test_in:
    ls = line.strip()
    if ls == pats[0]:
      pats.pop(0)
      if not pats:
        sys.exit(0)
sys.exit('Did not find "{p}" in "{t}"\n'.format(p=pats[0], t=testing))
