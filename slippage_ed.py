import sys

sites = {}

with open(sys.argv[1], 'U') as f:
  for line in f:
    l = line.split("\t")
    try:
      sites[int(l[1])] = int(l[0])
    except ValueError:
      1+1

id = 0
next_id = 0
results = {}

for s in sorted(sites):
  if s <= id:
    continue
  if sites.get(s, -1) == 0:
    continue

  local_max_id = s
  local_max = sites[s]
  sum = sites[s]
  id = s + 1
  while sites.get(id) > -1:
    sum = sum + sites.get(id)
    if sites.get(id) > local_max:
      local_max = sites.get(id)
      local_max_id = id

    id = id + 1

  results[local_max_id] = sum

for r in sorted(results):
  print(str(r) + "\t" + str(results[r]))
