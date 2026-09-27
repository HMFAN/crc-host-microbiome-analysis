import math
def rank(columns):
 basis=[]
 for col in columns:
  norm=math.sqrt(sum(v*v for v in col));assert norm>0
  x=[v/norm for v in col]
  for _ in range(2):
   for b in basis:
    dot=sum(v*w for v,w in zip(x,b));x=[v-dot*w for v,w in zip(x,b)]
  n=math.sqrt(sum(v*v for v in x))
  if n>1e-9:basis.append([v/n for v in x])
 return len(basis)
