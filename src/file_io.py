

def read_swap_trajectory ( fname ):

  inds,nswap,temps,enes = [],[],[],[]
  with open(fname, 'r') as f:
    for l in f.readlines():
      i,n,t,e = (float(v) for v in l.split())
      print(i)
      inds.append(i)
      nswap.append(n)
      temps.append(t)
      enes.append(e)

  return inds, nswap, temps, enes


def write_swap_accept ( fname,
                        write_mode,
                        step,
                        prev_accept,
                        temperature,
                        energy,
                        occupation_factors=None ):

  with open(fname, write_mode) as f:
    if occupation_factors is None:
      f.write(' '.join(map(str,[step, prev_accept, temperature, energy]))+'\n')
