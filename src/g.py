import time


class g:

  time_start = 0
  time_end = 0



  def set_time_start():
    g.start = time.time()

  def print(line="", fh=None):
    lines = line.strip().split("\n")
    for line in lines:
      t = time.time() - g.start
      t = g.pad("{:.3f}".format(t), 12)
      print(t, line)    
      if(fh != None):
        fh.write(t + line + '\n')

  def pad(inp="", l=5):
    inp = str(inp)
    while(len(inp)<l):
      inp = inp + " "
    return inp


