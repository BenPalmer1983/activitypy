import os
import sys
import numpy
import matplotlib.pyplot as plt
import time
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

from g import g
from read_input import read_input
from sim import sim

#from isotopes import isotopes




class activity:

  start = time.time()
  input_file = None

  def main():
    g.set_time_start()
    if(len(sys.argv) >= 2):
      if(sys.argv[1].lower().strip() == "test"):
        activity.test()
      elif(sys.argv[1].lower().strip() == "help"):
        activity.help()
      else:
        activity.input_file = sys.argv[1]

    if(activity.input_file == None):
      print("No input file")
      exit()
    else:
      if(os.path.isfile(activity.input_file)):
        activity.run()
      else: 
        print("File does not exist")
        exit()


  def run():
    print("###############################################")
    print("                  Activity                     ")
    print("###############################################")
    print()

    fhs = open('summary.csv', 'w')
    fhs.write("Sim,Flux,Target Depth,Projectile,Beam Duration, Sim End Time,EoB Activity,EoB Gamma Activity,EoB Gamma Power/mW,EoS Activity,EoS Gamma Activity,EoS Gamma Power/mW")
    fhs.close()
    
    g.print("Read input file")
    config = read_input.run(activity.input_file)


    sim.set(config['isotopes'], config['xs'])
    n = 0
    for this_sim in config['sims']:
      n = n + 1
      sim.run(this_sim, n)


    print("Run time: ", activity.time())

  def help():
    print("Help")

  def test():
    print("Test")

  def time():
    return time.time() - activity.start






if __name__ == "__main__":
  activity.main()    

