import numpy
import matplotlib.pyplot as plt

class energy:

  dist = 'flat'

  # Flat
  min = 0.0
  max = 1.0

  # 
  sigma = 0.0
  mean = 0.0
 
  @staticmethod  
  def set(reset=None, dist=None, min=None, max=None, sigma=None, mean=None):
    if(reset == True):
      energy.dist = 'flat'
      energy.min = 0.0
      energy.max = 1.0
      return True
    if(dist != None):
      energy.dist = dist
    if(min != None):
      energy.min = min
    if(max != None):
      energy.max = max
    if(sigma != None):
      energy.sigma = sigma
    if(mean != None):
      energy.mean = mean

  @staticmethod 
  def rand():
    if(energy.dist == 'gaussian'):
      return numpy.random.normal(energy.mean, energy.sigma, 1)[0]
    elif(energy.dist == 'maxwell'):
      return energy.rand_maxwell()
    # Flat (default)
    return numpy.random.uniform(energy.min, energy.max)


  @staticmethod 
  def rand_maxwell():
    a = energy.mean / (2.0 * 0.797884561)
    while(True):
      r1 = numpy.random.uniform(0.0, 5.0 * a)
      r2 = numpy.random.uniform(0.0, 1.0)
      if(r2 <= energy.maxwell(r1, a)):
        return r1
    #  return numpy.random.uniform(energy.min, energy.max)


  @staticmethod 
  def maxwell(x, a):
    return 0.797884561 * ((x * x * numpy.exp((-x**2)/(2*a**2))) / a**3)
    

  @staticmethod 
  def test():
    print("Test")
    
    energy.set(min=0.9)
    x = energy.rand()
    print(x)

    energy.set(reset=True)
    x = energy.rand()


    energy.set(dist='gaussian', mean=0.2, sigma=1.0)
    x = []
    for n in range(100000):
      x.append(energy.rand())
    plt.hist(x, density=True, bins=30)
    plt.show()



    energy.set(dist='maxwell', mean=12.0)
    x = []
    for n in range(100000):
      x.append(energy.rand())
    plt.hist(x, density=True, bins=30)
    plt.show()

    """
    x = numpy.linspace(0.0, 10.0, 101)
    a = energy.mean  / (2.0 * 0.797884561)
    y = energy.maxwell(x, a)
    plt.plot(x, y)
    plt.show()
    """

    """
    x = []
    for n in range(1000000):
      x.append(energy.rand())
    plt.hist(x, density=True, bins=100)
    plt.show()
    """





def main():
  energy.test()

if __name__ == "__main__":
    main()    
