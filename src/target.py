class target:

  avogadro = 6.0221409E23

  def __init__(self):
    self.composition = {}
    self.avgmass = None
    self.depth = None
    self.density = None
    self.loaded = False
    

  # example:   Fe,0.8,Cr,0.2 
  def set_composition(self, elements_in):  
    num = "0123456789"
    
    # String 
    elements = []
    masses = []
    for i in range(len(elements_in)//2):
      elements.append(elements_in[2*i])
      masses.append(float(elements_in[2*i+1]))
      
    # Mass to a percentage
    for m in range(len(masses)):
      masses[m] = float(masses[m])
    tm = sum(masses)    
    for m in range(len(masses)):
      masses[m] = masses[m] * (100 / tm)
      
     
    # Capitalise elements
    for e in range(len(elements)):
      elements[e] = elements[e].capitalize().strip()
      
    # Create target
    self.composition = {}
      
    self.avgmass = 0.0
    for n in range(len(elements)):
      element = ''
      isotope = ''
      for c in elements[n]:
        if(c in num):
          isotope = isotope + c
        else:
          element = element + c

      # Get isotope
      s = isotopes.get_isotopes(element, isotope)
      
      for sn in s:        
        k = int(1000 * int(sn['protons']) + int(sn['nucleons']))
       
        self.composition[k] = {
                   'element': isotopes.get_element(sn['protons']),
                   'protons': sn['protons'],
                   'neutrons': sn['neutrons'],
                   'nucleons': sn['nucleons'],
                   'mass': sn['mass'],
                   'percentage_by_mass': (masses[n] / 100) * sn['percentage'],
                   'percentage_by_number': 0.0,
                   'isotope_density': 0.0,
                   'isotope_number_density': 0.0,
                   }
        self.avgmass = self.avgmass + (self.composition[k]['percentage_by_mass'] / 100.0) * sn['mass']       
    
    # Calculate number percentage of each isotope
    s = 0.0
    for k in self.composition.keys():
      m = self.composition[k]['mass']
      s = s + self.composition[k]['percentage_by_mass'] / m
    for k in self.composition.keys():
      m = self.composition[k]['mass']
      self.composition[k]['percentage_by_number'] = ((self.composition[k]['percentage_by_mass'] / m) / s) * 100  
           
        
        
  def set_depth(self, depth, depth_unit):
    self.depth = units.convert(depth_unit, 'M', depth)
    
  def set_density(self, density, density_unit):
    self.density = units.convert(density_unit, 'KGM3', density)
    
    
  def calc_nd(self):
    for k in self.composition.keys():
      self.composition[k]['isotope_density'] = self.density * (self.composition[k]['percentage_by_mass'] / 100)
      self.composition[k]['isotope_number_density'] = ((self.composition[k]['isotope_density'] * 1E3) / self.composition[k]['mass']) * target.avogadro

    
  def get_isotope_list(self):
    i_list = []
    for k in self.composition.keys():
      i_list.append(k)
    return i_list  
    
  def get_depth(self):
    return self.depth

  def get_depth_ang(self):
    return 1.0E10 * self.depth
    
  def get_target_number_density(self, target):
    if(target in self.composition.keys()):
      return self.composition[target]['percentage_by_number']
    return 0.0
    
  def get_isotope_number_density(self, target):
    if(target in self.composition.keys()):
      return self.composition[target]['isotope_number_density']
    return 0.0
          
  def display(self):
    print("Target")
    print("========================================================")
    print("Depth:   ", self.depth, "m")
    print("Density: ", self.density, "kgm3")  
    print("Mass:    ", self.avgmass, "amu")  
    print("Isotope  Mass         PBM          PBN")
    print("========================================================")
    for k in self.composition.keys():
      isotope = '{:6}'.format(self.composition[k]['element'] + str(self.composition[k]['nucleons']))
      mass = str("{:.4E}".format(self.composition[k]['mass']))
      pbm = str("{:.4E}".format(self.composition[k]['percentage_by_mass']))
      pbn = str("{:.4E}".format(self.composition[k]['percentage_by_number']))
      id = str("{:.4E}".format(self.composition[k]['isotope_density']))
      ind = str("{:.4E}".format(self.composition[k]['isotope_number_density']))
      print(isotope, mass, pbm, pbn, id, ind, sep="   ")
  

    print("========================================================")
    print()
