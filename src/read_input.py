import os
from read_config import read_config
from units import units
from isotopes import isotopes

class read_input:

  
  def run(input_file):
  
    # Read input file
    inp = read_config.read_file(input_file)
    

    config = {}
   
    config['isotopes'] = str(inp['data']['isotopes']).strip()
    config['xs'] = str(inp['data']['xs']).strip()
    config['sims_in'] = []
    config['sims'] = []


    if(not os.path.isfile(config['isotopes'])):
      print(
"""
################################
Program terminated with an error
################################

Isotope file does not exist.
""")
      exit()
    isotopes.set(config['isotopes'])

   
      
    # Sim input
    loop = True
    n = 1
    m = 0
    while(loop):
      k = 'sim' + str(n)
      if(k in inp.keys()): 
        read_input.sim(config, k, m, inp)
        m = m + 1
      elif(n > 100 and k not in inp.keys()): 
        loop = False    
      n = n + 1   

    for i in range(len(config['sims_in'])):
      read_input.process_sim(config, i)

    return config

               
  def sim(config, k, n, inp): 
    config['sims_in'].append({
        'run': False, 
        'exyz': 'exyz.exyz',
        'neutron': None,
        'target_composition': 'FE,100',
        'target_depth': 1.0,
        'target_depth_unit': 'mm',
        'target_density': 10000,
        'target_density_unit': 'kgm3',
        'beam_projectile': 'proton',
        'beam_energy': 10,
        'beam_energy_unit': 'MeV',
        'beam_area': 100,
        'beam_area_unit': 'mm2',
        'beam_duration': 100,
        'beam_duration_unit': 's',
        'beam_current': 10,
        'beam_current_unit': 'uA',
        'beam_flux': None,
        'end_time': 100000,
        'end_time_unit': 's',
        'inbeampoints': 101,
        'outbeampoints': 101,
        'plot_w': 8,
        'plot_h': 6,
        'plot_dpi': 90,
        })

    # Load data             
    try:
      config['sims_in'][n]['run'] = False
      if(inp[k]['run'][0] == 'T'):
        config['sims_in'][n]['run'] = True
    except:  
      pass       
    try:
      config['sims_in'][n]['exyz'] = inp[k]['exyz']
    except:  
      pass
    try:
      config['sims_in'][n]['neutron'] = inp[k]['neutron']
    except:  
      pass
    try:
      config['sims_in'][n]['target_composition'] = inp[k]['target_composition']
    except:  
      pass
    try:
      config['sims_in'][n]['target_depth'] = float(inp[k]['target_depth'][0])
    except:  
      pass
    try:
      config['sims_in'][n]['target_depth_unit'] = inp[k]['target_depth'][1]
    except:  
      pass
    try:
      config['sims_in'][n]['target_density'] = float(inp[k]['target_density'][0])
    except:  
      pass
    try:
      config['sims_in'][n]['target_density_unit'] = inp[k]['target_density'][1]
    except:  
      pass
    try:
      config['sims_in'][n]['beam_projectile'] = inp[k]['beam_projectile']
    except:  
      pass
    try:
      config['sims_in'][n]['beam_energy'] = float(inp[k]['beam_energy'][0])
    except:  
      pass
    try:
      config['sims_in'][n]['beam_energy_unit'] = inp[k]['beam_energy'][1]
    except:  
      pass
    try:
      config['sims_in'][n]['beam_area'] = inp[k]['beam_area'][0]
    except:  
      pass
    try:
      config['sims_in'][n]['beam_area_unit'] = inp[k]['beam_area'][1]
    except:  
      pass
    try:
      config['sims_in'][n]['beam_duration'] = inp[k]['beam_duration'][0]
    except:  
      pass
    try:
      config['sims_in'][n]['beam_duration_unit'] = inp[ k]['beam_duration'][1]
    except:  
      pass
    try:
      config['sims_in'][n]['beam_current'] = inp[k]['beam_current'][0]
    except:  
      pass
    try:
      config['sims_in'][n]['beam_current_unit'] = inp[k]['beam_current'][1]
    except:  
      pass
    try:
      config['sims_in'][n]['beam_flux'] = inp[k]['beam_flux']
    except:  
      pass
    try:
      config['sims_in'][n]['end_time'] = inp[k]['end_time'][0]
    except:  
      pass
    try:
      config['sims_in'][n]['end_time_unit'] = inp[k]['end_time'][1]
    except:  
      pass
    try:
      config['sims_in'][n]['inbeampoints'] = inp[k]['inbeampoints']
    except:  
      pass
    try:
      config['sims_in'][n]['outbeampoints'] = inp[k]['outbeampoints']
    except:  
      pass
    try:
      config['sims_in'][n]['plot_w'] = float(inp[k]['plot_w'])
    except:  
      pass
    try:
      config['sims_in'][n]['plot_h'] = float(inp[k]['plot_h'])
    except:  
      pass
    try:
      config['sims_in'][n]['plot_dpi'] = int(inp[k]['plot_dpi'])
    except:  
      pass


               
  def process_sim(config, k): 

    config['sims'].append({
        'run': False,
        'exyz': None,
        'neutron': None,
        'target_composition': None,
        'target_composition_hr': '',
        'target_depth': 0.0,             # m
        'target_density': 0.0,           # kgm3
        'beam_projectile': None,         #
        'beam_energy': 0.0,              # MeV
        'beam_area': 0.0,                # m^2
        'beam_duration': 0.0,            # s
        'beam_flux': 0.0,                # projectiles/second
        'end_time': 0.0,                 # s
        'inbeampoints': 101,
        'outbeampoints': 101,
        'outbeampoints': 101,
        'plot_w': 8,
        'plot_h': 6,
        'plot_dpi': 90,
        })

    config['sims'][k]['run'] = config['sims_in'][k]['run']

    if(os.path.isfile(config['sims_in'][k]['exyz'])):
      config['sims'][k]['exyz'] = config['sims_in'][k]['exyz']

    if(config['sims_in'][k]['neutron'] != None):
      count = int(config['sims_in'][k]['neutron'][0])
      dist = config['sims_in'][k]['neutron'][1]

      config['sims'][k]['neutron'] = {}
      config['sims'][k]['neutron']['count'] = count
      config['sims'][k]['neutron']['dist'] = dist

      if(dist == 'mono'):
        config['sims'][k]['neutron']['energy'] = float(config['sims_in'][k]['neutron'][2])
      elif(dist == 'flat'):
        config['sims'][k]['neutron']['min'] = float(config['sims_in'][k]['neutron'][2])
        config['sims'][k]['neutron']['max'] = float(config['sims_in'][k]['neutron'][3])
      elif(dist == 'maxwell'):
        config['sims'][k]['neutron']['mean'] = float(config['sims_in'][k]['neutron'][2])
      
    
    # CONVERT UNITS
    # M, M2, S, A, MEV, KGM-3

    config['sims'][k]['target_depth'] = units.convert(
           config['sims_in'][k]['target_depth_unit'], 'M', config['sims_in'][k]['target_depth'])
    config['sims'][k]['target_density'] = units.convert(
           config['sims_in'][k]['target_density_unit'], 'KGM-3', config['sims_in'][k]['target_density'])
    config['sims'][k]['beam_energy'] = units.convert(
           config['sims_in'][k]['beam_energy_unit'], 'MEV', config['sims_in'][k]['beam_energy'])
    config['sims'][k]['beam_area'] = units.convert(
           config['sims_in'][k]['beam_area_unit'], 'M2', config['sims_in'][k]['beam_area'])
    config['sims'][k]['beam_duration'] = units.convert(
           config['sims_in'][k]['beam_duration_unit'], 'S', config['sims_in'][k]['beam_duration'])
    config['sims'][k]['end_time'] = units.convert(
           config['sims_in'][k]['end_time_unit'], 'S', config['sims_in'][k]['end_time'])

    config['sims'][k]['inbeampoints'] = int(config['sims_in'][k]['inbeampoints'])
    config['sims'][k]['outbeampoints'] = int(config['sims_in'][k]['outbeampoints'])


    config['sims'][k]['plot_w'] = int(config['sims_in'][k]['plot_w'])
    config['sims'][k]['plot_h'] = int(config['sims_in'][k]['plot_h'])
    config['sims'][k]['plot_dpi'] = int(config['sims_in'][k]['plot_dpi'])



    # BEAM PROJECTILE
    beam_projectile = str(config['sims_in'][k]['beam_projectile']).upper()
    if(beam_projectile[0] == 'P'):
      config['sims'][k]['beam_projectile'] = 'p'
      projectile_charge = 1
    elif(beam_projectile[0] == 'D'):
      config['sims'][k]['beam_projectile'] = 'd'
      projectile_charge = 1
    elif(beam_projectile[0] == 'T'):
      config['sims'][k]['beam_projectile'] = 't'
      projectile_charge = 1
    elif(beam_projectile[0:3] == 'HE3'):
      config['sims'][k]['beam_projectile'] = 'he3'
      projectile_charge = 2
    elif(beam_projectile[0] == 'A'):
      config['sims'][k]['beam_projectile'] = 'a'
      projectile_charge = 2
    elif(beam_projectile[0] == 'N'):
      config['sims'][k]['beam_projectile'] = 'n'
      projectile_charge = 0
    else:      
      print(
"""
################################
Program terminated with an error
################################

Invalid projectile code selected.
Accepted codes:
p = proton
d = Deuteron
t = Triton
he3 = Helium-3
a = alpha particle
""")
      exit()
      
    # BEAM FLUX
    beam_current = units.convert(config['sims_in'][k]['beam_current_unit'], 'A', config['sims_in'][k]['beam_current'])
    if(projectile_charge > 0):
      config['sims'][k]['beam_flux'] = beam_current / (projectile_charge * 1.60217662E-19)
    if(config['sims_in'][k]['beam_flux'] != None):
      config['sims'][k]['beam_flux'] = float(config['sims_in'][k]['beam_flux'])


    # TARGET
    target_composition = config['sims_in'][k]['target_composition']
    if(len(target_composition) % 2 != 0):   
      print(
"""
################################
Program terminated with an error
################################

Target list length must be even - isotope + relative weight
e.g. Fe,50,Ni,50
""")
      exit()

    # Valid numbers
    num = "0123456789"
    
    # String 
    elements = []
    masses = []
    for i in range(len(target_composition)//2):
      elements.append(target_composition[2*i])
      masses.append(float(target_composition[2*i+1]))
      
    # Mass to a percentage
    for m in range(len(masses)):
      masses[m] = float(masses[m])
    tm = sum(masses)    
    for m in range(len(masses)):
      masses[m] = masses[m] * (100 / tm)  
    
    # Capitalise elements
    for e in range(len(elements)):
      elements[e] = elements[e].capitalize().strip()

    target = {}

    for n in range(len(elements)):
      element = ''
      isotope = ''
      for c in elements[n]:
        if(c in num):
          isotope = isotope + c
        else:
          element = element + c
      if(isotope == ""):
        stable_isotopes = isotopes.get_stable_abundances(element)
        for icode in stable_isotopes:
          isotope = isotopes.get(icode)
          if(icode not in target.keys()):   
            target[icode] = {'pbm': 0.0, 'pbn': 0.0, 'nd': 0.0, 'amu': isotope['mass_amu'],}
          target[icode]['pbm'] = target[icode]['pbm'] + isotope['natural_abundance'] * masses[n]

      else:
        icode = isotopes.get_code(element, isotope)
        isotope = isotopes.get(icode)
        if(icode not in target.keys()):       
          target[icode] = {'pbm': 0.0, 'pbn': 0.0, 'nd': 0.0, 'amu': isotope['mass_amu'],}
        target[icode]['pbm'] = target[icode]['pbm'] + masses[n]

    # Renormalise to 100 to account for rounding errors in natural abundance data
    t = 0.0
    for icode in target.keys():
      t = t + target[icode]['pbm']
    for icode in target.keys():
      target[icode]['pbm'] = target[icode]['pbm'] * (100.0 / t)
    

    t = 0.0
    for icode in target.keys():
      t = t + target[icode]['pbm']
  
    
    # Calculate number percentage of each isotope
    s = 0.0
    for icode in target.keys():
      s = s + target[icode]['pbm'] / target[icode]['amu']
    for icode in target.keys():
      target[icode]['pbn'] = (100 / s) * (target[icode]['pbm'] / target[icode]['amu'])
    
    # Renormalise to 100 to account for rounding errors
    t = 0.0
    for icode in target.keys():
      t = t + target[icode]['pbn']
    for icode in target.keys():
      target[icode]['pbn'] = target[icode]['pbn'] * (100.0 / t)


    # Composition density
    avogadro = 6.0221409E23
    for icode in target.keys():
      isotope_density = config['sims'][k]['target_density'] * (target[icode]['pbm'] / 100) * 1000 # g per M3
      target[icode]['nd'] = (isotope_density / target[icode]['amu']) * avogadro
    config['sims'][k]['target_composition'] = target

    i_list = []
    pbm_list = []
    for icode in target.keys():
      i_list.append(icode)
      pbm_list.append(target[icode]['pbm'])

    config['sims'][k]['target_composition_hr'] = isotopes.isotope_table(i_list, pbm_list)



