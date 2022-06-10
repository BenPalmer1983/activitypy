class units:
  
  @staticmethod
  def convert(conv_from, conv_to, value_in):
    try:
      value_in = float(value_in)
    except:
      return None
    conv_from = conv_from.upper()
    conv_to = conv_to.upper()

    # LENGTH METERS
    length = {
    'M': 1.0,
    'CM': 100,
    'MM': 1E3,
    'UM': 1E6,
    'NM': 1E9,
    'ANG': 1E10,
    'BOHR': 1.89E10,
    }
    
    # AREA METERS SQUARED
    area = {
    'M2': 1.0,
    'CM2': 1E4,
    'MM2': 1E6,
    'UM2': 1E12,
    'NM2': 1E18,
    'ANG2': 1E20,
    }
    
    # VOLUME METERS CUBED
    volume = {
    'M2': 1.0,
    'CM2': 1E6,
    'MM2': 1E9,
    'UM2': 1E18,
    'NM2': 1E27,
    'ANG2': 1E30,
    }

    # ENERGY J
    energy = {
    'J': 1.0,
    'EV': 6.2415E18,
    'RY': 4.5874E17,
    'KEV': 6.2415E15,
    'MEV': 6.2415E12,
    }

    # FORCE N
    force = {
    'N': 1.0,
    'RY/BOHR': 2.4276e7,
    'EV/ANG':6.2414E8,    
    }
    
    # VELOCITY
    velocity = {
    'M/S': 1.0,
    'MPH': 2.25,    
    }
    
    # PRESSURE
    pressure = {
    'PA': 1.0,
    'GPA': 1.0E-9,    
    'BAR': 1.0E-5,    
    'ATMOSPHERE': 9.8692E-6,    
    'PSI': 1.45038E-4, 
    'KBAR': 1.0E-8,   
    'RY/BOHR3': 6.857E-14,   
    'EV/ANG3': 6.241E-12
    }
    
    # CHARGE DENSITY (UNIT CHARGE PER VOLUME - ANG^3)
    charge_density = {
    'ANG-3': 1.0,
    'BOHR-3': 0.14812,    
    }

    # TIME T
    time = {
    'S': 1.0,
    'MINUTE': 0.1666666666E-2,
    'HOUR': 2.7777777777E-4,
    }
    
    # CURRENT
    current = {
    'A': 1.0,
    'MA': 1.0E3,
    'UA': 1.0E6,
    }
    
    # DENSITY
    density = {
    'KGM3': 1.0,
    'GCM3': 1.0E-3,
    'KGM-3': 1.0,
    'GCM-3': 1.0E-3,
    }
    
    
    
    # TEMPERATURE
    
    
    unit_list = [length, area, volume, energy, force, velocity, pressure, charge_density, time, current, density]
    
    for l in unit_list:
      if(conv_from in l.keys() and conv_to in l.keys()):
        return round((l[conv_to] / l[conv_from]) * float(value_in),9)
  

