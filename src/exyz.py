import numpy


class exyz:

  @staticmethod
  def load(file_name):
      
    fh = open(file_name, 'r')
    in_data = False
    
    temp = {}
    
    for line in fh:
      line = line.strip()
      if(line[0:7] == "0000001"):
        in_data = True
      if(in_data and line != ""):
        ion = int(line[0:7])
        if(ion not in temp.keys()):
          temp[ion] = []
        energy = float(line[8:18])
        x = float(line[19:30])
        y = float(line[31:42])
        z = float(line[43:54])
        temp[ion].append([energy,x,y,z])
    
    
    data = {}
    for k in temp.keys():
      l = len(temp[k]) - 1
      temp_dr = []
      for i in range(len(temp[k]) - 1):
        e = 0.5 * (temp[k][i][0] + temp[k][i+1][0])
        dr = numpy.sqrt((temp[k][i][1] - temp[k][i+1][1])**2 + (temp[k][i][2] - temp[k][i+1][2])**2 + (temp[k][i][3] - temp[k][i+1][3])**2)
        if(dr > 0.0):
          temp_dr.append([e, dr, temp[k][i][0], temp[k][i+1][0], temp[k][i][1], temp[k][i+1][1]])
          
      data[k] = numpy.zeros((len(temp_dr), 6),)
      for i in range(len(temp_dr)):
        data[k][i,0] = temp_dr[i][0] * 1.0E-3   # Energy          Convert to MeV
        data[k][i,1] = temp_dr[i][1] * 1.0E-10  # m dr
        data[k][i,2] = temp_dr[i][2] * 1.0E-3   # Start Energy    Convert to MeV
        data[k][i,3] = temp_dr[i][3] * 1.0E-3   # End Energy      Convert to MeV
        data[k][i,4] = temp_dr[i][4] * 1.0E-10  # m start depth
        data[k][i,5] = temp_dr[i][5] * 1.0E-10  # m end depth

    fh.close()

    return data  