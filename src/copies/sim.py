import numpy
import matplotlib.pyplot as plt
import os
import shutil


from g import g
from isotopes import isotopes
from talys import talys
from exyz import exyz
from decay import decay
from std import std
from energy import energy


class sim:

  isotope_file = None
  talys_dir = None

  plot_w = 8
  plot_h = 6
  plot_dpi = 90

  piechart_threshold = 0.0001

  colours = []
  isotope_colour = {}
  
  @staticmethod
  def set(isotope_file, talys_dir):
    sim.isotope_file = isotope_file
    sim.talys_dir = talys_dir

    isotopes.set(sim.isotope_file)
    decay.set(sim.isotope_file)

    g.print("Loading XS data")
    talys.set(sim.talys_dir, sim.isotope_file) 
    talys.load()
    g.print("Load complete")


  @staticmethod
  def run(this_sim, counter):
    sim.make_colours()

    dir = "results/sim_" + str(counter)

    data_dir = dir + "/data"
    plot_dir = dir + "/plots"
    video_dir = dir + "/videos"
    logs_dir = dir + "/logs"


    if(os.path.isdir(data_dir)):
      shutil.rmtree(data_dir)
    if(os.path.isdir(plot_dir)):
      shutil.rmtree(plot_dir)
    if(os.path.isdir(video_dir)):
      shutil.rmtree(video_dir)
    if(os.path.isdir(logs_dir)):
      shutil.rmtree(logs_dir)

    std.make_dir(dir)
    std.make_dir(data_dir)
    std.make_dir(plot_dir)
    std.make_dir(video_dir)
    std.make_dir(logs_dir)

    # Open files
    fhr = open(dir + '/results.txt', 'w')

    # Log files
    sim.gamma_log = open(logs_dir + '/piecharts.txt', 'a+')
    sim.activity_log = open(logs_dir + '/activity.txt', 'a+')

    g.print('######################################################\n', fhr)
    g.print('Sim ' + str(counter) + '\n', fhr)
    g.print('######################################################\n', fhr)
    g.print('', fhr)

    beam_flux = this_sim['beam_flux']
    target_depth = this_sim['target_depth']
    projectile = this_sim['beam_projectile']
    beam_duration = this_sim['beam_duration']
    end_time = this_sim['end_time']


    g.print('Beam flux:          ' + str(beam_flux), fhr)
    g.print('Target depth:       ' + str(target_depth), fhr)
    g.print('Projectile:         ' + str(projectile), fhr)
    g.print('Beam duration/s:    ' + str(beam_duration), fhr)
    g.print('Sim end time/s:     ' + str(end_time), fhr)
    g.print('', fhr)


    if(this_sim['run'] == False):
      g.print('Run == False', fhr)
      g.print('', fhr)
      fhr.close()
      return None

    g.print('Starting Target Material:', fhr)
    g.print('', fhr)
    g.print(this_sim['target_composition_hr'], fhr)
    g.print('', fhr)
    fh = open(data_dir + "/target_material.txt", "w")
    fh.write(this_sim['target_composition_hr'])
    fh.close()

    ##############################################
    # Load trajectory data
    ##############################################

    procede = False
    if(this_sim['exyz'] != None):
      g.print('Loading exyz ' + this_sim['exyz'], fhr)
      exyz_data = exyz.load(this_sim['exyz'])
      g.print('Load complete', fhr)
      g.print('', fhr)
      procede = True

    ##############################################
    # Calculate reaction rates for residuals
    ##############################################

      g.print('Calculate residual reaction rates', fhr)
      residual_rrs = {}    
      particle_rrs = {}
      for icode in this_sim['target_composition'].keys():
        nd = this_sim['target_composition'][icode]['nd']
        sim.run_t(icode, projectile, nd, beam_flux, target_depth, exyz_data, residual_rrs, particle_rrs)

    elif(this_sim['neutron'] != None):
      procede = True
      
      dist = this_sim['neutron']['dist']
      if(this_sim['neutron']['dist'] == 'flat'):
        energy.set(dist='flat', min=this_sim['neutron']['min'], max=this_sim['neutron']['max'])

      neutron_energy = []
      for n in range(this_sim['neutron']['count']):
        neutron_energy.append(energy.rand())
      
      g.print('Calculate residual reaction rates', fhr)
      residual_rrs = {}    
      particle_rrs = {}
      for icode in this_sim['target_composition'].keys():
        nd = this_sim['target_composition'][icode]['nd']
        sim.run_tn(icode, projectile, nd, beam_flux, target_depth, neutron_energy, residual_rrs, particle_rrs)


    g.print('Save reaction rates', fhr)
    g.print('', fhr)
    g.print('', fhr)

    fh = open(data_dir + "/reaction_rates.txt", "w")
    fh.write("#################################################################\n")
    fh.write("Particle & Residual Reaction Rates\n")
    fh.write("#################################################################\n\n")
    g.print('#################################################################\n', fhr)
    g.print("                Particle & Residual Reaction Rates               \n", fhr)
    g.print('#################################################################\n\n', fhr)
    for target in residual_rrs.keys():
      thestr = "Source isotope: " + isotopes.get_hr(target) + "  (" + str(target) + ")" + "\n"
      fh.write(thestr + '\n')
      g.print(thestr, fhr)
      if(target in particle_rrs.keys()):
        for particle in particle_rrs[target].keys():
          thestr = std.pad_right(particle + ":", 20) + str(particle_rrs[target][particle])
          fh.write(thestr + "\n")
          g.print(thestr, fhr)
      for residual in residual_rrs[target].keys():
        thestr = std.pad_right(isotopes.get_hr(residual) + ":", 20) + str(residual_rrs[target][residual])
        fh.write(thestr + "\n")
        g.print(thestr, fhr)
      fh.write("\n\n")   
      g.print('', fhr)   
      g.print('', fhr)      
    fh.close()


    g.print("Calculate residual amounts over time", fhr)
    g.print('', fhr)

    steps_eob = this_sim['inbeampoints']    # From start of beam to end of beam
    steps_eos = this_sim['outbeampoints']     # From end of beam to end of sim
    steps_total = steps_eob + steps_eos - 1
    ta = 0
    tb = steps_eob
    tc = steps_eob - 1
    td = steps_total

    g.print('In beam time steps: ' + str(steps_eob), fhr)
    g.print('Out of beam time steps: ' + str(steps_eos), fhr)
    g.print('', fhr)

    # Time array
    t = numpy.zeros((steps_total,),)
    t[ta:tb] = numpy.linspace(0.0, beam_duration, steps_eob)
    t[tc:td] = numpy.linspace(beam_duration, end_time, steps_eos)

    """
    Isotope amounts:
    n_tree                    target -> residual -> amount array at each time step
    n_by_residual_isotope     
    n_by_target_isotope   

    Isotope activities:
    a_tree                     target -> residual -> amount(time/array)
    a_by_residual_isotope      residual (summed over target isotopes) -> amount(time/array)
    a_by_target_isotope        target (summed over residual isotopes) -> amount(time/array)

    Isotope activities over time:
    a_tree_over_time           time -> target -> residual -> amount
    a_tree_over_time_residual  time -> residual (summed over target isotopes) -> amount
    a_tree_over_time_target    time -> target (summed over residual isotopes) -> amount
    a_tree_over_time_total     time -> amount (summed over all targets and residuals)
    
    Gammas
    g_tree                     target -> residual -> gamma energy ->  
                                               ['activity'] and ['energy'] arrays
    g_by_gk_residual_isotope   residual -> gamma energy -> 
                                               ['activity'] and ['energy'] arrays   
    g_by_gk_target_isotope     target -> gamma energy -> 
                                               ['activity'] and ['energy'] arrays    
    g_by_gk_all                gamma energy -> 
                                               ['activity'] and ['energy'] arrays   
    g_by_residual_isotope      residual -> ['activity'] and ['energy'] arrays  
    g_by_target_isotope        target -> ['activity'] and ['energy'] arrays  
    g_by_all                   ['activity'] and ['energy'] arrays   
    
    Gammas over time
    g_tree_over_time           time -> target -> residual -> gamma energy -> amount 
    g_gk_residual_over_time    time -> residual -> gamma energy -> amount 
    g_gk_target_over_time
    g_gk_total_over_time
    
    """


    results = {
                   # Isotope Amounts
                   'n_tree': {}, 
                   'n_by_residual_isotope': {}, 
                   'n_by_target_isotope': {}, 
                   # Isotope Activities 
                   'a_tree': {}, 
                   'a_by_residual_isotope': {}, 
                   'a_by_target_isotope': {}, 
                   'a_total': None, 
                   # Isotope Activities over time
                   'a_tree_over_time': [], 
                   'a_tree_over_time_residual': [], 
                   'a_tree_over_time_target': [], 
                   'a_tree_over_time_total': [],
                   # Gammas by gamma energy
                   'g_tree': {},
                   'g_by_gk_residual_isotope': {}, 
                   'g_by_gk_target_isotope': {},   
                   'g_by_gk_all': {},   
                   'g_by_residual_isotope': {}, 
                   'g_by_target_isotope': {},   
                   'g_by_all': {},
                   # Gammas over time
                   'g_tree_over_time': [],  
                   'g_gk_residual_over_time': [],
                   'g_gk_target_over_time': [],
                   'g_gk_total_over_time': [],    
                   'g_residual_over_time': [],    
                   'g_target_over_time': [],     
                   'g_all_over_time': [],    
                   # Gammas over time with full details
                   'g_full_details': [],            
                   # For pie charts
                   'a_tree_over_time_residual_pie': [], 
                   'a_tree_over_time_target_pie': [], 
                   'g_residual_energy_pie': [], 
                   'g_target_energy_pie': [],           
                   # Particles
                   'particle_production': {}, 
                   'particle_production_by_particle': {}, 
                   # Totals for charts
                   'total_activity': [],
                   'total_gamma_activity': [],
                   'total_gamma_energy': [],
                  }    
    """
 
    """


    eob_idata = {}     # End of beam data

    # Create dictionary entries
    for target in residual_rrs.keys():
      if(target not in results['n_tree'].keys()):
        results['n_tree'][target] = {}  
      if(target not in eob_idata.keys()):
        eob_idata[target] = {}  

    # In beam activity
    g.print('In beam activity over ' + str(steps_eob) + ' time steps', fhr)
    for target in residual_rrs.keys():
      for residual in residual_rrs[target].keys():
        if(not isotopes.is_stable(residual)):
          idata = {}
          idata[residual] = {'w': residual_rrs[target][residual], 'n0': 0.0}
          parent = residual  
          # In beam calculations
          for i in range(steps_eob):
            time = t[i]

            sim.activity_log.write(str(i) + "   ")
            sim.activity_log.write(str(time) + "   ")
            sim.activity_log.write(str(target) + "   ")
            sim.activity_log.write(str(residual) + "\n")

            r = decay.calculate(parent, time, idata)
            for k in r['tally'].keys():
              if(k not in results['n_tree'][target].keys()):
                results['n_tree'][target][k] = numpy.zeros((steps_total,),)  
              results['n_tree'][target][k][i] = results['n_tree'][target][k][i] + r['tally'][k]['nend']   
              sim.activity_log.write("   ")
              sim.activity_log.write(str(k) + "   ")
              sim.activity_log.write(str(r['tally'][k]['nend']) + "   ")
              sim.activity_log.write(str(results['n_tree'][target][k][i]))
              sim.activity_log.write("\n")
            if(i==(steps_eob-1)):
              eob_idata[target][residual] = {}
              for k in r['tally'].keys():
                eob_idata[target][residual][k] = {}
                eob_idata[target][residual][k]['w'] = 0.0
                eob_idata[target][residual][k]['n0'] = results['n_tree'][target][k][i]

            sim.activity_log.write("\n")

 
    # End of beam activity      
    g.print('In beam activity over ' + str(steps_eos) + ' time steps', fhr)
    for target in residual_rrs.keys():
      for residual in residual_rrs[target].keys():
        if(not isotopes.is_stable(residual)):
          idata = eob_idata[target][residual]
          #print(idata)
          parent = residual  
          # Out of beam calculations
          for i in range(steps_eos-1):        
            time = t[i+steps_eob] - beam_duration

            sim.activity_log.write(str(i) + "   ")
            sim.activity_log.write(str(time) + "   ")
            sim.activity_log.write(str(target) + "   ")
            sim.activity_log.write(str(residual) + "\n")

            r = decay.calculate(parent, time, idata)
            for k in r['tally'].keys():
              if(k not in results['n_tree'][target].keys()):
                results['n_tree'][target][k] = numpy.zeros((steps_total,),)  
              results['n_tree'][target][k][i+steps_eob] = results['n_tree'][target][k][i+steps_eob] + r['tally'][k]['nend']   

              sim.activity_log.write("   ")
              sim.activity_log.write(str(k) + "   ")
              sim.activity_log.write(str(r['tally'][k]['nend']) + "   ")
              sim.activity_log.write(str(results['n_tree'][target][k][i]))
              sim.activity_log.write("\n") 

            sim.activity_log.write("\n")


    n = 0
    for target in results['n_tree'].keys():
      if(target not in sim.isotope_colour.keys()):
        sim.isotope_colour[target] = sim.colours[n]
        sim.isotope_colour[isotopes.get_hr(target)] = sim.colours[n]
        n = n + 1
      for residual in results['n_tree'][target].keys():
        if(residual not in sim.isotope_colour.keys()):
          sim.isotope_colour[residual] = sim.colours[n]
          sim.isotope_colour[isotopes.get_hr(residual)] = sim.colours[n]
          n = n + 1


    # By radioactive isotope
    g.print('Tally amount by residual isotope', fhr)
    for target in results['n_tree'].keys():
      for residual in results['n_tree'][target].keys():
        if(residual not in results['n_by_residual_isotope'].keys()):
          results['n_by_residual_isotope'][residual] = numpy.zeros((steps_total,),)  
        results['n_by_residual_isotope'][residual][:] = results['n_by_residual_isotope'][residual][:] + results['n_tree'][target][residual][:]


    # By target isotope
    g.print("Tally amount of residual isotopes by target isotope", fhr)
    for target in results['n_tree'].keys():
      if(target not in results['n_by_target_isotope'].keys()):
        results['n_by_target_isotope'][target] = numpy.zeros((steps_total,),)  
      for residual in results['n_tree'][target].keys():
        results['n_by_target_isotope'][target][:] = results['n_by_target_isotope'][target][:] + results['n_tree'][target][residual][:]




    ###################
    # Activities
    ###################

    # Activity n_tree
    g.print("Activity per target isotope by residual radioactive isotope over time", fhr)
    for target in results['n_tree'].keys():
      if(target not in results['a_tree'].keys()):
        results['a_tree'][target] = {}
      for residual in results['n_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          l = isotopes.get_decay_constant(residual)
          results['a_tree'][target][residual] = l * results['n_tree'][target][residual]

    # By radioactive isotope
    g.print("Activity by residual radioactive isotope over time", fhr)
    for residual in results['n_by_residual_isotope'].keys():
      if(target not in results['a_by_residual_isotope'].keys()):
        results['a_by_residual_isotope'][target] = numpy.zeros((steps_total,),)  
      if(isotopes.is_unstable(residual)):
        l = isotopes.get_decay_constant(residual)
        results['a_by_residual_isotope'][residual] = l * results['n_by_residual_isotope'][residual]
      
    # By target isotope
    g.print("Activity resulting from each target isotope over time", fhr)
    for target in results['n_tree'].keys():
      if(target not in results['a_by_target_isotope'].keys()):
        results['a_by_target_isotope'][target] = numpy.zeros((steps_total,),)  
      for residual in results['n_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          l = isotopes.get_decay_constant(residual)
          results['a_by_target_isotope'][target][:] = results['a_by_target_isotope'][target][:] + l * results['n_tree'][target][residual][:]
 
    # By total activity
    g.print("Total activity from all targets and residuals over time", fhr)
    results['a_total'] = numpy.zeros((steps_total,),)  
    for target in results['n_tree'].keys():
      for residual in results['n_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          l = isotopes.get_decay_constant(residual)
          results['a_total'] = results['a_total'] + l * results['n_tree'][target][residual][:]


    # Prep lists
    ar_max = []
    at_max = []
    ar_total = []
    at_total = []
    for tn in range(len(t)):
      results['a_tree_over_time'].append({})
      results['a_tree_over_time_residual'].append({})
      results['a_tree_over_time_target'].append({})
      results['a_tree_over_time_residual_pie'].append([[],[]])
      results['a_tree_over_time_target_pie'].append([[],[]])
      ar_max.append(0.0)
      at_max.append(0.0)
      ar_total.append(0.0)
      at_total.append(0.0)

    # Activity over time
    g.print("Activity per target and residual over time", fhr)
    for tn in range(len(t)):
      for target in results['a_tree'].keys():
        results['a_tree_over_time'][tn][target] = {}
        for residual in results['a_tree'][target].keys():
          if(isotopes.is_unstable(residual)):
            results['a_tree_over_time'][tn][target][residual] = results['a_tree'][target][residual][tn]


    g.print("Activity per residual over time", fhr)
    for tn in range(len(t)):
      for target in results['a_tree'].keys():
        for residual in results['a_tree'][target].keys():
          if(isotopes.is_unstable(residual)):
            if(residual not in results['a_tree_over_time_residual'][tn].keys()):
              results['a_tree_over_time_residual'][tn][residual] = 0.0
            results['a_tree_over_time_residual'][tn][residual] = results['a_tree_over_time_residual'][tn][residual] + results['a_tree'][target][residual][tn]
            ar_max[tn] = max(ar_max[tn], results['a_tree_over_time_residual'][tn][residual])
            ar_total[tn] = ar_total[tn] + results['a_tree'][target][residual][tn]

    g.print("Activity per target over time", fhr)
    for tn in range(len(t)):
      for target in results['a_tree'].keys():
        if(target not in results['a_tree_over_time_target'][tn].keys()):
          results['a_tree_over_time_target'][tn][target] = 0.0
        for residual in results['a_tree'][target].keys():
          if(isotopes.is_unstable(residual)):
            results['a_tree_over_time_target'][tn][target] = results['a_tree_over_time_target'][tn][target] + results['a_tree'][target][residual][tn]
            at_max[tn] = max(at_max[tn], results['a_tree_over_time_target'][tn][target])
            at_total[tn] = at_total[tn] + results['a_tree'][target][residual][tn]
    
    
    g.print("Activity data for pie charts", fhr)
    for tn in range(len(t)):
      for residual in results['a_tree_over_time_residual'][tn].keys():
        if(results['a_tree_over_time_residual'][tn][residual] > sim.piechart_threshold * ar_max[tn]):
          results['a_tree_over_time_residual_pie'][tn][0].append(isotopes.get_hr(residual))
          results['a_tree_over_time_residual_pie'][tn][1].append(results['a_tree_over_time_residual'][tn][residual])
    
    for tn in range(len(t)):
      for target in results['a_tree_over_time_target'][tn].keys():
        #if(results['a_tree_over_time_target'][tn][target] > 0.05 * at_max[tn]):
        results['a_tree_over_time_target_pie'][tn][0].append(isotopes.get_hr(target))
        results['a_tree_over_time_target_pie'][tn][1].append(results['a_tree_over_time_target'][tn][target])
    



    ###################
    # Gammas
    ###################

    # Tree

    # Target Residual Gamma tree 
    g.print("Tally gammas by target -> residual -> energy", fhr)
    for target in results['a_tree'].keys():
      if(target not in results['g_tree'].keys()):
        results['g_tree'][target] = {}
      for residual in results['a_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          results['g_tree'][target][residual] = {}          
          # Get gamma list
          ga = isotopes.get_gammas_array(residual)
          if(not isinstance(ga, type(None))):
            for n in range(len(ga)):
              gk = ga[n,0]
              g_energy_ev = ga[n,0]
              g_intensity = ga[n,1]
              if(gk not in results['g_tree'][target][residual].keys()):
                results['g_tree'][target][residual][gk] = {'activity': None, 'energy': None,}
                results['g_tree'][target][residual][gk]['activity'] = results['a_tree'][target][residual][:] * g_intensity
                results['g_tree'][target][residual][gk]['energy'] = g_energy_ev * results['g_tree'][target][residual][gk]['activity'][:]


    # Individual gamma energies

    # By residual isotope
    g.print("Tally gammas by residual -> energy", fhr)
    for target in results['g_tree'].keys():
      for residual in results['g_tree'][target].keys():
        if(residual not in results['g_by_gk_residual_isotope']):
          results['g_by_gk_residual_isotope'][residual] = {}
        if(isotopes.is_unstable(residual)):
          for gk in results['g_tree'][target][residual].keys():
            if(gk not in results['g_by_gk_residual_isotope'][residual].keys()):
              results['g_by_gk_residual_isotope'][residual][gk] = {'activity': numpy.zeros((steps_total,),), 'energy': numpy.zeros((steps_total,),),}
            # Activity
            results['g_by_gk_residual_isotope'][residual][gk]['activity'][:] = results['g_by_gk_residual_isotope'][residual][gk]['activity'][:] + results['g_tree'][target][residual][gk]['activity'][:]
            # Energy 
            results['g_by_gk_residual_isotope'][residual][gk]['energy'][:] = results['g_by_gk_residual_isotope'][residual][gk]['energy'][:] + results['g_tree'][target][residual][gk]['energy'][:]


    # By target isotope
    g.print("Tally gammas by target -> energy", fhr)
    for target in results['g_tree'].keys():
      if(target not in results['g_by_gk_target_isotope']):
        results['g_by_gk_target_isotope'][target] = {}
        for residual in results['g_tree'][target].keys():
          if(isotopes.is_unstable(residual)):
            for gk in results['g_tree'][target][residual].keys():
              if(gk not in results['g_by_gk_target_isotope'][target].keys()):
                results['g_by_gk_target_isotope'][target][gk] = {'activity': numpy.zeros((steps_total,),), 'energy': numpy.zeros((steps_total,),),}
              # Activity
              results['g_by_gk_target_isotope'][target][gk]['activity'][:] = results['g_by_gk_target_isotope'][target][gk]['activity'][:] + results['g_tree'][target][residual][gk]['activity'][:]
              # Energy 
              results['g_by_gk_target_isotope'][target][gk]['energy'][:] = results['g_by_gk_target_isotope'][target][gk]['energy'][:] + results['g_tree'][target][residual][gk]['energy'][:]

   
    # Gammas - tally of all
    g.print("Tally gammas -> by gamma energy", fhr)
    for target in results['g_tree'].keys():
      for residual in results['g_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          for gk in results['g_tree'][target][residual].keys():
            if(gk not in results['g_by_gk_all'].keys()):
              results['g_by_gk_all'][gk] = {'activity': numpy.zeros((steps_total,),), 'energy': numpy.zeros((steps_total,),),}
            # Activity
            results['g_by_gk_all'][gk]['activity'][:] = results['g_by_gk_all'][gk]['activity'][:] + results['g_tree'][target][residual][gk]['activity'][:]
            # Energy 
            results['g_by_gk_all'][gk]['energy'][:] = results['g_by_gk_all'][gk]['energy'][:] + results['g_tree'][target][residual][gk]['energy'][:]
    


    # Tallied over gamma energies

    # By residual isotope
    g.print("Tally gammas by residual and gamma energies", fhr)
    for target in results['g_tree'].keys():
      for residual in results['g_tree'][target].keys():
        if(residual not in results['g_by_residual_isotope']):
          results['g_by_residual_isotope'][residual] = {'activity': numpy.zeros((steps_total,),), 'energy': numpy.zeros((steps_total,),),}
        if(isotopes.is_unstable(residual)):
          for gk in results['g_tree'][target][residual].keys():
            # Activity
            results['g_by_residual_isotope'][residual]['activity'][:] = results['g_by_residual_isotope'][residual]['activity'][:] + results['g_tree'][target][residual][gk]['activity'][:]
            # Energy 
            results['g_by_residual_isotope'][residual]['energy'][:] = results['g_by_residual_isotope'][residual]['energy'][:] + results['g_tree'][target][residual][gk]['energy'][:]


    # By target isotope
    g.print("Tally gammas by target and gamma energies", fhr)
    for target in results['g_tree'].keys():
      if(target not in results['g_by_target_isotope']):
        results['g_by_target_isotope'][target] = {'activity': numpy.zeros((steps_total,),), 'energy': numpy.zeros((steps_total,),),}
      for residual in results['g_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          for gk in results['g_tree'][target][residual].keys():
            # Activity
            results['g_by_target_isotope'][target]['activity'][:] = results['g_by_target_isotope'][target]['activity'][:] + results['g_tree'][target][residual][gk]['activity'][:]
            # Energy 
            results['g_by_target_isotope'][target]['energy'][:] = results['g_by_target_isotope'][target]['energy'][:] + results['g_tree'][target][residual][gk]['energy'][:]


    # By total
    g.print("Tally gammas by target, residual and gamma energies", fhr)
    results['g_by_all'] = {'activity': numpy.zeros((steps_total,),), 'energy': numpy.zeros((steps_total,),),}
    for target in results['g_tree'].keys():
      for residual in results['g_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          for gk in results['g_tree'][target][residual].keys():
            # Activity
            results['g_by_all']['activity'][:] = results['g_by_all']['activity'][:] + results['g_tree'][target][residual][gk]['activity'][:]
            # Energy 
            results['g_by_all']['energy'][:] = results['g_by_all']['energy'][:] + results['g_tree'][target][residual][gk]['energy'][:]





    # Gammas at each time step

    # Gammas lines over time
    g.print("Per target and residual gamma lines over time", fhr)
    for tn in range(len(t)):
      results['g_tree_over_time'].append({})
    for tn in range(len(t)):
      for target in results['g_tree'].keys():
        results['g_tree_over_time'][tn][target] = {}
        for residual in results['g_tree'][target].keys():
          results['g_tree_over_time'][tn][target][residual] = {}
          if(isotopes.is_unstable(residual)):
            for gk in results['g_tree'][target][residual].keys():
              results['g_tree_over_time'][tn][target][residual][gk] = {}
              results['g_tree_over_time'][tn][target][residual][gk]['activity'] = results['g_tree'][target][residual][gk]['activity'][tn]
              results['g_tree_over_time'][tn][target][residual][gk]['energy'] = results['g_tree'][target][residual][gk]['energy'][tn]


    # By residual over time
    g.print("Per residual gamma lines over time", fhr)
    for tn in range(len(t)):
      results['g_gk_residual_over_time'].append({})
    for tn in range(len(t)):
      for target in results['g_tree_over_time'][tn].keys():
        for residual in results['g_tree_over_time'][tn][target].keys():
          if(residual not in results['g_gk_residual_over_time'][tn].keys()):
            results['g_gk_residual_over_time'][tn][residual] = {}
          for gk in results['g_tree_over_time'][tn][target][residual].keys():
            if(gk not in results['g_gk_residual_over_time'][tn][residual].keys()):
              results['g_gk_residual_over_time'][tn][residual][gk] = {'activity': 0.0, 'energy': 0.0}
            results['g_gk_residual_over_time'][tn][residual][gk]['activity'] = results['g_gk_residual_over_time'][tn][residual][gk]['activity'] + results['g_tree'][target][residual][gk]['activity'][tn]
            results['g_gk_residual_over_time'][tn][residual][gk]['energy'] = results['g_gk_residual_over_time'][tn][residual][gk]['energy'] + results['g_tree'][target][residual][gk]['energy'][tn]
    

    # By target over time
    g.print("Per target gamma lines over time", fhr)
    for tn in range(len(t)):
      results['g_gk_target_over_time'].append({})
    for tn in range(len(t)):
      for target in results['g_tree_over_time'][tn].keys():
        if(target not in results['g_gk_target_over_time'][tn].keys()):
          results['g_gk_target_over_time'][tn][target] = {}
        for residual in results['g_tree_over_time'][tn][target].keys():
          for gk in results['g_tree_over_time'][tn][target][residual].keys():
            if(gk not in results['g_gk_target_over_time'][tn][target].keys()):
              results['g_gk_target_over_time'][tn][target][gk] = {'activity': 0.0, 'energy': 0.0}
            results['g_gk_target_over_time'][tn][target][gk]['activity'] = results['g_gk_target_over_time'][tn][target][gk]['activity'] + results['g_tree'][target][residual][gk]['activity'][tn]
            results['g_gk_target_over_time'][tn][target][gk]['energy'] = results['g_gk_target_over_time'][tn][target][gk]['energy'] + results['g_tree'][target][residual][gk]['energy'][tn]



    # Total over time
    g.print("Tallied gamma lines over time", fhr)
    for tn in range(len(results['g_tree_over_time'])):
      results['g_gk_total_over_time'].append({})
      for target in results['g_tree_over_time'][tn].keys():
        for residual in results['g_tree_over_time'][tn][target].keys():
          for gk in results['g_tree_over_time'][tn][target][residual].keys():
            if(gk not in results['g_gk_total_over_time'][tn].keys()):
              results['g_gk_total_over_time'][tn][gk] = {'activity': 0.0, 'energy': 0.0}
            results['g_gk_total_over_time'][tn][gk]['activity'] = results['g_gk_total_over_time'][tn][gk]['activity'] + results['g_tree'][target][residual][gk]['activity'][tn]
            results['g_gk_total_over_time'][tn][gk]['energy'] = results['g_gk_total_over_time'][tn][gk]['energy'] + results['g_tree'][target][residual][gk]['energy'][tn]


    # By residual over time
    g.print("Per residual gamma lines over time", fhr)
    for tn in range(len(t)):
      results['g_residual_over_time'].append({})
    for tn in range(len(t)):
      for target in results['g_tree_over_time'][tn].keys():
        for residual in results['g_tree_over_time'][tn][target].keys():
          if(residual not in results['g_residual_over_time'][tn].keys()):
            results['g_residual_over_time'][tn][residual] = {'activity': 0.0, 'energy': 0.0}
          for gk in results['g_tree_over_time'][tn][target][residual].keys():
            results['g_residual_over_time'][tn][residual]['activity'] = results['g_residual_over_time'][tn][residual]['activity'] + results['g_tree'][target][residual][gk]['activity'][tn]
            results['g_residual_over_time'][tn][residual]['energy'] = results['g_residual_over_time'][tn][residual]['energy'] + results['g_tree'][target][residual][gk]['energy'][tn]


    # By target over time
    g.print("Per target gamma lines over time", fhr)
    for tn in range(len(t)):
      results['g_target_over_time'].append({})
    for tn in range(len(t)):
      for target in results['g_tree_over_time'][tn].keys():
        if(target not in results['g_target_over_time'][tn].keys()):
          results['g_target_over_time'][tn][target] = {'activity': 0.0, 'energy': 0.0}
        for residual in results['g_tree_over_time'][tn][target].keys():
          for gk in results['g_tree_over_time'][tn][target][residual].keys():
            results['g_target_over_time'][tn][target]['activity'] = results['g_target_over_time'][tn][target]['activity'] + results['g_tree'][target][residual][gk]['activity'][tn]
            results['g_target_over_time'][tn][target]['energy'] = results['g_target_over_time'][tn][target]['energy'] + results['g_tree'][target][residual][gk]['energy'][tn]


    # All over time
    g.print("Total activity and energy of gammas over time", fhr)
    for tn in range(len(t)):
      results['g_all_over_time'].append({'activity': 0.0, 'energy': 0.0})
    for tn in range(len(t)):
      for target in results['g_tree_over_time'][tn].keys():
        for residual in results['g_tree_over_time'][tn][target].keys():
          for gk in results['g_tree_over_time'][tn][target][residual].keys():
            results['g_all_over_time'][tn]['activity'] = results['g_all_over_time'][tn]['activity'] + results['g_tree'][target][residual][gk]['activity'][tn]
            results['g_all_over_time'][tn]['energy'] = results['g_all_over_time'][tn]['energy'] + results['g_tree'][target][residual][gk]['energy'][tn]



    # Gamma full details over time
    g.print("Gamma full details over time", fhr)
    results['g_full_details'] = []
    results['g_full_details_key'] = []
    for tn in range(len(t)):
      results['g_full_details'].append({})
      results['g_full_details_key'].append([])
    for tn in range(len(t)):
      for target in results['a_tree_over_time'][tn].keys():
        for residual in results['a_tree_over_time'][tn][target].keys():
          if(isotopes.is_unstable(residual)):          
            # Get data
            thalf = isotopes.get_half_life(residual)
            ga = isotopes.get_gammas_array(residual)
            if(not isinstance(ga, type(None))):
              for n in range(len(ga)):
                gk = str(ga[n,0]) 
                while(len(gk)<=10):
                  gk = "0" + gk
                gk = gk + "_" + str(residual)

                g_energy_ev = ga[n,0]
                g_intensity = ga[n,1]
                
                if(gk not in results['g_full_details'][tn].keys()):
                  results['g_full_details'][tn][gk] = {}
                  results['g_full_details'][tn][gk]['energy'] = g_energy_ev
                  results['g_full_details'][tn][gk]['intensity'] = g_intensity
                  results['g_full_details'][tn][gk]['residual'] = str(residual) + ',' + isotopes.get_hr(residual)
                  results['g_full_details'][tn][gk]['residual_half_life'] = thalf
                  results['g_full_details'][tn][gk]['activity'] = 0.0
                  results['g_full_details'][tn][gk]['energy_per_second'] = 0.0
                  results['g_full_details'][tn][gk]['target'] = []
                  results['g_full_details'][tn][gk]['target_activity'] = []

                a = results['a_tree_over_time'][tn][target][residual]
                results['g_full_details'][tn][gk]['activity'] = results['g_full_details'][tn][gk]['activity'] + a * g_intensity
                results['g_full_details'][tn][gk]['energy_per_second'] = results['g_full_details'][tn][gk]['energy_per_second'] + g_energy_ev * a * g_intensity
                results['g_full_details'][tn][gk]['target'].append(target)
                results['g_full_details'][tn][gk]['target_activity'].append(a * g_intensity)

 
    # 
    #results['g_residual_over_time'][tn][residual]['energy']
    g.print("Gamma data for pie charts", fhr)
    r_max = []
    for tn in range(len(t)):
      r_max.append(0.0)
      results['g_residual_energy_pie'].append([[],[]])
    for tn in range(len(t)):
      for residual in results['g_residual_over_time'][tn].keys():
        r_max[tn] = max(r_max[tn], results['g_residual_over_time'][tn][residual]['energy'])
    for tn in range(len(t)):
      for residual in results['g_residual_over_time'][tn].keys():
        if(results['g_residual_over_time'][tn][residual]['energy'] > sim.piechart_threshold * r_max[tn]):
          results['g_residual_energy_pie'][tn][0].append(isotopes.get_hr(residual))
          results['g_residual_energy_pie'][tn][1].append( results['g_residual_over_time'][tn][residual]['energy'])

    r_max = []
    for tn in range(len(t)):
      r_max.append(0.0)
      results['g_target_energy_pie'].append([[],[]])
    for tn in range(len(t)):
      for target in results['g_target_over_time'][tn].keys():
        r_max[tn] = max(r_max[tn], results['g_target_over_time'][tn][target]['energy'])
    for tn in range(len(t)):
      for target in results['g_target_over_time'][tn].keys():
        #if(results['g_target_over_time'][tn][target]['energy'] > 0.05 * r_max[tn]):
        results['g_target_energy_pie'][tn][0].append(isotopes.get_hr(target))
        results['g_target_energy_pie'][tn][1].append( results['g_target_over_time'][tn][target]['energy'])




    









    #########################
    # Particle Production
    #########################

    for target in particle_rrs.keys():
      results['particle_production'][target] = {}
      for particle in particle_rrs[target].keys():
        results['particle_production'][target][particle] = numpy.zeros((steps_total,),)  
        for tn in range(len(t)):
          if(t[tn] > beam_duration):
            results['particle_production'][target][particle][tn] = 0.0
          else:
            results['particle_production'][target][particle][tn] = particle_rrs[target][particle]
            
    for target in particle_rrs.keys():
      for particle in particle_rrs[target].keys():
        if(particle not in results['particle_production_by_particle'].keys()):
          results['particle_production_by_particle'][particle] = numpy.zeros((steps_total,),)  
          for tn in range(len(t)):
            if(t[tn] <= beam_duration):
              results['particle_production_by_particle'][particle][tn] = results['particle_production_by_particle'][particle][tn] + particle_rrs[target][particle]
                






    #########################
    # Totals
    #########################
    
    g.print("Activity and Gamma totals over time for charts", fhr)
    for tn in range(len(t)):
      results['total_activity'].append(results['a_total'][tn])
      results['total_gamma_activity'].append(results['g_all_over_time'][tn]['activity'])
      results['total_gamma_energy'].append(results['g_all_over_time'][tn]['energy'])

    g.print('', fhr)
    g.print('', fhr)
    time_list = []

    for i in range(10):
      tn = int(1 + i * ((steps_total - 2) / 9))
      time_list.append(tn)
    if(steps_eob-1 not in time_list):
      time_list.append(steps_eob-1)
    if(steps_total-1 not in time_list):
      time_list.append(steps_total-1)
    time_list.sort()
    
    fhr.write('\n')
    fhr.write('#######################################################\n')
    fhr.write('                      Results                          \n')
    fhr.write('#######################################################\n')
    fhr.write('\n')

    for tn in time_list:
      fhr.write('\n')
      fhr.write('##########################################################')
      fhr.write('##########################################################\n')
      fhr.write(std.pad_right("Time", 23))
      fhr.write(std.pad_right("Target", 23))
      fhr.write(std.pad_right("Residual", 23))
      fhr.write(std.pad_right("Activity", 23))
      fhr.write(std.pad_right("Half-life/s", 23))
      fhr.write('\n')
      fhr.write('##########################################################')
      fhr.write('##########################################################\n')
      for target in results['a_tree_over_time'][tn].keys():
        fhr.write('\n')
        for residual in results['a_tree_over_time'][tn][target].keys():
          fhr.write(std.pad_right(str(t[tn]) + 's  (' + str(round(t[tn] / 3600,3)) + 'hr)', 22))
          fhr.write(' ')
          fhr.write(std.pad_right(isotopes.get_hr(target), 22))
          fhr.write(' ')
          fhr.write(std.pad_right(isotopes.get_hr(residual), 22))
          fhr.write(' ')
          fhr.write(std.pad_right(results['a_tree_over_time'][tn][target][residual], 22))
          fhr.write(' ')
          fhr.write(std.pad_right(isotopes.get_half_life(residual), 22))
          fhr.write('\n')
      fhr.write('\n\n')
    fhr.write('\n')



    fhs = open('summary.csv', 'a')
    fhs.write('\n')
    fhs.write(str(counter))
    fhs.write(',')
    fhs.write(str(beam_flux))
    fhs.write(',')
    fhs.write(str(target_depth))
    fhs.write(',')
    fhs.write(projectile)
    fhs.write(',')
    fhs.write(str(beam_duration))
    fhs.write(',')
    fhs.write(str(end_time))
    fhs.write(',')
    fhs.write(str(results['total_activity'][tc]))
    fhs.write(',')
    fhs.write(str(results['total_gamma_activity'][tc]))
    fhs.write(',')
    fhs.write(str(1.60218E-16 * results['total_gamma_energy'][tc]))
    fhs.write(',')
    fhs.write(str(results['total_activity'][td-1]))
    fhs.write(',')
    fhs.write(str(results['total_gamma_activity'][td-1]))
    fhs.write(',')
    fhs.write(str(1.60218E-16 * results['total_gamma_energy'][td-1]))
    fhs.close()

    
    ############################################################################
    ############################################################################
    # Save Data
    ############################################################################
    ############################################################################

    g.print("Save data files - isotope amounts")
    dir_n_tree = data_dir + "/n_tree"
    std.make_dir(dir_n_tree)
    for target in results['n_tree'].keys():
      idir = dir_n_tree + "/" + isotopes.get_dir_hr(target)
      std.make_dir(idir)
      for residual in results['n_tree'][target].keys():
        out_array = numpy.zeros((steps_total,2,),) 
        out_array[:,0] = t[:]
        out_array[:,1] = results['n_tree'][target][residual][:]
        file_name = idir + "/" + isotopes.get_dir_hr(residual) + ".txt"
        numpy.savetxt(file_name, out_array, delimiter=",")    

    dir_n_tree = data_dir + "/n_by_residual_isotope"
    std.make_dir(dir_n_tree)
    for residual in results['n_by_residual_isotope'].keys():
      out_array = numpy.zeros((steps_total,2,),) 
      out_array[:,0] = t[:]
      out_array[:,1] = results['n_by_residual_isotope'][residual][:]
      file_name = dir_n_tree + "/" + isotopes.get_dir_hr(residual) + ".txt"
      numpy.savetxt(file_name, out_array, delimiter=",")

    dir_n_tree = data_dir + "/n_by_target_isotope"
    std.make_dir(dir_n_tree)
    for target in results['n_by_target_isotope'].keys():
      out_array = numpy.zeros((steps_total,2,),) 
      out_array[:,0] = t[:]
      out_array[:,1] = results['n_by_target_isotope'][target][:]
      file_name = dir_n_tree + "/" + isotopes.get_dir_hr(target) + ".txt"
      numpy.savetxt(file_name, out_array, delimiter=",")



    g.print("Save data files - isotope activities")
    dir_a_tree = data_dir + "/a_tree"
    std.make_dir(dir_a_tree)
    for target in results['a_tree'].keys():
      idir = dir_a_tree + "/" + isotopes.get_dir_hr(target)
      std.make_dir(idir)
      for residual in results['a_tree'][target].keys():
        out_array = numpy.zeros((steps_total,2,),) 
        out_array[:,0] = t[:]
        out_array[:,1] = results['a_tree'][target][residual][:]
        file_name = idir + "/" + isotopes.get_dir_hr(residual) + ".txt"
        numpy.savetxt(file_name, out_array, delimiter=",")

    g.print("Save data files - isotope activities")
    dir_n_tree = data_dir + "/a_by_residual_isotope"
    std.make_dir(dir_n_tree)
    for residual in results['a_by_residual_isotope'].keys():
      out_array = numpy.zeros((steps_total,2,),) 
      out_array[:,0] = t[:]
      out_array[:,1] = results['a_by_residual_isotope'][residual][:]
      file_name = dir_n_tree + "/" + isotopes.get_dir_hr(residual) + ".txt"
      numpy.savetxt(file_name, out_array, delimiter=",")

    dir_n_tree = data_dir + "/a_by_target_isotope"
    std.make_dir(dir_n_tree)
    for target in results['a_by_target_isotope'].keys():
      out_array = numpy.zeros((steps_total,2,),) 
      out_array[:,0] = t[:]
      out_array[:,1] = results['a_by_target_isotope'][target][:]
      file_name = dir_n_tree + "/" + isotopes.get_dir_hr(target) + ".txt"
      numpy.savetxt(file_name, out_array, delimiter=",")
    

    g.print("Save data files - tally of gamma lines by gamma energy over time")    
    gammas_data_dir = data_dir + "/gammas/over_time"
    std.make_dir(gammas_data_dir)

    for tn in range(len(t)):
      tn_str = str(tn)
      while(len(tn_str)<4):
        tn_str = "0" + tn_str
      tn_str = tn_str + "_" + str(int(t[tn])) + ".csv"

      g_out = numpy.zeros((len(results['g_gk_total_over_time'][tn]), 2),)
      gn = 0
      for gk in results['g_gk_total_over_time'][tn].keys():
        g_out[gn, 0] = gk / 1000.0  # Convert to KeV
        g_out[gn, 1] = results['g_gk_total_over_time'][tn][gk]['activity']
        gn = gn + 1

      g_out = g_out[g_out[:, 0].argsort(), :]
      file_name = gammas_data_dir + "/" + tn_str
      numpy.savetxt(file_name, g_out, delimiter=",")

    g.print("Save data files - tally of gamma lines by gamma energy over time")    
    gammas_data_dir = data_dir + "/gammas/over_time_full_details"
    std.make_dir(gammas_data_dir)

    for tn in range(len(t)):
      tn_str = str(tn)
      while(len(tn_str)<4):
        tn_str = "0" + tn_str
      tn_str = tn_str + "_" + str(int(t[tn])) + ".csv"

      fh = open(gammas_data_dir + "/" + tn_str, 'w')
      fh.write("Energy,")
      fh.write("Intensity,")
      fh.write("Residual Code,")
      fh.write("Residual,")
      fh.write("Residual Half Life,")
      fh.write("Activity,")
      fh.write("Energy per Second,")
      fh.write("Energy per Second Mw,")
      fh.write("Targets,")
      for target in results['a_tree'].keys():
        fh.write(",,")
      fh.write("\n")
 
      for key in sorted(results['g_full_details'][tn]):
        fh.write(str(results['g_full_details'][tn][key]['energy']) + ",")
        fh.write(str(results['g_full_details'][tn][key]['intensity']) + ",")
        fh.write(str(results['g_full_details'][tn][key]['residual']) + ",")
        fh.write(str(results['g_full_details'][tn][key]['residual_half_life']) + ",")
        fh.write(str(results['g_full_details'][tn][key]['activity']) + ",")
        fh.write(str(results['g_full_details'][tn][key]['energy_per_second']) + ",")
        fh.write(str(results['g_full_details'][tn][key]['energy_per_second'] * 1.60218e-16) + ",")
        for n in range(len(results['g_full_details'][tn][key]['target'])): 
          fh.write(str(results['g_full_details'][tn][key]['target'][n]) + ",")
          fh.write(str(results['g_full_details'][tn][key]['target_activity'][n]) + ",")
        fh.write('\n')
      fh.close()


      """
      g_out = numpy.zeros((len(results['g_gk_total_over_time'][tn]), 2),)
      gn = 0
      for gk in results['g_gk_total_over_time'][tn].keys():
        g_out[gn, 0] = gk / 1000.0  # Convert to KeV
        g_out[gn, 1] = results['g_gk_total_over_time'][tn][gk]['activity']
        gn = gn + 1

      g_out = g_out[g_out[:, 0].argsort(), :]
      file_name = gammas_data_dir + "/" + tn_str
      numpy.savetxt(file_name, g_out, delimiter=",")
      """



    
    dir_gamma_intensities = data_dir + "/gamma_intensities"
    std.make_dir(dir_gamma_intensities)
    for residual in results['a_by_residual_isotope'].keys():
      g_list = []
      ga = isotopes.get_gammas_array(residual)
      if(not isinstance(ga, type(None))):
        for n in range(len(ga)):
          gk = ga[n,0]
          g_energy_ev = ga[n,0] 
          g_intensity = ga[n,1]      
          g_list.append([gk, g_intensity])

        g_arr = numpy.zeros((len(g_list),len(g_list[0]),),)
        for n in range(len(g_list)):
          g_arr[n,:] = g_list[n][:]
        if(len(g_list) > 1):
          g_arr = g_arr[g_arr[:, 0].argsort(), :]
        file_name = dir_gamma_intensities + "/" + isotopes.get_dir_hr(residual) + ".csv"
        numpy.savetxt(file_name, g_arr, delimiter=",")    


    
    dir_data_totals = data_dir + "/totals"
    std.make_dir(dir_data_totals) 
    out_arr = numpy.zeros((len(t),4,),)
    out_arr[:,0] = t[:]
    out_arr[:,1] = results['a_total'][:]
    out_arr[:,2] = results['g_by_all']['energy'][:]
    out_arr[:,3] = results['g_by_all']['activity'][:]
    numpy.savetxt(dir_data_totals + "/totals.csv", out_arr, delimiter=",") 





    activity_a = results['total_activity'][tc]     
    gamma_energy_a = results['total_gamma_energy'][tc]
    gamma_power_a = 1.60218E-16 * results['total_gamma_energy'][tc]
    dose_s_a = (1.60218E-19 * results['total_gamma_energy'][tc] * (1 / (12.57 * 80))) / 1000
    dose_hr_a = (1.60218E-19 * results['total_gamma_energy'][tc] * 3600 * (1 / (12.57 * 80))) / 1000


    activity_b = results['total_activity'][td-1] 
    gamma_energy_b = results['total_gamma_energy'][td-1]
    gamma_power_b = 1.60218E-16 * results['total_gamma_energy'][td-1]
    dose_s_b = (1.60218E-19 * results['total_gamma_energy'][td-1] * (1 / (12.57 * 80))) / 1000
    dose_hr_b = (1.60218E-19 * results['total_gamma_energy'][td-1] * 3600 * (1 / (12.57 * 80))) / 1000

    
    output = """
Gamma Dose - Beam End  
===================================================
Activity (Bq)              {0:.4e}  
Gamma Energy (eV)          {1:.4e}    
Gamma Energy (mW)          {2:.4e}   
Absorbed Dose (mGy/s)      {3:.4e}    
Absorbed Dose (mGy/hr)     {4:.4e}                


Gamma Dose - Sim End
===================================================
Activity (Bq)              {5:.4e}
Gamma Energy (eV)          {6:.4e}    
Gamma Energy (mW)          {7:.4e}     
Absorbed Dose (mGy/s)      {8:.4e}     
Absorbed Dose (mGy/hr)     {9:.4e}   


Absorbed Dose Calculations
===================================================
Absorbed dose assumptions:
1. radiation from point, emitted isotropically
2. 80Kg human
3. 1m from point source
4. 1m squared surface area
5. all energy absorbed


Dose Limits
===================================================
employees 18+             20 millisieverts/year
trainees 18+              6 millisieverts/year
public and under 18s      1 millisievert/year
public and under 18s      1.140771128E-04 millisieverts/hour
public and under 18s      3.2E-08 millisieverts/second

Dose averaged over area of skin not exceeding 1cm2
Source: http://www.hse.gov.uk/radiation/ionising/doses/

1mW gammas completely absorbed by 80kg human will give public limit for entire year over 80 seconds.

Dose under 1W gammas absorbed for 1s no imediate changes (under 100 Rads)
Dose 1W-2W gammas absorbed for 1s acute radiation syndrome (100-200 Rads)
Dose 10W+ gammas absorbed for 1s mostly fatal (1000+ Rads)
    """

    output = output.format(activity_a, gamma_energy_a, gamma_power_a, dose_s_a, dose_hr_a, 
                           activity_b, gamma_energy_b, gamma_power_b, dose_s_b, dose_hr_b)
    g.print(output, fhr)
    fh = open(data_dir + '/absorbed_dose.txt', 'w')  
    fh.write(output)
    fh.close()

    exit()


    ############################################################################
    ############################################################################
    # Save Plots
    ############################################################################
    ############################################################################

    ###################
    # Plots - Activities only
    ###################
    
    g.print("Make plots")

    # Individual isotopes
    this_plot_dir = plot_dir + "/totals"
    std.make_dir(this_plot_dir)


    # Activity total over time
    g.print("Plots - activity total over time")
    plot_file = this_plot_dir + "/total_activity.eps"
    plot_list = []
    plot_array = numpy.zeros((steps_total,2,),) 
    plot_array[:,0] = t[:]
    plot_array[:,1] = results['total_activity'][:]
    plot_list.append(["Total gamma activity", plot_array])
    title = "Total activity over time"
    sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Activity (Bq)", 'log')


    # Activity total over time
    g.print("Plots - activity total over time")
    plot_file = this_plot_dir + "/total_activity_beam_on.eps"
    plot_list = []
    plot_array = numpy.zeros((steps_eob,2,),) 
    plot_array[:,0] = t[ta:tb]
    plot_array[:,1] = results['total_activity'][ta:tb]
    plot_list.append(["Total gamma activity", plot_array])
    title = "Total activity over time"
    sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Activity (Bq)", 'log')


    # Activity total over time
    g.print("Plots - activity total over time")
    plot_file = this_plot_dir + "/total_activity_beam_off.eps"
    plot_list = []
    plot_array = numpy.zeros((steps_eos,2,),) 
    plot_array[:,0] = t[tc:td]
    plot_array[:,1] = results['total_activity'][tc:td]
    plot_list.append(["Total gamma activity", plot_array])
    title = "Total activity over time"
    sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Activity (Bq)", 'log')


    # Activity n_tree
    g.print("Plots - activity by target and residual isotope")
    plot_file = this_plot_dir + "/total_gamma_activity.eps"
    plot_list = []
    plot_array = numpy.zeros((steps_total,2,),) 
    plot_array[:,0] = t[:]
    plot_array[:,1] = results['total_gamma_activity'][:]
    plot_list.append(["Total gamma activity", plot_array])
    title = "Total activity over time"
    sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Activity (Bq)", 'log')


    # Activity n_tree
    g.print("Plots - activity by target and residual isotope")
    plot_file = this_plot_dir + "/total_gamma_power.eps"
    plot_list = []
    plot_array = numpy.zeros((steps_total,2,),) 
    plot_array[:,0] = t[:]
    plot_array[:,1] = results['total_gamma_energy'][:]
    plot_list.append(["Total gamma output", plot_array])
    title = "Total gamma energy per second over time"
    sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Gamma Power (eV/s)", 'log')
   
    # Activity n_tree
    g.print("Plots - activity by target and residual isotope")
    plot_file = this_plot_dir + "/total_gamma_power_mw.eps"
    plot_list = []
    plot_array = numpy.zeros((steps_total,2,),) 
    plot_array[:,0] = t[:]
    plot_array[:,1] = results['total_gamma_energy'][:] 
    for n in range(len(plot_array[:,1])):
      plot_array[n,1] = plot_array[n,1] * 1.60218e-16
    plot_list.append(["Total gamma output", plot_array])
    title = "Total gamma power over time"
    sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Gamma Power (mW)", 'log')
   


    # Individual isotopes
    this_plot_dir = plot_dir + "/by_isotope_1"
    std.make_dir(this_plot_dir)

    # Activity n_tree
    g.print("Plots - activity by target and residual isotope")
    for target in results['a_tree'].keys():
      dir_a = this_plot_dir + "/"  + isotopes.get_dir_hr(target)
      std.make_dir(dir_a)
      for residual in results['a_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          plot_file = dir_a + "/" + isotopes.get_dir_hr(residual) + ".eps"
          plot_list = []
          plot_array = numpy.zeros((steps_total,2,),) 
          plot_array[:,0] = t[:]
          plot_array[:,1] = results['a_tree'][target][residual][:]
          plot_list.append([isotopes.get_hr(residual) + " from " + isotopes.get_hr(target), plot_array])
          title = "Activity - target: " + isotopes.get_hr(target) + "  residual: " + isotopes.get_hr(residual) 
          sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Activity (Bq)")


    # isotopes
    this_plot_dir = plot_dir + "/by_isotope_2"
    std.make_dir(this_plot_dir)

    # Activity n_tree
    g.print("Plots - residual activity by target")
    for target in results['a_tree'].keys():
      plot_list = []
      plot_file = this_plot_dir + "/" + isotopes.get_dir_hr(target) + ".eps"
      for residual in results['a_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          plot_array = numpy.zeros((steps_total,2,),) 
          plot_array[:,0] = t[:]
          plot_array[:,1] = results['a_tree'][target][residual][:]
          plot_list.append([isotopes.get_hr(residual), plot_array])
      title = "Activity - target: " + isotopes.get_hr(target) + "  residuals: all"
      sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Activity (Bq)")


    # isotopes
    this_plot_dir = plot_dir + "/by_isotope_3"
    std.make_dir(this_plot_dir)

    # Activity n_tree
    g.print("Plots - activity by residual tallied over all target isotopes")
    for residual in results['a_by_residual_isotope'].keys():
      plot_list = []
      plot_file = this_plot_dir + "/" + isotopes.get_dir_hr(residual) + ".eps"
      if(isotopes.is_unstable(residual) and results['a_by_residual_isotope'][residual] is not None):
        plot_array = numpy.zeros((steps_total,2,),) 
        plot_array[:,0] = t[:]
        plot_array[:,1] = results['a_by_residual_isotope'][residual][:]
        plot_list.append([isotopes.get_hr(residual), plot_array])
      title = "Activity - target: all  residual: " + isotopes.get_hr(residual) 
      sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Activity (Bq)")


    # isotopes
    this_plot_dir = plot_dir + "/by_isotope_4"
    std.make_dir(this_plot_dir)

    # Activity a_tree
    plot_list = []
    plot_file = this_plot_dir + "/all_radioactive_isotopes.eps"
    for target in results['a_tree'].keys():
      for residual in results['a_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          plot_array = numpy.zeros((steps_total,2,),) 
          plot_array[:,0] = t[:]
          plot_array[:,1] = results['a_tree'][target][residual][:]
          plot_list.append([isotopes.get_hr(target) + " -> " + isotopes.get_hr(residual), plot_array])
    title = "Activity - all residual isotopes" 
    sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Activity (Bq)")


    # Individual isotopes
    this_plot_dir = plot_dir + "/by_target_isotope"
    std.make_dir(this_plot_dir)

    # Activity n_tree
    g.print("Plots - activity by target isotope")
    plot_list = []
    plot_file = this_plot_dir + "/by_target_isotope.eps"
    for target in results['a_by_target_isotope'].keys():
      plot_array = numpy.zeros((steps_total,2,),) 
      plot_array[:,0] = t[:]
      plot_array[:,1] = results['a_by_target_isotope'][target][:]
      plot_list.append([isotopes.get_hr(target), plot_array])
    title = "Activity by Target Isotope"
    sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Activity (Bq)")



    # Individual isotopes
    this_plot_dir = plot_dir + "/particle_production"
    std.make_dir(this_plot_dir)

    # Activity n_tree
    g.print("Plots - particle production")
    plot_list = []
    plot_file = this_plot_dir + "/particle_production.eps"
    for particle in results['particle_production_by_particle'].keys():
      plot_array = numpy.zeros((steps_total,2,),) 
      plot_array[:,0] = t[:]
      plot_array[:,1] = results['particle_production_by_particle'][particle][:]
      plot_list.append([particle, plot_array])
    title = "Particle Production"
    sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Production Rate", 'log', 0.001)


    # Activity n_tree
    g.print("Plots - particle production in beam")
    plot_list = []
    plot_file = this_plot_dir + "/particle_production_in_beam.eps"
    for particle in results['particle_production_by_particle'].keys():
      plot_array = numpy.zeros((steps_eob,2,),) 
      plot_array[ta:tb,0] = t[ta:tb]
      plot_array[ta:tb,1] = results['particle_production_by_particle'][particle][ta:tb]
      plot_list.append([particle, plot_array])
    title = "Particle Production"
    sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Production Rate", 'log', 0.001)


    # Individual isotopes
    g.print("Plots - gamma lines over time")
    this_plot_dir = plot_dir + "/gammas/gamma_lines"
    this_plot_dir_eps = plot_dir + "/gammas/gamma_lines/eps"
    this_plot_dir_png = plot_dir + "/gammas/gamma_lines/png"
    std.make_dir(this_plot_dir)
    std.make_dir(this_plot_dir_eps)
    std.make_dir(this_plot_dir_png)
    for tn in range(1, len(t)):
      tn_str = str(tn)
      while(len(tn_str)<4):
        tn_str = "0" + tn_str
      tn_str = tn_str + "_" + str(int(t[tn]))
      file_path = this_plot_dir_eps + "/" + tn_str + ".eps"
      file_path_png = this_plot_dir_png + "/img" + str(tn - 1) + ".png"

      plot_data = [[],[]]
      for gk in results['g_gk_total_over_time'][tn].keys():
        plot_data[0].append(gk / 1000.0)    # convert to KeV
        plot_data[1].append(results['g_gk_total_over_time'][tn][gk]['activity'])

      plt.clf()
      plt.figure(figsize=(sim.plot_w, sim.plot_h))    
      plt.title('Gamma Lines at time ' + str(int(t[tn])) + 's')
      plt.xlabel('Energy (keV)')
      plt.ylabel('Activity (Bq)')
      if(tn < tb):
        plt.stem(plot_data[0][:], plot_data[1][:], 'r', markerfmt='ro')
      else:
        plt.stem(plot_data[0][:], plot_data[1][:], 'g', markerfmt='go')
      plt.savefig(file_path, format='eps', transparent=False)
      plt.savefig(file_path_png, format='png', dpi=sim.plot_dpi, transparent=False)
      plt.close('all') 



    
    ############################################
    # Pie Charts
    ############################################

    # Individual isotopes
    g.print("Plots - residual activity piecharts over time")
    this_plot_dir = plot_dir + "/activity_residual_piechart"
    this_plot_dir_eps = plot_dir + "/activity_residual_piechart/eps"
    this_plot_dir_png = plot_dir + "/activity_residual_piechart/png"
    std.make_dir(this_plot_dir)
    std.make_dir(this_plot_dir_eps)
    std.make_dir(this_plot_dir_png)
    for tn in range(1, len(t)):
      tn_str = str(tn)
      while(len(tn_str)<4):
        tn_str = "0" + tn_str
      tn_str = tn_str + "_" + str(int(t[tn]))
      file_path = this_plot_dir_eps + "/" + tn_str + ".eps"
      file_path_png = this_plot_dir_png + "/img" + str(tn - 1) + ".png"

      #colours = ["#009933","#CC9900","#0033CC","#CC3300","#6600CC","#000066","#333300","#FF0000","#FF9966","#009999","#993366","#66FF33","#FF0066","#00FFFF","#FFFF00","#990000","#663300"]
      labels = results['a_tree_over_time_residual_pie'][tn][0][:]
      values = results['a_tree_over_time_residual_pie'][tn][1][:]
      colours = []
      for label in labels:
        colours.append(sim.isotope_colour[label])
      sim.plot_pie(file_path, 'eps', labels, values, colours, 'Residual activity at time t=' + str(t[tn]) + '\nTotal activity: ' + str(results['total_activity'][tn]))
      sim.plot_pie(file_path_png, 'png', labels, values, colours, 'Residual activity at time t=' + str(t[tn]) + '\nTotal activity: ' + str(results['total_activity'][tn]),16,12)

      """
      plt.clf()
      plt.figure(figsize=(sim.plot_w,sim.plot_h))  
      plt.pie(values, labels=labels, autopct='%1.1f%%',
        shadow=True, startangle=90, colors=colours)
      plt.title('Residual activity at time t=' + str(t[tn]) + '\nTotal activity: ' + str(results['total_activity'][tn]))
      plt.legend(loc='upper right', fontsize='xx-small')
      plt.savefig(file_path, format='eps')
      plt.savefig(file_path_png, format='png', dpi=sim.plot_dpi)
      plt.close('all') 
      """

    # Individual isotopes
    g.print("Plots - target activity piecharts over time")
    this_plot_dir = plot_dir + "/activity_target_piechart"
    this_plot_dir_eps = plot_dir + "/activity_target_piechart/eps"
    this_plot_dir_png = plot_dir + "/activity_target_piechart/png"
    std.make_dir(this_plot_dir)
    std.make_dir(this_plot_dir_eps)
    std.make_dir(this_plot_dir_png)
    for tn in range(1, len(t)):
      tn_str = str(tn)
      while(len(tn_str)<4):
        tn_str = "0" + tn_str
      tn_str = tn_str + "_" + str(int(t[tn]))
      file_path = this_plot_dir_eps + "/" + tn_str + ".eps"
      file_path_png = this_plot_dir_png + "/img" + str(tn - 1) + ".png"

      #colours = ["#009933","#CC9900","#0033CC","#CC3300","#6600CC","#000066","#333300","#FF0000","#FF9966","#009999","#993366","#66FF33","#FF0066","#00FFFF","#FFFF00","#990000","#663300"]
      labels = results['a_tree_over_time_target_pie'][tn][0][:]
      values = results['a_tree_over_time_target_pie'][tn][1][:]
      colours = []
      for label in labels:
        colours.append(sim.isotope_colour[label])

      
      sim.plot_pie(file_path, 'eps', labels, values, colours, 'Target activity at time t=' + str(t[tn]) + '\nTotal activity: ' + str(results['total_activity'][tn]))
      sim.plot_pie(file_path_png, 'png', labels, values, colours, 'Target activity at time t=' + str(t[tn]) + '\nTotal activity: ' + str(results['total_activity'][tn]),16,12)

      """
      plt.clf()
      plt.figure(figsize=(sim.plot_w,sim.plot_h))  
      plt.pie(values, labels=labels, autopct='%1.1f%%',
        shadow=True, startangle=90, colors=colours)
      plt.title('Target activity at time t=' + str(t[tn]) + '\nTotal activity: ' + str(results['total_activity'][tn]))
      plt.savefig(file_path, format='eps', transparent=False)
      plt.savefig(file_path_png, format='png', dpi=sim.plot_dpi, transparent=False)
      plt.close('all') 
      """

    # Individual isotopes
    g.print("Plots - residual gamma energy piecharts over time")
    this_plot_dir = plot_dir + "/gamma_energy_by_residual_piechart"
    this_plot_dir_eps = plot_dir + "/gamma_energy_by_residual_piechart/eps"
    this_plot_dir_png = plot_dir + "/gamma_energy_by_residual_piechart/png"
    std.make_dir(this_plot_dir)
    std.make_dir(this_plot_dir_eps)
    std.make_dir(this_plot_dir_png)
    for tn in range(1, len(t)):
      tn_str = str(tn)
      while(len(tn_str)<4):
        tn_str = "0" + tn_str
      tn_str = tn_str + "_" + str(int(t[tn]))
      file_path = this_plot_dir_eps + "/" + tn_str + ".eps"
      file_path_png = this_plot_dir_png + "/img" + str(tn - 1) + ".png"

      #colours = ["#009933","#CC9900","#0033CC","#CC3300","#6600CC","#000066","#333300","#FF0000","#FF9966","#009999","#993366","#66FF33","#FF0066","#00FFFF","#FFFF00","#990000","#663300"]
      labels = results['g_residual_energy_pie'][tn][0][:]
      values = results['g_residual_energy_pie'][tn][1][:]
      colours = []
      for label in labels:
        colours.append(sim.isotope_colour[label])
      
      sim.plot_pie(file_path, 'eps', labels, values, colours, 'Residual gamma power at time t=' + str(t[tn]) + '\nTotal gamma power: ' + str(results['total_gamma_energy'][tn] * 1.60218e-16) + ' mW',8,6,1.60218e-16)
      sim.plot_pie(file_path_png, 'png', labels, values, colours, 'Residual gamma power at time t=' + str(t[tn]) + '\nTotal gamma power: ' + str(results['total_gamma_energy'][tn] * 1.60218e-16) + ' mW',16,12,1.60218e-16)

      """
      plt.clf()
      plt.figure(figsize=(sim.plot_w,sim.plot_h))  
      plt.pie(values, labels=labels, autopct='%1.1f%%',
        shadow=True, startangle=90, colors=colours)
      plt.title('Residual gamma activity at time t=' + str(t[tn]) + '\nTotal gamma power: ' + str(results['total_gamma_energy'][tn] * 1.60218e-16) + ' mW')
      plt.savefig(file_path, format='eps', transparent=False)
      plt.savefig(file_path_png, format='png', dpi=sim.plot_dpi, transparent=False)
      plt.close('all') 
      """

    # Individual isotopes
    g.print("Plots - target gamma energy piecharts over time")
    this_plot_dir = plot_dir + "/gamma_energy_by_target_piechart"
    this_plot_dir_eps = plot_dir + "/gamma_energy_by_target_piechart/eps"
    this_plot_dir_png = plot_dir + "/gamma_energy_by_target_piechart/png"
    std.make_dir(this_plot_dir)
    std.make_dir(this_plot_dir_eps)
    std.make_dir(this_plot_dir_png)
    for tn in range(1, len(t)):
      tn_str = str(tn)
      while(len(tn_str)<4):
        tn_str = "0" + tn_str
      tn_str = tn_str + "_" + str(int(t[tn]))
      file_path = this_plot_dir_eps + "/" + tn_str + ".eps"
      file_path_png = this_plot_dir_png + "/img" + str(tn - 1) + ".png"

      #colours = ["#009933","#CC9900","#0033CC","#CC3300","#6600CC","#000066","#333300","#FF0000","#FF9966","#009999","#993366","#66FF33","#FF0066","#00FFFF","#FFFF00","#990000","#663300"]
      labels = results['g_target_energy_pie'][tn][0][:]
      values = results['g_target_energy_pie'][tn][1][:]
      colours = []
      for label in labels:
        colours.append(sim.isotope_colour[label])

      sim.plot_pie(file_path, 'eps', labels, values, colours, 'Target gamma power at time t=' + str(t[tn]) + '\nTotal gamma power: ' + str(results['total_gamma_energy'][tn] * 1.60218e-16) + ' mW',8,6,1.60218e-16)
      sim.plot_pie(file_path_png, 'png', labels, values, colours, 'Target gamma power at time t=' + str(t[tn]) + '\nTotal gamma power: ' + str(results['total_gamma_energy'][tn] * 1.60218e-16) + ' mW',16,12,1.60218e-16)

      """
      plt.clf()
      plt.figure(figsize=(sim.plot_w,sim.plot_h))  
      plt.pie(values, labels=labels, autopct='%1.1f%%',
        shadow=True, startangle=90, colors=colours)
      plt.title('Target gamma activity at time t=' + str(t[tn]) + '\nTotal gamma power: ' + str(results['total_gamma_energy'][tn] * 1.60218e-16) + ' mW')
      plt.savefig(file_path, format='eps', transparent=False)
      plt.savefig(file_path_png, format='png', dpi=sim.plot_dpi, transparent=False)
      plt.close('all') 
      """

    g.print("Attempt to make videos.")

    # Make videos
    try: 
      this_plot_dir_png = plot_dir + "/activity_residual_piechart/png"
      os.system("ffmpeg -r 15 -i " + this_plot_dir_png + "/img%01d.png -codec:v libx264 -crf 18 -y " + video_dir + "/activity_residual_piechart.mp4 > " + logs_dir + "/video1.txt 2>&1")
    except:
      pass

    try: 
      this_plot_dir_png = plot_dir + "/activity_target_piechart/png"
      os.system("ffmpeg -r 15 -i " + this_plot_dir_png + "/img%01d.png -codec:v libx264 -crf 18 -y " + video_dir + "/activity_target_piechart.mp4 > " + logs_dir + "/video2.txt 2>&1")
    except:
      pass

    try: 
      this_plot_dir_png = plot_dir + "/gamma_energy_by_residual_piechart/png"
      os.system("ffmpeg -r 15 -i " + this_plot_dir_png + "/img%01d.png -codec:v libx264 -crf 18 -y " + video_dir + "/gamma_energy_by_residual_piechart.mp4 > " + logs_dir + "/video3.txt 2>&1")
    except:
      pass

    try: 
      this_plot_dir_png = plot_dir + "/gamma_energy_by_target_piechart/png"
      os.system("ffmpeg -r 15 -i " + this_plot_dir_png + "/img%01d.png -codec:v libx264 -crf 18 -y " + video_dir + "/gamma_energy_by_target_piechart.mp4 > " + logs_dir + "/video4.txt 2>&1")
    except:
      pass

    try: 
      this_plot_dir_png = plot_dir + "/gammas/gamma_lines/png"
      os.system("ffmpeg -r 15 -i " + this_plot_dir_png + "/img%01d.png -codec:v libx264 -crf 18 -y " + video_dir + "/gamma_lines.mp4 > " + logs_dir + "/video5.txt 2>&1")
    except:
      pass





    fhr.close()
    sim.gamma_log.close()
    sim.activity_log.close()


    fh = open(dir + '/readme.txt', 'w')  
    fh.write("""
##############################################################
                            Read Me                           
##############################################################

Each sim has its own directory.


data
########################


a_tree

A directory for each isotope in the target material.  Contains a data file for the activity of each residual isotope created from that target isotope.  Time in s, activity in Bq.


a_by_radioactive_isotope

A data file for each radioactive isotope created from all target isotopes.  Time in s, activity in Bq.


a_by_target_isotope

A data file for each target isotope with a tally of radioactivity resulting from each target isotope.  Time in s, activity in Bq.


gammas eob_gammas.txt

The energy (keV) and activity (Bq) of gammas from all radioactive isotopes at the time the beam stops.


gammas eos_gammas.txt

The energy (keV) and activity (Bq) of gammas from all radioactive isotopes at the time the simulation stops.




plots
########################


by_isotope_1

A directory for each target isotope, each containing a time-activity plot for the residual radioactive isotopes created.


by_isotope_2

A time-activity plot for each target isotope showing the activity over time of all radioactive residual isotopes created from that target isotope.


by_isotope_3

A time-activity plot for each radioactive residual isotope collated across all target isotopes.


by_isotope_4

A time-activity plot for each radioactive residual isotope.


by_target_isotope

A time-activity plot for each target isotope collated over all radioactive residual isotopes for each target isotope.



gammas - all

A time-activity plot for each gamma across all radioactive residual isotopes created by all target isotopes.




videos
########################

activity_residual_piechart.mp4
The activity grouped by residual isotopes tallied over all target isotopes as it changes over time.


activity_target_piechart.mp4
The activity grouped by target isotopes as it changes over time.


gamma_energy_by_residual_piechart.mp4
The gamma energy grouped by residual isotopes tallied over all target isotopes as it changes over time.


gamma_energy_by_target_piechart.mp4
The gamma energy grouped by target isotopes as it changes over time.


gamma_lines.mp4
Predicted gamma lines from the target over time (due to the decay of residual isotopes).



""")
    fh.write('\n')
    fh.close()






  # rr = f * nd * dt * 1.0E-28  
  @staticmethod
  def run_t(icode, projectile, nd, beam_flux, target_depth, exyz_data, residual_rrs, particle_rrs):
    residual_rr = {}
    particle_rr = {}
    dnisotope = 0.0
    ion_count = len(exyz_data)

    for ion in exyz_data.keys():
      for i in range(len(exyz_data[ion])): 
        e = exyz_data[ion][i,0]             # Energy in MeV
        dt = exyz_data[ion][i,1]
        start_depth = exyz_data[ion][i,4]
        end_depth = exyz_data[ion][i,5]
        if(start_depth > target_depth):
          break       
        else:
          m = 1.0
          if(end_depth > target_depth):
            m = (target_depth - start_depth) / (end_depth - start_depth)
          dt = m * dt  
 
          # Residuals
          rxs = talys.search_rxs(e, projectile, icode)
          for rcode in rxs.keys():
            if(rcode not in residual_rr.keys()):
              residual_rr[rcode] = 0.0
            dn = (1 / ion_count) * rxs[rcode] * beam_flux * nd * dt * 1.0E-28
            residual_rr[rcode] = residual_rr[rcode] + dn
            dnisotope = dnisotope - dn

          # Particles
          pxs = talys.search_pxs(e, projectile, icode)
          for pcode in pxs.keys():
            if(pcode not in particle_rr.keys()):
              particle_rr[pcode] = 0.0
            dn = (1 / ion_count) * pxs[pcode] * beam_flux * nd * dt * 1.0E-28
            particle_rr[pcode] = particle_rr[pcode] + dn
            dnisotope = dnisotope - dn

    residual_rrs[icode] = residual_rr
    particle_rrs[icode] = particle_rr


  # rr = f * nd * dt * 1.0E-28  
  @staticmethod
  def run_tn(icode, projectile, nd, beam_flux, target_depth, neutron_energy, residual_rrs, particle_rrs):
    residual_rr = {}
    particle_rr = {}
    dnisotope = 0.0
    neutron_count = len(neutron_energy)
    for e in neutron_energy:
      dt = target_depth

      # Residuals
      rxs = talys.search_rxs(e, projectile, icode)
      for rcode in rxs.keys():
        if(rcode not in residual_rr.keys()):
          residual_rr[rcode] = 0.0
        dn = (1 / neutron_count) * rxs[rcode] * beam_flux * nd * dt * 1.0E-28
        residual_rr[rcode] = residual_rr[rcode] + dn
        dnisotope = dnisotope - dn

      # Particles
      pxs = talys.search_pxs(e, projectile, icode)
      for pcode in pxs.keys():
        if(pcode not in particle_rr.keys()):
          particle_rr[pcode] = 0.0
        dn = (1 / neutron_count) * pxs[pcode] * beam_flux * nd * dt * 1.0E-28
        particle_rr[pcode] = particle_rr[pcode] + dn
        dnisotope = dnisotope - dn

    residual_rrs[icode] = residual_rr
    particle_rrs[icode] = particle_rr

  

  @staticmethod
  def plot(file_path, plot_data, title="", x_axis=None, y_axis=None):
    plt.clf()
    plt.figure(figsize=(sim.plot_w, sim.plot_h))    
    plt.title(title)
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.stem(plot_data[:,0], plot_data[:,1])
    plt.savefig(file_path, format='eps', transparent=False)
    plt.close('all') 


  @staticmethod
  def plot_activity(file_path, plot_list, title="", x_axis=None, y_axis=None, scale="linear", ymin=None, ymax=None):

    if(x_axis == None):
      x_axis = 'Time (s)'
    if(y_axis == None):
      y_axis = 'Activity (Bq)'

    linestyle = ["solid", "dashed", "dashdot", "dotted"]
    colours = ["#009933","#CC9900","#0033CC","#CC3300","#6600CC","#000066","#333300","#FF0000","#FF9966","#009999","#993366","#66FF33","#FF0066","#00FFFF","#FFFF00","#990000","#663300"]

    
    plt.figure(figsize=(sim.plot_w, sim.plot_h))
    n = 0
    for p in plot_list:
      label = p[0]
      points = p[1]
      lc = n % len(colours)
      ls = int(numpy.floor(n / len(colours))) % len(linestyle)

      plt.plot(points[:,0], points[:,1], label=label, color=colours[lc], linestyle=linestyle[ls])
      n = n + 1
    cols = (max(1,int(numpy.floor(len(plot_list)/10))))
    if(scale == 'log'):
      plt.yscale('log')
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    if(ymin != None and ymax == None):
      plt.ylim(bottom=ymin)
    elif(ymin == None and ymax != None):
      plt.ylim(top=ymax)
    elif(ymin != None and ymax != None):
      plt.ylim(bottom=ymin, top=ymax)
    plt.title(title)
    plt.legend(fontsize="xx-small", ncol=cols)
    plt.grid(True)
    plt.savefig(file_path, format='eps', transparent=False)
    plt.close('all') 


  @staticmethod
  def pfunc(pct, allvals, show_value):
    if(show_value == 0.0):
      if(pct < 5.0):
        return ""
      elif(pct >= 5.0):
        return "{:.1f}%".format(pct)
    else:
      if(pct < 10.0):
        return ""
      else:
        absolute = int(pct/100.*numpy.sum(allvals))
        return "{:.3e}".format(show_value * absolute)


  @staticmethod
  def plot_pie(file_path, file_format, labels, values, colours, title="", w=None, h=None, show_value=1.0):
    if(w == None):
      w = sim.plot_w
    if(h == None):
      h = sim.plot_h

    sim.gamma_log.write(file_path + '\n')

    s = numpy.sum(values)
    p = (values[:] / s) * 100

    l_labels = []
    for ln in range(len(labels)):
      l_labels.append(labels[ln] + " {:.1f}%".format(100 * (values[ln] / s)))
      sim.gamma_log.write(str(labels[ln]) + " " + str(p[ln]) + " " + str(values[ln]) + "\n")
    sim.gamma_log.write('\n')

    fig, ax = plt.subplots(figsize=(w,h), subplot_kw=dict(aspect="equal"))
    wedges, texts, autotexts = ax.pie(p, 
                                      autopct='%1.1f%%',
                                      textprops=dict(color="w"), colors=colours)
  
    ax.legend(wedges, l_labels,
          title='',
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1),
          fontsize=7)
    plt.setp(autotexts, size=7, weight="bold")
    ax.set_title(title)
    plt.savefig(file_path, format=file_format, transparent=False)
    plt.close('all') 

    """
    fig, ax = plt.subplots(figsize=(w,h))
    s = numpy.sum(values)
    l_labels = []
    for ln in range(len(labels)):
      l_labels.append(labels[ln] + " {:.1f}%".format(100 * (values[ln] / s)))
    
    wedges, texts, autotexts = ax.pie(values, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90, textprops=dict(color="w"), colors=colours)
    ax.axis('equal') 
    ax.legend(wedges, l_labels,
      title='',
      loc="center left",
      bbox_to_anchor=(1, 0, 0.5, 1),
      fontsize=7)
    plt.savefig(file_path, format=file_format, transparent=False)
    plt.close('all') 
    """

  def make_colours():
    sim.colours = []
    for i in range(1,2048):
      r = (i * 769927) % 256
      g = (i * 1209139) % 256
      b = (i * 490951) % 256
      h = hex(256*256*r + 256*g + b)
      if(len(h) == 7):
        h = "0" + h[2:9]
      else:
        h = h[2:10]
      sim.colours.append("#"+h)

