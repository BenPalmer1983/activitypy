import numpy
import matplotlib.pyplot as plt
import os


from g import g
from isotopes import isotopes
from talys import talys
from exyz import exyz
from decay import decay
from std import std


class sim:

  isotope_file = None
  talys_dir = None

  
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
    dir = "results/sim_" + str(counter)
    data_dir = dir + "/data"
    plot_dir = dir + "/plots"
    video_dir = dir + "/videos"
    std.make_dir(dir)
    std.make_dir(data_dir)
    std.make_dir(plot_dir)
    std.make_dir(video_dir)


    beam_flux = this_sim['beam_flux']
    target_depth = this_sim['target_depth']
    projectile = this_sim['beam_projectile']
    beam_duration = this_sim['beam_duration']
    end_time = this_sim['end_time']


    g.print()
    g.print("Sim " + str(counter))
    g.print("###################################")
    g.print()
    if(this_sim['run'] == False):
      g.print("Run == False")
      g.print()
      return None

    g.print("Starting Target Material:")
    g.print()
    g.print(this_sim['target_composition_hr'])
    g.print()
    fh = open(data_dir + "/target_material.txt", "w")
    fh.write(this_sim['target_composition_hr'])
    fh.close()

    ##############################################
    # Load trajectory data
    ##############################################

    g.print("Loading exyz " + this_sim['exyz'])
    exyz_data = exyz.load(this_sim['exyz'])
    g.print("Load complete")
    g.print()



    ##############################################
    # Calculate reaction rates for residuals
    ##############################################

    g.print("Residual reaction rates")
    residual_rrs = {}    
    for icode in this_sim['target_composition'].keys():
      nd = this_sim['target_composition'][icode]['nd']
      sim.run_t(icode, projectile, nd, beam_flux, target_depth, exyz_data, residual_rrs)
 
    g.print("Save reaction rates")
    fh = open(data_dir + "/rr_residuals.txt", "w")
    fh.write("Residual Reaction Rates\n\n")
    for target in residual_rrs.keys():
      fh.write("Source isotope: " + isotopes.get_hr(target) + "  (" + str(target) + ")" + "\n\n")
      for residual in residual_rrs[target].keys():
        fh.write(std.pad_right(isotopes.get_hr(residual) + ":", 20) + str(residual_rrs[target][residual]) + "\n")
      fh.write("\n\n")        
    fh.close()


    g.print("Calculate residual amounts over time")

    steps_eob = 151    # From start of beam to end of beam
    steps_eos = 181    # From end of beam to end of sim
    steps_total = steps_eob + steps_eos - 1
    ta = 0
    tb = steps_eob
    tc = steps_eob - 1
    td = steps_total

    # Time array
    t = numpy.zeros((steps_total,),)
    t[ta:tb] = numpy.linspace(0.0, beam_duration, steps_eob)
    t[tc:td] = numpy.linspace(beam_duration, end_time, steps_eos)

    """
    Isotope amounts:
    n_tree                    target -> residual -> amount array at each time step
    n_by_residual_isotope     
    n_by_target_isotope       
    """


    results = {
                   'n_tree': {}, 
                   'n_by_residual_isotope': {}, 
                   'n_by_target_isotope': {}, 
                   'a_tree': {}, 
                   'a_by_radioactive_isotope': {}, 
                   'a_by_target_isotope': {}, 
                   'a_total': None, 
                   'a_tree_over_time': [], 
                   'a_tree_over_time_residual': [], 
                   'a_tree_over_time_target': [], 
                   'a_tree_over_time_residual_pie': [], 
                   'a_tree_over_time_target_pie': [], 
                   'a_tree_over_time_total': [],
                   'g_tree': {},
                   'g_by_gk_residual_isotope': {}, 
                   'g_by_gk_target_isotope': {},   
                   'g_by_gk_all': {},   
                   'g_by_residual_isotope': {}, 
                   'g_by_target_isotope': {},   
                   'g_by_all': {},
                   'g_tree_over_time': [],  
                   'g_tree_over_time_residual': [],
                   'g_tree_over_time_target': [],
                   'g_tree_over_time_total': [],
                   'g_total': [],      
                   'g_tree_eob': {}, 
                   'g_eob': None, 
                   'g_eos': None, 
                   'g_activity': None, 
                   'g_energy': None, 
                   'g_energy_by_target_isotope': {},
                   'g_energy_by_residual_isotope': {},
                  }    
    eob_idata = {}     # End of beam data

    # Create dictionary entries
    for target in residual_rrs.keys():
      if(target not in results['n_tree'].keys()):
        results['n_tree'][target] = {}  
      if(target not in eob_idata.keys()):
        eob_idata[target] = {}  

    # In beam activity
    g.print("In beam activity over 101 time steps")
    for target in residual_rrs.keys():
      for residual in residual_rrs[target].keys():
        if(not isotopes.is_stable(residual)):
          idata = {}
          idata[residual] = {'w': residual_rrs[target][residual], 'n0': 0.0}
          parent = residual  
          # In beam calculations
          for i in range(steps_eob):
            time = t[i]
            r = decay.calculate(parent, time, idata)
            for k in r['tally'].keys():
              if(k not in results['n_tree'][target].keys()):
                results['n_tree'][target][k] = numpy.zeros((steps_total,),)  
              results['n_tree'][target][k][i] = results['n_tree'][target][k][i] + r['tally'][k]['nend']   
            if(i==(steps_eob-1)):
              eob_idata[target][residual] = {}
              for k in r['tally'].keys():
                eob_idata[target][residual][k] = {}
                eob_idata[target][residual][k]['w'] = 0.0
                eob_idata[target][residual][k]['n0'] = results['n_tree'][target][k][i]

    # End of beam activity      
    g.print("End of beam activity over 101 time steps")
    for target in residual_rrs.keys():
      for residual in residual_rrs[target].keys():
        if(not isotopes.is_stable(residual)):
          idata = eob_idata[target][residual]
          #print(idata)
          parent = residual  
          # Out of beam calculations
          for i in range(steps_eos-1):        
            time = t[i+steps_eob] - beam_duration
            r = decay.calculate(parent, time, idata)
            for k in r['tally'].keys():
              if(k not in results['n_tree'][target].keys()):
                results['n_tree'][target][k] = numpy.zeros((steps_total,),)  
              results['n_tree'][target][k][i+steps_eob] = results['n_tree'][target][k][i+steps_eob] + r['tally'][k]['nend']    


    # By radioactive isotope
    g.print("Tally amount by residual isotope")
    for target in results['n_tree'].keys():
      for residual in results['n_tree'][target].keys():
        if(residual not in results['n_by_residual_isotope'].keys()):
          results['n_by_residual_isotope'][residual] = numpy.zeros((steps_total,),)  
        results['n_by_residual_isotope'][residual][:] = results['n_by_residual_isotope'][residual][:] + results['n_tree'][target][residual][:]


    # By target isotope
    g.print("Tally amount of residual isotopes by target isotope")
    for target in results['n_tree'].keys():
      if(target not in results['n_by_target_isotope'].keys()):
        results['n_by_target_isotope'][target] = numpy.zeros((steps_total,),)  
      for residual in results['n_tree'][target].keys():
        results['n_by_target_isotope'][target][:] = results['n_by_target_isotope'][target][:] + results['n_tree'][target][residual][:]




    ###################
    # Activities
    ###################

    # Activity n_tree
    g.print("Activity per target isotope by residual radioactive isotope over time")
    for target in results['n_tree'].keys():
      if(target not in results['a_tree'].keys()):
        results['a_tree'][target] = {}
      for residual in results['n_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          l = isotopes.get_decay_constant(residual)
          results['a_tree'][target][residual] = l * results['n_tree'][target][residual]

    # By radioactive isotope
    g.print("Activity by residual radioactive isotope over time")
    for residual in results['n_by_residual_isotope'].keys():
      if(target not in results['a_by_radioactive_isotope'].keys()):
        results['a_by_radioactive_isotope'][target] = numpy.zeros((steps_total,),)  
      if(isotopes.is_unstable(residual)):
        l = isotopes.get_decay_constant(residual)
        results['a_by_radioactive_isotope'][residual] = l * results['n_by_residual_isotope'][residual]
      
    # By target isotope
    g.print("Activity resulting from each target isotope over time")
    for target in results['n_tree'].keys():
      if(target not in results['a_by_target_isotope'].keys()):
        results['a_by_target_isotope'][target] = numpy.zeros((steps_total,),)  
      for residual in results['n_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          l = isotopes.get_decay_constant(residual)
          results['a_by_target_isotope'][target][:] = results['a_by_target_isotope'][target][:] + l * results['n_tree'][target][residual][:]
 
    # By total activity
    g.print("Total activity from all targets and residuals over time")
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
    g.print("Activity per target and residual over time")
    for tn in range(len(t)):
      for target in results['a_tree'].keys():
        results['a_tree_over_time'][tn][target] = {}
        for residual in results['a_tree'][target].keys():
          if(isotopes.is_unstable(residual)):
            results['a_tree_over_time'][tn][target][residual] = results['a_tree'][target][residual][tn]


    for tn in range(len(t)):
      for target in results['a_tree'].keys():
        for residual in results['a_tree'][target].keys():
          if(isotopes.is_unstable(residual)):
            if(residual not in results['a_tree_over_time_residual'][tn].keys()):
              results['a_tree_over_time_residual'][tn][residual] = 0.0
            results['a_tree_over_time_residual'][tn][residual] = results['a_tree_over_time_residual'][tn][residual] + results['a_tree'][target][residual][tn]
            ar_max[tn] = max(ar_max[tn], results['a_tree_over_time_residual'][tn][residual])
            ar_total[tn] = ar_total[tn] + results['a_tree'][target][residual][tn]

    for tn in range(len(t)):
      for target in results['a_tree'].keys():
        if(target not in results['a_tree_over_time_target'][tn].keys()):
          results['a_tree_over_time_target'][tn][target] = 0.0
        for residual in results['a_tree'][target].keys():
          if(isotopes.is_unstable(residual)):
            results['a_tree_over_time_target'][tn][target] = results['a_tree_over_time_target'][tn][target] + results['a_tree'][target][residual][tn]
            at_max[tn] = max(at_max[tn], results['a_tree_over_time_target'][tn][target])
            at_total[tn] = at_total[tn] + results['a_tree'][target][residual][tn]
      
    for tn in range(len(t)):
      for residual in results['a_tree_over_time_residual'][tn].keys():
        if(results['a_tree_over_time_residual'][tn][residual] > 0.05 * ar_max[tn]):
          results['a_tree_over_time_residual_pie'][tn][0].append(isotopes.get_hr(residual))
          results['a_tree_over_time_residual_pie'][tn][1].append(results['a_tree_over_time_residual'][tn][residual])
    
    for tn in range(len(t)):
      for target in results['a_tree_over_time_target'][tn].keys():
        if(results['a_tree_over_time_target'][tn][target] > 0.05 * at_max[tn]):
          results['a_tree_over_time_target_pie'][tn][0].append(isotopes.get_hr(target))
          results['a_tree_over_time_target_pie'][tn][1].append(results['a_tree_over_time_target'][tn][target])
    





    ###################
    # Gammas
    ###################

    # Tree

    # Target Residual Gamma tree 
    g.print("Tally gammas by target -> residual -> energy")
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
              g_energy_ev = ga[n,0] * 1000.0
              g_intensity = ga[n,1]
              if(gk not in results['g_tree'][target][residual].keys()):
                results['g_tree'][target][residual][gk] = {'activity': None, 'energy': None,}
                results['g_tree'][target][residual][gk]['activity'] = results['a_tree'][target][residual][:] * g_intensity
                results['g_tree'][target][residual][gk]['energy'] = g_energy_ev * results['g_tree'][target][residual][gk]['activity'][:]



    # Individual gamma energies

    # By residual isotope
    g.print("Tally gammas by residual -> energy")
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
    g.print("Tally gammas by target -> energy")
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
    g.print("Tally gammas -> by gamma energy")
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
    g.print("Tally gammas by residual and gamma energies")
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
    g.print("Tally gammas by target and gamma energies")
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
    g.print("Tally gammas by target, residual and gamma energies")
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
    g.print("Per target and residual gamma lines over time")
    for tn in range(len(t)):
      results['g_tree_over_time'].append({})
    for tn in range(len(t)):
      for target in results['g_tree'].keys():
        results['g_tree_over_time'][tn][target] = {}
        for residual in results['g_tree'][target].keys():
          results['g_tree_over_time'][tn][target][residual] = {}
          if(isotopes.is_unstable(residual)):
            for gk in results['g_tree'][target][residual].keys():
              results['g_tree_over_time'][tn][target][residual][gk] = results['g_tree'][target][residual][gk]['activity'][tn]
    

    # By target over time
    g.print("Per target gamma lines over time")
    for tn in range(len(t)):
      results['g_tree_over_time_target'].append({})
    for tn in range(len(t)):
      for target in results['g_tree_over_time'][tn].keys():
        if(target not in results['g_tree_over_time_target'][tn].keys()):
          results['g_tree_over_time_target'][tn][target] = {}
        for residual in results['g_tree_over_time'][tn][target].keys():
          for gk in results['g_tree_over_time'][tn][target][residual].keys():
            if(gk not in results['g_tree_over_time_target'][tn][target].keys()):
              results['g_tree_over_time_target'][tn][target][gk] = 0.0
            results['g_tree_over_time_target'][tn][target][gk] = results['g_tree_over_time_target'][tn][target][gk] + results['g_tree'][target][residual][gk]['activity'][tn]


    # By residual over time
    g.print("Per residual gamma lines over time")
    for tn in range(len(t)):
      results['g_tree_over_time_residual'].append({})
    for tn in range(len(t)):
      for target in results['g_tree_over_time'][tn].keys():
        for residual in results['g_tree_over_time'][tn][target].keys():
          if(residual not in results['g_tree_over_time_residual'][tn].keys()):
            results['g_tree_over_time_residual'][tn][residual] = {}
          for gk in results['g_tree_over_time'][tn][target][residual].keys():
            if(gk not in results['g_tree_over_time_residual'][tn][residual].keys()):
              results['g_tree_over_time_residual'][tn][residual][gk] = 0.0
            results['g_tree_over_time_residual'][tn][residual][gk] = results['g_tree_over_time_residual'][tn][residual][gk] + results['g_tree'][target][residual][gk]['activity'][tn]


    # Total over time
    g.print("Tallied gamma lines over time")
    for tn in range(len(results['g_tree_over_time'])):
      results['g_tree_over_time_total'].append({})
      for target in results['g_tree_over_time'][tn].keys():
        for residual in results['g_tree_over_time'][tn][target].keys():
          for gk in results['g_tree_over_time'][tn][target][residual].keys():
            if(gk not in results['g_tree_over_time_total'][tn].keys()):
              results['g_tree_over_time_total'][tn][gk] = 0.0
            results['g_tree_over_time_total'][tn][gk] = results['g_tree_over_time_total'][tn][gk] + results['g_tree'][target][residual][gk]['activity'][tn]






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

    dir_n_tree = data_dir + "/a_by_radioactive_isotope"
    std.make_dir(dir_n_tree)
    for residual in results['a_by_radioactive_isotope'].keys():
      out_array = numpy.zeros((steps_total,2,),) 
      out_array[:,0] = t[:]
      out_array[:,1] = results['a_by_radioactive_isotope'][residual][:]
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
      tn_str = tn_str + "_" + str(int(t[tn])) + ".txt"

      g_out = numpy.zeros((len(results['g_tree_over_time_total'][tn]), 2),)
      gn = 0
      for gk in results['g_tree_over_time_total'][tn].keys():
        g_out[gn, 0] = gk / 1000.0  # Convert to KeV
        g_out[gn, 1] = results['g_tree_over_time_total'][tn][gk]
        gn = gn + 1

      g_out = g_out[g_out[:, 0].argsort(), :]
      file_name = gammas_data_dir + "/" + tn_str
      numpy.savetxt(file_name, g_out, delimiter=",")

      





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
    this_plot_dir = plot_dir + "/total_activity"
    std.make_dir(this_plot_dir)

    # Activity n_tree
    g.print("Plots - activity by target and residual isotope")
    plot_file = this_plot_dir + "/total_activity.eps"
    plot_list = []
    plot_array = numpy.zeros((steps_total,2,),) 
    plot_array[:,0] = t[:]
    plot_array[:,1] = results['a_total'][:]
    plot_list.append([isotopes.get_hr(residual) + " from " + isotopes.get_hr(target), plot_array])
    title = "Total activity over time"
    sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Activity (Bq)")
   

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
    for residual in results['a_by_radioactive_isotope'].keys():
      plot_list = []
      plot_file = this_plot_dir + "/" + isotopes.get_dir_hr(residual) + ".eps"
      if(isotopes.is_unstable(residual) and results['a_by_radioactive_isotope'][residual] is not None):
        plot_array = numpy.zeros((steps_total,2,),) 
        plot_array[:,0] = t[:]
        plot_array[:,1] = results['a_by_radioactive_isotope'][residual][:]
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
    this_plot_dir = plot_dir + "/residual_activity_piechart"
    this_plot_dir_eps = plot_dir + "/residual_activity_piechart/eps"
    this_plot_dir_png = plot_dir + "/residual_activity_piechart/png"
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

      colours = ["#009933","#CC9900","#0033CC","#CC3300","#6600CC","#000066","#333300","#FF0000","#FF9966","#009999","#993366","#66FF33","#FF0066","#00FFFF","#FFFF00","#990000","#663300"]

      plt.clf()
      plt.figure(figsize=(16,12))  
      plt.pie(results['a_tree_over_time_residual_pie'][tn][1][:], labels=results['a_tree_over_time_residual_pie'][tn][0][:], autopct='%1.1f%%',
        shadow=True, startangle=90, colors=colours)
      plt.title('Residual activity at time t=' + str(t[tn]), bbox={'facecolor':'0.8', 'pad':5})
      plt.savefig(file_path, type='eps')
      plt.savefig(file_path_png, type='png', dpi=200)
      plt.close('all') 



    # Individual isotopes
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
      for gk in results['g_tree_over_time_total'][tn].keys():
        plot_data[0].append(gk)
        plot_data[1].append(results['g_tree_over_time_total'][tn][gk] / 1000.0)   # convert to KeV

      plt.clf()
      plt.figure(figsize=(16,12))    
      plt.title('Gamma Lines at time ' + str(int(t[tn])) + 's')
      plt.xlabel('Energy (keV)')
      plt.ylabel('Activity (Bq)')
      plt.stem(plot_data[0][:], plot_data[1][:])
      plt.savefig(file_path, type='eps')
      plt.savefig(file_path_png, type='png', dpi=200)
      plt.close('all') 





    # Make videos
    try: 
      this_plot_dir_png = plot_dir + "/residual_activity_piechart/png"
      os.system("ffmpeg -r 15 -i " + this_plot_dir_png + "/img%01d.png -codec:v libx264 -crf 18 -y " + video_dir + "/residual_isotopes.mp4")
    except:
      pass

    try: 
      this_plot_dir_png = plot_dir + "/gammas/gamma_lines/png"
      os.system("ffmpeg -r 15 -i " + this_plot_dir_png + "/img%01d.png -codec:v libx264 -crf 18 -y " + video_dir + "/gamma_lines.mp4")
    except:
      pass







    exit()

    # end of beam
    g_eob = {}
    for target in results['a_tree'].keys():
      for residual in results['a_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          ga = isotopes.get_gammas_array(residual)
          if(not isinstance(ga, type(None))):
            for n in range(len(ga)):
              if(ga[n,0] not in g_eob.keys()):
                g_eob[ga[n,0]] = 0.0
              g_eob[ga[n,0]] = g_eob[ga[n,0]] + results['a_tree'][target][residual][100] * ga[n,1]
          
    results['g_eob'] = numpy.zeros((len(g_eob), 2),)
    n = 0
    for k in g_eob.keys():
      results['g_eob'][n,0] = float(k) / 1000.0   # Convert from eV to keV
      results['g_eob'][n,1] = float(g_eob[k])
      n = n + 1
    results['g_eob'] = results['g_eob'][results['g_eob'][:, 0].argsort(), :]


    # end of sim
    g_eos = {}
    for target in results['a_tree'].keys():
      for residual in results['a_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          ga = isotopes.get_gammas_array(residual)
          if(not isinstance(ga, type(None))):
            for n in range(len(ga)):
              if(ga[n,0] not in g_eos.keys()):
                g_eos[ga[n,0]] = 0.0
              g_eos[ga[n,0]] = g_eos[ga[n,0]] + results['a_tree'][target][residual][-1] * ga[n,1]

    results['g_eos'] = numpy.zeros((len(g_eos), 2),)
    n = 0
    for k in g_eos.keys():
      results['g_eos'][n,0] = float(k) / 1000.0   # Convert from eV to keV
      results['g_eos'][n,1] = float(g_eos[k])
      n = n + 1
    results['g_eos'] = results['g_eos'][results['g_eos'][:, 0].argsort(), :]

    
    # Gammas over time    
    for target in results['a_tree'].keys():
      for residual in results['a_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          ga = isotopes.get_gammas_array(residual)
          if(not isinstance(ga, type(None))):
            for n in range(len(ga)):
              k = ga[n,0] / 1000.0
              if(k not in results['g_all'].keys()):
                results['g_all'][k] = numpy.zeros((201,),) 
              results['g_all'][k][:] = results['g_all'][k][:] + results['a_tree'][target][residual][:] * ga[n,1]    


  



    ###################
    # Plots - Activities only
    ###################
    
    g.print("Make plots")

    # Individual isotopes
    this_plot_dir = plot_dir + "/total_activity"
    std.make_dir(this_plot_dir)

    # Activity n_tree
    for target in results['a_tree'].keys():
      dir_a = this_plot_dir + "/"  + isotopes.get_dir_hr(target)
      std.make_dir(dir_a)
      for residual in results['a_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          plot_file = dir_a + "/" + isotopes.get_dir_hr(residual) + ".eps"
          plot_list = []
          plot_array = numpy.zeros((201,2,),) 
          plot_array[:,0] = t[:]
          plot_array[:,1] = results['a_tree'][target][residual][:]
          plot_list.append([isotopes.get_hr(residual) + " from " + isotopes.get_hr(target), plot_array])
          title = "Activity - target: " + isotopes.get_hr(target) + "  residual: " + isotopes.get_hr(residual) 
          sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Activity (Bq)")

   



    # Individual isotopes
    this_plot_dir = plot_dir + "/by_isotope_1"
    std.make_dir(this_plot_dir)

    # Activity n_tree
    for target in results['a_tree'].keys():
      dir_a = this_plot_dir + "/"  + isotopes.get_dir_hr(target)
      std.make_dir(dir_a)
      for residual in results['a_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          plot_file = dir_a + "/" + isotopes.get_dir_hr(residual) + ".eps"
          plot_list = []
          plot_array = numpy.zeros((201,2,),) 
          plot_array[:,0] = t[:]
          plot_array[:,1] = results['a_tree'][target][residual][:]
          plot_list.append([isotopes.get_hr(residual) + " from " + isotopes.get_hr(target), plot_array])
          title = "Activity - target: " + isotopes.get_hr(target) + "  residual: " + isotopes.get_hr(residual) 
          sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Activity (Bq)")


    # isotopes
    this_plot_dir = plot_dir + "/by_isotope_2"
    std.make_dir(this_plot_dir)

    # Activity n_tree
    for target in results['a_tree'].keys():
      plot_list = []
      plot_file = this_plot_dir + "/" + isotopes.get_dir_hr(target) + ".eps"
      for residual in results['a_tree'][target].keys():
        if(isotopes.is_unstable(residual)):
          plot_array = numpy.zeros((201,2,),) 
          plot_array[:,0] = t[:]
          plot_array[:,1] = results['a_tree'][target][residual][:]
          plot_list.append([isotopes.get_hr(residual), plot_array])
      title = "Activity - target: " + isotopes.get_hr(target) + "  residuals: all"
      sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Activity (Bq)")


    # isotopes
    this_plot_dir = plot_dir + "/by_isotope_3"
    std.make_dir(this_plot_dir)

    # Activity n_tree
    for residual in results['a_by_radioactive_isotope'].keys():
      plot_list = []
      plot_file = this_plot_dir + "/" + isotopes.get_dir_hr(residual) + ".eps"
      if(isotopes.is_unstable(residual) and results['a_by_radioactive_isotope'][residual] is not None):
        plot_array = numpy.zeros((201,2,),) 
        plot_array[:,0] = t[:]
        plot_array[:,1] = results['a_by_radioactive_isotope'][residual][:]
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
          plot_array = numpy.zeros((201,2,),) 
          plot_array[:,0] = t[:]
          plot_array[:,1] = results['a_tree'][target][residual][:]
          plot_list.append([isotopes.get_hr(target) + " -> " + isotopes.get_hr(residual), plot_array])
    title = "Activity - all residual isotopes" 
    sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Activity (Bq)")





    # Individual isotopes
    this_plot_dir = plot_dir + "/by_target_isotope"
    std.make_dir(this_plot_dir)

    # Activity n_tree
    plot_list = []
    plot_file = this_plot_dir + "/by_target_isotope.eps"
    for target in results['a_by_target_isotope'].keys():
      plot_array = numpy.zeros((201,2,),) 
      plot_array[:,0] = t[:]
      plot_array[:,1] = results['a_by_target_isotope'][target][:]
      plot_list.append([isotopes.get_hr(target), plot_array])
    title = "Activity by Target Isotope"
    sim.plot_activity(plot_file, plot_list, title, "Time (s)", "Activity (Bq)")




    
    gammas_data_dir = data_dir + "/gammas"
    std.make_dir(gammas_data_dir)
    
    file_name = gammas_data_dir + "/eob_gammas.txt"
    numpy.savetxt(file_name, results['g_eob'], delimiter=",")

    file_name = gammas_data_dir + "/eos_gammas.txt"
    numpy.savetxt(file_name, results['g_eos'], delimiter=",")




    gammas_plot_dir = plot_dir + "/gammas"
    std.make_dir(gammas_plot_dir)

    plt.clf()
    plt.figure(figsize=(12,8))    
    plt.title('Predicted Gamma Lines End of Beam')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Activity (Bq)')
    plt.stem(results['g_eob'][:,0], results['g_eob'][:,1])
    plt.savefig(gammas_plot_dir + '/eob_gammas.eps', format='eps')
    plt.close('all') 



    gammas_plot_dir = plot_dir + "/gammas"
    std.make_dir(gammas_plot_dir)

    plt.clf()
    plt.figure(figsize=(12,8))    
    plt.title('Predicted Gamma Lines End of Sim')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Activity (Bq)')
    plt.stem(results['g_eos'][:,0], results['g_eos'][:,1])
    plt.savefig(gammas_plot_dir + '/eos_gammas.eps', format='eps')
    plt.close('all') 


    
    gammas_plot_dir = plot_dir + "/gammas/all"
    std.make_dir(gammas_plot_dir)
   
    for gk in results['g_all'].keys():
      plot_list  = []
      plot_array = numpy.zeros((201,2,),) 
      plot_array[:,0] = t[:]
      plot_array[:,1] = results['g_all'][gk][:]
      name = str(gk) + "_keV"
      plot_file = gammas_plot_dir + "/" + name + ".eps"
      plot_list.append([name, plot_array])
      sim.plot_activity(plot_file, plot_list, '')

    
    gammas_plot_dir = plot_dir + "/gammas/over_time"
    std.make_dir(gammas_plot_dir)

    per_plot = 20
    n = 0
    m = 0
    plot_list  = []
    for gk in sorted(results['g_all'].keys()):
      plot_array = numpy.zeros((201,2,),) 
      plot_array[:,0] = t[:]
      plot_array[:,1] = results['g_all'][gk][:]
      name = str(gk) + "_keV"
      plot_list.append([name, plot_array])
      n = n + 1
      if(n == per_plot):
        n = 0
        m = m + 1

        plot_file = gammas_plot_dir + "/gammas_" + str(m) + ".eps"
        sim.plot_activity(plot_file, plot_list, '')
        plot_list  = []
 
    if(len(plot_list)>0):  
      m = m + 1  
      plot_file = gammas_plot_dir + "/gammas_" + str(m) + ".eps"
      sim.plot_activity(plot_file, plot_list, '')



    # Total gamma activity and energy
    results['g_activity'] = numpy.zeros((201,),) 
    results['g_energy'] = numpy.zeros((201,),) 
    for target in results['a_tree'].keys():
      if(target not in results['g_energy_by_target_isotope'].keys()):
        results['g_energy_by_target_isotope'][target] = numpy.zeros((201,),) 
      for residual in results['a_tree'][target].keys():
        if(residual not in results['g_energy_by_residual_isotope'].keys()):
          results['g_energy_by_residual_isotope'][residual] = numpy.zeros((201,),) 
        if(isotopes.is_unstable(residual)):
          # Get array of gammas for the residual
          ga = isotopes.get_gammas_array(residual)
          if(not isinstance(ga, type(None))):

            # Loop through gammas
            for n in range(len(ga)):
              energy = ga[n,0]
              intensity = ga[n,1]

              # Save overall activity and energy from gammas
              results['g_activity'] = results['g_activity'] + results['a_tree'][target][residual][:] * intensity
              results['g_energy'][:] = results['g_energy'][:] + results['a_tree'][target][residual][:] * energy * intensity

              results['g_energy_by_target_isotope'][target][:]  = results['g_energy_by_target_isotope'][target][:] + results['a_tree'][target][residual][:] * energy * intensity

              results['g_energy_by_residual_isotope'][residual][:]  = results['g_energy_by_residual_isotope'][residual][:] + results['a_tree'][target][residual][:] * energy * intensity


    # Save gamma data
    gammas_data_dir_target = data_dir + "/gammas/target_isotopes"
    gammas_plot_dir_target = plot_dir + "/gammas/target_isotopes"

    std.make_dir(gammas_data_dir_target)
    std.make_dir(gammas_plot_dir_target)
    for target in results['g_energy_by_target_isotope'].keys():      
      out_array = numpy.zeros((201,2,),) 
      out_array[:,0] = t[:]
      out_array[:,1] = results['g_energy_by_target_isotope'][target]
      file_name = gammas_data_dir_target + "/" + isotopes.get_dir_hr(target) + ".txt"
      numpy.savetxt(file_name, out_array, delimiter=",")

      file_name = gammas_plot_dir_target + "/" + isotopes.get_dir_hr(target) + ".eps"
      sim.plot(file_name, out_array, title="", x_axis=None, y_axis=None)






    gammas_data_dir_residual = data_dir + "/gammas/residual_isotopes"
    std.make_dir(gammas_data_dir_residual)
    for residual in results['g_energy_by_residual_isotope'].keys():      
      out_array = numpy.zeros((201,2,),) 
      out_array[:,0] = t[:]
      out_array[:,1] = results['g_energy_by_residual_isotope'][residual]
      file_name = gammas_data_dir_residual + "/" + isotopes.get_dir_hr(residual) + ".txt"
      numpy.savetxt(file_name, out_array, delimiter=",")



    gammas_plot_dir = plot_dir + "/gammas"

    plot_list = []
    plot_array = numpy.zeros((201,2,),)
    plot_array[:,0] = t[:]
    plot_array[:,1] = results['g_activity'][:]
    plot_file = gammas_plot_dir + "/gamma_activity.eps"
    plot_list.append(['Gamma activity', plot_array])
    sim.plot_activity(plot_file, plot_list, '', 'Time (s)', 'Activity (Bq)')

    plot_list = []
    plot_array = numpy.zeros((201,2,),)
    plot_array[:,0] = t[:]
    plot_array[:,1] = results['g_energy'][:]
    plot_file = gammas_plot_dir + "/gamma_energy_ev.eps"
    plot_list.append(['Gamma activity', plot_array])
    sim.plot_activity(plot_file, plot_list, '', 'Time (s)', 'Energy (eV)')
    
    plot_list = []
    plot_array = numpy.zeros((201,2,),)
    plot_array[:,0] = t[:]
    plot_array[:,1] = 1.60218E-19 * results['g_energy'][:]
    plot_file = gammas_plot_dir + "/gamma_energy_j.eps"
    plot_list.append(['Gamma activity', plot_array])
    sim.plot_activity(plot_file, plot_list, '', 'Time (s)', 'Energy (J)')
    
    plot_list = []
    plot_array = numpy.zeros((201,2,),)
    plot_array[:,0] = t[:]
    plot_array[:,1] = (1.60218E-19 / (12.57 * 80)) * results['g_energy'][:]
    plot_file = gammas_plot_dir + "/gamma_energy_dose.eps"
    plot_list.append(['Gamma activity', plot_array])
    sim.plot_activity(plot_file, plot_list, '', 'Time (s)', 'Dose in Si')    
    








    plot_list = []
    plot_array = numpy.zeros((201,2,),)
    plot_array[:,0] = t[:]
    plot_array[:,1] = (((1.60218E-19 / (12.57 * 80)) * 3600) / 1.140771128E-04) * 100 * results['g_energy'][:]
    plot_file = gammas_plot_dir + "/gamma_energy_annual_dose.eps"
    plot_list.append(['Gamma activity', plot_array])
    sim.plot_activity(plot_file, plot_list, 'Percentage of Annual Dose per Hour', 'Time (s)', '% annual dose per hour')


    

    
    fh = open(data_dir + '/absorbed_dose.txt', 'w')  
    fh.write('\n')
    fh.write('\n')  
    fh.write('Gamma Dose - Beam End\n')
    fh.write('===================================================\n')
    fh.write('Activity/Bq                                    ' + str(results['g_activity'][100]) + '\n')
    fh.write('Power eV/s                                     ' + str(results['g_energy'][100]) + '\n')
    fh.write('Power J/s                                      ' + str(1.60218E-19 * results['g_energy'][100]) + '\n')
    fh.write('Dose Gy/s                                      ' + str((1.60218E-19 / (12.57 * 80)) * results['g_energy'][100]) + '\n')
    fh.write('Dose Gy/hr                                     ' + str(((1.60218E-19 / (12.57 * 80)) * 3600) * results['g_energy'][100]) + '\n')
    fh.write('Fraction of annual public dose per second      ' + str(((1.60218E-19 / (12.57 * 80)) * results['g_energy'][100]) / 3.2E-11) + '\n')
    fh.write('\n')
    fh.write('\n')
    fh.write('Gamma Dose - Sim End\n')
    fh.write('===================================================\n')
    fh.write('Activity/Bq                                    ' + str(results['g_activity'][-1]) + '\n')
    fh.write('Power eV/s                                     ' + str(results['g_energy'][-1]) + '\n')
    fh.write('Power J/s                                      ' + str(1.60218E-19 * results['g_energy'][-1]) + '\n')
    fh.write('Dose Gy/s                                      ' + str((1.60218E-19 / (12.57 * 80)) * results['g_energy'][-1]) + '\n')
    fh.write('Dose Gy/hr                                     ' + str(((1.60218E-19 / (12.57 * 80)) * 3600) * results['g_energy'][-1]) + '\n')
    fh.write('Fraction of annual public dose per second      ' + str(((1.60218E-19 / (12.57 * 80)) * results['g_energy'][-1]) / 3.2E-11) + '\n')
    fh.write('\n')
    fh.write('\n')
    fh.write('Absorbed Dose Calculations\n')
    fh.write('===================================================\n')
    fh.write('Absorbed dose assumptions:\n')
    fh.write('1. radiation from point, emitted isotropically\n')
    fh.write('2. 80Kg human\n')
    fh.write('3. 1m from point source\n')
    fh.write('4. 1m squared surface area\n')
    fh.write('5. all energy absorbed\n')
    fh.write('\n')
    fh.write('\n')
    fh.write('Dose Limits\n')
    fh.write('===================================================\n')
    fh.write('employees 18+             20 millisieverts/year\n')
    fh.write('trainees 18+              6 millisieverts/year\n')
    fh.write('public and under 18s      1 millisievert/year\n')
    fh.write('public and under 18s      1.140771128E-04 millisieverts/hour\n')
    fh.write('public and under 18s      3.2E-08 millisieverts/second\n')
    fh.write('public and under 18s      3.2E-11 sieverts/second\n')
    fh.write('\n')
    fh.write('Dose averaged over area of skin not exceeding 1cm2\n')
    fh.write('Source: http://www.hse.gov.uk/radiation/ionising/doses/\n')
    fh.write('\n')
    fh.write('\n')
    fh.close()
  

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


""")
    fh.write('\n')
    fh.close()






  # rr = f * nd * dt * 1.0E-28  
  @staticmethod
  def run_t(icode, projectile, nd, beam_flux, target_depth, exyz_data, residual_rrs):
    residual_rr = {}
    dnisotope = 0.0
    ion_count = len(exyz_data)
    for ion in exyz_data.keys():
      for i in range(len(exyz_data[ion])): 
        e = exyz_data[ion][i,0]
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
          rxs = talys.search_rxs(e, projectile, icode)

          for rcode in rxs.keys():
            if(rcode not in residual_rr.keys()):
              residual_rr[rcode] = 0.0
            dn = (1 / ion_count) * rxs[rcode] * beam_flux * nd * dt * 1.0E-28
            residual_rr[rcode] = residual_rr[rcode] + dn
            dnisotope = dnisotope - dn


    residual_rrs[icode] = residual_rr

  
  @staticmethod
  def plot(file_path, plot_data, title="", x_axis=None, y_axis=None):
    plt.clf()
    plt.figure(figsize=(12,8))    
    plt.title(title)
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.stem(plot_data[:,0], plot_data[:,1])
    plt.savefig(file_path, format='eps')
    plt.close('all') 


  @staticmethod
  def plot_activity(file_path, plot_list, title="", x_axis=None, y_axis=None, scale="linear"):

    if(x_axis == None):
      x_axis = 'Time (s)'
    if(y_axis == None):
      y_axis = 'Activity (Bq)'

    linestyle = ["solid", "dashed", "dashdot", "dotted"]
    colours = ["#009933","#CC9900","#0033CC","#CC3300","#6600CC","#000066","#333300","#FF0000","#FF9966","#009999","#993366","#66FF33","#FF0066","#00FFFF","#FFFF00","#990000","#663300"]

    
    plt.figure(figsize=(12,8))
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
    plt.yscale('log')
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.ylim(bottom=1.0E-2, top=1.0E15)
    plt.title(title)
    plt.legend(fontsize="xx-small", ncol=cols)
    plt.grid(True)
    plt.savefig(file_path, type='eps')
    plt.close('all') 


    """
    
    plt.figure(figsize=(12,8))
    n = 0
    for particle in pxs.keys():     
      pp_xs = pxs[particle]
      t_hr = isotopes.get_hr(target_code) 
      r_hr = particle
      label = t_hr.strip() + " - " + r_hr.strip()
      if(pp_xs is not None):
        lc = n % len(colours)
        ls = int(numpy.floor(n / len(colours))) % len(linestyle)

        plt.plot(pp_xs[:,0], pp_xs[:,1], label=label, color=colours[lc], linestyle=linestyle[ls])
        n = n + 1
    cols = (max(1,int(numpy.floor(len(pxs)/10))))

    talys.make_dir(path)
    if(len(pxs)>1):
      file_name = "pxs_" + str(target_code) + "_all.eps"
    else:
      file_name = "pxs_" + str(target_code) + "_" + str(particle) + ".eps"
    file_path = path + "/" + file_name

    plt.xlabel('Energy (MeV)')
    plt.ylabel('Cross-section (Barns)')
    plt.title('Target-Residual Cross Section')
    plt.legend(fontsize="xx-small", ncol=cols)
    plt.grid(True)
    plt.savefig(file_path, type='eps')
    plt.close('all') 

    """