import numpy as np 

ANGULAR_SYMBOLS = ['s','p','d','f','g','h','i','k','l','m','o','p']

def display_states(states,header,ground):
    
    #displays the eigenvec
    
    print(header)
    
    for state in states:
        shift = states[ground].energy_ryd
        state.display_state(shift)
        
    return 0 

def readConfigFile(path_to_config_dot_dat):
    '''
    This routine will read the CONFIG.DAT file,
    which is perhaps a more reliable way of obtaining the data than in das,
    which may or may not be formatted in a particular way.
    
    The code is (probably) similar to that of autoDarc 
    https://github.com/LeoMul/autoDarc/blob/main/write_dstg2inp.py
    
    '''
    
    f = open('CONFIG.DAT','r')

    firstline = f.readline().split()


    second_line = f.readline().split()
    
    num_orbs = int(len(second_line) / 2 )

    third_line = f.readline().split()
    
    num_csfs = int(third_line[0])
    
    princN = []
    orbital_L = [] #np.zeros(num_orbs,dtype=int)
    orbital_strings = [] 
    
    for ii in range(0,len(second_line),2):
        princN.append(int(second_line[ii]))
        orbital_L.append(int(second_line[ii+1]))
        string = second_line[ii]+translate_angular(int(second_line[ii+1]))
        #print(string)
        orbital_strings.append(string)

        #print(orbL[ii])
        #string = str(princN[ii])+angular_symbols[orbL[ii]]
        #orbs.append(string)

#   orbital_strings = []
#   orbital_L = []
#
#   for kk in range(0,len(orbital_data),2):
#       orbital_strings.append(orbital_data[kk]+translate_angular(int(orbital_data[kk+1])))
#       orbital_L.append(int(orbital_data[kk+1]))

    #num_csfs = int(f.readline().split()[0])

    f.close()


    #easiest way for me to do this:
    das_file_numpy = np.loadtxt(path_to_config_dot_dat,skiprows=5,dtype=int)[:,0:-1]
    #print(das_file_numpy)
    #temperory
    lambda_array = np.ones(num_orbs)
    
    #print(orbital_strings)
    
    return das_file_numpy,orbital_strings,num_csfs,lambda_array,orbital_L





def read_das(path_to_das):

    das_file_opened = open(path_to_das,'r')
    das_file_read = das_file_opened.readlines()

    orbital_data = das_file_read[2].split()

    total_num_orbs = int( len(orbital_data) / 2 ) 
    
    orbital_strings = []
    orbital_L = []

    for kk in range(0,len(orbital_data),2):
        orbital_strings.append(orbital_data[kk]+translate_angular(int(orbital_data[kk+1])))
        orbital_L.append(int(orbital_data[kk+1]))

    #print(orbital_strings)

    for jj in range(0,len(das_file_read)): 
        current_line = das_file_read[jj]

        if current_line.split()[0] == '&SMINIM':
            #print(jj-3,'csfs found')
            break    
    
    num_csfs = jj - 3

    das_file_opened.close() 

    #assumes that the namelist is one line.. porbably needs to he changed.
    das_file_numpy = np.loadtxt(path_to_das,skiprows=3,max_rows=num_csfs,dtype=int)
    #print(np.shape(das_file_numpy))
    
    #if there's only 1 csf, numpy breaks so just rehsape
    if np.shape(das_file_numpy) == (total_num_orbs,):
        das_file_numpy = das_file_numpy.reshape((1,total_num_orbs))

    lambda_array = np.zeros(total_num_orbs)
    for jj in range(0,len(das_file_read)):
        if das_file_read[jj].split()[0] == '&SMINIM':
            #print('lambda namelist found')
            break 
        
    #need to make this a bit better
    
    #right now as assume includ and nvar = 0
    lambda_collection = []
    for ii in range(jj+1,len(das_file_read)):
        if len(das_file_read[ii].split()[0])==0:
            print('skipping lambdas - couldnt find them')
        else:
            if das_file_read[ii].split()[0] != '&SRADWIN':
                lambda_collection.extend(das_file_read[ii].split())
            else:
                break
    try:
        if len(lambda_array) >0:     
            for ii in range(0, len(lambda_array)):
                lambda_array[ii] = float(lambda_collection[ii])
    except:
        lambda_array = np.zeros(100)
    #print('Lambda array found:')
    #print(lambda_array)

    return das_file_numpy,orbital_strings,num_csfs,lambda_array,orbital_L

def findCore(das_file_numpy,orbital_L):
    ncsfs,num_orbitals = np.shape(das_file_numpy)
    print(num_orbitals,ncsfs)
    
    print(das_file_numpy)
    
    orbital_L = np.array(orbital_L)  
    print(orbital_L)  
    full_L = (2*(2*orbital_L+1))
    print(full_L)
    core = 0
    exit_outer = False
    for j in range(0,num_orbitals):
        for i in range(0,ncsfs):
            if das_file_numpy[i,j] != full_L[j]:
                print(orbital_L[j],das_file_numpy[i,j])
                exit_outer = True
        if(exit_outer):
            break
        core+= 1 
    
    return core

def make_csf_strings(das_file_numpy,orbital_strings,num_csfs,num_requested_orbitals,lambda_array,orbital_L,core):
    csf_strings = []
    pseudo_array = []
    
    #print(findCore(das_file_numpy,orbital_L))
    

    for jj in range(0, num_csfs):
        current_csf_ints = das_file_numpy[jj][core:]
        concerned_occupations_locations = np.argwhere(current_csf_ints > 0 )

        #print(concerned_occupations_locations)
        concerned_occupations_numbers = current_csf_ints[concerned_occupations_locations]
        concerned_lambdas = lambda_array[concerned_occupations_locations]
        concerned_occupations_orbitals = []

        for kk in range(0,len(concerned_occupations_locations)):
            concerned_occupations_orbitals.append(orbital_strings[concerned_occupations_locations[kk][0]+core])

        pseudo = 0

        if any(concerned_lambdas < 0.0):
            pseudo = 1
        pseudo_array.append(pseudo)
        
        num_orbitals = len(concerned_occupations_orbitals)
        current_csf_string = ''

        min = max(0,len(concerned_occupations_locations) - num_requested_orbitals)
        cf_format = '{:>3}{:3} '
        for kk in range(min,num_orbitals):
            current_csf_string += cf_format.format(concerned_occupations_orbitals[kk], str(concerned_occupations_numbers[kk][0]))# + " "
        #  print(current_csf_string)

        csf_strings.append(current_csf_string)
    return csf_strings,pseudo_array

def read_terms_and_output(terms,csf_strings,num_levels,pseudo_array):

    terms_data = np.loadtxt(terms,skiprows=1) 
    #print(np.shape(terms_data))
    if num_levels > np.shape(terms_data)[0]:
        num_levels = np.shape(terms_data)[0]

    header = 'index, level (ryd), config, psuedo indicator'
    print(header)
            
    output_string = '{:5},   {:10.6f},     {},     {:5},     {}'


    #print(num_levels)

    for kk in range(0,num_levels-1):
        line = terms_data[kk]
        multiplicity = line[0]
        angular = line[1]
        parity = line[2]
        cf_number = int(line[3]-1)
        energy_ryd = line[5]

        term_string = '('+str(int(multiplicity)) + translate_angular(int(angular))

        if int(parity) == 0: 
            term_string += ' '
        elif int(parity) == 1:
            term_string += '*'
        else:
            print('failure in parity idenficiation at level',kk+1,'parity found',parity)

        term_string += ')'
        print(output_string.format(kk+1,energy_ryd,csf_strings[cf_number] + term_string.upper(),kk+1,pseudo_array[cf_number]))



    return 0

def read_levels_and_output(levels,csf_strings,num_levels):

    terms_data = np.loadtxt(levels,skiprows=1) 
    #print(np.shape(terms_data))
    if num_levels > np.shape(terms_data)[0]:
        num_levels = np.shape(terms_data)[0]


            
    output_string = '{:5},   {:12.8f},     {}  {}'


    #print(num_levels)

    header = 'Index       Energy(Ry)     CSF(TERM)        J'
    print(header)
    for kk in range(0,num_levels-1):
        line = terms_data[kk]

        j = int(line[0]) / 2
        parity = line[1]

        multiplicity = abs(line[2])
        angular = line[3]
        cf_number = int(line[4]-1)
        energy_ryd = line[-1]

        term_string = '('+str(int(multiplicity)) + translate_angular(int(angular))

        if int(parity) == 0: 
            term_string += ' '
        elif int(parity) == 1:
            term_string += '*'
        else:
            print('failure in parity idenficiation at level',kk+1,'parity found',parity)

        term_string += ')'
        print(output_string.format(kk+1,energy_ryd,csf_strings[cf_number] + term_string.upper(),j))



    return 0

import linecache

def read_oic_and_output(oic,csf_strings,num_levels,unit):
    #format string in as for oic output is
    #10240 FORMAT(7I5,F15.6,I10)
    #for this routine i modified it to: 
    #10240 FORMAT(7(I5,1X),F15.6,I10)
    #better solution is a fortran format read using the 
    #fortranformat package.
    
    x = linecache.getline(oic, len(csf_strings) + 6).split()
    level_data = np.loadtxt(oic,skiprows=len(csf_strings) + 7,max_rows=int(x[1])) 

    #print(np.shape(terms_data))
    if num_levels > np.shape(level_data)[0]:
        num_levels = np.shape(level_data)[0]

    output_string = '{:5},   {:12.8f},     {}  {}   {:4}'


    #print(num_levels)

    header = 'Index,       Energy(Ry),     CSF(TERM),        J,     LV'
    print(header)
    for kk in range(0,num_levels):
        line = level_data[kk]

        j = int(line[5]) / 2
        lv = int(line[1])
        multiplicity = line[3]
        angular = line[4]
        cf_number = int(line[6]-1)
        energy_ryd = line[-1]

        term_string = '('+str(int(abs(multiplicity))) + translate_angular(int(angular))

        if int(multiplicity) > 0: 
            term_string += ' '
        elif int(multiplicity) < 0:
            term_string += '*'
        else:
            print('failure in parity idenficiation at level',kk+1,'parity found',multiplicity)

        term_string += ')'
        print(output_string.format(kk+1,energy_ryd,csf_strings[cf_number] + term_string.upper(),j,lv))

    return 0

def read_oic_into_list_of_eigenstates(oic,csf_strings,num_levels,factor,unit,core,override):

    #skip override
    
    skip = len(csf_strings)

    x = linecache.getline(oic, 4).split()
    #print(x)
    #for each config, as outputs a hex code. this can somestimes extend to two lines
    #this detectst this.
    if len(x) == 1:
        skip *= 2

    print('skip=',skip)
    x = linecache.getline(oic, skip + 6).split()
    #print(x)
    #evil hack, idk why it works this way. 
    if len(x) >4:
        skip = skip - 1
        x = linecache.getline(oic, skip + 6).split()
        #print(x)

    if override != 0:
        print('override=',override,'skiprows=',skip + 7 + override)
    #print(int(x[1]))
    level_data = np.loadtxt(oic,skiprows=skip + 7 + override,max_rows=int(x[1])) 
    #print(level_data)
    #print(np.shape(terms_data))
    if num_levels > np.shape(level_data)[0]:
        num_levels = np.shape(level_data)[0]

    #header = 'Index       Energy(Ry)     CSF(TERM)        J     LV'
    #print(header)

    states = []
    fracJ = False
    for kk in range(0,num_levels):
        line = level_data[kk]
        #print(line)
        if int(line[5]) % 2 ==1:
            fracJ = True
            
        j = int(line[5]) / 2
        #j ="1/2"
        if fracJ:
            j_string = str(int(line[5]))+'/2'
        else:
            j_string = str(j)
        
        lv = int(line[1])
        multiplicity = line[3]
        angular = line[4]
        cf_number = int(line[6]-1)
        energy_ryd = line[-1]*factor

        parity = 0
        if multiplicity < 0:
            parity = 1 

        state = energy_eigenstate_as_ic(
            energy_ryd,
            kk+1,
            j,
            angular,
            abs(multiplicity),
            parity,
            cf_number,
            csf_strings,
            lv,
            unit,
            factor,
            j_string
        )
        states.append(state)

    lengths = []

    for state in states:
        lengths.append(len(state.label_string))

    if (core == 0):
        max_length = max(lengths)
        #print(max_length)
        for state in states:
            length = len(state.label_string)
            if len(state.label_string) < max_length:
                #print('hello')
                state.label_string = state.label_string + (max_length-length)*' '

    return states

def translate_angular(angular_number):
    if angular_number < len(ANGULAR_SYMBOLS):
        return ANGULAR_SYMBOLS[angular_number]
    else:
        return str(angular_number)

RYDBERG_CM = 109737.316 


class energy_eigenstate_as_ic:

    def __init__(self,
                 energy_ryd,
                 level_index,
                 angular_momentum_j,
                 angular_momentum_L,
                 angular_momentum_S,
                 parity,
                 cf_number,
                 csf_strings,
                 lv_number,
                 unit,
                 factor,
                 j_string,
                 numorbsshown=[]
                 ):

        self.energy_ryd         = energy_ryd
        self.energy_wavenumber  = energy_ryd * RYDBERG_CM
        self.angular_momentum_j = angular_momentum_j
        self.angular_momentum_L = angular_momentum_L
        self.angular_momentum_S = angular_momentum_S
        self.parity             = parity
        self.cf_number          = cf_number
        self.level_index        = level_index
        self.lv_number = lv_number
        self.stat_weight = 2.0 * angular_momentum_j + 1.0
        self.j_string = j_string
        
        #output stuff
        self.output_string = '{:5},'
        if(unit==2):
            self.output_string+= '{:12.3f}'
        else:
            self.output_string+= '{:12.8f}'    
        self.output_string += ',     {},  {:>4},     {}'
        
        
        term_string = '('+str(int(angular_momentum_S)) + translate_angular(int(angular_momentum_L))

        if int(parity) == 0: 
            term_string += ' '
        elif int(parity) == 1:
            term_string += '*'
        else:
            print('failure in parity idenficiation at level',level_index,'parity found',parity)
        
        term_string += ')'

        self.label_string = csf_strings[cf_number]+term_string.upper()



    def display_state(self,shift):

            self.output_string = self.output_string.format(
            self.level_index,
            self.energy_ryd-shift,
            self.label_string,
            self.j_string,
            self.lv_number
        )
            print(self.output_string)

limit = 3

def get_transitions(path,states):
    f = open(path,'r')
    limit = 2**63-1
    numStates = len(states)
    avalues = np.zeros([numStates,numStates])
    
    
    #The formatting for OIC E1 data is: 10420, line 46400 is the call.

    for ii in range(limit):
        line = f.readline().split()
        if len(line) > 1:
            if line[1] == '1-DATA':
                break 
    transitions = []
    runlimit = 1
    add = False
    for jj in range(0,limit):
        line = f.readline().split()
        if runlimit>=5:
            break

        if len(line) > 0:
            #todo: check if it runs over.
            #print(line) 

            upper_level = int(line[1])
            lower_level = int(line[2])
            a_value_len = float(line[3])
            
            avalues[upper_level-1,lower_level-1] += abs(a_value_len)
            avalues[lower_level-1,upper_level-1] += abs(a_value_len)

            if (add):
                avalues[upper_level-1,lower_level-1] += abs(float(line[4]))
                avalues[lower_level-1,upper_level-1] += abs(float(line[4]))
            #print(a_value_len)
            gf_len = float(line[5])
            if '*' not in line[8]:
                lambda_ang_vac = float(line[8])
            else:
                lambda_ang_vac = np.inf

            gf_vel = float(line[9])

            transition = transition_as_ic(
                upper_level,
                lower_level,
                a_value_len,
                gf_len,
                gf_vel,
                lambda_ang_vac,
                states
            )
            transitions.append(transition)
        else:
            line = f.readline().split()
            
            print('breaking? ',line)
    
            runlimit = runlimit + 1
            add = True 
            #line = f.readline().split()
        
        
        
        
    return transitions,avalues

class transition_as_ic:
    def __init__(self,
                 upper_level,
                 lower_level,
                 a_value_len,
                 gf_len,
                 gf_vel,
                 lambda_ang_vac,
                 energy_eigenstates_array:list[energy_eigenstate_as_ic]):
        
        self.upper_level_int = upper_level
        self.lower_level_int = lower_level



        self.upper_level_csf_string = energy_eigenstates_array[upper_level-1].label_string
        self.lower_level_csf_string = energy_eigenstates_array[lower_level-1].label_string

        self.lower_ang = energy_eigenstates_array[lower_level-1].angular_momentum_j
        self.upper_ang = energy_eigenstates_array[upper_level-1].angular_momentum_j

        self.a_value_len = abs(a_value_len)

        if gf_len != 0:
            self.vel_len_ratio = abs(gf_vel / gf_len)
        else:
            self.vel_len_ratio = np.inf

        self.gf_len = abs(gf_len)
        self.f_len = abs(gf_len / energy_eigenstates_array[lower_level-1].stat_weight)
        self.ga = abs(a_value_len * energy_eigenstates_array[upper_level-1].stat_weight)
        self.log_gf = np.log10(abs(gf_len))
        self.gf_vel = abs(gf_vel)

        
        self.upper_level_ryd = energy_eigenstates_array[upper_level-1].energy_ryd
        self.lower_level_ryd = energy_eigenstates_array[lower_level-1].energy_ryd

        self.upper_level_wavenumber = energy_eigenstates_array[upper_level-1].energy_wavenumber
        self.lower_level_wavenumber = energy_eigenstates_array[lower_level-1].energy_wavenumber

        self.wavelength_ang_vac = lambda_ang_vac 
        self.wavelength_nnm_vac = lambda_ang_vac * 0.1 
        self.wavelength_ang_air = vac_to_air(lambda_ang_vac)
        self.wavelength_nnm_air = self.wavelength_ang_air * 0.1 


    
    def display_transition(self):
        string = '{:5} {:5} {:14.4F}  {:3.1F}   {:14.4F}  {:3.1F}    {:14.4F}    {:14.4F}  {:10.2E} {:10.2E}    {:10.2E}  {:10.2E}  {:10.4f} {:10.4f}    {}    {}'
        
        #if (self.lower_eigenstate.shifted_energy_cm != -1.0) and (self.upper_eigenstate.shifted_energy_cm != -1.0):
        #    lower_energy = self.lower_eigenstate.shifted_energy_cm
        #    upper_energy = self.upper_eigenstate.shifted_energy_cm
        #else:
        #    lower_energy = self.lower_eigenstate.wavenumber
        #    upper_energy = self.upper_eigenstate.wavenumber           
        
        print(string.format(self.lower_level_int,\
                            self.upper_level_int,\
                            self.lower_level_wavenumber,\
                            self.lower_ang,\
                            self.upper_level_wavenumber,\
                            self.upper_ang,\
                            self.wavelength_nnm_vac,\
                            self.wavelength_nnm_air,\
                            self.a_value_len,\
                            self.ga,\
                            self.f_len,\
                            self.gf_len,\
                            self.log_gf,\
                            self.vel_len_ratio,\
                            self.lower_level_csf_string,\
                            self.upper_level_csf_string))

def print_out_a_values(total_data_class_array:list[transition_as_ic]):
    
    
    #string = '{:5} {:5} {:14.4F}  {:3.1F}   {:14.4F}  {:3.1F}    {:14.4F}  {:10.2E} {:10.2E}    {:10.2E}  {:10.2E}  {:10.2E}    {}    {}'

    header = '#{:>4s} {:>5s}   {:14} {:3}    {:14} {:3}   {:>14s}    {:>14s}  {:>10s} {:>10s}     {:>10s} {:>10s}  {:>10s} {:>10s}   {:>20s}  {:>20s}'
    print(header.format(
            'Low','Upp'
            ,'Lower (cm-1)','JL'
            ,'Upper (cm-1)','JU'
            ,'λ (vac,nm)'
            ,'λ (air,nm)'
            ,'A (s-1)','gA (s-1)'
            ,'f ','gf','log(gf)'
            ,'vel\len'
            ,'LowDesignation','UppDesignation'
    ))
    for trans in total_data_class_array:
        trans.display_transition()

    return 0

def vac_to_air(wl_ang_vac):

    #converts a vaccumm wl (ang) to air. there are supposedly many conversion.
    #the commented one is from vald, the other is from https://www.sdss3.org/dr8/spectro/spectra.php
    
    #they give identical results as far as I can see. The other one might be cheaper tbf.

    if (wl_ang_vac > 2000.0 and wl_ang_vac < 100000.0):
        ssq = 1e8 / (wl_ang_vac * wl_ang_vac)
        #n =  1.0 + 0.0000834254 + 0.02406147 / (130.0 - ssq) + 0.00015998  / (38.9 - ssq)
        n = 1.0 + 2.735182E-4 + 131.4182 / wl_ang_vac**2 + 2.76249E8 / wl_ang_vac**4
        return wl_ang_vac / n
    else:
        return wl_ang_vac
    
    
def read_olg_for_groups():
    f = open('olg','r')
    checker = True 
    ii = 0
    #find data
    while checker: 
        ii+=1 
        line = f.readline().split()
        if len(line) > 3: 
            if line[-3] == 'LEVEL':
                checker = False
    print('jgroup data found at olg line',ii)
    
    done = False 
    
    current_group = []
    groups = []
    group_sizes = []
    group_2j = []
    num_levels = 0
    num_blanks = 0
    num_in_this_group = 0
    while not done:
        line = f.readline().split()
        #we are in a block
        if   len(line) == 6:
            if (line[0]!= 'LV'):
                num_blanks = 0
                #print(line)
                current_group.append(line)
                num_levels = num_levels+1
                num_in_this_group+= 1
 
        
        #end of a group
        elif len(line) == 7 or len(line) == 8:
            #print(line)
            #print(current_group[0])
            data = np.array(current_group,dtype=int)
            #print(data)
            groups.append(data)
            group_sizes.append(num_in_this_group)
            current_group = [] 
            num_in_this_group = 0

            
        elif len(line) == 0:
            num_blanks +=1 
        if num_blanks == 3:
            done = True
    #print(groups[0][0,1])
    f.close()
    print('total levels = ',num_levels)
    print('jpi block sizes are',group_sizes)
    return groups,group_sizes



    



def read_olg_for_ci_matrix(groups,group_sizes):
    total = 0
    f = open('olg','r')
    def read_eigenvector(blocksize,line):
        #blocksize = 41
        vector = np.zeros(blocksize)
        veciter = 0
        vector_isolated = line[10:]
        #print(line[10:])
        limit = 6 
        for ii in range(0,min(len(vector_isolated),blocksize)):
            veciter+=1 
            if line[ii][0] == '-':
                limit = 8 
            else:
                limit = 7
            #vector[ii]
            vector[ii] = float(vector_isolated[ii][0:limit].replace("*",""))
        #print(vector,veciter)

        while (veciter < blocksize):
            line = f.readline().split()
            for ii in range(0,len(line)):
                if veciter>=blocksize:
                    break
                #print(line[ii])
                if line[ii][0] == '-':
                    limit = 7 
                else:
                    limit = 6
                if ((veciter == blocksize-1) and len(line[ii] )>limit ):
                    vector[veciter] = float(line[ii][0:limit])
                else:
                    vector[veciter] = line[ii]
                veciter+=1 
                #print(veciter,blocksize)


            #print(vector,veciter)
        norm = np.dot(vector,vector)
        if (abs(norm-1.0) >1e-3):
            print('potential norm problem')
        return vector

    checker = True 
    ii=0
    while checker: 
        ii+=1 
        #TODO: this relies on there being a space.. 
        #this was gauranteed by making the global repaclement
        # 10050 FORMAT(I5,I4,4I3,F14.5,F14.6,I3,I4,10F8.4)

        # 10050 FORMAT(I5,1X,I4,4I3,F14.5,F14.6,I3,I4,10F8.4)
                #but - probably better to count character numbers and read that way. Or!
                #reading with fortran format.
        line = f.readline().split()
        if len(line) > 9: 
            if line[9] == 'CI-MATRIX':
                checker = False 
                #print(ii)
                
    level_counter = 0
    
    strings_in_current_vector = []
    
    num_blocks = len(groups)
    num_levels = sum(group_sizes)
    block_lv_map = np.zeros(num_levels,dtype=int)

    requested_level = 0
    vector_list = []
    for jj in range(0,num_blocks):
        blocksize = group_sizes[jj]
        #print('new block',blocksize)
        for kk in range(0,blocksize):
            total+=1

            #print(jj,blocksize,kk,total)
            checker = True 
            #print(ii)
            block_lv_map[requested_level] = int(jj)
            requested_level +=1
            #print(requested_level)

            while checker:
                line = f.readline().split()
                if len(line) > 1:
                    if float(line[0]) == requested_level:
                        #we are in the correct level
                        checker = False
            #print('reading level',requested_level,'block size = ',blocksize)
            vector = read_eigenvector(blocksize,line)
            vector_list.append(vector)

    f.close()
    

    return vector_list,block_lv_map