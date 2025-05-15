from as_lib import * 
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-d', '--das',  help='Specify path of das')
parser.add_argument('-t', '--terms',  help='Specify path of TERMS')
parser.add_argument('-l', '--oic',  help='Specify path of oic')

parser.add_argument('-v', '--eigenvector',  help='Print out eigenvector expansion',action='store_true')

parser.add_argument('-n', '--num',  help='Requested number of states (all states by default)',type=int)
parser.add_argument('-o', '--num_orbitals',  help='number of orbitals to show in NRCSF (2 default)',type=int)

parser.add_argument('-u', '--unit',  help='Unit: Ryd=0,eV=1,cm=2. Default: 0',type=int)

args = parser.parse_args()


if not args.das:
    print('no das given. stopping')
else:

    if (not args.terms) and (not args.oic):
        print('no terms or oic given. stopping')
    else:

        if args.num:
            user_num_levels = args.num 
        else:
            user_num_levels = 2**63

        num_orbitals_shown = 2
        if args.num_orbitals:
            num_orbitals_shown = max(num_orbitals_shown,args.num_orbitals)


        das_file_numpy,orbital_strings,num_csfs,lambda_array = read_das(args.das)
        csf_strings,pseudo_array = make_csf_strings(das_file_numpy,orbital_strings,num_csfs,num_orbitals_shown,lambda_array)

        #todo: case statement for this 
        unit = 0
        factor=1
        unitstring = 'Ryd'
        if args.unit:
            if int(args.unit) == 0:
                unit = 0
                unitstring = ' Ryd '
                factor=1
            elif int(args.unit) == 1:
                unit = 1
                unitstring = 'eV '
                factor=13.605693122990

            elif int(args.unit) == 2:
                unit = 2
                unitstring = 'cm^-1'
                factor=109737.31568157
            else:
                unit = 0
                unitstring = 'Ryd'
                factor = 1.0
                print("Invalid unit specification. Defaulting to Ryd.")
    
    
        if args.eigenvector:
            #print('eigenvector argument given')
            g,gs = read_olg_for_groups()

            vector_list,block_lv_map = read_olg_for_ci_matrix(g,gs)
        
        
        if args.terms:
            read_terms_and_output(args.terms,csf_strings,user_num_levels,pseudo_array)
        
        if args.oic:
            #to do: put these into a subroutine
            states = read_oic_into_list_of_eigenstates(args.oic,csf_strings,user_num_levels,factor,unit)
            
            #header = 'Index,  Energy({:5}),     CSF(TERM),        J,     LV'.format(unitstring)
            
            
            
            #print(header)
            for (level_index,state) in enumerate(states):
                if args.eigenvector:
                    lv_m1 = states[level_index].lv_number - 1
                    block_index = block_lv_map[lv_m1]
                    group = g[block_index]
                    csfs_list = group[:,-2]
                    #print()
                    #print(csfs_list)
                    #print(vector_list[lv_m1])
                    orbLmap = ['S','P','D','F','G','H','I','K','L','M','N','O']
                    formatee = '{:6.2f}%'
                    arguments = np.argsort(abs(vector_list[lv_m1]))[::-1]
                    #print(17*' '+45*'*')
                    #print(17*' '+'Expansion:')
                    summ = 0 
                    counter = 0
                    string = ""
                    for ii in range(0,4):
                        vec_ind = arguments[ii]
                        percent = abs(vector_list[lv_m1][vec_ind]) * vector_list[lv_m1][vec_ind]*100
                        summ+= abs(percent)
                        counter += 1 
                        
                
                        
                        multiplicity = group[vec_ind,1]
                        orbL = orbLmap[group[vec_ind,2]]
                        term = '('+str(multiplicity)+orbL
                        if(int(state.parity) ==1):
                            term = term+'*'
                        else:
                            term = term+' '
                        term = term + ')'
                        string_print = '['+csf_strings[csfs_list[vec_ind]-1]+term+']'
                        string = string + formatee.format(percent) + string_print
                        #print(17*' ',formatee.format(percent)
                        #      ,string_print
                        #      )
                        
                        if((counter > 20)or(summ > 99)):
                            break
                    print('{:2}, {:2}, {:12}, {:5}, '.format(level_index+1,state.angular_momentum_j,state.energy_ryd,state.lv_number),string)
                    #print("Total % comp: ",summ)
                    #print(17*' '+45*'*')
                else:
                    state.display_state()




