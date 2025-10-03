#!/usr/bin/python3
from library import * 
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--file',        help='Specify path of GRASP.OUT')
parser.add_argument('-n', '--num',         help='Requested number of states (all states by default)',type=int)
parser.add_argument('-i', '--inner_terms', help='show inner terms',action='store_true')
parser.add_argument('-nc', '--ncorinc',    help='print csfs needed for energy less than ncorinc',type=float)
parser.add_argument('-c', '--core',        help='Core   override, orbital index.',type=float)
parser.add_argument('-u', '--unit',        help='Unit: Ryd=0,eV=1,cm=2. Default: 0',type=int)
parser.add_argument('-a', '--adas',        help='Print an adasexj.in',action='store_true')

args = parser.parse_args()

def main(
         grasp_out_path,
         unit,
         user_num_levels,
         display_inner_terms
         ):
    
    csf_array,orb_labels,num_nrcsf,num_orbs,mode = create_csf_array_from_output(grasp_out_path)
    
    sorted_orbital_array,csf_sorted_by_orbital = sort_orbitals(orb_labels,orbitals_order,csf_array)

    if mode == 0:
        print('-------------')
        print('your GRASP printing is in mode 0 - pressing on but it is not likely your grasp.out contains the jj2ls data i want.')
        print('-------------')

    nelec = int(np.average(np.sum(csf_array,axis=1)))
    #print(nelec)
    
    core = 0 
    user_override_core = False
    try:
        if (args.core+1):
            user_override_core = True  
    except:
        user_override_core = False
        
    if user_override_core:
        core  = int(args.core)
    else:
        core  = find_core(csf_sorted_by_orbital,sorted_orbital_array)

    csf_strings_prepared,adas_strings = make_csf_strings(csf_sorted_by_orbital,sorted_orbital_array,num_nrcsf,core)
    #print(csf_strings_prepared)
    
    map,num_rcsf = find_relativistic_csfs(grasp_out_path,num_nrcsf)
    #print(csf_strings_prepared)
    if args.inner_terms:
        inner_terms = find_inner_occupation_terms(grasp_out_path,mode,num_rcsf)
    else:
        inner_terms = []
    #print(inner_terms)
    states,charge = find_levels(grasp_out_path,inner_terms,csf_strings_prepared,map,adas_strings)

    if display_inner_terms:
        states = add_inner_occupation_strings_to_eigenclass(grasp_out_path,mode,states,display_inner_terms,csf_strings_prepared,map)

    output_table(states,user_num_levels,unit)
    
    adas = False 
    if (args.adas):
        adas = True
    if adas:
        write_adasexjin(states,user_num_levels,charge,nelec)        
    
    if args.ncorinc: 
        cut = float(args.ncorinc)
        
        print("Finding all csfs with energy less than ",cut)
        csfs = []
        for (index,state) in enumerate(states):
            check = False 

            if(state.eigenenergy >  cut):
                print("exiting at state",index)
                break 
            stringIwant = state.leading_csf
            #print(stringIwant)
            for string in csfs:
                if (string == stringIwant): 
                    check = True 
                    break
            if not(check):
                csfs.append(stringIwant)
        for string in csfs:
            print(string)
        #print(csf_strings_prepared)

    return 0

try:
    if (args.unit+1):
        unit = int(args.unit) 
except:
    unit = 0 

if not args.file:
    print('no grasp out given. stopping')
else:
    if args.num:
        user_num_levels = args.num 
    else:
        user_num_levels = 0
        

    #num_orbitals_shown = 2
    #if args.num_orbitals:
    #    num_orbitals_shown = max(num_orbitals_shown,args.num_orbitals)
    display_inner_terms = False
    if args.inner_terms:
        display_inner_terms = True
        
    main(args.file,
         unit,
         user_num_levels,
         display_inner_terms)



