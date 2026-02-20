from jj_library import * 
import library

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--file',  help='Specify path of GRASP.OUT')
parser.add_argument('-n', '--num',  help='Requested number of states (all states by default)',type=int)
parser.add_argument('-o', '--num_orbitals',  help='number of orbitals to show in NRCSF (2 default)',type=int)

args = parser.parse_args()



if not args.file:
    print('no grasp out given. stopping')
else:

    print('NOTE THAT YOU NEED MCDF 32 AND NOT ANG 7 AS THIS SUPPRESSES THE JJCSFS')

    num_orbitals_shown = 2
    if args.num_orbitals:
        num_orbitals_shown = max(num_orbitals_shown,args.num_orbitals)

    #grasp_out_path = 'legacy_testing/GRASP.OUT'
    grasp_out_path = args.file
    csf_array,orb_labels,num_nrcsf,num_orbs,mode = library.create_csf_array_from_output(grasp_out_path)
    sorted_orbital_array,csf_sorted_by_orbital = library.sort_orbitals(orb_labels,library.orbitals_order,csf_array)

    csf_strings_prepared,ads = library.make_csf_strings(csf_sorted_by_orbital,sorted_orbital_array,num_nrcsf,num_orbitals_shown)

    #print(csf_strings_prepared)

    map,num_rcsf = library.find_relativistic_csfs(grasp_out_path,num_nrcsf)

    print(num_rcsf,'rcsfs found')
    num_levels_to_be_shown = 1
    if args.num:
        if args.num < num_rcsf:
            num_levels_to_be_shown = args.num 
        else:
            num_levels_to_be_shown = num_rcsf
    else:
        num_levels_to_be_shown = num_rcsf 



    index = locate_eigenvectors(grasp_out_path)

    #print(index)
    num_components_to_be_shown = 3
    eigenvectors = read_in_eigenvectors(index,num_rcsf,grasp_out_path)

    x = find_levels_jj(grasp_out_path)

    levels = read_in_eigenlevels(grasp_out_path,x,num_rcsf)

    pos = locate_cfout_jj(grasp_out_path)
    jj_terms_raw = read_cfout_jj(grasp_out_path,pos,num_rcsf)
    jj_terms = decode_jj_terms(jj_terms_raw)
    #print(jj_terms)
    eigenstate_array = produce_eigenstate_array(levels,eigenvectors)
    #print(csf_strings_prepared[0])
    eigenstate_array = make_many_eigenstate_strings(eigenstate_array,num_components_to_be_shown,jj_terms,map,csf_strings_prepared)
    output_table(eigenstate_array,num_levels_to_be_shown)

