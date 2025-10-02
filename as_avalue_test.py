from as_lib import * 

path = 'olg'

das_file_numpy,orbital_strings,num_csfs,lambda_array = read_das('/Users/leomulholland/CeIII_AS/opt3_icr/shift/das_test')
csf_strings,pseudo_array = make_csf_strings(das_file_numpy,orbital_strings,num_csfs,5,lambda_array)

states = read_oic_into_list_of_eigenstates('oic',csf_strings,2**63,1,'Ryd')


transitions,avalues = get_transitions(path,states)
numStates = len(states)
avalues = avalues + 1e-30
for ii in range(0,37):
    for jj in range(ii+1,37):
        print('{:4} {:4} {:10.2E}'.format(jj+1,ii+1,avalues[ii,jj]))

#print_out_a_values(transitions)
