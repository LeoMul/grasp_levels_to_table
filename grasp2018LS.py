#!/usr/bin/python3
import numpy as np  
RYDBERG_CM = 109737.316 
DIMENSION = 1000 
import argparse
import sys 
orbital_angular_momentum_dictionary = np.array(['S','P','D','F','G','H','I','K','L','M','N','O'])
full_orbital_dic = np.zeros_like(orbital_angular_momentum_dictionary,dtype=int)
for ii in range(0,len(full_orbital_dic)):
    full_orbital_dic[ii] = 2 * (2*ii+1)


parser = argparse.ArgumentParser()
parser.add_argument('-p', '--prefix',help='Specify prefix of jj2lsj output',type=str)
parser.add_argument('-n', '--num',   help='Requested number of states (all states by default)',type=int)
parser.add_argument('-s', '--shift',   help='path of shifts (assumes in the same order as unshifted levels)',type=str)

args = parser.parse_args()

prefix = args.prefix 

if not(args.prefix):
    print(' no prefix specified, e.g the name of your vector.')
    sys.exit()    


adas_format = '{:>3}{:1}' 
occupation_adas_character = [None,'1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H']
elements=['H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg'
        ,'Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn'
        ,'Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr'
        ,'Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb'
        ,'Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd'
        ,'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir'
        ,'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th'
        ,'Pa','U ']

energies  = []
vectors   = [] 
jvalues   = [] 

f = open(args.prefix+'.lsj.lbl')
data = f.readlines()
#skip header
numblocks  = 0 
numvectors = 0 
thisvector = []
for ii in range(1, len(data)):
    string = data[ii]
    stringsplit = string.split()
    length = len(stringsplit)
    match length:
        case 5:
            #print('Vector head')
            if len(thisvector) > 0 :
                vectors.append(thisvector)
            energies.append(float ( stringsplit[3]))
            jvalues.append(stringsplit[1])
            numvectors += 1 
            if (numvectors > DIMENSION):
                print('dimension problem ?')
            thisvector = [] 
        case 3:
            #print('Vector element')
            thisvector.append(stringsplit[-1])
        case 0:
            numblocks += 1
        case _: 
            print('Invalid, error? idk',length,string)
            

f.close()
if len(thisvector) > 0 :
    vectors.append(thisvector)
energies = 2.0 * np.array(energies) #(Rydberg)
ground = min(energies)
energies = energies - ground 
#
numrequest = args.num 
if (not args.num):
    numrequest = numvectors
numrequest = min(numrequest,numvectors )
#
sortargs = np.argsort(energies)

if args.shift:
    shift = np.loadtxt(args.shift)
    energies_copy = energies[sortargs]
    energies_copy[0:len(shift)] = shift 
    sortargs_two = np.argsort(energies_copy)
    energies[sortargs] = shift
    sortargs = sortargs[sortargs_two]

def getFullOcc(shell:str):
    angL = shell[-1].upper()
    try:
        occ = np.argwhere(orbital_angular_momentum_dictionary == angL)
        occ = full_orbital_dic[occ[0][0]]
    except:
        print(angL)
        print('ailure in getFullOcc')
        sys.exit()
    return int(occ)

def isShellFull(shell:str,thisocc:int):
    occ = getFullOcc(shell)
    #print(occ,angL,thisocc,int(occ) == int(thisocc))
    return int(occ) == int(thisocc)


def processCSFStringIntoShells(csf:str):
    #print(csf)
    occupations = [] 
    shells = [] 
    term = ''
    string = ''
    numshells = 0
    bb = False 
    for ii in range(0,len(csf)):
        char = csf[ii]
        #print(char)
        match char:
            case '.':
                #print('single occupied shell, or full shell!')
                if len(string) > 0:
                    occupations.append(1)
                    shells.append(string)
                    numshells += 1
                string = ''
                if isShellFull(shells[numshells-1],occupations[numshells-1]):
                    string = csf[ii+1:ii+2]  #hack fixe
                bb = False
            case '(':
                bb = True 
                #print('more than 1 electron in shell')
                occupations.append(int(csf[ii+1]))
                shells.append(string)
                numshells += 1
                string = ''
            case '_':
                #print('overall term')
                term = csf[ii+1:ii+3]
            case _:
                if ((csf[ii-1] != '(' and csf[ii-2] != '(') and csf[ii] != '('):
                    if ((csf[ii-1] != ')' and csf[ii-2] != ')') and csf[ii] != ')'):
                        if ((csf[ii-1] != '_' and csf[ii-2] != '_') and csf[ii] != '_'):
                            #print('     ',csf[ii-1] ,csf[ii-2])
                            string = string + csf[ii]
    if len(string) > 1:
        shells.append(string)
        if (not bb):
            occupations.append(1)
    #print(occupations,shells,term)
    adasstring = ''

    for ii in range(0, len(occupations)):
        adasstring += adas_format.format(shells[ii].upper(),occupation_adas_character[occupations[ii]])

    if (sum(occupations) != 3):
        print(occupations,shells)
        #sys.exit()
    
    occ= sum(occupations)
    return adasstring,term,occ

def proceessTermStringIntoInts(term: str):
    try:
        L = np.argwhere(orbital_angular_momentum_dictionary == term[-1])[0][0]
    except:
        L = -1
        print('failure in term L determination')
    twoSP1= int(term[0])
    return L,twoSP1 

def processJStringIntoFloat(j:str):
    l = len(j)
    pos = j.find('/')
    if pos > 0: 
        jfloat = 0.5 * float(j[0:pos]) 
    else:
        jfloat = float(j)
        
    return  jfloat

from pathlib import Path
#
adasin_format = '{:>5}{:<19}({:1}){:1}({:4.1f}){:>21.4f}'
ROOT_DIR = Path(__file__).parent
TEXT_FILE = ROOT_DIR / 'ion_energy.dat'
ip = np.loadtxt(TEXT_FILE)
TEXT_FILE_terms = ROOT_DIR / 'ion_terms.dat'
tt = np.loadtxt(TEXT_FILE_terms,dtype='>U2')

#get charge
g = open('isodata','r')
g.readline()
x = g.readline().split()
charge = int(float(x[0]))
g.close()
#get num elec:
#g = open()
g = open(args.prefix+'.c')
g.readline()
ff = g.readline().split() 

gg = g.readline().split()
if gg[0] != 'Peel':
    ff.extend(gg)
g.close()
nelec = 0 
for shell in ff:
    if not ('-' in shell):
        nelec += getFullOcc(shell)
#
adas_string,term,occ = processCSFStringIntoShells(vectors[0][0])
nelec += occ
ion_stage_index = charge - nelec 
atomic_number_index = charge - 1 
ion_pot = ip[ion_stage_index,atomic_number_index]
ion_term = tt[ion_stage_index,atomic_number_index]
#
print("&ADASEX NLEVS= {} NUMTMP=19 IRDTMP=1 ITCC=1 IBORN=-2 IRMPS=-1  IEL='{:2}' FIPOT={:11.1f} IONTRM='{:2}'/".format(
    numrequest,
    elements[charge-1],
    ion_pot,
    ion_term))
#
print("1.00+03 1.50+03 1.80+03 2.00+03 2.50+03 5.00+03 7.50+03 1.00+04 1.50+04 1.80+04 2.00+04 3.00+04 4.00+04 5.00+04 6.00+04 7.00+04 8.00+04 9.00+04 1.00+05")
#
old_occ = 0
for ii in range(0,numrequest):
    index_ii = sortargs[ii]
    csf = vectors[index_ii][0]
    adas_string,term,occ = processCSFStringIntoShells(csf)

    L,twoSP1 = proceessTermStringIntoInts(term)
    j   = processJStringIntoFloat(jvalues[index_ii])
    #print(energies_print[ii],s,L,twoSP1)
    
    print(adasin_format.format(
        ii+1,
        adas_string,
        twoSP1,
        L,
        j,
        energies[index_ii]  * RYDBERG_CM
    ))
    if (ii >= 1 and old_occ != occ):
        print('inconsistent occupation ',energies[index_ii]*0.5 + ground*0.5)
        sys.exit()
    old_occ = occ


    

print('NAME:')
print('DATE:')
print('                      ENTER DETAILS OF CALCULATION')
print('.')
