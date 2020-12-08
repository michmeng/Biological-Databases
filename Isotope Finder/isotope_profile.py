from itertools import product
from itertools import  combinations
import math

def atom_composition(peptide):
    comp = {'C':0, 'H':0, 'N':0, 'O':0, 'S':0}
    r = {'A': [1,3,0,0],'V':[3,7,0,0],'L':[4,9,0,0],'G':[0,1,0,0],'I':[4,9,0,0],'M':[3,7,0,0,1],'W':[9,8,1,0],'F':[7,7,0,0],'P':[3,5,0,0],
        'S':[1,3,0,1],'C':[1,3,0,0,1],'N':[2,4,1,1],'Q':[3,6,1,1],'T':[2,5,0,1],'Y':[7,7,0,1],
        'D':[2,3,0,2],'E':[3,5,0,2],'K':[4,10,1,0],'R':[4,10,3,0],'H':[4,5,2,0]}

    for each in peptide:
        comp['C'] += (r[each][0] + 2)
        comp['H'] += (r[each][1] + 4)
        comp['N'] += (r[each][2] + 1)
        comp['O'] += (r[each][3] + 2)
        if len(r[each]) == 5:
            comp['S'] += r[each][4]
    comp['H'] -= ((len(peptide)-1) * 2)
    comp['O'] -= (len(peptide)-1)
    return comp

def find_prob(total, iso, prob):
    nCr = math.factorial(total) / math.factorial(iso) / math.factorial(total-iso)
    return nCr * (prob**iso) * ((1-prob)**(total-iso))

def isotope_profile(composition, p_C13, p_N15, p_O18):
    profile = []

    high = composition['C'] + composition['N'] + (2 * composition['O'])
    
    combinations = [p for p in product(list(range(high + 1)), repeat=3)]

    for each in range(high + 1):
        prob = 0.0
        for combo in combinations:
            if combo[0] + combo[1] + (2 * combo[2]) != each or combo[0] > composition['C'] or combo[1] > composition['N'] or combo[2] > composition['O']:   
                continue
            prob += find_prob(composition['C'], combo[0], p_C13) * find_prob(composition['N'], combo[1], p_N15) * find_prob(composition['O'], combo[2], p_O18)

        profile.append(('P' + str(each),prob))

    return profile

def isotope_profile_2(composition, p_C13, p_N15, p_O18, p_2H):
    profile = []

    # what is high 
    high = composition['C'] + composition['N'] + (2 * composition['O']) + composition['H']

    for each in range(high + 1):
        if each % 25 == 0:
            print ('HELLO')
        prob = 0.0
        for C in range(each + 1):
            for N in range(each + 1):
                for O in range(each + 1):
                    for H in range(each + 1):
                        if C + N + (2 * O) + H != each or C > composition['C'] or N > composition['N'] or O > composition['O'] or H > composition['H']:   
                            continue
                        prob += find_prob(composition['C'], C, p_C13) * find_prob(composition['N'], N, p_N15) * find_prob(composition['O'], O, p_O18) * find_prob(composition['H'], H, p_2H)
        profile.append(('P' + str(each),prob))

    return profile

if __name__ == "__main__":
    peptide = 'SLAMMER'

    # find atomic composition
    comp = atom_composition(peptide)

    # abundance probabilities 
    p_C13 = 0.0111
    p_N15 = 0.0037
    p_O18 = 0.0020
    p_2H = 0.6000

    # find isotopic composition 
    iso_prof_1 = isotope_profile(comp, p_C13, p_N15, p_O18)
    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 5/iso_profile.txt", 'w') as f:
        for each in iso_prof_1:
            f.write(str(each[0]) + ": " + str(each[1]) + "\n")
    f.close()

    # for 60% of H are deurterium
    iso_prof_2 = isotope_profile_2(comp, p_C13, p_N15, p_O18, p_2H)
    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 5/iso_profile_60.txt", 'w') as f:
        for each in iso_prof_2:
            f.write(str(each[0]) + ": " + str(each[1]) + "\n")
    f.close()

 
    