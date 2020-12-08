from HMM import HMM, HMM_Node
import math

def find_transmission(seq, in_islands, islands):
    trans = {'II': 0, 'IO': 0, 'OI': 0, 'OO': 0}
    vals = {'nII': 0, 'nIO': 0, 'nOI': 0, 'nOO': 0}

    vals['nIO'] = len(in_islands)
    vals['nOI'] = len(in_islands)
    vals['nII'] = len(islands) - len(in_islands)
    vals['nOO'] = len(seq) - 1 - vals['nIO'] - vals['nOI'] - vals['nII']

    trans['II'] = (float(vals['nII'])/float(vals['nII'] + vals['nIO']))
    trans['IO'] = (float(vals['nIO'])/float(vals['nIO'] + vals['nII']))
    trans['OI'] = (float(vals['nOI'])/float(vals['nOI'] + vals['nOO']))
    trans['OO'] = (float(vals['nOO'])/float(vals['nOO'] + vals['nOI']))

    return trans


if __name__ == "__main__":

    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 4/chrA.fasta") as f:
        sequence = [line.rstrip() for line in f]
    f.close()
    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 4/chrA.islands") as f:
        in_islands = [line.rstrip() for line in f]
    f.close()

    sequence = sequence[1:]
    seq = ''.join(sequence)

    islands = []
    for each in in_islands:
        curr = (tuple(each.split(' ')))
        for i in range(int(curr[0]), int(curr[1])+1):
            islands.append(i)

    trans = find_transmission(seq, in_islands, islands)

    # find freq of CG in CPG islands and out of CPG islands
    # for each in the tuples of segments find count of CG and do same for not CPG

    in_other_count = 0
    cg_c, ca_c, cc_c, ct_c = 0,0,0,0
    for each in in_islands:
        curr = (tuple(each.split(' ')))
        for i in range(int(curr[0]), int(curr[1]) + 1):
            dinuc = seq[i : i + 2]
            if dinuc == "CA":
                ca_c += 1
            elif dinuc == 'CC':
                cc_c += 1
            elif dinuc == 'CG':
                cg_c += 1
            elif dinuc == 'CT':
                ct_c += 1
    in_other_count = ca_c + cc_c + cg_c + ct_c
    in_freq = float(cg_c/in_other_count)

    out_other_count = 0
    cg_c, ca_c, cc_c, ct_c = 0,0,0,0
    starts = []
    starts.append(0)
    ends = []
    for each in in_islands:
        curr = (tuple(each.split(' ')))
        starts.append(curr[1])
        ends.append(curr[0])
    starts = starts[:-1]
    for i in range(len(starts)):
        for j in range(int(starts[i]), int(ends[i]) + 1):
            dinuc = seq[j : j + 2]
            if dinuc == "CA":
                ca_c += 1
            elif dinuc == 'CC':
                cc_c += 1
            elif dinuc == 'CG':
                cg_c += 1
            elif dinuc == 'CT':
                ct_c += 1
    out_other_count =  ca_c + cc_c + cg_c + ct_c
    out_freq = float(cg_c/out_other_count)

    # print (in_freq)
    # print (out_freq)

    # create transitions
    N_isl_N_reg = float(trans['IO']/4)
    N_reg_N_isl = float(trans['OI']/4)
    N_isl_N_isl = float((1-trans['IO'])/4)
    N_reg_N_reg = float((1-trans['OI'])/4)
    C_isl_G_isl = float(in_freq)
    C_isl_N_isl = float((1 - C_isl_G_isl - N_isl_N_reg)/3) 
    C_reg_G_reg = float(out_freq)
    C_reg_N_reg = float((1 - C_reg_G_reg - N_reg_N_isl)/3)

    # print (N_isl_N_reg, N_reg_N_isl, N_isl_N_isl, N_reg_N_reg, C_isl_G_isl, C_isl_N_isl, C_reg_G_reg, C_reg_N_reg)


    A_CPG = HMM_Node('A', True, {'A_CPG': N_isl_N_isl, 'C_CPG': N_isl_N_isl, 'G_CPG': N_isl_N_isl,'T_CPG': N_isl_N_isl, 'A_REG': N_isl_N_reg, 'C_REG': N_isl_N_reg, 'G_REG': N_isl_N_reg, 'T_REG': N_isl_N_reg})
    C_CPG = HMM_Node('C', True, {'A_CPG': C_isl_N_isl, 'C_CPG': C_isl_N_isl, 'G_CPG': C_isl_G_isl,'T_CPG': C_isl_N_isl, 'A_REG': N_isl_N_reg, 'C_REG': N_isl_N_reg, 'G_REG': N_isl_N_reg, 'T_REG': N_isl_N_reg})
    G_CPG = HMM_Node('G', True, {'A_CPG': N_isl_N_isl, 'C_CPG': N_isl_N_isl, 'G_CPG': N_isl_N_isl,'T_CPG': N_isl_N_isl, 'A_REG': N_isl_N_reg, 'C_REG': N_isl_N_reg, 'G_REG': N_isl_N_reg, 'T_REG': N_isl_N_reg})
    T_CPG = HMM_Node('T', True, {'A_CPG': N_isl_N_isl, 'C_CPG': N_isl_N_isl, 'G_CPG': N_isl_N_isl,'T_CPG': N_isl_N_isl, 'A_REG': N_isl_N_reg, 'C_REG': N_isl_N_reg, 'G_REG': N_isl_N_reg, 'T_REG': N_isl_N_reg})
    A_REG = HMM_Node('A', False, {'A_CPG': N_reg_N_isl, 'C_CPG': N_reg_N_isl, 'G_CPG': N_reg_N_isl,'T_CPG': N_reg_N_isl, 'A_REG': N_reg_N_reg, 'C_REG': N_reg_N_reg, 'G_REG': N_reg_N_reg, 'T_REG': N_reg_N_reg})
    C_REG = HMM_Node('C', False, {'A_CPG': N_reg_N_isl, 'C_CPG': N_reg_N_isl, 'G_CPG': N_reg_N_isl,'T_CPG': N_reg_N_isl, 'A_REG': C_reg_N_reg, 'C_REG': C_reg_N_reg, 'G_REG': C_reg_G_reg, 'T_REG': C_reg_N_reg})
    G_REG = HMM_Node('G', False, {'A_CPG': N_reg_N_isl, 'C_CPG': N_reg_N_isl, 'G_CPG': N_reg_N_isl,'T_CPG': N_reg_N_isl, 'A_REG': N_reg_N_reg, 'C_REG': N_reg_N_reg, 'G_REG': N_reg_N_reg, 'T_REG': N_reg_N_reg})
    T_REG = HMM_Node('T', False, {'A_CPG': N_reg_N_isl, 'C_CPG': N_reg_N_isl, 'G_CPG': N_reg_N_isl,'T_CPG': N_reg_N_isl, 'A_REG': N_reg_N_reg, 'C_REG': N_reg_N_reg, 'G_REG': N_reg_N_reg, 'T_REG': N_reg_N_reg})
    hmm = HMM([A_CPG, C_CPG, G_CPG, T_CPG, A_REG, C_REG, G_REG, T_REG])


    for each in hmm.model:
        print (each.succ)

    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 4/chrB.fasta") as f:
        sequence = [line.rstrip() for line in f]
    f.close()
    sequence = sequence[1:]
    seq = ''.join(sequence)

    # find CPG islands
    # make matrices to save scores and backtracking
    score = [[float('-inf') for x in range(len(seq)+1)] for y in range(8)] 
    back = [[0 for x in range(len(seq)+1)] for y in range(8)] 
     
    # initialize the backtrack matrix 
    for i in range(8):
        back[i][0] = 'START'

    indices = ['A_CPG', 'C_CPG', 'G_CPG', 'T_CPG', 'A_REG', 'C_REG', 'G_REG', 'T_REG']
    first = seq[0] + '_REG'
    score[indices.index(first)][0] = 0

    # fill in score and back matrix 
    for i in range(1, len(seq)+1):
        curr = seq[i-1]
        prev = []
        for k in range(8):
            prev.append(score[k][i-1])
        for j in range(8):
            em = float('-inf')
            if indices[j][0] == curr:
                em = 0
            score_max = []
            for l in range(len(prev)):
                score_max.append(prev[l] + math.log(hmm.model[l].succ.get(indices[j])) + em)
            score[j][i] = max(score_max)
            index_score = score_max.index(max(score_max))
            back[j][i] = index_score

    # backtrack by getting max of last column and finding where that came from on back and backtrack all the way and keep track of cpg or not 
    max_last = float('-inf')
    max_ind = -1
    for i in range(8): 
        if score[i][len(seq)-1] > max_last:
            max_last = score[i][len(seq)-1] 
            max_ind = i
    
    # get states 
    states = ''
    if max_ind <= 3: 
        states = states + ('I')
    else:
        states = states + ('O')
    next = back[max_ind][len(seq)]
    i = len(seq)-1
    while i > 0:
        if next <= 3:
            states = states + ('I')
        else:
            states = states + ('O')
        next = back[next][i]
        i -= 1
    states = states[::-1]

    # find positions of cpg islands 
    pos = []
    start = 0
    end = 0
    for i in range(len(states)-1): 
        if states[i] == 'O' and states[i+1] == 'I':
            start = i+1
        elif states[i] =='I' and states[i+1] == 'O':
            end = i
            pos.append((start,end))
        elif states[i] == 'I' and states[i+1] == 'I' and i == len(states)-2:
            end = i + 1
            pos.append((start,end))

    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 4/q5_results.txt", 'w') as f:
        for each in pos:
            f.write(str(each[0]) + " " + str(each[1]) + "\n")
    f.close()
    print (len(pos))
