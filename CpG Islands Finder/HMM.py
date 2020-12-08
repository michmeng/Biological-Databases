import math

class HMM_Node():

    def __init__(self, nuc, cpg, succ):
        self.succ = succ
        self.nuc = nuc
        self.cpg = cpg

class HMM():

    def __init__(self, model):
        self.model = model

if __name__ == "__main__":
    # read in file 
    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 4/chrB.fasta") as f:
        sequence = [line.rstrip() for line in f]
    f.close()
    sequence = sequence[1:]
    seq = ''.join(sequence)

    # calculate transition probabilities for number 3
    N_isl_N_reg = float(0.005/4)
    N_reg_N_isl = float(0.001/4)
    N_isl_N_isl = float(0.995/4)
    N_reg_N_reg = float(0.999/4)
    C_isl_G_isl = float(0.995 * (0.06/(0.06 + (3 * (0.94/15)))))
    C_isl_N_isl = float((1 - C_isl_G_isl - N_isl_N_reg)/3)
    C_reg_G_reg = float(0.999 * (0.01/(0.01 + (3 * (0.99/15)))))
    C_reg_N_reg = float((1 - C_reg_G_reg - N_reg_N_isl)/3)

    print (N_isl_N_reg, N_reg_N_isl, N_isl_N_isl, N_reg_N_reg, C_isl_G_isl, C_isl_N_isl, C_reg_G_reg, C_reg_N_reg)


    # create 8 states(nodes) and add them to the model 
    A_CPG = HMM_Node('A', True, {'A_CPG':N_isl_N_isl,'C_CPG':N_isl_N_isl,'G_CPG':N_isl_N_isl,'T_CPG':N_isl_N_isl,'A_REG':N_isl_N_reg,'C_REG':N_isl_N_reg,'G_REG':N_isl_N_reg,'T_REG':N_isl_N_reg})
    C_CPG = HMM_Node('C', True, {'A_CPG':C_isl_N_isl,'C_CPG':C_isl_N_isl,'G_CPG':C_isl_G_isl,'T_CPG':C_isl_N_isl,'A_REG':N_isl_N_reg,'C_REG':N_isl_N_reg,'G_REG':N_isl_N_reg,'T_REG':N_isl_N_reg})
    G_CPG = HMM_Node('G', True, {'A_CPG':N_isl_N_isl,'C_CPG':N_isl_N_isl,'G_CPG':N_isl_N_isl,'T_CPG':N_isl_N_isl,'A_REG':N_isl_N_reg,'C_REG':N_isl_N_reg,'G_REG':N_isl_N_reg,'T_REG':N_isl_N_reg})
    T_CPG = HMM_Node('T', True, {'A_CPG':N_isl_N_isl,'C_CPG':N_isl_N_isl,'G_CPG':N_isl_N_isl,'T_CPG':N_isl_N_isl,'A_REG':N_isl_N_reg,'C_REG':N_isl_N_reg,'G_REG':N_isl_N_reg,'T_REG':N_isl_N_reg})
    A_REG = HMM_Node('A', False, {'A_CPG':N_reg_N_isl,'C_CPG':N_reg_N_isl,'G_CPG':N_reg_N_isl,'T_CPG':N_reg_N_isl,'A_REG':N_reg_N_reg,'C_REG':N_reg_N_reg,'G_REG':N_reg_N_reg,'T_REG':N_reg_N_reg})
    C_REG = HMM_Node('C', False, {'A_CPG':N_reg_N_isl,'C_CPG':N_reg_N_isl,'G_CPG':N_reg_N_isl,'T_CPG':N_reg_N_isl,'A_REG':C_reg_N_reg,'C_REG':C_reg_N_reg,'G_REG':C_reg_G_reg,'T_REG':C_reg_N_reg})
    G_REG = HMM_Node('G', False, {'A_CPG':N_reg_N_isl,'C_CPG':N_reg_N_isl,'G_CPG':N_reg_N_isl,'T_CPG':N_reg_N_isl,'A_REG':N_reg_N_reg,'C_REG':N_reg_N_reg,'G_REG':N_reg_N_reg,'T_REG':N_reg_N_reg})
    T_REG = HMM_Node('T', False, {'A_CPG':N_reg_N_isl,'C_CPG':N_reg_N_isl,'G_CPG':N_reg_N_isl,'T_CPG':N_reg_N_isl,'A_REG':N_reg_N_reg,'C_REG':N_reg_N_reg,'G_REG':N_reg_N_reg,'T_REG':N_reg_N_reg})
    hmm = HMM([A_CPG,C_CPG,G_CPG,T_CPG,A_REG,C_REG,G_REG,T_REG])


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

    #fill in score and back matrix 
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

    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 4/q3_results.txt", 'w') as f:
        for each in pos:
            f.write(str(each[0]) + " " + str(each[1]) + "\n")
    f.close()
    print (len(pos))