def find_peptide(parent, curr, peaks, masses, amino_acids, mass_aa, peptide):
    while abs(parent - curr) > 0.1:
        next = (None, None)
        trash = []
        for i in range(len(peaks)):
            b_ion = peaks[i][0] - curr
            b_poss = [amino_acids[each] for each in amino_acids if abs(amino_acids[each] - b_ion) < 0.1]
            y_ion = peaks[i][1] - curr
            y_poss = [amino_acids[each] for each in amino_acids if abs(amino_acids[each] - y_ion) < 0.1]
            if len(b_poss) > 1:
                b_poss = [b_poss[1]]
            if len(y_poss) > 1:
                y_poss = [y_poss[0]]
            if len(b_poss) == 1:
                next = (mass_aa[b_poss[0]], b_poss[0])
                trash.append(i)
            if len(y_poss) == 1:
                next = (mass_aa[y_poss[0]], y_poss[0])
                trash.append(i)
        peptide = peptide + next[0]
        curr += next[1]
        next = (None, None)
        for each in trash:
            peaks.pop(each)

    return peptide

if __name__ == '__main__':
    # read in peak data 
    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 5/unknown.mgf") as f:
        lines = [line.rstrip().split('\t') for line in f]
    lines = lines[2:]
    f.close()
    # manually looked through peak data to find parent 
    parent = (565.2968 * 2) - 18 - 2
    peaks = []
    for each in lines: 
        if float(each[0]) != 565.2968:
            peaks.append((float(each[0]) - 1, parent - float(each[0]) + 19))
    peaks.append((parent, 0))

    # read in amino acid data 
    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 5/aa_masses.csv") as f:
        lines = [line.rstrip().split(',') for line in f]
    lines = lines[3:]
    f.close()
    amino_acids = {}
    mass_aa = {}
    masses = []
    for each in lines:
        amino_acids[each[1]] = float(each[2])
        mass_aa[float(each[2])] = each[1]
        masses.append(each[2])

    curr_pep = 'L'
    curr = amino_acids['L']

    # read in database

    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 5/randomacid.txt") as f:
        text = f.read()
    result = text.split("\n\n")

    sequences = []
    for each in result: 
        lines = each.split('\n')
        sequences.append((lines[0], ''.join(lines[1:])))
    sequences = sequences[:-1]
    
    peptide = find_peptide(parent, curr, peaks, masses, amino_acids, mass_aa, curr_pep)

    # search database 
    for each in sequences: 
        if peptide in each[1]:
            print (peptide + " is in " + str(each[0][1:]))
