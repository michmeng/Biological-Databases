def dinucleotide_freq(seq, sig, chr):
    dinuc_counts = {'AA':0,'AC':0,'AG':0,'AT':0,'CA':0,'CC':0,'CG':0,'CT':0,'GA':0,'GC':0,'GG':0,'GT':0,'TA':0,'TC':0,'TG':0,'TT':0}
    nuc_counts = {'A':0,'C':0,'G':0,'T':0}

    for j in range(0, len(seq) - 1):
        dinuc = seq[j : j + 2]
        dinuc_counts[dinuc] = dinuc_counts[dinuc] + 1
    for i in range(0, len(seq)):
        nuc_counts[seq[i]] += 1

    for each in dinuc_counts:
        dinuc_counts[each] = float(dinuc_counts[each]/(len(seq) - 1))
    for each in nuc_counts:
        nuc_counts[each] = float(nuc_counts[each]/(len(seq)))

    if chr == False:
        nuc_counts = {'A':0.29,'C':0.21,'G':0.21,'T':0.29}

    diff = []
    for each in dinuc_counts: 
        # print (each, (nuc_counts[each[0]] * nuc_counts[each[1]]), dinuc_counts[each])
        if (abs((nuc_counts[each[0]] * nuc_counts[each[1]]) - dinuc_counts[each])) > sig:
            diff.append(each)
    return diff

if __name__ == "__main__":
    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 4/chrB.fasta") as f:
        sequence = [line.rstrip() for line in f]
    f.close()
    sequence = sequence[1:]
    seq = ''.join(sequence)
    sig = 0.01
    print ("Using the expected nucleotide frequencies for chrB: " + str(dinucleotide_freq(seq, sig, True)))
    print ("Using the expected nucleotide frequencies for the human genome" + str(dinucleotide_freq(seq, sig, False)))
