import random
import sys

num = int(sys.argv[1])
size = int(sys.argv[2])

def generate_random(l):
    #A: 0-0.25 / C: 0.26-0.5 / G: 0.51-0.25 / T: 0.75-1
    seq = ""
    target = [0, 0.25, 0.5, 0.75, 1]
    for i in range(0, l):
        num = random.random()
        for i in range(len(target)-1):
            if num >= target[i] and num < target[i+1]:
                index = i
                break
        if index == 0:
            seq = seq + "A"
        elif index == 1:
            seq = seq + "C"
        elif index == 2:
            seq = seq + "G"
        else:
            seq = seq + "T"
    return seq

if __name__== "__main__":
    seqs = []
    for i in range(0, num):
        seqs.append(generate_random(size))
    output_seq = open("/Users/michm/Documents/BENG 182/Assignment 2/seq_pairs.txt","w")
    for each in seqs:
        output_seq.write(each + "\n")
