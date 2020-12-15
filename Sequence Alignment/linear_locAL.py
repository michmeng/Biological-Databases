import argparse
import sys
parser = argparse.ArgumentParser(description='hi michelle is cool')

OUTDIR = "/Users/michm/Documents/BENG 182/Assignment 2/"
parser.add_argument('infile', action="store")
parser.add_argument('-m', action="store", dest = "match", required = True, type = int)
parser.add_argument('-s', action="store", dest="mismatch", required = True, type = int)
parser.add_argument('-d', action="store", dest="indel", required = True, type=int)
parser.add_argument('-a', action="store_true", dest = "a")
args = parser.parse_args()

in_file = OUTDIR + "/" + str(args.infile)

def linear_local(str1, str2, match, mismatch, indel, type):
    max_score = -1000000000
    max_position = (0,0)

    # create and initialize window
    window = [[0 for i in range(2)] for j in range(len(str2) + 1)]
    j = 1
    for i in range(len(str1)+1):
        if i == 0:
            window[i][1] = 0
        else:
            check1 = window[i-1][1] - indel
            check2 = window[i][0] - indel
            dec = -1 * mismatch
            if str1[i-1] == str2[j-1]:
                dec = match
            check3 = window[i-1][0] + dec
            window[i][1] = max(0, check1, check2, check3)
            if type == "end":
                if window[i][1] >= max_score:
                    max_score = window[i][1]
                    max_position = (i,1)
            else:
                if window[i][1] > max_score:
                    max_score = window[i][1]
                    max_position = (i,1)

    i = 0
    #run through matrix moving window one column over
    j = 2
    while i < len(str1) - 1:
        #move window over one
        for k in range(len(str2)+1):
            window[k][0] = window[k][1]

        #fill in next column
        for k in range(len(str2)+1):
            if k == 0:
                window[k][1] = 0
            else:
                check1 = window[k-1][1] - indel
                check2 = window[k][0] - indel
                dec = -1 * mismatch
                if str1[j-1] == str2[k-1]:
                    dec = match
                check3 = window[k-1][0] + dec
                window[k][1] = max(0, check1, check2, check3)
                if window[k][1] > max_score:
                    max_score = window[k][1]
                    max_position = (k,j)
        j = j + 1
        i = i + 1

    return max_score, max_position

def backtracker(str1,str2, match, mismatch, indel):
    backtrack = [[0 for i in range(len(str2) + 1)] for j in range(len(str1) + 1)]
    #initialize matrix and 0 row and column
    score = [[0 for i in range(len(str2) + 1)] for j in range(len(str1) + 1)]
    max_score = 0
    max_location = []

    for i in range(1, len(str1)+1):
        for j in range(1, len(str2)+1):
            check1 = score[i-1][j] - indel
            check2 = score[i][j-1] - indel
            dec = (-1) * mismatch
            if str1[i-1] == str2[j-1]:
                dec = match
            check3 = score[i-1][j-1] + dec
            score[i][j] = max(0, check1, check2, check3)

            if score[i][j] >= max_score:
                max_score = score[i][j]
                max_location = []
                max_location.append(i)
                max_location.append(j)
            if score[i][j] == score[i-1][j] - indel:
                backtrack[i][j] = "D"
            elif score[i][j] == score[i][j-1] - indel:
                backtrack[i][j] = "R"
            elif score[i][j] == score[i-1][j-1] + dec:
                backtrack[i][j] = "DIAG"

    return backtrack, max_score, max_location

def local_allignment(str1,str2,match,mismatch,indel):
    l_allign = []
    l_allign_two = []
    if len(str1) > len(str2):
        l_string = str1
        s_string = str2
    else:
        l_string = str2
        s_string = str1

    l_string = str1
    s_string = str2

    backtrack, max_score, max_location = backtracker(l_string,s_string,match,mismatch,indel)
    i = max_location[0]
    j = max_location[1]
    allignment = ""

    check = False
    while i > 0 and j > 0:
        if backtrack[i][j] == "DIAG":
            l_allign.append(l_string[i-1])
            l_allign_two.append(s_string[j-1])
            i = i-1
            j = j-1
        elif backtrack[i][j] == "D":
            l_allign.append(l_string[i-1])
            l_allign_two.append("-")
            i = i-1
        elif backtrack[i][j] == "R" :
            l_allign.append("-")
            l_allign_two.append(s_string[j-1])
            j = j-1
        elif backtrack[i][j] == 0:
            break


    l_allign = l_allign[::-1]
    l_allign_two = l_allign_two[::-1]


    allign_one = "".join(l_allign)
    allign_two = "".join(l_allign_two)

    if len(str1) > len(str2):
        str1 = allign_one
        str2 = allign_two
    else:
        str2 = allign_one
        str1 = allign_two

    return max_score, str1, str2

if __name__== "__main__":
    with open(in_file) as f:
        lines = [line.strip() for line in f]
    in_one = lines[1]
    in_two = lines[3]
    in_one_inv = in_one[::-1]
    in_two_inv = in_two[::-1]
    score, end_position =  linear_local(in_one, in_two,args.match,args.mismatch,args.indel, "end")
    score, s = linear_local(in_one_inv, in_two_inv,args.match,args.mismatch,args.indel, "start")
    start_position = (len(in_two)-s[0], len(in_one)-s[1])
    in_one = in_one[start_position[1]:end_position[1]]
    in_two = in_two[start_position[0]:end_position[0]]
    score, str1, allignment = local_allignment(in_one, in_two, args.match, args.mismatch, args.indel)

    output = open("/Users/michm/Documents/BENG 182/Assignment 2/linear_locAL_output.txt", "w")
    output.write("Score: " + str(score) + "\n")
    output.write("Length: " + str(len(str1)) + "\n")
    output.write(str1 + "\n")
    output.write(allignment)
    output.close()

    print ("Score: " + str(score))
    print ("Length: " + str(len(str1)))
    if args.a == True:
        print (str1)
        print (allignment)
