if __name__ == '__main__':
    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 5/iso_profile.txt") as f:
        og = [line.rstrip().split('\t') for line in f]
    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 5/iso_profile_60.txt") as f:
        mod = [line.rstrip().split('\t') for line in f]

    original = []
    for each in og: 
        original.append((each[0].split(' ')[0], each[0].split(' ')[1]))

    modified = []
    for each in mod: 
        modified.append((each[0].split(' ')[0], each[0].split(' ')[1]))

    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 5/compare_iso_profile.txt", 'w') as f:
        for i in range(len(modified)):
            if i < len(original):
                f.write(str(modified[i][0]) + ": " + str(float(modified[i][1]) - float(original[i][1])) + "\n")
            else:
                f.write(str(modified[i][0]) + ": " + str(modified[i][1]) + "\n")
    f.close()