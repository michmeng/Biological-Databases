with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 4/q3_results.txt") as f:
   sequence = [line.rstrip() for line in f]
f.close()
with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/Assignment 4/q5_results.txt") as f:
   sequence_two = [line.rstrip() for line in f]
f.close()

count = 0
for each in sequence:
    curr = (tuple(each.split(' ')))
    for check in sequence_two:
        curr_check = (tuple(check.split(' ')))
        if curr_check[0] <= curr[0] and curr_check[1] >= curr[1]:
            count += 1

print (count)

