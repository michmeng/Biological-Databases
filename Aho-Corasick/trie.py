class TrieNode():

    def __init__(self, data, count):
        self.children = {}
        self.data = data
        self.count = count
        self.terminate = False
        self.matches = 0

class Trie():

    def __init__(self):
        self.root = TrieNode("*", 0)

    def insert(self, pattern, count):
        # count = 1
        curr = self.root
        for each in pattern:
            if each not in curr.children:
                curr.children[each] = TrieNode(each, count)
                count = count + 1
                #print (str(curr.count) + '->' + str(curr.children[each].count) + ":" + str(curr.children[each].data))
            curr = curr.children[each]
        curr.terminate = True
        return count

    def search(self, pattern):
        curr = self.root
        while len(pattern) > 0:
            each = pattern[0]
            if each not in curr.children:
                return False
            curr = curr.children[each]
            pattern = pattern[1:]
            if len(pattern) == 0: 
                if curr.terminate == True:
                    return True 
                else:
                    return False
        return True

    def pattern_match(self, sequence, file): 
        matches = {}
        l = 0
        c = 0
        v = self.root
        # n = 1
        while c < len(sequence):
            # if c > 0:
                # if c == (len(sequence)//10) * n:
                #     print (str(c) + "/" + str(len(sequence)))
                #     for each in matches:
                #         file.write(str(each) + " " + str(matches[each]) + " ; ")
                #     n += 1
                #     file.write("\n")
            if sequence[c] in v.children: 
                v = v.children[sequence[c]]
                c += 1
                if v.terminate == True:
                    if sequence[l:c] not in matches:
                        matches[sequence[l:c]] = 1
                    else:
                        matches[sequence[l:c]] = matches[sequence[l:c]] + 1
            else:
                c = l + 1
                l = c
                v = self.root
        for each in matches:
            file.write(str(each) + " " + str(matches[each]) + "\n")
        return matches 