from trie import TrieNode
from trie import Trie

if __name__ == "__main__":
    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/queries.txt") as f:
        q1 = [line.rstrip() for line in f]
    f.close()
    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/queries2.txt") as f:
        q2 = [line.rstrip() for line in f]
    f.close()
    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/DNA.txt") as f:
        seq = [line.rstrip() for line in f]
    seq = seq[1:]
    sequence = ''.join(seq)
    f.close()
    # query = q2

    trie1 = Trie()
    count = 1 
    for q in q1:
        count = trie1.insert(q, count)
    
    trie2 = Trie()
    count = 1
    for q in q2:
        count = trie2.insert(q, count)

    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/results_query1.txt", 'w') as f:
        # print(trie.pattern_match(sequence, f))
        trie1.pattern_match(sequence, f)
    f.close()

    with open("/Users/michmeng/Documents/UCSD 2020/BENG 182/results_query2.txt", 'w') as f:
        # print(trie.pattern_match(sequence, f))
        trie2.pattern_match(sequence, f)
    f.close()

