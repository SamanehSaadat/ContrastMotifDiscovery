import numpy as np


class SequenceDistance:
    def hamming_distance(self, s1, s2):
        n = len(s1)
        dist = 0
        for k in range(n):
            if s1[k] != s2[k]:
                dist += 1
        return dist

    def LCS_distance(self, s1, s2):
        l1 = len(s1)
        l2 = len(s2)
        dist = (l1 + l2) // 2 - self.LCS(s1, s2)
        return dist

    def LCS(self, s1, s2):
        s1 = "-" + s1
        s2 = "-" + s2
        l1 = len(s1)
        l2 = len(s2)
        d = np.zeros((l1, l2))
        for i in range(1, l1):
            for j in range(1, l2):
                if s1[i] == s2[j]:
                    d[i, j] = d[i - 1, j - 1] + 1
                else:
                    d[i, j] = max(d[i - 1, j], d[i, j - 1])
        return d[l1 - 1, l2 - 1]