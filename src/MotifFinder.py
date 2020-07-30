import pandas as pd
import numpy as np
from src.Distance import SequenceDistance


class CandidMotifFinder:
    def window_profile_min(self, seq, window, distance_function):
        window_size = len(window)
        msf = np.inf
        for i in range(len(seq) - window_size + 1):
            subseq = seq[i: i + window_size]
            if distance_function == 'hamming':
                dist = SequenceDistance.hamming_distance(window, subseq)
            elif distance_function == 'LCS':
                dist = SequenceDistance.LCS_distance(window, subseq)
            else:
                print("ERROR: INVALID DISTANCE FUNCTION")
                return None
            msf = min(msf, dist)
        return msf

    def distance_matrix(self, sequences, window_size, distance_function):
        seq_lens = [len(s) for s in sequences]
        seqs_count = len(sequences)
        windows_count = sum(seq_lens) - (seqs_count * (window_size - 1))
        print(windows_count, seqs_count, sum(seq_lens), (seqs_count * (window_size - 1)))
        print("Total number of windows is", windows_count)
        dist_mat = np.full((windows_count, seqs_count), np.inf)
        windows = []
        
        for m in range(seqs_count):
            for i in range(len(sequences[m]) - window_size + 1):
                w = sequences[m][i: i + window_size]
                windows.append(w)
                for j in range(seqs_count):
                    dist = self.window_profile_min(sequences[j], w, distance_function)
                    dist_mat[len(windows) - 1, j] = dist
        return dist_mat, windows

    def motif_finder(self, sequences, motifs_count, window_size, dist_func, dist_mat=None, windows=None):
        motifs, motifs_farc, motifs_loc = [], [], []
        motif_profiles = []
        sequences = [s for s in sequences if str(s) != 'nan' and len(s) >= window_size]
        
        if dist_mat is None:
            dist_mat, windows = self.distance_matrix(sequences, window_size, distance_function=dist_func)
        profile_len = len(sequences)
        Q = np.full(profile_len, np.inf)

        while True:
            min_area, idx, snpt = np.inf, -1, None
            for i in range(len(windows)):
                motif = windows[i]
                element_wise_min = np.minimum(np.array(dist_mat[i]), Q)
                profile_area = sum(element_wise_min)
                if profile_area < min_area:
                    min_area = profile_area
                    idx = i
                    snpt = motif
            Q = np.minimum(np.array(dist_mat[idx]), Q)
                 
            motifs.append(snpt)
            motifs_loc.append(idx)
            motif_profiles.append(dist_mat[idx])
            if len(motifs_loc) == motifs_count:
                break
      
        total_min = np.full((profile_len), np.inf)
        for i in range(motifs_count):
            total_min = np.minimum(total_min, motif_profiles[i])

        for i in range(motifs_count):
            f = 0
            for j in range(profile_len):
                if motif_profiles[i][j] == total_min[j]:
                    f += 1
            frac = f / profile_len
            motifs_farc.append(frac)
        res = dict(zip(motifs, motifs_farc))
        res_df = pd.DataFrame.from_dict(res, orient='index', columns=['Fractions'])
        return res_df
