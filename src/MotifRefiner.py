import pandas as pd
import numpy as np
from scipy import stats
from src.Distance import SequenceDistance


class ContrastMotifRefiner:
    def window_profile_min(self, seq, window, distance_function='LCS'):
        window_size = len(window)
        msf, dist = np.inf, np.inf
        for i in range(len(seq) - window_size + 1):
            subseq = seq[i: i + window_size]
            if distance_function == 'LCS':
                dist = SequenceDistance.LCS_distance(window, subseq)
            else:
                print("ERROR: INVALID DISTANCE FUNCTION")
            msf = min(msf, dist)
        return msf

    def refine_motifs(self, seqs, window_size, motif, motif_category):
        dist_dict = {}
        avg_dist = {}
        for seq_category, seq_list in seqs.items():
            seq_list = [s for s in seq_list if len(s) >= window_size]
            dists = []
            for seq in seq_list:
                dists.append(self.window_profile_min(seq, motif))
            dist_dict[seq_category] = dists
            avg_dist[seq_category] = np.mean(dists)

        significant = True
        pvalue_threshold = 0.05
        for seq_category, seq_list in seqs.items():
            if seq_category != motif_category:
                statistics, pvalue = stats.mannwhitneyu(dist_dict[motif_category], dist_dict[seq_category])
                if pvalue > pvalue_threshold:
                    significant = False
                    break

        has_min_mean = False
        temp = min(avg_dist.values())
        min_mean = [key for key in avg_dist if avg_dist[key] == temp]
        if len(min_mean) == 1 and min_mean[0] == motif_category:
            has_min_mean = True

        if significant and has_min_mean:
            df = pd.DataFrame.from_dict(avg_dist, orient='index').T
            df['motif'] = motif
            df['motif category'] = motif_category
            return df
        return None