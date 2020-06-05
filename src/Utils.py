import pandas as pd
import numpy as np
import numpy as np
from scipy import stats


def remove_repeats(s):
    if s is np.nan:
        return np.nan
    if len(s) == 0:
        return np.nan
    res = s[0]
    for c in s[1:]:
        if res[-1] == c:
            continue
        else:
            res += c
    return res


def clean_sequences(df):
    df = df[df['category'] != 'o']
    df['seq_no_move'] = df['seq'].apply(lambda x: x.replace('C', ''))
    df['seq_no_move'] = df['seq_no_move'].apply(lambda x: x.replace('d', ''))
    df['seq_no_move'] = df['seq_no_move'].apply(lambda x: x.replace('j', ''))
    df['seq_no_move'] = df['seq_no_move'].apply(lambda x: x.replace('D', ''))
    df['seq_no_move'] = df['seq_no_move'].apply(lambda x: x.replace('E', ''))

    df['seq_rep_rem'] = df['seq_no_move'].apply(remove_repeats)
    return df


def categorize_sequences(df, categories, seq_column, category_column):
    seqs = {}
    for c in categories:
        seqs[c] = df[df[category_column] == c][seq_column].tolist()
    return  seqs


def LCS( s1, s2):
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


def LCS_distance(s1, s2):
    l1 = len(s1)
    l2 = len(s2)
    dist = (l1 + l2) // 2 - LCS(s1, s2)
    return dist


def window_profile_min(seq, window, distance_function='LCS'):
    window_size = len(window)
    msf, dist = np.inf, np.inf
    for i in range(len(seq) - window_size + 1):
        subseq = seq[i: i + window_size]
        if distance_function == 'LCS':
            dist = LCS_distance(window, subseq)
        else:
            print("ERROR: INVALID DISTANCE FUNCTION")
        msf = min(msf, dist)
    return msf


def find_motifs_avg_dist(seqs, window_size, motif, motif_category):
    dist_dict = {}
    avg_dist = {}
    for seq_category, seq_list in seqs.items():
        seq_list = [s for s in seq_list if len(s) >= window_size]
        dists = []
        for seq in seq_list:
            dists.append(window_profile_min(seq, motif))
        dist_dict[seq_category] = dists
        avg_dist[seq_category] = np.mean(dists)

    significant = True
    pvalue_threshold = 0.05
    for seq_category, seq_list in seqs.items():
        if seq_category != motif_category:
            #             print(motif_category_seq_list)
            #             print(seq_list)
            statistics, pvalue = stats.ttest_ind(dist_dict[motif_category], dist_dict[seq_category])
            if pvalue > pvalue_threshold:
                significant = False
                break

    has_min_mean = False
    temp = min(avg_dist.values())
    min_mean = [key for key in avg_dist if avg_dist[key] == temp]
    if len(min_mean) == 1 and min_mean[0] == motif_category:
        has_min_mean = True

    if significant == True and has_min_mean == True:
        df = pd.DataFrame.from_dict(avg_dist, orient='index').T
        df['motif'] = motif
        df['motif category'] = motif_category
        return df
    return None


def translate_seq(event_map, s):
    events = []
    for c in s:
        events.append(event_map[c][:-5])
    return events