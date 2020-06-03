import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math


df = pd.read_csv("./labeled_seqs_scaled_move_subev.csv")
df = df[df['action'] != 'o']
df['seq_no_move'] = df['seq'].apply(lambda x: x.replace('C', ''))
df['seq_no_move'] = df['seq_no_move'].apply(lambda x: x.replace('d', ''))
df['seq_no_move'] = df['seq_no_move'].apply(lambda x: x.replace('j', ''))
df['seq_no_move'] = df['seq_no_move'].apply(lambda x: x.replace('D', ''))
df['seq_no_move'] = df['seq_no_move'].apply(lambda x: x.replace('E', ''))

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

df['seq_rep_rem'] = df['seq_no_move'].apply(remove_repeats)

seqs = {}
seqs['f'] = df[df['action'] == 'f']['seq_rep_rem'].tolist()
seqs['e'] = df[df['action'] == 'e']['seq_rep_rem'].tolist()
seqs['m'] = df[df['action'] == 'm']['seq_rep_rem'].tolist()
seqs['b'] = df[df['action'] == 'b']['seq_rep_rem'].tolist()

snippets_dfs = []
for w in range(3, 10):
    for action, action_seqs in seqs.items():
        print(w, action)
        mat_prof = MatrixProfile()
        snippets = mat_prof.snippet_finder(action_seqs, snippets_count=10, window_size=w, dist_func='LCS')
        snippets['action'] = action
        snippets['window_size'] = w
        print(snippets.head())
        snippets_dfs.append(snippets)
concated_snippet_df = pd.concat(snippets_dfs)

concated_snippet_df.to_csv("snippets_various_window_sizes_final_w-3-10_new.csv")


from scipy import stats
def find_snippets_avg_dist(seqs, window_size, snippet, snippet_action):
    dist_dict = {}
    avg_dist = {}
    for seq_action, seq_list in seqs.items():
        seq_list = [s for s in seq_list if len(s) >= window_size]
        dists = []
        for seq in seq_list:
            dists.append(window_profile_min(seq, snippet))
        dist_dict[seq_action] = dists
        avg_dist[seq_action] = np.mean(dists)

    significant = True
#     snippet_action_seq_list = seqs[snippet_action]
    pvalue_threshold = 0.05
    for seq_action, seq_list in seqs.items():
        if seq_action != snippet_action:
#             print(snippet_action_seq_list)
#             print(seq_list)
            statistics, pvalue = stats.ttest_ind(dist_dict[snippet_action], dist_dict[seq_action])
            if pvalue > pvalue_threshold:
                significant = False
                break
    
    has_min_mean = False
    temp = min(avg_dist.values()) 
    min_mean = [key for key in avg_dist if avg_dist[key] == temp]
    if len(min_mean) == 1 and min_mean[0] == snippet_action:
        has_min_mean = True
    
    if significant == True and has_min_mean == True:
        df = pd.DataFrame.from_dict(avg_dist, orient='index').T
        df['snippet'] = snippet
        df['Snippet action'] = snippet_action
        return df
    return None
    
    
w = 5
grp_snpts = snippets_df[snippets_df['window_size'] == w]
window_snippets, window_actions = grp_snpts['Snippet'].tolist(), grp_snpts['action'].tolist()
window_dfs = []
for i in range(len(window_snippets)):
    res = find_snippets_avg_dist(seqs, w, window_snippets[i], window_actions[i])
    if res is not None:
        window_dfs.append(res)
selected_snippets = pd.concat(window_dfs)
selected_snippets.reset_index(drop=True, inplace=True)
selected_snippets

event_map_df = pd.read_csv("../event_map_subev.csv")
event_map = dict(zip(event_map_df['symbol'], event_map_df['event']))
def translate_seq(s):
    events = []
    for c in s:
        events.append(event_map[c][:-5])
    return events