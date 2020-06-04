import pandas as pd
import src.Utils as utils
import src.MotifFinder as MotifFinder

window_size = 5

df = pd.read_csv("../data/labeled_seqs_scaled_move_subev.csv")
df = utils.clean_sequences(df)
seqs = utils.categorize_sequences(df,
                                  categories=['f', 'e', 'm', 'b'],
                                  seq_column='seq_rep_rem',
                                  category_column='action')

motifs_dfs = []
for action, action_seqs in seqs.items():
    print(action)
    motif_finder = MotifFinder.CandidMotifFinder()
    motifs = motif_finder.motif_finder(action_seqs, motifs_count=10, window_size=window_size, dist_func='LCS')
    motifs['action'] = action
    motifs['window_size'] = window_size
    motifs_dfs.append(motifs)
concated_motif_df = pd.concat(motifs_dfs)
concated_motif_df.to_csv("./motifs_window_sizes_%d.csv" % window_size)

    
w = 5
grp_snpts = motifs_df[motifs_df['window_size'] == w]
window_motifs, window_actions = grp_snpts['motif'].tolist(), grp_snpts['action'].tolist()
window_dfs = []
for i in range(len(window_motifs)):
    res = find_motifs_avg_dist(seqs, w, window_motifs[i], window_actions[i])
    if res is not None:
        window_dfs.append(res)
selected_motifs = pd.concat(window_dfs)
selected_motifs.reset_index(drop=True, inplace=True)
selected_motifs

event_map_df = pd.read_csv("../event_map_subev.csv")
event_map = dict(zip(event_map_df['symbol'], event_map_df['event']))
def translate_seq(s):
    events = []
    for c in s:
        events.append(event_map[c][:-5])
    return events