import pandas as pd
import src.Utils as utils
import src.MotifFinder as MotifFinder

window_size = 5

print("Reading data")
df = pd.read_csv("../data/labeled_seqs.csv")
cleaned_df = utils.clean_sequences(df)
seqs = utils.categorize_sequences(cleaned_df,
                                  categories=['f', 'e', 'm', 'b'],
                                  seq_column='seq_rep_rem',
                                  category_column='category')
event_map_df = pd.read_csv("../data/event_map.csv")
event_map = dict(zip(event_map_df['symbol'], event_map_df['event']))


# print()
motifs_dfs = []
for category, category_seqs in seqs.items():
    print(category)
    motif_finder = MotifFinder.CandidMotifFinder()
    motifs = motif_finder.motif_finder(category_seqs, motifs_count=10, window_size=window_size, dist_func='LCS')
    motifs['category'] = category
    motifs['window_size'] = window_size
    motifs_dfs.append(motifs)
concated_motif_df = pd.concat(motifs_dfs)
# concated_motif_df.to_csv("./motifs_window_sizes_%d.csv" % window_size)


window_motifs, window_categories = concated_motif_df['motif'].tolist(), concated_motif_df['category'].tolist()
window_dfs = []
for i in range(len(window_motifs)):
    res = utils.find_motifs_avg_dist(seqs, window_size, window_motifs[i], window_categories[i])
    if res is not None:
        window_dfs.append(res)
selected_motifs = pd.concat(window_dfs)
selected_motifs.reset_index(drop=True, inplace=True)
selected_motifs.to_csv("./contrast_motifs.csv")

