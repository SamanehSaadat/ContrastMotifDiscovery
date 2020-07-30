import numpy as np

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
    return seqs


def translate_seq(event_map, s):
    events = []
    for c in s:
        events.append(event_map[c][:-5])
    return events