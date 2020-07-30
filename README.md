# Contrast Motif Discovery in Minecraft
## Abstract
Understanding event sequences is an important aspect of game analytics, since it is relevant to many player modeling questions. This paper introduces a method for analyzing event sequences by detecting contrasting motifs; the aim is to discover subsequences that are significantly more similar to one set of sequences vs. other sets. Compared to existing methods, our technique is scalable, capable of handling long event sequences. We applied our proposed sequence mining approach to analyze player behavior in Minecraft. Minecraft is a massively multiplayer online game that supports many forms of player collaboration. As a sandbox game, it provides players a large amount of flexibility in deciding how to complete tasks; this lack of goal-orientation makes the problem of analyzing Minecraft event sequences more challenging than event sequences from more structured games.  Using our approach, we were able to discover contrast motifs for many player actions, despite variability in how different players accomplished the same tasks. Furthermore, we explored how the level of player collaboration affects the contrast motifs. Although this paper focuses on applications within Minecraft, our tool, which we have made publicly available along with our dataset, can be used on any set of game event sequences.


_Samaneh Saadat and Gita Sukthankar, “Contrast Motif Discovery in Minecraft”, To appear in the Proceedings of the AAAI Conference on Artificial Intelligence and Interactive Digital Entertainment (AIIDE), Oct 2020_

## Data
All data used in our study is stored in the data folder of this repository and the description of the data files is as follows:
* `labeled_seqs.csv`: Minecraft action sequences extracted from the HeapCraft data set.
* `player_seqs.csv`: Player sequences extracted from the HeapCraft data set.
* `event_map.csv`: Symbols used to represent different Minecraft events.
* `collaboration_indicies.csv`: Collaboration index of players.

## Code
Our algorithm first finds the candidate motifs (`MotifFinder.py`) and then refine the motifs to select the motifs that are significantly more similar to the sequences of their own group compared to sequences of other groups (`MotifRefiner.py`).

Running the `main.py` script generate all the contrast motifs for Minecraft action sequences.


## Contact
If you have any question, don't hesitate to reach out to us at `ssaadat[AT]cs.ucf.edu`