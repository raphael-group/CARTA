### Quantitative fate maps specifications

'fate_map_panel_metadata.csv' lists the fate map categories, number of terminal types and Colless index

Fate maps are saved in .json format, and have the following attributes:

1. **'node_id'**: progenitor states. Progenitor states are represented by positive integers.
2. **'root-id'**: root state. Root state is one of the progenitor states represented by the largest positive integer.
2. **'tip_id'**: terminal types. Terminal types are represented by negative integers.
3. **'merge'**: list of daughter state(s)/type(s) for each progenitor state.
4. **'diff_time'**: time of progenitor state commitments. Terminal types have infinite time of commitment.
5. **'diff_mode_probs'**: commitment bias of progenitor states. Values corresponds to three modes of commitments: symmetric to downstream state X (first element in 'merge'), symmetric to downstream states Y (second element in 'merge') and assymetric commitment.
6. **'lifetime'**: the doubling time of cells of each state/type.
7. **'founder_size'**: the number of cell(s) of the root state at time zero.
8. **'target_time'**: time of sample collection.
9. **'prob_loss'**: cell death probabilities for each state/type.
10. **'prob_pm'**: non-doubling proabilities for each state/type.
11. **'edges'**: edge list. This can be derived from the information above.
12. **'root_time'**: time until first commitment event of the root state.
