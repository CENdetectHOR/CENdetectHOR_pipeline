from Bio.SeqFeature import SeqFeature
from featureUtils import SeqFeaturesByContiguity
from loops import LoopInSeq, find_loops

def find_inversion_loops(
    seq_features: list[SeqFeature],
    max_allowed_gap: int = 10,
    min_loop_size: int = 2,
    max_loop_size: int = 30,
    min_loops: int = 5,
    allowed_mismatch_rate: float = 0.0,
    allow_hor_overlap: bool = False
) -> list[LoopInSeq]:
    sfbc = SeqFeaturesByContiguity(
            seq_features=seq_features,
            max_allowed_gap=max_allowed_gap,
            independent_strands=False
    )
    split_limits = [0] + sfbc.gap_indices + [len(seq_features)]
    num_splits = len(sfbc.gap_indices) + 1
    return find_loops(
        whole_seq=[(1 - seq_feature.strand) // 2 for seq_feature in sfbc.sorted_seq_features],
        seq_spans=[
            (split_limits[split_index], split_limits[split_index + 1])
            for split_index in range(num_splits)
        ],
        min_loop_size=min_loop_size,
        max_loop_size=max_loop_size,
        min_loops=min_loops,
        allowed_mismatch_rate=allowed_mismatch_rate,
        allow_overlap=allow_hor_overlap
    )    
