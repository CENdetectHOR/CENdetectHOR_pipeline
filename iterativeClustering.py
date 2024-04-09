from cluster import build_string_distance_matrix, get_seq_as_txt, merge_clusters, min_distance
from featureUtils import feature_to_seq, label_to_phyloxml_sequence, location_to_feature, order_by_indices, order_matrix_by_indeces, sorted_locations_indeces
from hor import hor_tree_to_phylogeny, loops_to_HORs, name_hor_tree
from loops import find_loops
import numpy as np
from Bio.Phylo.PhyloXML import Phyloxml, Other

from treeFromClusters import merge_clades, new_leaves, new_phylogeny, features_to_leaves

class ClusteredSeq:
    def __init__(self, clusters_expansion, loops = [], clades = None, gap_indeces = [], seq_locations=None):
        self.clades = clades
        self.clusters_expansion = clusters_expansion
        self.loops = []
        self.add_loops(loops)
        whole_seq_as_clusters = list(np.arange(len(clusters_expansion)) @ clusters_expansion)
        seq_split_indeces = [0] + gap_indeces + [len(whole_seq_as_clusters)]
        self.seqs_as_clusters = [whole_seq_as_clusters[seq_split_indeces[i]:seq_split_indeces[i+1]] for i in range(len(gap_indeces) + 1)]
        # self.seq_as_clusters = list(np.arange(len(clusters_expansion)) @ clusters_expansion)

    def add_loops(self, loops, seq_locations=None):
        self.loops.extend(loops)
        if self.clades is not None:
            self.hors = loops_to_HORs(self.loops, self.clades, seq_locations=seq_locations)

    def __str__(self):
        return (
            f'Num clusters: {len(self.clusters_expansion)}, ' +
            f'Seqs: {[get_seq_as_txt(seq_as_clusters) for seq_as_clusters in self.seqs_as_clusters]}, ' +
            f'Loops: {[str(loop) for loop in self.loops]}'
        )
    
    def to_dict(self):
        return {
            "num_clusters": len(self.clusters_expansion),
            "seqs": [get_seq_as_txt(seq_as_clusters) for seq_as_clusters in self.seqs_as_clusters],
            "loops": [{
                "loop_seq": get_seq_as_txt(l.loop.loop_seq),
                "spans": [{
                    "span_start": span_inseq.span_start,
                    "span_length": span_inseq.span_length,
                    "num_of_laps": span_inseq.num_of_laps,
                    "in_loop_start": span_inseq.in_loop_start
                } for span_inseq in l.spans_in_seq]
            } for l in self.loops]
        }        

def matrix_sparsity(matrix):
    return 1.0 - np.count_nonzero(matrix) / matrix.size

# Compares two collections of spans in a sequence, called "coverages"
# Each coverage implicitly represents the set of items of the sequence included in some span
# Returns true iff the set of items associated with coverage_a contains the one associated with coverage_b
def coverage_includes(coverage_a, coverage_b):
    return all([
        any([
            span_a.span_start <= span_b.span_start and span_a.span_start + span_a.span_length >= span_b.span_start + span_b.span_length
            for span_a in coverage_a
        ]) for span_b in coverage_b
    ])

def coverage_diff(coverage_a, loops_b):
    return [
        loop_b
        for loop_b in loops_b
        if any([
            all([
                span_a.span_start > span_b.span_start or span_a.span_start + span_a.span_length <span_b.span_start + span_b.span_length
                for span_a in coverage_a
            ])
            for span_b in loop_b.spans_in_seq
        ])
    ]

def clusterings_with_hors(
        seqs=None,
        seqs_as_features=None,
        seq_locations=None,
        references=None,
        sorted=False,
        sorted_by_positive_strand_location=False,
        gap_indeces=[],
        max_allowed_bases_gap_in_hor = 10,
        distance_matrix=None,
        seq_labels_prefix='',
        starting_distance=0,
        min_num_clusters=2, max_num_clusters=None, order_clusters=False,
        require_loop=True, min_len_loop=2, max_len_loop=30, min_loop_reps=3,
        require_increasing_loop_coverage=True,
        require_relevant_loop_coverage=True,
        incremental_loops=False,
        closure_sparsity_threshold = 0.97,
        build_tree = True):
    print(f'Start of clusterings_with_hors')

    if seqs_as_features is not None:
        seq_locations = [feature.location for feature in seqs_as_features]

    if seq_locations is not None:

        if not sorted:
            print(f'Reorder')
            if sorted_by_positive_strand_location:
                print(f'Reorder negative strand')
                indexed_locations = list(enumerate(seq_locations))
                positive_location_indexes = [
                    indexed_location[0] for indexed_location in indexed_locations
                    if indexed_location[1].strand is None or indexed_location[1].strand == 1
                ]
                negative_location_indexes = list(reversed([
                    indexed_location[0] for indexed_location in indexed_locations
                    if indexed_location[1].strand is not None and indexed_location[1].strand == -1
                ]))
                reordered_indeces = positive_location_indexes + negative_location_indexes
            else:
                print(f'Reorder all')
                reordered_indeces = sorted_locations_indeces(seq_locations)
            print(f'Indexes built, now reordering lists...')
            seq_locations = order_by_indices(seq_locations, reordered_indeces)
            seqs_as_features = order_by_indices(seqs_as_features, reordered_indeces)
            seqs = order_by_indices(seqs, reordered_indeces)
            print(f'Lists reordered, now reordering matrix...')
            distance_matrix = order_matrix_by_indeces(distance_matrix, reordered_indeces)
            print(f'Matrix reordered')

        if seqs_as_features is None:
            seqs_as_features = [location_to_feature(location) for location in seq_locations]

        if seqs is None:
            seqs = [feature_to_seq(feature, references) for feature in seqs_as_features]

        gap_indeces = [
            monomer_index + 1
            for monomer_index in range(len(seq_locations) - 1)
            if seq_locations[monomer_index + 1].ref != seq_locations[monomer_index].ref
                or seq_locations[monomer_index + 1].strand != seq_locations[monomer_index].strand
                or seq_locations[monomer_index + 1].start - seq_locations[monomer_index].end > max_allowed_bases_gap_in_hor
        ]

    if distance_matrix is None:
        plain_seqs = [str(seq.seq) for seq in seqs]
        print(f'Computing distances...')
        distance_matrix = build_string_distance_matrix(plain_seqs)
        print(f'Distance matrix computed!')

    if max_num_clusters is None:
        max_num_clusters = distance_matrix.shape[0] / 4

    if require_relevant_loop_coverage:
        require_increasing_loop_coverage=True

    if build_tree:
        if seqs_as_features is not None:
            curr_clades = features_to_leaves(seqs_as_features)
        else:
            seq_labels = [seq_labels_prefix + str(i) for i in range(distance_matrix.shape[0])]
            curr_clades = new_leaves(seq_labels)
    else:
        curr_clades = None

    clusterings = []
    curr_distance_matrix = distance_matrix
    curr_clusters_expansion = None
    curr_clusters_max_internal_distance = 0
    last_loops_coverage = []
    while curr_clusters_expansion is None or len(curr_clusters_expansion) >= min_num_clusters:
        curr_distance_matrix, curr_clusters_expansion, merged_clusters_expansion, merged_clusters_distance = merge_clusters(
            curr_distance_matrix,
            clusters_expansion=curr_clusters_expansion,
            sparsity_threshold=closure_sparsity_threshold,
            max_distance=max(starting_distance, min_distance(curr_distance_matrix)))
        if build_tree:
            curr_clades = merge_clades(
                clades=curr_clades, new_clusters_matrix=merged_clusters_expansion,
                branch_length=(merged_clusters_distance-curr_clusters_max_internal_distance)/2)
        if len(curr_clusters_expansion) <= max_num_clusters:
            if order_clusters:
                clusters_size = np.sum(curr_clusters_expansion, axis=1)
                ordered_clusters_expansion = curr_clusters_expansion[clusters_size.argsort()[::-1]]
                clustered_seq = ClusteredSeq(ordered_clusters_expansion, clades=curr_clades, gap_indeces=gap_indeces)
            else:
                clustered_seq = ClusteredSeq(curr_clusters_expansion, clades=curr_clades, gap_indeces=gap_indeces)
            print(f'Looking for loops in {str(clustered_seq)}...')
            loops = find_loops(
                clustered_seq.seqs_as_clusters,
                min_loop_size=min_len_loop, max_loop_size=max_len_loop,
                min_loops=min_loop_reps
            )
            print(f'Loops found: {len(loops)}')
            if incremental_loops:
                clustered_seq.add_loops(coverage_diff(last_loops_coverage, loops), seq_locations=seq_locations)
            else:
                clustered_seq.add_loops(loops, seq_locations=seq_locations)


            if len(loops) > 0 or not require_loop:
                loops_coverage = [span_in_seq for loop in loops for span_in_seq in loop.spans_in_seq]
                if not require_increasing_loop_coverage or not coverage_includes(last_loops_coverage, loops_coverage):
                    clusterings.append(clustered_seq)
                    last_loops_coverage = loops_coverage

        curr_clusters_max_internal_distance = merged_clusters_distance

    if require_relevant_loop_coverage:
        new_clusterings_reversed = [clusterings[-1]]
        curr_preserved_loops = clusterings[-1].loops
        curr_hor_tree_nodes = clusterings[-1].hors
        for clustering_level in reversed(range(len(clusterings) - 1)):
            new_preserved_loops = []
            new_hor_tree_nodes = []
            new_loops = []
            clustered_seq = clusterings[clustering_level]
            for existing_loop_index in range(len(curr_preserved_loops)):
                existing_loop = curr_preserved_loops[existing_loop_index]
                existing_hor = curr_hor_tree_nodes[existing_loop_index]
                corresponding_candidate_loops = []
                corresponding_candidate_hors = []
                specificity_increased = False
                for loop_index in range(len(clustered_seq.loops)):
                    loop = clustered_seq.loops[loop_index]
                    hor = clustered_seq.hors[loop_index]
                    if coverage_includes(existing_loop.spans_in_seq, loop.spans_in_seq):
                        corresponding_candidate_loops.append(loop)
                        corresponding_candidate_hors.append(hor)
                        if (len(loop.loop.loop_seq) > len(existing_loop.loop.loop_seq) or
                            len(set(loop.loop.loop_seq)) > len(set(existing_loop.loop.loop_seq))):
                            specificity_increased = True
                if len(corresponding_candidate_loops) > 1 or specificity_increased:
                    new_preserved_loops.extend(corresponding_candidate_loops)
                    new_hor_tree_nodes.extend(corresponding_candidate_hors)
                    new_loops.extend(corresponding_candidate_loops)
                    existing_hor.sub_hors = corresponding_candidate_hors
                    for sub_hor in corresponding_candidate_hors:
                        sub_hor.super_hor = existing_hor
                        sub_hor.sub_hors = []
                else:
                    new_preserved_loops.append(existing_loop)
                    new_hor_tree_nodes.append(existing_hor)
            clustered_seq.loops = new_loops
            if len(new_loops) > 0 or not require_loop:
                new_clusterings_reversed.append(clustered_seq)
            curr_preserved_loops = new_preserved_loops
            curr_hor_tree_nodes = new_hor_tree_nodes

        hor_tree_roots = new_clusterings_reversed[0].hors

        clusterings = list(reversed(new_clusterings_reversed))

    if build_tree:
        seq_tree_root = curr_clades[0]
        hor_tree_root = hor_tree_roots[0]
        name_hor_tree(hor_tree_root)
        reference_seqs_element = (
            Other(
                tag='reference-sequences',
                children=[label_to_phyloxml_sequence(ref_id) for ref_id in references.keys()]
            )
            if references is not None else None
        )
        return (
            Phyloxml(
                phylogenies=[new_phylogeny(seq_tree_root), hor_tree_to_phylogeny(hor_tree_root)],
                attributes={'xsd':'http://www.w3.org/2001/XMLSchema'},
                other=reference_seqs_element
            ),
            hor_tree_root, clusterings
        )
    else:
        return clusterings

def bfs_merge_rec(hor_tree_level):
    if not any([hor_in_level.sub_hors for hor_in_level in hor_tree_level]):
        return [hor_tree_level]
    sub_hors = [sub_hor for hor_in_level in hor_tree_level for sub_hor in (hor_in_level.sub_hors or [hor_in_level])]
    return [sub_hors] + bfs_merge_rec(sub_hors)

def bfs_merge(hor_tree_root):
    return bfs_merge_rec([hor_tree_root])

def bfs(hor_in_seq, tree):
    sub_tree = tree.common_ancestor(hor_in_seq.hor.clade_seq)
    hor_in_seq.sub_hors