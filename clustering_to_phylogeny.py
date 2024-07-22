from dataclasses import dataclass
import json
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from treeFromClusters import feature_to_leave, new_phylogeny
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.PhyloXML import Clade, Phylogeny
from featureUtils import feature_to_seq

class SimplePhylogeny:
    num_leaves: int
    children: list[list[int]]
    root_clade_index: int
    num_clades: int
    
    def __init__(
        self,
        num_leaves: int,
        children: list[list[int]]
    ):
        self.num_leaves = num_leaves
        self.children = children
        self.num_clades = num_leaves + len(children)
        self.root_clade_index = self.num_clades - 1
    
    def get_clade_children(self, clade_index: int) -> list[int]:
        return (
            [] if clade_index < self.num_leaves
            else self.children[clade_index - self.num_leaves]
        )
    
    # def set_clade_children(self, clade_index: int, children: list[int]):
    #     return (
    #         [] if clade_index < self.num_leaves
    #         else self.children[clade_index - self.num_leaves]
    #     )
    
    def get_clade_height(self, clade_index: int) -> int:
        return (
            0 if clade_index < self.num_leaves
            else 1 + max([
                self.get_clade_height(subclade_index)
                for subclade_index in self.get_clade_children(clade_index)
            ])
        )
        
    def get_leaves_in_order(
        self, clade_index: int = None
    ) -> list[int]:
        if clade_index is None:
            clade_index = self.root_clade_index
        subclades = self.get_clade_children(clade_index)
        if len(subclades) == 0:
            return [clade_index]
        else:
            return [
                leave
                for subclade in subclades
                for leave in self.get_leaves_in_order(subclade)
            ]

class SimplePhylogenyWithDistances(SimplePhylogeny):
    max_distances: list[int]
    
    def __init__(
        self,
        num_leaves: int,
        children: list[list[int]],
        max_distances: list[int]
    ):
        super().__init__(num_leaves=num_leaves, children=children)
        if max_distances is None:
            max_distances = [
                self.get_clade_height(num_leaves + internal_clade_index)
                for internal_clade_index in enumerate(children)
            ]
        self.max_distances = max_distances
    
    def get_clade_distance(self, clade_index: int):
        return (
            0 if clade_index < self.num_leaves
            else (self.max_distances[clade_index - self.num_leaves])
        )
    
class SimplePhylogenyWithBranchLengths(SimplePhylogeny):
    branch_lengths: list[int] = None
    
    def __init__(
        self,
        num_leaves: int,
        children: list[list[int]],
        branch_lengths: list[int]
    ):
        super().__init__(num_leaves=num_leaves, children=children)
        self.branch_lengths = branch_lengths
    
    def get_branch_length(self, clade_index: int):
        return self.branch_lengths[clade_index]
    
def save_phylogeny_with_branch_lengths(
    phylogeny: SimplePhylogenyWithBranchLengths,
    filename: str
):
    with open(filename + '.json', 'w') as f:
        json.dump({
            'num_leaves': phylogeny.num_leaves,
            'children': phylogeny.children,
            'branch_lengths': phylogeny.branch_lengths
        }, f)
        
def save_phylogeny_with_distances(
    phylogeny: SimplePhylogenyWithDistances,
    filename: str
):
    with open(filename + '.json', 'w') as f:
        json.dump({
            'num_leaves': phylogeny.num_leaves,
            'children': phylogeny.children,
            'max_distances': phylogeny.max_distances
        }, f)
        
def load_phylogeny_with_branch_lengths(
    filename: str
) -> SimplePhylogenyWithBranchLengths:
    with open(filename + '.json') as f:
        dict = json.load(f)
        return SimplePhylogenyWithBranchLengths(
            num_leaves=dict['num_leaves'],
            children=dict['children'],
            branch_lengths=dict['branch_lengths']
        )
        
def load_phylogeny_with_distances(
    filename: str
) -> SimplePhylogenyWithDistances:
    with open(filename + '.json') as f:
        dict = json.load(f)
        return SimplePhylogenyWithDistances(
            num_leaves=dict['num_leaves'],
            children=dict['children'],
            max_distances=dict['max_distances']
        )
        
def compact_phylogeny(
    input_phylogeny: SimplePhylogenyWithDistances
) -> SimplePhylogenyWithDistances:
    def clades_to_descendants_at_distance(clade_indeces: list[int], distance: int) -> list[int]:
        return [
            descendant_clade
            for child_clade_index in clade_indeces
            for descendant_clade in clade_to_descendants_at_distance(child_clade_index, distance)
        ]

    def clade_to_descendants_at_distance(clade_index: int, distance: int) -> list[int]:
        if input_phylogeny.get_clade_distance(clade_index) < distance:
            return [clade_index]
        descendant_clades_at_distance = clades_to_descendants_at_distance(
            input_phylogeny.get_clade_children(clade_index), distance
        )
        if len(descendant_clades_at_distance) == 0:
            return [clade_index]
        return descendant_clades_at_distance
    
    new_internal_clades_distances = []
    new_internal_clades_children = []
    
    def compact_clade(clade_index: int) -> int:
        subclades = input_phylogeny.get_clade_children(clade_index)
        max_distance = input_phylogeny.get_clade_distance(clade_index)
        if len(subclades) == 0:
            return clade_index
        descendants = clades_to_descendants_at_distance(
            subclades,
            max_distance
        )
        new_internal_clades_children.append([compact_clade(clade) for clade in descendants])
        new_internal_clades_distances.append(max_distance)
        new_clade_index = input_phylogeny.num_leaves + len(new_internal_clades_distances) - 1
        return new_clade_index 

    new_root_clade_index = compact_clade(input_phylogeny.root_clade_index)
    
    return SimplePhylogenyWithDistances(
        num_leaves=input_phylogeny.num_leaves,
        children=new_internal_clades_children,
        max_distances=new_internal_clades_distances
    )
    
def distances_to_branch_lengths(
    input_phylogeny: SimplePhylogenyWithDistances
) -> SimplePhylogenyWithBranchLengths:
    branch_lengths = [0] * (input_phylogeny.num_clades)
    def set_clade_branch_length(clade_index: int):
        for subclade_index in input_phylogeny.get_clade_children(clade_index):
            branch_lengths[subclade_index] = (
                input_phylogeny.get_clade_distance(clade_index) -
                input_phylogeny.get_clade_distance(subclade_index)
            ) / 2
            set_clade_branch_length(subclade_index)
    
    set_clade_branch_length(input_phylogeny.root_clade_index)
    
    return SimplePhylogenyWithBranchLengths(
        num_leaves=input_phylogeny.num_leaves,
        children=input_phylogeny.children,
        branch_lengths=branch_lengths
    )
    

def build_phylogeny(
    simple_phylogeny: SimplePhylogenyWithBranchLengths,
    items_as_seq_features: list[SeqFeature]
) -> Phylogeny:
    
    def build_clade(clade_index: int):
        subclade_indices = simple_phylogeny.get_clade_children(clade_index)
        branch_length = simple_phylogeny.get_branch_length(clade_index)
        if len(subclade_indices) == 0:
            return feature_to_leave(items_as_seq_features[clade_index], branch_length=branch_length)
        subclades = [build_clade(subclade_index) for subclade_index in subclade_indices]
        return Clade(clades=subclades, branch_length=branch_length)
        
    root_clade = build_clade(simple_phylogeny.root_clade_index)
    phylogeny = new_phylogeny(root_clade)

    return phylogeny

@dataclass
class CladeSortResult:
    subclade_index: int
    min_index: int    

def sort_by_leaf_indexes(phylogeny: SimplePhylogeny):
    
    def sort_clade_by_leaf_indices(clade_index: int):
        if len(phylogeny.get_clade_children(clade_index)) == 0:
            return clade_index
        subclade_results = [
            CladeSortResult(
                subclade_index,
                sort_clade_by_leaf_indices(subclade_index)
            )
            for subclade_index in phylogeny.get_clade_children(clade_index)
        ]
        subclade_results.sort(
            key=lambda c: c.min_index
        )
        phylogeny.children[clade_index - phylogeny.num_leaves] = [
            subclade_result.subclade_index
            for subclade_result in subclade_results
        ]
        return subclade_results[0].min_index
    
    sort_clade_by_leaf_indices(phylogeny.root_clade_index)
        
def get_clades_by_level(
    phylogeny: SimplePhylogeny
) -> tuple[list[list[int]],list[list[int]]]:
    
    clades_by_level = []
    children_num_by_level = []

    def set_by_height(clade_height: int, clade_index: int, children_num: int):
        if clade_height >= len(clades_by_level):
            missing_levels = clade_height - len(clades_by_level) + 1
            clades_by_level.extend([[]] * missing_levels)
            children_num_by_level.extend([[]] * missing_levels)
        clades_by_level[clade_height].append(clade_index)
        children_num_by_level[clade_height].append(children_num)

    def extract_clades_by_height(clade_index: int):
        subclade_indices = phylogeny.get_clade_children(clade_index)
        subclade_heights = [
            phylogeny.get_clade_height(subclade_index)
            for subclade_index in subclade_indices
        ]
        clade_height = (
            0 if len(subclade_indices) == 0
            else 1 + max(subclade_heights)
        )
        for subclade_position, subclade_height in enumerate(subclade_heights):
            subclade_index = subclade_indices[subclade_position]
            extract_clades_by_height(subclade_index)
            for height in range(subclade_height + 1, clade_height):
                set_by_height(clade_height=height, clade_index=subclade_indices[subclade_position], children_num=1)            
        set_by_height(clade_height=clade_height, clade_index=clade_index, children_num=len(subclade_indices))

    extract_clades_by_height(phylogeny.root_clade_index)
    return clades_by_level, children_num_by_level

@dataclass
class ClusteringToPhylogenyResult:
    phylogeny: Phylogeny
    item_position_to_leaf_index: list[int]

def clustering_to_phylogeny(
    clustering: AgglomerativeClustering = None,
    items_vs_feature_array = None,
    items_as_seq_records: list[SeqRecord] = None,
    items_as_seq_features: list[SeqFeature] = None,
    seq_references: dict[SeqRecord] = None,
    dist_matrix: np.ndarray = None,
    compute_distances: bool = True,
    linkage='single',
    metric='euclidean',
    sort: bool = True
) -> ClusteringToPhylogenyResult:
    
    if clustering is None:
        clustering = AgglomerativeClustering(
            metric=metric if dist_matrix is None else 'precomputed',
            compute_full_tree=True,
            linkage=linkage,
            compute_distances = compute_distances,
            n_clusters=1)


    if (items_as_seq_features is not None and
        seq_references is not None and
        items_as_seq_records is None
    ):
        items_as_seq_records = [
            feature_to_seq(seq_feature, references=seq_references)
            for seq_feature in items_as_seq_features
        ]
    
    if not hasattr(clustering, 'children_'):
        
        if dist_matrix is not None:
            clustering.fit(dist_matrix)
        else:
            if (items_as_seq_records is not None and 
                items_vs_feature_array is None
            ):
                items_vs_feature_array = [
                    str(seq.seq) for seq in items_as_seq_records
                ]
            if items_vs_feature_array is not None:
                clustering.fit(items_vs_feature_array)
            else:
                raise Exception('No data available to perform clustering')
    
    aggregation_result = SimplePhylogenyWithDistances(
        num_leaves=int(dist_matrix.shape[0]),
        children=[[int(subclade) for subclade in subclades] for subclades in clustering.children_],
        max_distances=[int(distance) for distance in clustering.distances_] if hasattr(clustering, 'distances_') else None
    )

    compacted_result = compact_phylogeny(aggregation_result)
    simple_phylogeny = distances_to_branch_lengths(compacted_result)
    
    if sort:
        sort_by_leaf_indexes(simple_phylogeny)
        
    phylogeny = build_phylogeny(
        simple_phylogeny,
        items_as_seq_features = items_as_seq_features
    )
    
    return ClusteringToPhylogenyResult(
        phylogeny=phylogeny,
        item_position_to_leaf_index=simple_phylogeny.get_leaves_in_order()
    )

