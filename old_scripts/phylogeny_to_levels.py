from dataclasses import dataclass
import re
from Bio.Phylo.BaseTree import Tree, Clade

from featureUtils import SeqFeaturesByContiguity, label_to_feature, label_to_location, location_to_feature
from Bio.SeqFeature import SimpleLocation

@dataclass
class CladesByLevelResult:
    clades_by_level: list[list[Clade]]
    children_num_by_level: list[list[int]]
    
class CladeHeights:
    __clade_heights_dict: dict[Clade, int]
    tree_height: int
    
    def calc_clade_height(self, clade: Clade) -> int:
        clade_height = (
            0 if len(clade.clades) == 0
            else 1 + max([
                self.calc_clade_height(subclade)
                for subclade in clade.clades
            ])
        )
        self.__clade_heights_dict[clade] = clade_height
        return clade_height
    
    def __init__(self, phylogeny: Tree) -> None:
        self.__clade_heights_dict = {}
        self.tree_height = self.calc_clade_height(phylogeny.root)

    def get_clade_height(self, clade: Clade) -> int:
        return self.__clade_heights_dict[clade]

def get_clades_by_level(
    phylogeny: Tree
) -> CladesByLevelResult:
    
    clades_by_level: list[list[Clade]] = []
    children_num_by_level: list[list[int]] = []

    # def get_clade_height(clade: Clade) -> int:
    #     return (
    #         0 if len(clade.clades) == 0
    #         else 1 + max([
    #             get_clade_height(subclade)
    #             for subclade in clade.clades
    #         ])
    #     )
    
    clade_heights = CladeHeights(phylogeny=phylogeny)

    def set_by_height(clade_height: int, clade: Clade, children_num: int):
        if clade_height >= len(clades_by_level):
            missing_levels = clade_height - len(clades_by_level) + 1
            clades_by_level.extend([[]] * missing_levels)
            children_num_by_level.extend([[]] * missing_levels)
        clades_by_level[clade_height].append(clade)
        children_num_by_level[clade_height].append(children_num)

    def extract_clades_by_height(clade: Clade):
        subclades = clade.clades
        subclade_heights = [
            clade_heights.get_clade_height(subclade)
            for subclade in subclades
        ]
        clade_height = (
            0 if len(subclades) == 0
            else 1 + max(subclade_heights)
        )
        for subclade_position, subclade_height in enumerate(subclade_heights):
            subclade = subclades[subclade_position]
            extract_clades_by_height(subclade)
            for height in range(subclade_height + 1, clade_height):
                set_by_height(clade_height=height, clade=subclades[subclade_position], children_num=1)            
        set_by_height(clade_height=clade_height, clade=clade, children_num=len(subclades))

    extract_clades_by_height(phylogeny.root)
    return CladesByLevelResult(
        clades_by_level=clades_by_level,
        children_num_by_level=children_num_by_level[1:]
    )

def get_clade_contraction_by_level(
    children_num_by_level: list[list[int]]
) -> list[list[int]]:
    contraction_by_level = []
    for children_num_by_clade in children_num_by_level:
        contraction = []
        for clade_index, children_num in enumerate(children_num_by_clade):
            contraction.extend([clade_index] * children_num)
        contraction_by_level.append(contraction)
    return contraction_by_level
        
def get_labelled_items_by_level(
    clade_contraction_by_level: list[list[int]],
    item_position_to_leaf_index: list[int]
) -> list[list[int]]:
    labelled_items_by_level = [item_position_to_leaf_index]
    for level, contraction in enumerate(clade_contraction_by_level):
        labelled_items_by_level.append([
            contraction[label]
            for label in labelled_items_by_level[level]
        ])
    return labelled_items_by_level

location_pattern = re.compile("(.*)\[([0-9]+):([0-9]+)\]\((\+|-)\)")

def parse_location(location_str: str) -> SimpleLocation:
    match = location_pattern.match(location_str)
    return SimpleLocation(
        ref=match.group(1),
        start=int(match.group(2)),
        end=int(match.group(3)),
        strand=int(match.group(4) + '1'))

def extract_features_from_leaves(phylogeny: Tree):
    return [
        location_to_feature(
            parse_location(leave.sequences[0].location)
            if isinstance(leave.sequences[0].location, str)
            else leave.sequences[0].location
        )
        if (
            hasattr(leave, 'sequences') and
            len(leave.sequences) == 1 and
            hasattr(leave.sequences[0], 'location')
        )
        else label_to_feature(leave.name)
        for leave in phylogeny.get_terminals()
    ]

@dataclass
class PhylogenyToLevelsResult:
    clades_by_level: list[list[Clade]]
    clade_contraction_by_level: list[list[int]]
    labelled_items_by_level: list[list[int]]
    
def phylogeny_to_levels(
    phylogeny: Tree,
    item_position_to_leaf_index: list[int] = None,
    max_allowed_gap: int = 10
) -> PhylogenyToLevelsResult:
    if item_position_to_leaf_index is None:
        item_position_to_leaf_index = SeqFeaturesByContiguity(
            seq_features=extract_features_from_leaves(phylogeny),
            max_allowed_gap=max_allowed_gap
        ).reordered_indices
    clades_by_level_res = get_clades_by_level(phylogeny)
    clade_contraction_by_level = get_clade_contraction_by_level(
        clades_by_level_res.children_num_by_level
    )
    labelled_items_by_level = get_labelled_items_by_level(
        clade_contraction_by_level=clade_contraction_by_level,
        item_position_to_leaf_index=item_position_to_leaf_index
    )
    return PhylogenyToLevelsResult(
        clades_by_level=clades_by_level_res.clades_by_level,
        clade_contraction_by_level=clade_contraction_by_level,
        labelled_items_by_level=labelled_items_by_level
    )