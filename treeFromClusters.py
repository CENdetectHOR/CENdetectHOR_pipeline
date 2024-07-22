from collections.abc import Iterable
import numpy as np
from Bio.Phylo.PhyloXML import Clade, Phylogeny, Sequence, Phyloxml
from Bio.SeqFeature import SeqFeature


def new_clade(
    label: str = None,
    branch_length: float = None,
    clades: list[Clade] = None
) -> Clade:
    if clades is not None:
        if branch_length is not None:
            for clade in clades:
                clade.branch_length += branch_length
        if len(clades) == 1:
            return clades[0]
    return Clade(name=label, clades=clades, branch_length=0)


def merge_clades(
    clades: Iterable[Clade],
    new_clusters_matrix: np.ndarray,
    branch_length:float = None
) -> list[Clade]:
    return [
        new_clade(
            clades=[
                clade
                for index, clade in enumerate(clades)
                if new_cluster_row[index] > 0
            ],
            branch_length=branch_length
        )
        for new_cluster_row in new_clusters_matrix
    ]


def new_leaves(leaf_labels: Iterable[str]) -> list[Clade]:
    return [new_clade(label=leaf_label) for leaf_label in leaf_labels]


def feature_to_leave(feature: SeqFeature, branch_length: int = 0) -> Clade:
    return Clade(
        name=feature.id,
        branch_length=branch_length,
        sequences=[Sequence(type='dna', location=feature.location)])


def features_to_leaves(features: Iterable[SeqFeature]) -> list[Clade]:
    return [feature_to_leave(feature) for feature in features]


def new_phylogeny(root_clade: Clade, name: str = 'monomers') -> Phylogeny:
    return Phylogeny(root=root_clade, name=name)


def new_phyloXML(phylogenies: list[Phylogeny], attributes: dict={}) -> Phyloxml:
    return Phyloxml(attributes, phylogenies=phylogenies)
