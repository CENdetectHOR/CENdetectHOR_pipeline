import csv
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
import numpy as np
from Bio.Phylo.PhyloXML import Sequence


def BED_file_to_features(BED_filename):
    with open(BED_filename) as tsv:
        return [
            SeqFeature(
                SimpleLocation(
                    int(bedLine['chromStart']),
                    int(bedLine['chromEnd']),
                    strand=(1 if bedLine['strand'] == '+' else (-1 if bedLine['strand']
                            == '-' else None)) if 'strand' in bedLine else None,
                    ref=bedLine['chrom']
                ),
                id=bedLine['name'], type='repeat'
            )
            for bedLine in csv.DictReader(tsv, delimiter="\t", fieldnames=[
                'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
                'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'])
        ]


def feature_to_seq(
    feature: SeqFeature,
    references: dict[SeqRecord]
) -> SeqRecord:
    seq = feature.extract(None, references=references)
    seq.id = feature.id
    seq.description = ''
    return seq


def id_from_location(location: SimpleLocation) -> str:
    return (
        f'{location.ref}:{location.start}-{location.end}'
        + f'({"-" if location.strand == -1 else "+"})' if location.strand is not None else ''
    )


def location_to_feature(location: SimpleLocation) -> SeqFeature:
    return SeqFeature(location=location, id=id_from_location(location))


def location_to_seq(
    location: SimpleLocation,
    references: dict[SeqRecord]
) -> SeqRecord:
    return feature_to_seq(location_to_feature(location), references)


def label_to_location(seq_id: str) -> SimpleLocation:
    coordinates_parts = seq_id.split(':')
    seq_label = coordinates_parts[0]
    start_end = coordinates_parts[1].split('-')
    start = int(start_end[0])
    end = int(start_end[1])
    return SimpleLocation(start, end, strand=1, ref=seq_label)


def label_to_feature(seq_id: str) -> SeqFeature:
    return SeqFeature(location=label_to_location(seq_id), id=seq_id)


def label_to_phyloxml_sequence(seq_id: str) -> Sequence:
    return Sequence(name=seq_id, location=label_to_location(seq_id))


def extract_indices(indexed_list):
    return [indexed_item[0] for indexed_item in indexed_list]


def sorted_locations_indices(
    locations,
    independent_strands: bool = True
):
    indexed_locations = list(enumerate([location for location in locations]))
    sorted_indices = []
    for curr_ref in sorted(set([location.ref for location in locations])):
        curr_ref_indexed_locations = [
            indexed_location for indexed_location in indexed_locations if indexed_location[1].ref == curr_ref
        ]

        def indexed_locations_by_strand(strand):
            return [
                indexed_location
                for indexed_location in curr_ref_indexed_locations
                if indexed_location[1].strand == strand
            ]

        def get_start(indexed_location):
            return indexed_location[1].start
        
        if independent_strands:
            sorted_indices.extend(
                extract_indices(
                    sorted(
                        indexed_locations_by_strand(1),
                        key=get_start
                    )
                ) +
                extract_indices(
                    sorted(
                        indexed_locations_by_strand(-1),
                        key=get_start,
                        reverse=True
                    )
                )
            )
        else:
            sorted_indices.extend(
                extract_indices(
                    sorted(
                        curr_ref_indexed_locations,
                        key=get_start
                    )
                )

            )
    return sorted_indices


def order_by_indices(list_to_order, reordered_indeces):
    return [list_to_order[old_index] for old_index in reordered_indeces] if list_to_order is not None else None


def sorted_locations(locations):
    order_by_indices(locations, sorted_locations_indices(locations))


def sorted_features(features):
    order_by_indices(features, sorted_locations_indices(
        [feature.location for feature in features]))


def indeces_to_ordering_matrix(reordered_indeces):
    list_size = len(reordered_indeces)
    return np.array([[i == old_index for i in range(list_size)] for old_index in reordered_indeces]).T


def order_matrix_by_indeces(matrix_to_order, reordered_indeces, reorder_rows=True, reorder_cols=True):
    if matrix_to_order is None:
        return None
    return np.array([
        [
            matrix_to_order[
                reordered_indeces[row_index] if reorder_rows else row_index
            ][
                reordered_indeces[col_index] if reorder_cols else col_index
            ]
            for col_index in range(matrix_to_order.shape[1])
        ] for row_index in range(matrix_to_order.shape[0])]
    )

class SeqFeaturesByContiguity:
    sorted_seq_features: list[SeqFeature]
    gap_indices: list[int]
    max_allowed_gap : int
    reordered_indices: list[int]
    independent_strands: bool
    
    def __init__(
        self,
        seq_features: list[SeqFeature] = None,
        seq_locations: list[SimpleLocation] = None,
        seq_labels: list[str] = None,
        max_allowed_gap: int = 10,
        sorted_by_positive_strand_location: bool = False,
        sorted: bool = False,
        independent_strands: bool = True
    ) -> None:
        
        self.independent_strands = independent_strands
        
        if (
            seq_labels is not None and
            seq_features is None and seq_locations is None
        ):
            seq_features = [label_to_feature(label) for label in seq_labels]
        
        if seq_features is not None and seq_locations is None:
            seq_locations = [
                feature.location for feature in seq_features]
            
        self.max_allowed_gap = max_allowed_gap

        if seq_locations is None:
            raise Exception('Sequence location information is not available')

        if seq_features is None:
            seq_features = [
                location_to_feature(location)
                for location in seq_locations
        ]

        if not sorted and (independent_strands or not sorted_by_positive_strand_location):
            if sorted_by_positive_strand_location:
                indexed_locations = list(enumerate(self.seq_locations))
                positive_location_indices = [
                    indexed_location[0] for indexed_location in indexed_locations
                    if indexed_location[1].strand is None or indexed_location[1].strand == 1
                ]
                negative_location_indices = list(reversed([
                    indexed_location[0] for indexed_location in indexed_locations
                    if indexed_location[1].strand is not None and indexed_location[1].strand == -1
                ]))
                self.reordered_indices = positive_location_indices + negative_location_indices
            else:
                self.reordered_indices = sorted_locations_indices(
                    seq_locations,
                    independent_strands=independent_strands
                )
            self.sorted_seq_features = order_by_indices(
                seq_features,
                self.reordered_indices
            )
            
        sorted_seq_locations = [
            seq_feature.location
            for seq_feature in self.sorted_seq_features
        ]

        self.gap_indices = [
            i + 1
            for i in range(len(sorted_seq_locations) - 1)
            if sorted_seq_locations[i + 1].ref != sorted_seq_locations[i].ref
                or (
                    independent_strands and
                    sorted_seq_locations[i + 1].strand != sorted_seq_locations[i].strand
                )
                or sorted_seq_locations[i + 1].start - sorted_seq_locations[i].end > self.max_allowed_gap
        ]
        
def feature_len(feature: SeqFeature) -> int:
    return feature.location.end - feature.location.start   

def remove_overlapping_features(
    features: list[SeqFeature],
    max_allowed_overlap: int = None,
    max_allowed_overlap_fraction: float = None,
    expected_feature_size: int = None,
    keep_overlap_warning: bool = True,
    remove_overlap_warning: bool = True
) -> list[SeqFeature]:
    if (
        max_allowed_overlap is None and
        max_allowed_overlap_fraction is not None and
        expected_feature_size is not None
    ):
        max_allowed_overlap = expected_feature_size * max_allowed_overlap_fraction
    if max_allowed_overlap is None and not keep_overlap_warning:
        return
    filtered_features: list[SeqFeature] = []
    sorted_features = sorted(features, key=lambda ref: (ref.ref, ref.location.start))
    for curr_feature in sorted_features:
        prev_feature = filtered_features[-1] if len(filtered_features) > 0 else None
        if prev_feature is not None and curr_feature.ref == prev_feature.ref and curr_feature.location.start < prev_feature.location.end:
            if prev_feature.location.end - curr_feature.location.start > max_allowed_overlap:
                keep_previous = (
                    (
                        expected_feature_size is None and
                        feature_len(prev_feature) > feature_len(curr_feature)
                    ) or (
                        abs(feature_len(prev_feature) - expected_feature_size) <=
                        abs(feature_len(curr_feature) - expected_feature_size)
                    )
                )
                if not keep_previous:
                    filtered_features[-1] = curr_feature
                if remove_overlap_warning:
                    print(
                        f"Removed sequence feature {(curr_feature if keep_previous else prev_feature).location} " +
                        f"overlapping ({prev_feature.location.end - curr_feature.location.start} bps) " +
                        f"with sequence feature {(prev_feature if keep_previous else curr_feature).location}."
                    )
            else:
                filtered_features.append(curr_feature)
                if keep_overlap_warning:
                    print(
                        f"Kept sequence feature {curr_feature.location} " +
                        f"overlapping ({prev_feature.location.end - curr_feature.location.start} bps) " +
                        f"with sequence feature {prev_feature.location}."
                    )
        else:
            filtered_features.append(curr_feature)
    return filtered_features
