import csv
from Bio.SeqFeature import SeqFeature, FeatureLocation 
import numpy as np
from Bio.Phylo.PhyloXML import Sequence

def BED_file_to_features(BED_filename):
    with open(BED_filename) as tsv:
        return [
            SeqFeature(
                FeatureLocation(
                    int(bedLine['chromStart']),
                    int(bedLine['chromEnd']),
                    strand=(1 if bedLine['strand'] == '+' else (-1 if bedLine['strand'] == '-' else None)) if 'strand' in bedLine else None,
                    ref=bedLine['chrom']
                ),
                id=bedLine['name'], type='repeat'
            )
            for bedLine in csv.DictReader(tsv, delimiter="\t", fieldnames=[
            'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
            'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'])
        ]


def feature_to_seq(feature, references):
    seq = feature.extract(None, references=references)
    seq.id = feature.id
    return seq 

def id_from_location(location):
    return (
        f'{location.ref}:{location.start}-{location.end}'
        + f'({"-" if location.strand == -1 else "+"})' if location.strand is not None else ''
    )

def location_to_feature(location):
    return SeqFeature(location=location, id=id_from_location(location))

def location_to_seq(location, references):
    return feature_to_seq(location_to_feature(location), references)

def label_to_location(seq_id):
    coordinates_parts = seq_id.split(':')
    seq_label = coordinates_parts[0]
    start_end = coordinates_parts[1].split('-')
    start = int(start_end[0])
    end = int(start_end[1])
    return FeatureLocation(start, end, strand=1, ref=seq_label)

def label_to_feature(seq_id):
    return SeqFeature(location=label_to_location(seq_id), id=seq_id)

def label_to_phyloxml_sequence(seq_id):
    return Sequence(name=seq_id, location=label_to_location(seq_id))

def extract_indices(indexed_list):
    return [indexed_item[0] for indexed_item in indexed_list]

def sorted_locations_indeces(locations):
    indexed_locations = list(enumerate([location for location in locations]))
    sorted_indeces = []
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
        sorted_indeces.extend(
            extract_indices(sorted(indexed_locations_by_strand(1), key=get_start))
            + extract_indices(sorted(indexed_locations_by_strand(-1), key=get_start, reverse=True))
        )
    return sorted_indeces

def order_by_indices(list_to_order, reordered_indeces):
    return [list_to_order[old_index] for old_index in reordered_indeces] if list_to_order is not None else None

def sorted_locations(locations):
    order_by_indices(locations, sorted_locations_indeces(locations))

def sorted_features(features):
    order_by_indices(features, sorted_locations_indeces([feature.location for feature in features]))

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


