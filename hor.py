from typing import List
from loops import LoopSpanInSeq
from Bio.Phylo import BaseTree
from Bio.Phylo.PhyloXML import Clade, Phylogeny, Sequence, Phyloxml, Property
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, Location

class HOR:
    clade_seq: List[BaseTree.Clade]

    def __init__(self, clade_seq):
        self.clade_seq = clade_seq

    def __str__(self):
        return ''.join([clade.id if clade.id is not None else '*' for clade in self.clade_seq])

class HORInSeq:
    hor: HOR
    locations: List[Location]
    spans_in_seq: List[LoopSpanInSeq]
    super_hor: any
    sub_hors: List[any]

    def __init__(self, hor, spans_in_seq = [], locations = None):
        self.hor = hor
        self.spans_in_seq = spans_in_seq
        self.locations = locations

    def add_span(self, span_in_seq):
        self.spans_in_seq.append(span_in_seq)

    def __str__(self):
        return (
            f'{self.loop}' +
            (
                f' in {",".join([str(span) for span in self.spans_in_seq])}'
                    if len(self.spans_in_seq) > 0 else ''
            )
        )

def seq_span_to_location(span, seq_locations):
    start_location = seq_locations[span.span_start]
    end_location = seq_locations[span.span_start + span.span_length - 1]
    ref = start_location.ref
    strand = start_location.strand
    return FeatureLocation(
        ref=ref,
        strand=strand,
        start=start_location.start if strand is None or strand == 1 else end_location.start,
        end=end_location.end if strand is None or strand == 1 else start_location.end
    )

def seq_spans_to_compound_location(spans_in_seq, seq_locations):
    return CompoundLocation([seq_span_to_location(span, seq_locations) for span in spans_in_seq])

def loop_to_HOR(loop_in_seq, clades, seq_locations=None):
    hor = HOR([clades[clade_index] for clade_index in loop_in_seq.loop.loop_seq])
    return HORInSeq(
        hor,
        spans_in_seq=loop_in_seq.spans_in_seq,
        locations=(
            [seq_span_to_location(span, seq_locations) for span in loop_in_seq.spans_in_seq]
            if seq_locations is not None else None
        )
    )

def loops_to_HORs(loops_in_seq, clades, seq_locations=None):
    return [loop_to_HOR(loop_in_seq, clades, seq_locations=seq_locations) for loop_in_seq in loops_in_seq]

def name_hor_tree(hor, node_prefix='', clade_name_prefix='F', hor_name_prefix='H', level_separator='_'):
    hor.id = f'{hor_name_prefix}{node_prefix}'
    hor.feature = SeqFeature(
        id=hor.id,
        location=(
            None
                if hor.locations is None or len(hor.locations) == 0
            else hor.locations[0]
                if len(hor.locations) == 1
            else CompoundLocation(hor.locations)
        )
    )
    clade_count = 0
    if len(hor.hor.clade_seq) == 1:
        clade = hor.hor.clade_seq[0]
        if clade.name is None:
            clade.name = f'{clade_name_prefix}{node_prefix}'
    for clade in hor.hor.clade_seq:
        if clade.name is None:
            clade_count += 1
            clade.name = f'{clade_name_prefix}{node_prefix}#{clade_count}'
    common_prefix_for_sub_hors = f"{node_prefix}{level_separator if len(node_prefix) > 0 else ''}"
    for sub_hor_index, sub_hor in enumerate(hor.sub_hors):
        name_hor_tree(
            sub_hor,
            node_prefix=f'{common_prefix_for_sub_hors}{sub_hor_index + 1}')
        
def hor_to_clade(hor):
    clade_seq_str = ",".join([clade.name for clade in hor.hor.clade_seq])
    return Clade(
        # node_id=hor.id,
        name=hor.id,
        sequences=[Sequence(type='dna', location=location) for location in hor.locations],
        clades=[hor_to_clade(sub_hor) for sub_hor in hor.sub_hors],
        properties=[Property(value=clade_seq_str, ref='monomer_clade_seq', applies_to='clade', datatype='xsd:string')]
    )

def hor_tree_to_phylogeny(hor_tree_root, name='hors'):
    return Phylogeny(root=hor_to_clade(hor_tree_root), name=name)

