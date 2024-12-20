from collections.abc import Iterable
#from typing import Self
from loops import LoopSpanInSeq, LoopInSeq
from Bio.Phylo import BaseTree
from Bio.Phylo.PhyloXML import Clade, Phylogeny, Sequence, Phyloxml, Property
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation, Location

class HOR:
    clade_seq: list[BaseTree.Clade]

    def __init__(self, clade_seq: list[BaseTree.Clade]):
        self.clade_seq = clade_seq

    def __str__(self):
        return ''.join([clade.id if clade.id is not None else '*' for clade in self.clade_seq])

class HORInSeq:
    hor: HOR
    locations: list[Location]
    spans_in_seq: list[LoopSpanInSeq]
    super_hor: object
    sub_hors: list[object]

    def __init__(
        self,
        hor: HOR,
        spans_in_seq: list[LoopSpanInSeq] = [],
        locations: list[Location] = None
    ):
        self.hor = hor
        self.spans_in_seq = spans_in_seq
        self.locations = locations

    def add_span(self, span_in_seq: LoopSpanInSeq):
        self.spans_in_seq.append(span_in_seq)

    def __str__(self):
        return (
            f'{self.loop}' +
            (
                f' in {",".join([str(span) for span in self.spans_in_seq])}'
                    if len(self.spans_in_seq) > 0 else ''
            )
        )

def seq_span_to_location(
    span: LoopSpanInSeq,
    seq_locations: list[SimpleLocation]
) -> SimpleLocation:
    start_location = seq_locations[span.span_start]
    end_location = seq_locations[span.span_start + span.span_length - 1]
    ref = start_location.ref
    strand = start_location.strand
    return SimpleLocation(
        ref=ref,
        strand=strand,
        start=min(start_location.start,end_location.start),
        end=max(end_location.end,start_location.end)
    )

def seq_spans_to_compound_location(
    spans_in_seq: Iterable[LoopSpanInSeq],
    seq_locations: list[SimpleLocation]
) -> CompoundLocation:
    return CompoundLocation([seq_span_to_location(span, seq_locations) for span in spans_in_seq])

def loop_to_HOR(
    loop_in_seq: LoopInSeq,
    clades: list[Clade],
    seq_locations: list[SimpleLocation]=None
) -> HORInSeq:
    hor = HOR([clades[clade_index] for clade_index in loop_in_seq.loop.loop_seq])
    return HORInSeq(
        hor,
        spans_in_seq=loop_in_seq.spans_in_seq,
        locations=(
            [seq_span_to_location(span, seq_locations) for span in loop_in_seq.spans_in_seq]
            if seq_locations is not None else None
        )
    )

def loops_to_HORs(
    loops_in_seq: Iterable[LoopInSeq], clades: list[Clade], seq_locations=None
) -> list[HORInSeq]:
    return [loop_to_HOR(loop_in_seq, clades, seq_locations=seq_locations) for loop_in_seq in loops_in_seq]

def name_hor_tree(
    hor: HORInSeq, 
    node_prefix: str = '',
    clade_name_prefix: str = 'F',
    hor_name_prefix: str = 'H',
    level_separator: str = '_'
) -> None:
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
    if hasattr(hor, 'sub_hors'):
        for sub_hor_index, sub_hor in enumerate(hor.sub_hors):
            name_hor_tree(
                sub_hor,
                node_prefix=f'{common_prefix_for_sub_hors}{sub_hor_index + 1}')
        
def hor_to_clade(hor: HORInSeq) -> Clade:
    clade_seq_str = ",".join([clade.name for clade in hor.hor.clade_seq])
    return Clade(
        # node_id=hor.id,
        name=hor.id,
        sequences=[Sequence(type='dna', location=location) for location in hor.locations],
        clades=
            [hor_to_clade(sub_hor) for sub_hor in hor.sub_hors]
            if hasattr(hor, 'sub_hors') else [],
        properties=[Property(value=clade_seq_str, ref='monomer_clade_seq', applies_to='clade', datatype='xsd:string')]
    )

def hor_tree_to_phylogeny(
    hor_tree_root: HORInSeq,
    name: str = 'hors'
) -> Phylogeny:
    return Phylogeny(root=hor_to_clade(hor_tree_root), name=name)

