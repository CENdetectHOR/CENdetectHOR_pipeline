from itertools import accumulate
from cluster import get_seq_as_txt

def normalize_loop(loop_seq: list[int]) -> tuple[list[int], int]:
    def invert_pos(pos: int) -> int:
        return len(loop_seq) - pos if pos > 0 else 0
    options = []
    for start_index in range(len(loop_seq)):
        options.append(loop_seq[start_index:] + loop_seq[:start_index] + [start_index])
    options.sort()
    return (options[0][:-1],invert_pos(options[0][-1]))

def denormalize_loop(
    loop_seq: list[int],
    in_loop_start: int
) -> list[int]:
    return loop_seq[in_loop_start:] + loop_seq[:in_loop_start]

class Loop:
    loop_seq: list[int]

    def __init__(self, loop_seq: list[int]):
        self.loop_seq = loop_seq

    def __str__(self):
        return get_seq_as_txt(self.loop_seq)

class LoopSpanInSeq:
    span_start: int
    span_length: int
    num_of_laps: int
    in_loop_start: int
    num_mismatches: int

    def __init__(
        self,
        span_start: int, span_length: int,
        num_of_laps: int, in_loop_start: int,
        num_mismatches: int = 0
    ):
        if span_length <= 0:
            raise Exception(f"Span length {span_length} not permitted")
        self.span_start = span_start
        self.span_length = span_length
        self.num_of_laps = num_of_laps
        self.in_loop_start = in_loop_start
        self.num_mismatches = num_mismatches

    def __str__(self):
        return (
            f'[{self.span_start}:{self.span_start + self.span_length}]' +
            (f'#{self.in_loop_start}' if self.in_loop_start != 0 else '') +
            (f'mm({self.num_mismatches})' if self.num_mismatches != 0 else '')
        )

class LoopInSeq:
    loop: Loop
    spans_in_seq: list[LoopSpanInSeq]

    def __init__(self, loop: Loop, spans_in_seq: list[LoopSpanInSeq] = []):
        self.loop = loop
        self.spans_in_seq = spans_in_seq

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
        
def loop_to_spans(loop_in_seq: LoopInSeq) -> list[tuple[int,int]]:
    spans: list[tuple[int,int]] = []
    for loop_span in loop_in_seq.spans_in_seq:
        loop_span_end = loop_span.span_start + loop_span.span_length
        if len(spans) > 0 and loop_span.span_start <= spans[-1][1]:
            if loop_span_end > spans[-1][1]:
                spans.append((spans[-1][1], loop_span_end))
        else:
            spans.append((loop_span.span_start, loop_span_end))
    return spans

def extract_spans_in_seqs(
    spans_in_seq: list[LoopSpanInSeq],
    spans_to_extract: list[tuple[int,int]],
    loop_size: int
) -> list[LoopSpanInSeq]:
    new_spans_in_seq: list[LoopSpanInSeq] = []
    span_to_extract_index = 0
    for span_in_seq in spans_in_seq:
        span_in_seq_end = span_in_seq.span_start + span_in_seq.span_length
        while span_to_extract_index < len(spans_to_extract):
            span_to_extract_start, span_to_extract_end = spans_to_extract[span_to_extract_index]
            if (span_to_extract_start < span_in_seq_end and span_to_extract_end > span_in_seq.span_start):
                new_span_start = max(span_in_seq.span_start, span_to_extract_start)
                new_span_end = min(span_in_seq_end, span_to_extract_end)
                new_span_length = new_span_end - new_span_start
                new_in_loop_start = (
                    (span_in_seq.in_loop_start + new_span_start - span_in_seq.span_start) % loop_size
                )
                new_spans_in_seq.append(LoopSpanInSeq(
                    span_start=new_span_start,
                    span_length=new_span_length,
                    num_of_laps=new_span_length // loop_size,
                    in_loop_start=new_in_loop_start,
                    num_mismatches=span_in_seq.num_mismatches
                ))
            if span_to_extract_start < span_in_seq.span_start:
                spans_to_extract.insert(
                    span_to_extract_index,
                    (
                        span_to_extract_start,
                        min(span_to_extract_end, span_in_seq.span_start)
                    )
                )
                span_to_extract_index += 1
            if (span_to_extract_end > span_in_seq_end):
                spans_to_extract[span_to_extract_index] = (
                    max(span_in_seq_end, span_to_extract_start),
                    span_to_extract_end
                 )
                break
            else:
                del spans_to_extract[span_to_extract_index]
    return new_spans_in_seq

def intersect_spans(
    spans_a: list[tuple[int,int]],
    spans_b: list[tuple[int,int]]
) -> list[tuple[int,int]]:
    intersection_spans = []
    next_span_b_index = 0
    for start_pos_a, end_pos_a in spans_a:
        for span_b_index_offset, span_b_boundaries in enumerate(spans_b[next_span_b_index:]):
            start_pos_b, end_pos_b = span_b_boundaries
            if (start_pos_b < end_pos_a and end_pos_b > start_pos_a):
                new_span_start = max(start_pos_a, start_pos_b)
                new_span_end = min(end_pos_a, end_pos_b)
                intersection_spans.append((new_span_start, new_span_end))
            if (end_pos_b > end_pos_a):
                next_span_b_index += span_b_index_offset
                break
    return intersection_spans

def complement_of_spans(
    spans: list[tuple[int,int]],
    whole_range: int
) -> list[tuple[int,int]]:
    if len(spans) == 0:
        return [(0, whole_range)]
    complement = []
    if spans[0][0] > 0:
        complement.append((0, spans[0][0]))
    complement.extend([
        (spans[i][1], spans[i+1][0])
        for i in range(len(spans) - 1)
    ])
    if spans[-1][1] < whole_range:
        complement.append((spans[-1][1], whole_range))
    return complement
    

def spans_difference(
    spans_a: list[tuple[int,int]],
    spans_b: list[tuple[int,int]]
):
    if len(spans_a) == 0:
        return []
    return intersect_spans(
        spans_a,
        complement_of_spans(spans_b, whole_range=spans_a[-1][1]))
        

def find_loops(
    seqs: list[list[int]] = None,
    whole_seq: list[int] = None,
    gap_indices: list[int] = None,
    seq_spans: list[tuple[int,int]] = None,
    min_loop_size: int = 2, max_loop_size: int = 30, min_loops: int = 3,
    allowed_mismatch_rate: float = 0,
    allow_overlap: bool = False
) -> list[LoopInSeq]:
    
    if seqs is None:
        if seq_spans is None:
            if gap_indices is None:
                gap_indices = []
            seq_offsets = [0] + gap_indices + [len(whole_seq)]
            seq_spans = [
                (seq_offsets[i],seq_offsets[i+1])
                for i in range(len(gap_indices) + 1)
            ]
            seqs = [
                whole_seq[seq_offsets[i]:seq_offsets[i+1]]
                for i in range(len(gap_indices) + 1)
            ]
        else:
            seq_offsets = [seq_start for seq_start, seq_end in seq_spans]
            seqs = [
                whole_seq[seq_start:seq_end]
                for seq_start, seq_end in seq_spans
            ]
    else:
        seq_offsets = [0] + list(accumulate([len(seq) for seq in seqs]))
        seq_spans = [
            (seq_offsets[i],seq_offsets[i+1])
            for i in range(len(seqs))
        ]
        
    loops_found: dict[str,LoopInSeq] = {} #defaultdict(list)
    for seqIndex, seq in enumerate(seqs):

        curr_loops = {loop_size:0 for loop_size in range(1, max_loop_size + 1)}
        curr_loops_mismatches = {loop_size:0 for loop_size in range(1, max_loop_size + 1)}

        def last_of_size(curr_position, loop_size):
            if curr_loops[loop_size] >= (min_loops - 1) * loop_size:
                loop_start = curr_position - curr_loops[loop_size] - loop_size
                loop_length = curr_loops[loop_size] + loop_size
                loop_laps = loop_length // loop_size
                loop_items = seq[loop_start:loop_start + loop_size]
                normal_loop, in_loop_start_position = normalize_loop(loop_items)
                normal_loop_str = str(normal_loop)
                loop_span = LoopSpanInSeq(
                    seq_offsets[seqIndex] + loop_start,
                    loop_length, loop_laps, in_loop_start_position,
                    num_mismatches=curr_loops_mismatches[loop_size])
                if normal_loop_str not in loops_found:
                    loops_found[normal_loop_str] = LoopInSeq(
                        Loop(normal_loop),
                        [loop_span]
                    )
                else:
                    loops_found[normal_loop_str].add_span(loop_span)
        
        for curr_position, curr_symbol in enumerate(seq):
            max_loop_length_closed = 0
            for loop_size in range(1, min(max_loop_size, curr_position) + 1):
                if seq[curr_position - loop_size] == curr_symbol:
                    curr_loops[loop_size] += 1
                elif curr_loops_mismatches[loop_size] + 1 <= (curr_loops[loop_size] + 1) * allowed_mismatch_rate:
                    curr_loops[loop_size] += 1
                    curr_loops_mismatches[loop_size] += 1
                else:
                    if curr_loops[loop_size] > max_loop_length_closed:
                        last_of_size(curr_position, loop_size)
                        max_loop_length_closed = curr_loops[loop_size]
                    curr_loops[loop_size] = 0
                    curr_loops_mismatches[loop_size] = 0
        
        max_loop_length_closed = 0
        for loop_size in range(1, max_loop_size + 1):
            if curr_loops[loop_size] > max_loop_length_closed:
                last_of_size(len(seq), loop_size)
                max_loop_length_closed = curr_loops[loop_size]

    loops = list(loops_found.values())

    if min_loop_size > 1:
        loops = [loop for loop in loops if len(loop.loop.loop_seq) >= min_loop_size]
        
    if not allow_overlap:
        loops.sort(key=lambda loop: len(loop.loop.loop_seq), reverse=True)
        filtered_loops = []
        available_spans = seq_spans
        for loop in loops:
            extracted_spans = extract_spans_in_seqs(
                spans_in_seq=loop.spans_in_seq,
                spans_to_extract=available_spans,
                loop_size=len(loop.loop.loop_seq)
            )
            if len(extracted_spans) > 0:
                filtered_loop = LoopInSeq(
                    loop=loop.loop,
                    spans_in_seq=extracted_spans
                )
                filtered_loops.append(filtered_loop)
        loops = filtered_loops

    for loop in loops:
        spans = loop.spans_in_seq
        if all([span.in_loop_start == spans[0].in_loop_start for span in spans]):
            loop.loop.loop_seq = denormalize_loop(loop.loop.loop_seq, spans[0].in_loop_start)
            for span in spans:
                span.in_loop_start = 0

    return loops
