from typing import List
from itertools import accumulate
from cluster import get_seq_as_txt

def normalize_loop(loop_seq):
    def invert_pos(pos):
        return len(loop_seq) - pos if pos > 0 else 0
    options = []
    for start_index in range(len(loop_seq)):
        options.append(loop_seq[start_index:] + loop_seq[:start_index] + [start_index])
    options.sort()
    return (options[0][:-1],invert_pos(options[0][-1]))

def denormalize_loop(loop_seq, in_loop_start):
    return loop_seq[in_loop_start:] + loop_seq[:in_loop_start]

class Loop:
    loop_seq: List[int]

    def __init__(self, loop_seq):
        self.loop_seq = loop_seq

    def __str__(self):
        return get_seq_as_txt(self.loop_seq)

class LoopSpanInSeq:
    def __init__(self, span_start, span_length, num_of_laps, in_loop_start):
        self.span_start = span_start
        self.span_length = span_length
        self.num_of_laps = num_of_laps
        self.in_loop_start = in_loop_start

    def __str__(self):
        return (
            f'[{self.span_start}:{self.span_start + self.span_length}]' +
            (f'#{self.in_loop_start}' if self.in_loop_start != 0 else '')
        )

class LoopInSeq:
    loop: Loop
    spans_in_seq: List[LoopSpanInSeq]

    def __init__(self, loop, spans_in_seq = []):
        self.loop = loop
        self.spans_in_seq = spans_in_seq

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

def find_loops(seqs, min_loop_size = 2, max_loop_size = 30, min_loops = 3):
    seq_offsets = [0] + list(accumulate([len(seq) for seq in seqs]))
    loops_found = {} #defaultdict(list)
    for seqIndex, seq in enumerate(seqs):

        curr_loops = {loop_size:0 for loop_size in range(1, max_loop_size + 1)}

        def last_of_size(curr_position, loop_size):
            if curr_loops[loop_size] >= (min_loops - 1) * loop_size:
                loop_start = curr_position - curr_loops[loop_size] - loop_size
                loop_length = curr_loops[loop_size] + loop_size
                loop_laps = loop_length // loop_size
                loop_items = seq[loop_start:loop_start + loop_size]
                (normal_loop, in_loop_start_position) = normalize_loop(loop_items)
                normal_loop_str = str(normal_loop)
                loop_span = LoopSpanInSeq(seq_offsets[seqIndex] + loop_start, loop_length, loop_laps, in_loop_start_position)
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
                else:
                    if curr_loops[loop_size] > max_loop_length_closed:
                        last_of_size(curr_position, loop_size)
                        max_loop_length_closed = curr_loops[loop_size]
                    curr_loops[loop_size] = 0
        
        max_loop_length_closed = 0
        for loop_size in range(1, max_loop_size + 1):
            if curr_loops[loop_size] > max_loop_length_closed:
                last_of_size(len(seq), loop_size)
                max_loop_length_closed = curr_loops[loop_size]

    loops = list(loops_found.values())

    if min_loop_size > 1:
        loops = [loop for loop in loops if len(loop.loop.loop_seq) >= min_loop_size]

    for loop in loops:
        spans = loop.spans_in_seq
        if all([span.in_loop_start == spans[0].in_loop_start for span in spans]):
            loop.loop.loop_seq = denormalize_loop(loop.loop.loop_seq, spans[0].in_loop_start)
            for span in spans:
                span.in_loop_start = 0

    return loops
