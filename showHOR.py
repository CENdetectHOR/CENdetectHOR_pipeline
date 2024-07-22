import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from Bio import Phylo
from treeFromClusters import new_phylogeny


def show_hor(hor_in_seq, tree=None, seq_length=None,
             label='unnamed',
             color_palette=[mcolors.TABLEAU_COLORS[color_id]
                            for color_id in mcolors.TABLEAU_COLORS],
             seq_fg_color='black',
             seq_bg_color='#D3D3D3',
             tree_bg_color='#D3D3D3',
             hor_item_size_ratio=(0.04, 0.05),
             hor_item_height=1,
             tree_height=10,
             figure_width=7,
             #  seq_item_size=(0.01,0.2),
             seq_height_ratio=0.05,
             tree_height_ratio=0.9,
             h_space=0.1,
             monomer_size=107,
             label_func=(lambda clade_name: None)):

    spans = hor_in_seq.spans_in_seq

    bp_start = spans[0].span_start * monomer_size
    bp_end = (spans[-1].span_start + spans[-1].span_length) * monomer_size

    print(f'HOR: {label}')
    print(f'Coverage (bp): {bp_start}-{bp_end}')

    if seq_length is None:
        seq_length = tree.count_terminals()

    clade_seq = hor_in_seq.hor.clade_seq
    for clade_index, clade in enumerate(set(clade_seq)):
        clade.color = color_palette[clade_index]

    fig = plt.figure(constrained_layout=True, figsize=(
        figure_width, len(clade_seq) * hor_item_height + tree_height))
    fig.suptitle(f'HOR {label} ({bp_start}-{bp_end})')

    hor_size_ratio = hor_item_size_ratio[0]*len(clade_seq)

    gs = fig.add_gridspec(
        3, 2,
        width_ratios=[hor_size_ratio, 1-hor_size_ratio],
        height_ratios=[hor_item_size_ratio[0],
                       tree_height_ratio, seq_height_ratio],
        hspace=h_space
    )
    ax_hor = plt.subplot(gs[0, 0])
    ax_tree = plt.subplot(gs[1, :])
    ax_seq = plt.subplot(gs[2, :])

    ax_hor.set_xlim([0, len(clade_seq)])
    ax_hor.grid(True)
    ax_hor.set_xticks(np.arange(0, len(clade_seq) + 1, 1))
    ax_hor.get_yaxis().set_visible(False)
    for clade_pos, clade in enumerate(clade_seq):
        ax_hor.add_patch(patches.Rectangle((clade_pos, 0), 1,
                         1, facecolor=clade.color.to_hex()))

    bp_seq_length = seq_length * monomer_size

    ax_seq.set_xlim([0, bp_seq_length])
    ax_seq.get_yaxis().set_visible(False)

    ax_seq.add_patch(patches.Rectangle(
        (0, 0), bp_seq_length, 1, facecolor=seq_bg_color))
    for span in spans:
        span_bp_start = span.span_start * monomer_size
        span_bp_length = span.span_length * monomer_size
        ax_seq.add_patch(patches.Rectangle((span_bp_start, 0),
                         span_bp_length, 1, facecolor=seq_fg_color))

    ax_tree.get_yaxis().set_visible(False)
    ax_tree.get_xaxis().set_visible(False)
    # ax_tree.spines['top'].set_visible(False)
    # ax_tree.spines['right'].set_visible(False)
    # ax_tree.spines['bottom'].set_visible(False)
    # ax_tree.spines['left'].set_visible(False)
    ax_tree.get_yaxis().set_ticks([])
    # ax_tree.axis('off')
    if tree is not None:
        tree.root.color = tree_bg_color
        Phylo.draw(tree, axes=ax_tree, label_func=label_func)

    for clade_index, clade in enumerate(set(hor_in_seq.hor.clade_seq)):
        clade.color = None
    if tree is not None:
        tree.root.color = None


def show_hors(hors_in_seq, reference_seq,
              tree=None,
              label='unnamed',
              color_palette=[mcolors.TABLEAU_COLORS[color_id]
                             for color_id in mcolors.TABLEAU_COLORS],
              seq_fg_color='black',
              seq_bg_color='#D3D3D3',
              tree_bg_color='#D3D3D3',
              hor_item_size_ratio=(0.04, 0.5),
              seq_height_ratio=0.5,
              fig_width=7,
              tree_height=5,
              hor_height=1,
              h_space=0.1,
              h_space_in_hor=0.02,
              label_func=(lambda clade_name: None)):

    clade_set = set(
        [clade for hor_in_seq in hors_in_seq for clade in hor_in_seq.hor.clade_seq])
    for clade_index, clade in enumerate(clade_set):
        # this is bad, but at least avoid error
        clade.color = color_palette[clade_index % 10]

    hors_height = hor_height * len(hors_in_seq)
    fig = plt.figure(constrained_layout=True, figsize=(
        fig_width, tree_height + hors_height))
    # fig.suptitle(f'HOR {label} ({bp_start}-{bp_end})')

    print([hors_height for i in range(len(hors_in_seq))].append(tree_height))

    gs = fig.add_gridspec(
        len(hors_in_seq) + 1, 1,
        # width_ratios=[hor_size_ratio,1-hor_size_ratio],
        height_ratios=[hor_height for i in range(
            len(hors_in_seq))] + [tree_height],
        hspace=h_space
    )
    # ax_hors = [plt.subplot(gs[i, :]) for i in range(len(hors_in_seq))]
    ax_tree = plt.subplot(gs[len(hors_in_seq), :])

    for hor_index in range(len(hors_in_seq)):
        # ax_hor_and_seq = plt.subplot(gs[hor_index, 0])

        hor_in_seq = hors_in_seq[hor_index]
        clade_seq = hor_in_seq.hor.clade_seq
        hor_size_ratio = hor_item_size_ratio[0]*len(clade_seq)

        gs_hor = gridspec.GridSpecFromSubplotSpec(2, 2,
                                                  subplot_spec=gs[hor_index, 0],
                                                  width_ratios=[
                                                      hor_size_ratio, 1-hor_size_ratio],
                                                  height_ratios=[1, 1],
                                                  hspace=h_space_in_hor
                                                  )

        # gs_hor = ax_hor_and_seq.add_gridspec(
        #     2, 2,
        #     width_ratios=[hor_size_ratio,1-hor_size_ratio],
        #     height_ratios=[1,1],
        #     hspace=h_space
        # )

        ax_hor = plt.subplot(gs_hor[0, 0])
        ax_seq = plt.subplot(gs_hor[1, :])

        ax_hor.set_xlim([0, len(clade_seq)])
        ax_hor.grid(True)
        ax_hor.set_xticks(np.arange(0, len(clade_seq) + 1, 1))
        ax_hor.get_yaxis().set_visible(False)
        for clade_pos, clade in enumerate(clade_seq):
            ax_hor.add_patch(patches.Rectangle(
                (clade_pos, 0), 1, 1, facecolor=clade.color.to_hex()))

        bp_seq_length = len(reference_seq)

        ax_seq.set_xlim([0, bp_seq_length])
        ax_seq.get_yaxis().set_visible(False)

        # spans = hor_in_seq.spans_in_seq
        locations = hor_in_seq.locations

        # bp_start = spans[0].span_start * monomer_size
        # bp_end = (spans[-1].span_start + spans[-1].span_length) * monomer_size

        ax_seq.add_patch(patches.Rectangle(
            (0, 0), bp_seq_length, 1, facecolor=seq_bg_color))
        for location in locations:
            ax_seq.add_patch(patches.Rectangle(
                (location.start, 0),
                location.end - location.start,
                1,
                facecolor=seq_fg_color))

    ax_tree.get_yaxis().set_visible(False)
    ax_tree.get_xaxis().set_visible(False)
    # ax_tree.spines['top'].set_visible(False)
    # ax_tree.spines['right'].set_visible(False)
    # ax_tree.spines['bottom'].set_visible(False)
    # ax_tree.spines['left'].set_visible(False)
    ax_tree.get_yaxis().set_ticks([])
    # ax_tree.axis('off')
    if tree is not None:
        tree.root.color = tree_bg_color
        Phylo.draw(tree, axes=ax_tree, label_func=label_func)
        # Phylo.draw_graphviz(tree, axes=ax_tree)

    for clade in clade_set:
        clade.color = None
    if tree is not None:
        tree.root.color = None


def show_hor_tree(hor_tree_root, reference_seq, tree, path=[], level=0):
    if not hor_tree_root.sub_hors:
        return
    print(f'Subtree: {path}')
    show_hors(hor_tree_root.sub_hors, reference_seq, tree=tree)
    for branch_index, sub_hor in enumerate(hor_tree_root.sub_hors):
        show_hor_tree(sub_hor, reference_seq, tree=new_phylogeny(
            tree.common_ancestor(sub_hor.hor.clade_seq)), path=path + [branch_index + 1])
