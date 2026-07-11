#!/usr/bin/env python3
"""Generate styled code-box figures for the ProtOS chapter.

Each processor gets a minimal code example rendered as a publication-ready PNG.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.font_manager import FontProperties

THESIS_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(THESIS_DIR / "shared"))
from thesis_style import DARK_GRAY, WARM_WHITE, LIGHT_GRAY, SLATE, SAGE, OCHRE, RUST, apply_style

apply_style()

OUTPUT_DIR = Path(__file__).resolve().parent
MONO = FontProperties(family="monospace", size=10)
MONO_BOLD = FontProperties(family="monospace", size=10, weight="bold")
TITLE_FONT = FontProperties(family="sans-serif", size=9, weight="bold")

# Syntax colors — muted, thesis-consistent
C_KEYWORD = SLATE       # from, import, class, def, ...
C_STRING = SAGE          # "strings"
C_COMMENT = LIGHT_GRAY   # # comments
C_CALL = OCHRE           # function calls
C_DEFAULT = DARK_GRAY    # everything else

KEYWORDS = {
    "from", "import", "as", "class", "def", "return", "if", "else",
    "for", "in", "not", "and", "or", "True", "False", "None", "with",
}


def tokenize_line(line: str) -> list[tuple[str, str]]:
    """Split a line into (text, color) tokens for basic Python highlighting."""
    tokens = []
    if not line.strip():
        return [("", C_DEFAULT)]

    # Comment line
    stripped = line.lstrip()
    if stripped.startswith("#"):
        return [(line, C_COMMENT)]

    # Inline comment
    comment_idx = None
    in_str = None
    for i, ch in enumerate(line):
        if ch in ('"', "'") and not in_str:
            in_str = ch
        elif ch == in_str:
            in_str = None
        elif ch == "#" and not in_str:
            comment_idx = i
            break

    main_part = line[:comment_idx] if comment_idx is not None else line
    comment_part = line[comment_idx:] if comment_idx is not None else None

    # Tokenize main part
    pattern = r'("(?:[^"\\]|\\.)*"|\'(?:[^\'\\]|\\.)*\'|\w+|[^\w\s]+|\s+)'
    for match in re.finditer(pattern, main_part):
        tok = match.group()
        if tok.startswith(("'", '"')):
            tokens.append((tok, C_STRING))
        elif tok in KEYWORDS:
            tokens.append((tok, C_KEYWORD))
        elif re.match(r'\w+$', tok) and tokens and tokens[-1][0].rstrip().endswith('.'):
            # method call after dot
            tokens.append((tok, C_DEFAULT))
        elif re.match(r'\w+$', tok):
            # Check if followed by '(' -> function/class call
            end = match.end()
            rest = main_part[end:]
            if rest.lstrip().startswith('('):
                tokens.append((tok, C_CALL))
            else:
                tokens.append((tok, C_DEFAULT))
        else:
            tokens.append((tok, C_DEFAULT))

    if comment_part:
        tokens.append((comment_part, C_COMMENT))

    return tokens


PAD_X = 0.30        # horizontal padding (left = right)
PAD_Y = 0.20        # vertical padding above/below code
TITLE_H = 0.35      # title bar height
LINE_H = 0.24       # height per code line
FONT_SIZE = 9.5     # monospace code font size


def measure_global_width(figures: list[dict]) -> float:
    """Measure the widest code line across *all* figures and return box width."""
    tmp_fig = plt.figure(figsize=(14, 1))
    tmp_fig.canvas.draw()
    renderer = tmp_fig.canvas.get_renderer()
    max_px = 0
    for spec in figures:
        for line in spec["code"].rstrip("\n").split("\n"):
            if not line.strip():
                continue
            t = tmp_fig.text(
                0, 0, line,
                fontsize=FONT_SIZE, fontfamily="monospace",
            )
            bb = t.get_window_extent(renderer)
            max_px = max(max_px, bb.width)
    max_text_in = max_px / tmp_fig.dpi
    plt.close(tmp_fig)
    # Symmetric padding: pad_x on left + text + pad_x on right
    return max(max_text_in + 2 * PAD_X, 4.0)


def render_code_figure(
    code: str,
    title: str,
    filename: str,
    box_w_in: float,
) -> Path:
    """Render a code snippet as a styled figure with a fixed box width."""
    lines = code.rstrip("\n").split("\n")
    n_lines = len(lines)

    code_h = n_lines * LINE_H
    box_h_in = TITLE_H + PAD_Y + code_h + PAD_Y
    fig_w = box_w_in + 0.3
    fig_h = box_h_in + 0.2

    fig = plt.figure(figsize=(fig_w, fig_h))

    # Box in figure-fraction coordinates
    bx = 0.15 / fig_w
    by = 0.1 / fig_h
    bw = box_w_in / fig_w
    bh = box_h_in / fig_h

    box = mpatches.FancyBboxPatch(
        (bx, by), bw, bh,
        boxstyle="round,pad=0.012",
        facecolor=WARM_WHITE,
        edgecolor=LIGHT_GRAY,
        linewidth=0.8,
        transform=fig.transFigure,
    )
    fig.patches.append(box)

    # Title separator
    sep_y = by + bh - TITLE_H / fig_h
    fig.lines.append(plt.Line2D(
        [bx + 0.01, bx + bw - 0.01], [sep_y, sep_y],
        color=LIGHT_GRAY, linewidth=0.5,
        transform=fig.transFigure, figure=fig,
    ))

    # Title
    fig.text(
        bx + PAD_X / fig_w, sep_y + (TITLE_H / fig_h) * 0.5,
        title,
        fontsize=9, fontweight="bold", fontfamily="sans-serif",
        color=DARK_GRAY, va="center",
        transform=fig.transFigure,
    )

    # Render code lines with syntax highlighting (token-by-token)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    code_top = sep_y - PAD_Y / fig_h
    left = bx + PAD_X / fig_w

    for i, line in enumerate(lines):
        y = code_top - i * (LINE_H / fig_h)
        tokens = tokenize_line(line)
        x_px = fig.transFigure.transform((left, 0))[0]

        for text, color in tokens:
            if not text:
                continue
            x_fig = fig.transFigure.inverted().transform((x_px, 0))[0]
            t = fig.text(
                x_fig, y, text,
                fontsize=FONT_SIZE, fontfamily="monospace",
                color=color, va="top",
                transform=fig.transFigure,
            )
            bb = t.get_window_extent(renderer)
            x_px += bb.width

    out = OUTPUT_DIR / filename
    fig.savefig(str(out), dpi=300, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"  {out.name}")
    return out


# =============================================================================
# Code examples — one per processor, max 6-7 lines
# =============================================================================
FIGURES = [
    {
        "title": "Sequence Processor",
        "filename": "fig_2.2_opsin_diversity/fig_2.2_code.png",
        "code": """\
seq_proc = SequenceProcessor()
loader = SequenceLoader(processor=seq_proc)
loader.download_and_register("uniprot:P02699", name="RHO_BOVIN")

ncbi = NCBILoader(processor=seq_proc)
hits = ncbi.blast_search(
    sequence=seq_proc.load_entity("RHO_BOVIN"),
    database="swissprot", hitlist_size=5000
)
ncbi.download_batch(accessions, dataset_name="opsin_atlas")\
""",
    },
    {
        "title": "Structure Processor",
        "filename": "fig_2.3_br_binding_pocket/fig_2.3_code.png",
        "code": """\
struct_proc = StructureProcessor()
loader = StructureLoader(processor=struct_proc)
loader.download_batch(
    ["1C3W", "1U19"], dataset_name="opsins"
)

contacts = struct_proc.get_ligand_interactions(
    "1C3W", ligand_id="RET",
    chain_id="A", cutoff=4.0
)\
""",
    },
    {
        "title": "Graph Processor",
        "filename": "fig_2.4_pocket_graphs/fig_2.4_code.png",
        "code": """\
graph_proc = GraphProcessor(
    default_cutoff=7.0,
    default_level="residue"
)
graph_name = graph_proc.generate_graph(
    "1U19", chain="A",
    near_hetatm=("RET", 7.0)
)
graph_data = graph_proc.load_entity(graph_name)\
""",
    },
    {
        "title": "GRN Processor",
        "filename": "fig_2.5_grn_comparison/fig_2.5_code.png",
        "code": """\
grn_table = seq_proc.annotate_with_grn(
    dataset_name="opsin_atlas",
    reference_table="type II",
    protein_family="gpcr_a",
)
grn_proc = GRNProcessor()
table = grn_proc.load_table("atlas_grn")\
""",
    },
    {
        "title": "Embedding Processor",
        "filename": "fig_2.6_embeddings/fig_2.6_code.png",
        "code": """\
emb_proc = EmbeddingProcessor(
    model_name="ankh_large"
)
embeddings = emb_proc.embed_sequences(
    sequences,
    embedding_type="per_residue"
)\
""",
    },
    {
        "title": "Property Processor",
        "filename": "fig_2.7_properties/fig_2.7_code.png",
        "code": """\
prop_proc = PropertyProcessor()
prop_proc.record_properties(
    "atlas_properties", rows,
    allow_create=True
)
table = prop_proc.load_table(
    "atlas_properties"
)\
""",
    },
    {
        "title": "Model Manager",
        "filename": "fig_2.8_model_manager/fig_2.8_code.png",
        "code": """\
manager = ModelManager()
inv = manager.prepare(
    "boltz2",
    inputs={"sequence": seq,
            "name": "design_01"}
)
result = manager.run_and_ingest(inv)\
""",
    },
]


def main():
    print("Generating code figures...")
    box_w = measure_global_width(FIGURES)
    print(f"  Uniform box width: {box_w:.2f} in")
    for spec in FIGURES:
        render_code_figure(spec["code"], spec["title"], spec["filename"], box_w)
    print("Done.")


if __name__ == "__main__":
    main()
