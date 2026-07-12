# ProtOS GRN annotation: detailed pseudocode and behaviour

This document describes the implementation, not an aspirational replacement.
The primary path is `GRNProcessor.annotate_sequences`, called by
`SequenceProcessor.annotate_with_grn` and the structure annotation workflows.
ProtOS also retains an older numeric expansion path in
`assign_grns_to_sequences` / `expand_annotation`; that distinction is important.

## Data model

A reference table has:

- one protein or design per row;
- one GRN per column;
- cells such as `A123`, meaning residue `A` at one-based sequence position 123;
- `-` for a GRN absent from that reference.

Ordinary non-flexible coordinates use the general grammar:

```text
GRN := SEGMENT_ID "." POSITION
SEGMENT_ID := one or more non-empty dot-separated string components
POSITION := positive integer
```

Consequently all of these are valid:

```text
1.50                 numeric GPCR segment
TM1.50               named segment
G.HN.03              CGN segment
H.ha-hb.01           CGN connecting segment
N.S1.01              CAN segment
A.H1.50              class-specific hierarchical segment
```

Only the final component is the position. `G.HN` is one segment identifier; it
is never converted to a float. A segment is not necessarily a helix.

Flexible regions between two intrinsically ordered segments use a separate,
directional grammar:

```text
FLEXIBLE_GRN := NEARER_SEGMENT_INDEX FARTHER_SEGMENT_INDEX "." DISTANCE
12.003       := between segments 1 and 2, three residues from segment 1
21.003       := the same region, three residues from segment 2
```

The segment indices are one-based positions in protein sequence order, so this
notation is class-agnostic and also works when the real segment identifiers are
strings. The compact representation currently supports up to nine ordered
segments. Numeric GPCR utilities remain available for thesis-era arithmetic,
but they do not define the meaning of directional flexible coordinates.

An insertion inside a non-flexible segment adds `.001` to the position in the
GPCRdb sense. In string notation this means appending the insertion digit to the
two-digit final position: `1.41 -> 1.411`, `G.HN.03 -> G.HN.031`.

## Primary annotation algorithm

### 1. Validate inputs and load the table

```text
ANNOTATE_SEQUENCES(queries, reference_table_name, options):
    if queries is empty:
        raise "No sequences provided"

    reference_table = load bundled CSV by name
    normalize only recognized legacy numeric labels
    preserve valid named segment labels exactly

    (reference_sequences, residue_to_grn_maps) =
        PREPARE_REFERENCE_SEQUENCES(reference_table)

    if no row contains usable residue-position cells:
        raise "no usable reference sequences"

    if search mode is not pairwise:
        raise unsupported-search error

    aligner = global BLOSUM62 aligner configured with:
        internal gap-open penalty
        internal gap-extension penalty
        terminal gap penalty

    sequence_grn_order = INFER_REFERENCE_GRN_ORDER(reference_table)
```

Gap parameters are caller-visible. Defaults are `-10.0` to open an internal
gap, `-0.05` to extend one, and `0.0` for terminal overhangs. The earlier
wrapper accidentally passed its argument as a match score; tests now assert
the actual Biopython internal and terminal gap attributes.

### 2. Reconstruct each reference in sequence order

Physical CSV column order is not trusted. Some older GPCR tables place the
`0.*` segment after `9.*`, even though its residues occur at the N terminus.

```text
PREPARE_REFERENCE_SEQUENCES(table):
    for each reference row:
        mapped = empty list

        for each (grn, cell) in row:
            if cell matches RESIDUE + POSITIVE_INTEGER:
                append (sequence_position, residue, grn) to mapped

        sort mapped by sequence_position
        reject duplicate sequence_position values
        skip rows with no mapped cells

        reference_sequence[row] = concatenate mapped residues
        residue_to_grn_map[row] = mapped GRNs in the same order

    return both dictionaries
```

The reconstructed sequence is the annotated portion of the protein. Unannotated
N/C termini are not invented. When a full UniProt sequence is aligned, those
regions therefore become terminal overhangs rather than internal insertions.

### 3. Infer intrinsic segment and GRN order

Named segments cannot be ordered alphabetically, and numeric conversion would
destroy labels such as `G.HN`. Order is inferred from sequence positions across
the reference table.

```text
INFER_REFERENCE_GRN_ORDER(table):
    for each GRN column:
        split at final dot -> (segment_id, integer_position)
        collect all sequence positions present in that column
        column_sequence_location = median(collected positions), if any

    for each segment_id:
        segment_location = median(locations of its populated columns)
        segment_fallback = index of its first physical table column
        segment_direction =
            ascending if larger GRN positions occur later in sequence
            descending if larger GRN positions occur earlier in sequence
            ascending if insufficient evidence

    order segments by (segment_location, segment_fallback)

    within each segment:
        order every GRN by its final position in the inferred direction
        treat insertion 1.411 or G.HN.031 as locally following its
        two-digit base coordinate
        place 12.xxx / 21.xxx between ordinal segments 1 and 2;
        order the second direction in reverse distance
        therefore all-gap columns inherit the same local segment order

    retain unparseable legacy columns at stable fallback positions
    return ordered GRN labels
```

This order supplies insertion candidates and segment ordinals. Distributed
reference CSVs are not rewritten; newly required insertion/flexible columns
are added to the returned annotation table next to their preceding anchor.

### 4. Sanitize each query

```text
SANITIZE(query):
    require a string
    uppercase it
    remove whitespace, periods, and existing gap characters
    allow standard amino acids plus B, J, Z, X, and stop

    if another symbol remains:
        return an all-gap annotation with status "invalid_sequence"
    if the result is empty:
        return an all-gap annotation with status "empty_sequence"
```

Pre-aligned input is intentionally treated as an ordinary sequence; ProtOS
computes its own alignment using the configured scoring model.

### 5. Select a reference

```text
FIND_BEST_REFERENCE(query, references):
    if query exactly equals a reconstructed reference:
        align once and return that reference

    best_score = negative infinity
    for each non-empty reference:
        alignment = global_pairwise_align(query, reference, BLOSUM62, gap options)
        normalized_score = raw_alignment_score /
                           max(query_length, reference_length, 1)
        retain the largest normalized_score

    return reference id, alignment, normalized score
```

The normalized score is suitable for ranking within this implementation, but
is not a probability and can exceed 1 because BLOSUM62 match scores exceed 1.
No universal scientific cutoff is assumed. Callers may provide a minimum score
and/or coverage; rejected results retain diagnostics but their output row is
replaced by gaps.

### 6. Project GRNs and correct sequence positions through gaps

The formatted alignment contains query, match, and reference strings. Query and
reference indices advance independently.

```text
PROJECT_ALIGNMENT(query_alignment, reference_alignment, reference_grns):
    output every table GRN as "-"
    query_index = 0
    reference_index = 0

    for each aligned character pair (query_char, reference_char):
        if reference_char is not a gap:
            current_grn = reference_grns[reference_index]
        else:
            current_grn = none

        if query_char is not a gap:
            query_index += 1

        if both characters are residues and current_grn exists:
            output[current_grn] = query_char + query_index

        if query is residue and reference is gap:
            extend/create an insertion event

        if query is gap and reference is residue:
            extend/create a deletion event containing current_grn

        if reference_char is not a gap:
            reference_index += 1
```

This is the position correction: after a three-residue query insertion,
downstream GRNs use sequence positions shifted by +3; after a deletion they
shift by -1. No post-hoc renumbering is needed.

### 7. Distinguish alignment errors, segment insertions, and flexible regions

```text
FINALIZE_UNALIGNED_QUERY_RUN(event, intrinsic_grn_order, output, config):
    event.after_grn  = closest preceding reference GRN
    event.before_grn = closest following reference GRN

    if after is missing and before exists:
        kind = n_terminal_overhang
    else if after exists and before is missing:
        kind = c_terminal_overhang
    else if both are missing:
        kind = unaligned_query
    else:
        kind = internal_insertion

    candidate_grns = currently empty table GRNs strictly between anchors

    if kind is internal_insertion
       and caller enabled assign_unambiguous_insertions
       and candidate count exactly equals inserted residue count:
        assign residues to candidates in intrinsic sequence order
        mark "assigned_exact_candidate_count"

    else if kind is internal_insertion
            and run length > max_alignment_gap
            and both anchors lie in the same configured STRICT region:
        standard_extent = matching STANDARD region, else strict region
        shift each projected GRN from before_grn through standard_extent.end
              left by run length in the true query sequence
        replace both residue letter and sequence position from the query
        mark "alignment_gap_compressed"

    else if kind is internal_insertion
            and both anchors belong to the same segment:
        anchor = preceding two-digit ordinary coordinate
        generate anchor+.001, anchor+.002, ...
        assign inserted residues and add those columns to the output

    else if kind is internal_insertion
            and anchors belong to different ordered segments:
        n_count = ceil(inserted_count / 2)
        c_count = floor(inserted_count / 2)
        label the first n_count residues from the N-side segment:
            12.001, 12.002, ...
        label the last c_count residues from the C-side segment in sequence order:
            ..., 21.002, 21.001
        add those columns to the output

    else:
        leave residues unnumbered and report why
```

Insertion annotation is enabled by default. Existing exact-fit columns take
precedence because a family table may already define the coordinates. A
non-exact candidate set is retained in diagnostics but does not block creation
of the coordinate required by the rules above.

`strict` and `standard` are deliberately different. `strict` describes the
minimum trusted non-flexible core and is used to decide whether a long gap is an
alignment error. `standard` describes the usual full non-flexible extent and
sets the boundary through which corrected assignments may be shifted. A long
run without two strict anchors is never compressed: inside one segment it is a
biological insertion, and between segments it is a legitimate flexible region.
Runs of one residue are retained as insertions even between strict anchors.

### 8. Deletions, coverage, and result status

```text
FOR EACH DELETION RUN:
    record deleted GRNs, reference residues, and preceding query position
    leave every deleted GRN as "-"

coverage = assigned projectable GRNs / projectable GRNs

if score or coverage is below caller threshold:
    status = below_threshold
    return an all-gap row while retaining alignment/indel diagnostics
else if coverage > 0:
    status = ok
else:
    status = no_overlap
```

Per-query metadata includes the selected reference, normalized score, coverage,
assigned count, insertion/deletion event lists, true internal insertion count,
deletion count, and terminal-overhang count.

## Parameterized numeric expansion path

`assign_grns_to_sequences` is separate from the public processor path. For
numeric GPCR-like labels it can attempt to annotate positions absent from the
selected reference.

```text
ASSIGN_GRNS_WITH_EXPANSION(queries, table, family_config):
    sanitize queries
    reconstruct references by cell sequence position
    optionally use MMseqs to choose a candidate, then perform pairwise alignment
    otherwise score all references with Biopython

    project aligned reference GRNs into a seed row
    determine strict GRNs from caller input or family configuration
    filter seed row to strict anchors

    if every seed label is numeric thesis-era segment/flexible/tail notation
       and none uses a structure-corrected insertion/bulge encoding:
        EXPAND_NUMERIC_ANNOTATION(seed, query, family_config)
    else:
        bypass all float arithmetic
        return segment projection with expansion_method="segment_projection"
```

The numeric expansion routine then approximately does:

```text
EXPAND_NUMERIC_ANNOTATION(seed, query, config):
    invert seed mapping to residue-number -> GRN
    enumerate every query residue number
    generate allowed standard GRNs from configured numeric intervals

    fill N-terminal residues before first anchor
    fill C-terminal residues after last anchor
    find unassigned residues between anchors
    assign missing standard GRNs where flanking distances are consistent
    classify remaining residues as N-side flexible region, gap, or C-side region
    generate directional labels from neighboring numeric segments
    sort numeric segment/flexible/tail GRNs and report still-unassigned residues
```

That routine remains useful for thesis-era GPCR/opsin configurations, but it is
not applied to `G.HN.03`, CAN, arbitrary named segments, or GPCR bulge labels
such as `1.411`. Its float arithmetic cannot represent insertion encodings
without loss and cannot be generalized merely by relaxing a regular expression.

## Known exceptions and handling

| Case | Current handling |
| --- | --- |
| Exact reference or known wild type | Projects all mapped positions; exact match avoids all-vs-all scanning |
| Arbitrary string segment | Preserved; ordered from reference sequence evidence |
| Query insertion inside a segment | Reuses exact-fit columns or creates `position + .001` columns |
| Residues between ordered segments | Split between directional N-side and C-side coordinates such as `12.001` and `21.001` |
| Query deletion | Deleted GRNs remain `-`; event records the run and downstream positions are corrected |
| Long unannotated N/C terminus | Reported as terminal overhang, not an internal insertion |
| Gap longer than configured tolerance inside strict core | Treated as an alignment error and compressed through the standard extent |
| Long gap without two strict anchors | Never compressed; treated as an insertion or flexible-region run |
| Empty reference table | Explicit error; current human Class D1 table is empty because D1 is fungal |
| Invalid/empty query | All-gap row with `invalid_sequence` or `empty_sequence` status |
| Distant sequence | Optional score/coverage thresholds can reject it; benchmark-derived defaults are still needed |
| More than nine insertions after one coordinate | Left unassigned: compact `.001` notation is ambiguous beyond one insertion digit |
| More than nine ordered segments | Directional compact labels cannot encode the ordinal; explicit table coordinates are required |
| Table with duplicate sequence positions in one row | Rejected as corrupt rather than silently scrambling the reference |
| Custom/design tables without UniProt accessions | Tested by round-tripping their stored real member sequences |

## Test coverage

`tests/test_grn_annotation.py` covers notation compatibility, segment ordering,
scrambled columns, configurable gap penalties, insertions, deletions, terminal
overhangs, exact and ambiguous insertion correction, invalid inputs, thresholds,
the numeric-expansion safety boundary, manifest/checksum coverage, one stored
member from every bundled table, and accession-traceable frozen full-length
UniProt fixtures for every biological table.

`tests/test_grn_uniprot_integration.py` uses ProtOS's own `SequenceLoader` to
fetch, register, and annotate live UniProt entries. It covers GPCR A, B1, B2, C,
F, O1, O2, T2 and unclassified tables; aggregate/core GPCR tables; every human
G-alpha family and the aggregate CGN table; and arrestin CAN. Network tests are
opt-in with `PROTOS_RUN_NETWORK_TESTS=1` so ordinary CI remains deterministic.
