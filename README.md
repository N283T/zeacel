# zeasel

A Zig reimplementation of [Easel](https://github.com/EddyRivasLab/easel), the C library for biological sequence analysis from the Eddy/Rivas laboratory. Easel underpins [HMMER](http://hmmer.org) and [Infernal](http://eddylab.org/infernal/). zeasel serves as the foundation for a Zig reimplementation of HMMER.

## Features

- **Biological alphabets** — DNA, RNA, amino acid with comptime lookup tables, degeneracy resolution, scoring/counting
- **Sequence I/O** — FASTA, Stockholm, GenBank, EMBL, Clustal, aligned FASTA, PHYLIP, A2M, PSI-BLAST, SELEX
- **Multiple sequence alignment** — MSA operations, weighting (PB/BLOSUM), clustering, gap manipulation, metadata
- **SIMD vector operations** — comptime-generic `@Vector`-based ops with automatic ISA selection
- **Statistical distributions** — Gumbel, gamma, exponential, normal, Weibull, Dirichlet, mixture Dirichlet
- **Statistical functions** — logGamma, digamma, chi-squared test, G-test, linear regression
- **Numerical optimization** — minimizer (conjugate gradient), root finder (Brent's method)
- **Score matrices** — BLOSUM45/62/80, PAM30/70/120/250 with named lookup and probability analysis
- **Evolutionary models** — rate matrix Q, probability matrix P=exp(tQ), PAML model reader
- **Genetic code** — codon translation with NCBI genetic code tables
- **Sequence distance** — pairwise identity and Jukes-Cantor correction
- **Phylogenetic trees** — UPGMA, single/complete linkage, WPGMA, Newick read/write
- **SSI index** — sequence/subsequence index with aliases, multi-file, subseq positioning
- **dsqdata** — binary packed sequence database with threaded prefetch reader
- **Threading** — thread pool, work queue (producer-consumer), CPU detection
- **Data structures** — red-black tree, Huffman coding, varint encoding, bipartite matching
- **WUSS notation** — RNA secondary structure parsing and validation

## Requirements

- [Zig](https://ziglang.org/) 0.15.0 or later

## Build

```bash
# Build the library
zig build

# Run all tests
zig build test

# Build CLI tools
zig build tools
```

## CLI Tools

| Tool | Description | Easel equivalent |
|------|-------------|------------------|
| `zeasel-seqstat` | Report sequence file statistics | `esl-seqstat` |
| `zeasel-reformat` | Convert between sequence formats | `esl-reformat` |
| `zeasel-seqfetch` | Fetch sequences by name | `esl-seqfetch` |

```bash
# After `zig build tools`, binaries are in zig-out/bin/
./zig-out/bin/zeasel-seqstat sequences.fasta
./zig-out/bin/zeasel-reformat stockholm input.sto
./zig-out/bin/zeasel-seqfetch sequences.fasta seq1 seq2
```

## Usage as a Library

Add zeasel as a Zig dependency, then import:

```zig
const zeasel = @import("zeasel");

// Digitize a DNA sequence
const abc = &zeasel.alphabet.dna;
const dsq = try abc.digitize(allocator, "ACGTACGT");
defer allocator.free(dsq);

// Read sequences from a FASTA file
var reader = try zeasel.io.Reader.fromMemory(allocator, abc, file_data, .fasta);
defer reader.deinit();
const sequences = try reader.readAll();
```

## Project Structure

```
src/
├── root.zig             # Public API re-exports
├── alphabet.zig         # Biological alphabets + degeneracy tables
├── sequence.zig         # Sequence type + SequenceBlock
├── io.zig               # I/O module hub
│   └── io/              # Format-specific parsers/writers
│       ├── reader.zig, writer.zig
│       ├── fasta.zig, stockholm.zig, genbank.zig
│       ├── clustal.zig, afa.zig, phylip.zig
│       ├── a2m.zig, psiblast.zig, selex.zig
├── msa.zig              # Multiple sequence alignment
├── msa_ops.zig          # MSA clustering and manipulation
├── msa_weight.zig       # MSA sequence weighting
├── score_matrix.zig     # Substitution matrices (BLOSUM, PAM)
├── composition.zig      # Background residue frequencies
├── distance.zig         # Sequence distance calculations
├── tree.zig             # Phylogenetic trees (UPGMA, Newick)
├── matrix.zig           # Matrix operations + LU/inverse/exp
├── ratematrix.zig       # Evolutionary rate matrices
├── paml.zig             # PAML model reader
├── mixdchlet.zig        # Mixture Dirichlet priors
├── ssi.zig              # Sequence/Subsequence Index
├── dsqdata.zig          # Binary sequence database + prefetch
├── gencode.zig          # Genetic code / translation
├── wuss.zig             # RNA secondary structure (WUSS)
├── simd.zig             # SIMD vector operations
│   └── simd/vector_ops.zig
├── stats.zig            # Statistical distributions hub
│   └── stats/           # Distribution modules + functions
├── threads.zig          # Thread pool + work queue
├── cpu.zig              # CPU/SIMD feature detection
├── cluster.zig          # Generalized clustering
├── red_black.zig        # Red-black tree
├── huffman.zig          # Huffman coding
├── varint.zig           # Variable-length integer encoding
├── graph.zig            # Bipartite matching
├── iset.zig             # Independent set splitting
├── recorder.zig         # I/O recorder with rewind
├── util.zig             # Utilities (RNG, random sequences)
└── tools/               # CLI tools
    ├── seqstat.zig
    ├── reformat.zig
    └── seqfetch.zig
```

## Easel Module Coverage

Full comparison of Easel's C modules and their zeasel status.

### Implemented in zeasel

| Easel module | zeasel module | Notes |
|---|---|---|
| `esl_alphabet` | `alphabet.zig` | Comptime alphabets + degeneracy tables, scoring/counting functions |
| `esl_sq` | `sequence.zig` | Digital sequences + SequenceBlock for threaded pipelines |
| `esl_msa` | `msa.zig` | Full MSA operations: subset, gaps, RC, metadata, validation |
| `esl_msaweight` | `msa_weight.zig` | PB weights, BLOSUM weights |
| `esl_msacluster` | `msa_ops.zig` | Single-linkage clustering by % identity |
| `esl_msashuffle` | `msa_ops.zig` | Column shuffling (Fisher-Yates) |
| `esl_scorematrix` | `score_matrix.zig` | BLOSUM45/62/80, PAM30/70/120/250, named lookup, entropy |
| `esl_composition` | `composition.zig` | BL62, WAG, SW34, SW50 background frequencies |
| `esl_distance` | `distance.zig` | Pairwise identity, Jukes-Cantor correction |
| `esl_tree` | `tree.zig` | UPGMA, single/complete linkage, WPGMA, Newick read/write |
| `esl_dmatrix` | `matrix.zig` | Dense matrix + LU decomposition, inverse, exp(tQ) |
| `esl_ssi` | `ssi.zig` | SSI index with aliases, multi-file, subseq positioning |
| `esl_dsqdata` | `dsqdata.zig` | Binary sequence DB + threaded PrefetchReader |
| `esl_gencode` | `gencode.zig` | NCBI genetic code tables, codon translation |
| `esl_wuss` | `wuss.zig` | RNA secondary structure WUSS notation |
| `esl_cluster` | `cluster.zig` | Generalized single-linkage clustering |
| `esl_vectorops` | `simd/vector_ops.zig` | Comptime-generic SIMD via `@Vector` |
| `esl_gumbel` | `stats/gumbel.zig` | Gumbel distribution (E-value core) |
| `esl_exponential` | `stats/exponential.zig` | Exponential distribution |
| `esl_gamma` | `stats/gamma.zig` | Gamma distribution + incomplete gamma |
| `esl_normal` | `stats/normal.zig` | Normal distribution |
| `esl_weibull` | `stats/weibull.zig` | Weibull distribution |
| `esl_dirichlet` | `stats/dirichlet.zig` | Dirichlet distribution |
| `esl_histogram` | `stats/histogram.zig` | Histogram binning |
| `esl_minimizer` | `stats/minimizer.zig` | Conjugate gradient minimizer |
| `esl_rootfinder` | `stats/rootfinder.zig` | Brent's method root finder |
| `esl_stats` | `stats/functions.zig` | logGamma, psi, trigamma, erfc, chi-squared, G-test, linear regression |
| `esl_mixdchlet` | `mixdchlet.zig` | Mixture Dirichlet priors for hmmbuild |
| `esl_ratematrix` | `ratematrix.zig` | Rate matrix Q, P=exp(tQ), normalization |
| `esl_paml` | `paml.zig` | PAML exchangeability matrix reader |
| `esl_red_black` | `red_black.zig` | Red-black tree (HMMER sparse DP) |
| `esl_huffman` | `huffman.zig` | Huffman coding for compressed sequences |
| `esl_varint` | `varint.zig` | Google varint (base-128) encoding |
| `esl_graph` | `graph.zig` | Maximum bipartite matching |
| `esl_iset` | `iset.zig` | Independent set splitting (training/test) |
| `esl_recorder` | `recorder.zig` | Line-based I/O with rewind/replay |
| `esl_cpu` | `cpu.zig` | SIMD feature detection, CPU core count |
| `esl_threads` | `threads.zig` | Thread pool + WorkQueue |
| `esl_random` | `util/random.zig` | Mersenne Twister RNG |
| `esl_randomseq` | `util/random_seq.zig` | Random sequence generation |
| `esl_sqio_ascii` | `io/fasta.zig` | FASTA/EMBL sequence I/O |
| `esl_msafile_stockholm` | `io/stockholm.zig` | Stockholm format |
| `esl_msafile_clustal` | `io/clustal.zig` | Clustal format |
| `esl_msafile_afa` | `io/afa.zig` | Aligned FASTA |
| `esl_msafile_phylip` | `io/phylip.zig` | PHYLIP interleaved/sequential |
| `esl_msafile_a2m` | `io/a2m.zig` | A2M (UCSC SAM) format |
| `esl_msafile_psiblast` | `io/psiblast.zig` | PSI-BLAST format |
| `esl_msafile_selex` | `io/selex.zig` | SELEX format |

### Not needed (replaced by Zig standard library or language features)

| Easel module | Zig replacement | Why not needed |
|---|---|---|
| `esl_alloc` | `std.mem.Allocator` | Zig allocator interface supports alignment natively |
| `esl_arr2`, `esl_arr3` | Native slices | Zig multi-dimensional arrays are first-class |
| `esl_avx`, `esl_sse`, `esl_neon`, `esl_vmx`, `esl_avx512` | `@Vector` built-in | Zig comptime SIMD replaces per-ISA C files |
| `esl_bitfield` | `std.StaticBitSet` | Standard library bit set |
| `esl_buffer` | `[]const u8` slices | Zig slice-based parsers don't need buffered I/O abstraction |
| `esl_fileparser` | `std.mem.tokenizeAny` | Standard library tokenizer |
| `esl_getopts` | `std.process.argsAlloc` | Zig arg parsing or third-party libraries |
| `esl_heap` | `std.PriorityQueue` | Standard library priority queue |
| `esl_json` | `std.json` | Standard library JSON parser |
| `esl_keyhash` | `std.StringHashMap` | Standard library hash map |
| `esl_mem` | Zig slice operations | `std.mem.eql`, `std.mem.indexOf`, etc. |
| `esl_quicksort` | `std.mem.sort` | Standard library sort |
| `esl_rand64` | `std.Random` | Standard library RNG |
| `esl_regexp` | `std.mem` string ops | Pattern matching via standard library |
| `esl_stack` | `std.ArrayList` | Standard library dynamic array as stack |
| `esl_stopwatch` | `std.time.Timer` | Standard library timer |
| `esl_subcmd` | Direct dispatch | Zig switch-based subcommand dispatch |
| `esl_workqueue` | `threads.WorkQueue` | Integrated into threads module |

### Not applicable

| Easel module | Reason |
|---|---|
| `esl_mpi` | MPI cluster parallelism — out of scope |
| `esl_msafile2` | Legacy Pfam Stockholm reader — superseded |
| `esl_sqio_ncbi` | NCBI BLAST database format — specialized |
| `esl_swat` | Smith-Waterman — marked UNFINISHED in Easel |
| `esl_hmm` | General discrete HMMs — not profile HMMs (HMMER has its own) |
| `esl_gev`, `esl_mixgev` | Generalized extreme value — research only |
| `esl_hyperexp`, `esl_stretchexp`, `esl_lognormal` | Niche distributions — not used by HMMER core |
| `esl_matrixops` | C array helpers — replaced by `matrix.zig` |
| `interface_gsl`, `interface_lapack` | External C library wrappers — pure Zig implementations instead |

## Design Decisions

| Aspect | Easel (C) | zeasel (Zig) |
|--------|-----------|--------------|
| Digital sequences | `ESL_DSQ` (uint8_t) + sentinel bytes | `[]u8` slice, no sentinels |
| Format dispatch | Function pointer vtable | Tagged union + switch |
| Memory management | malloc/realloc + manual free | `std.mem.Allocator` interface |
| Error handling | Macros (`ESL_FAIL`) | Zig error unions (`!T`) |
| Type-specific functions | Separate `DSet/FSet/ISet` | `comptime` generics |
| SIMD | Separate files per ISA | `@Vector` + `suggestVectorLength` |
| Hash tables | Custom `ESL_KEYHASH` | `std.StringHashMap` |

## License

See [LICENSE](LICENSE).
