# How Discriminase works (step by step)

A plain-language, no-skipping walkthrough of how a guide goes from raw DNA to a
"keep it / drop it" decision. It uses a deliberately tiny toy size so every number
fits on screen; the **real sizes** are listed at the end. Nothing here is skipped.

Throughout the demo: **guide = 4 letters, seed = first 2 letters.**

---

## Step 0 — the pieces and their lengths

On the DNA, every candidate looks like this:

```
… TTT · [ guide, seed first ] …
```

| Piece | What it is | Real default | In this demo |
|---|---|---|---|
| PAM | a fixed motif the nuclease needs, sitting **5′** (before) the guide | `TTT` (3 letters) | `TTT` |
| **guide** | the part that matches the DNA you want to cut | **23 letters** (max 31) | **4 letters** |
| seed | the guide letters nearest the PAM — where CRISPR is most sensitive | first **10** | first **2** |

The guide sits **immediately next to** the PAM. If a nuclease needs slack between
the two, you express it inside the PAM with `N` (e.g. `TTTN` = `TTT` then any one
base), so there is no separate "gap" knob.

The PAM is Cas12a-style (5′ `TTT`), **not** SpCas9 (`NGG`, 3′). All of these are
configurable in `config.py`.

---

## Step 1 — walk the genome, find a PAM, cut out the guide

Slide along one strand, left to right. At each position ask: *"are the next 3 letters
`TTT`?"*

```
position:  0   1   2   3   4   5   6
letter:    T   T   T   C   A   G   T
           └── PAM ──┘   └──── guide ────┘
```

- `TTT` is found at positions 0–2 → this is a PAM site.
- Take the next **4** letters (positions 3–6) = **`C A G T`** → one guide.
- Its **seed** = first 2 letters = **`C A`**.

Do this at every position, then again on the **reverse-complement strand** (to catch
PAMs on the other strand). A ~5 Mbp genome yields hundreds of thousands of guides.

> For **commensals** we keep *all* extracted guides (they are the cut-sites to avoid).
> For the **target** we also remember each guide's position and strand, so the output
> can report where it is.

---

## Step 2 — turn each guide into ONE number

There are 4 letters, so each takes **2 bits** (2 bits = 4 combinations; 1 bit would
only give 2):

```
A = 00      C = 01      G = 10      T = 11
```

Pack the guide `C A G T` left to right, **seed first**:

```
C    A    G    T
01   00   10   11      →   0100 1011
```

That bit-string is just a number. Read it with powers-of-two place values:

```
place: 128  64  32  16   8   4   2   1
bit:     0   1   0   0   1   0   1   1
                64        + 8     + 2 + 1   =  75
```

So `CAGT` → **75**. (The value "75" doesn't matter on its own; what matters is that
every guide gets a number and we can always convert back.)

**The crucial choice: the seed goes on the LEFT (the high-value bits).** In any number
the left digits dominate its size — in `523` the `5` (hundreds) matters most. So **two
guides with the same seed become numbers in the same range.** For `CAGT`, seed `CA` =
`0100` = the top 4 bits, so *every* `CA`-seed guide is a number in `[64, 80)`.

This single idea is what lets us replace a trie + BK-tree with one sorted list.

---

## Step 3 — sort all the commensal numbers

Convert every commensal guide to a number, then **sort ascending**. Example panel of 8:

```
index:   0      1      2      3      4      5      6      7
guide:  AACC   CAAA   CAGT   CATT   CCCC   GGGG   TACA   TTTT
number:   5     64     75     79     85    170    196    255
seed:    AA     CA     CA     CA     CC     GG     TA     TT
                └──── seed "CA" ────┘
```

Sorting did the hard part for free: **all three `CA`-seed guides (indices 1–3) are now
contiguous**, because their numbers all fall in `[64, 80)`. Finding "everything with
this seed" (the old trie's job) is now just "everything in this number range."

This sorted array of numbers **is** the database. It's saved as a plain NumPy file and
memory-mapped, so loading is instant and uses almost no RAM.

---

## Step 4 — binary search (lo / mid / hi, in full)

A **target** guide arrives: `CAGA`.
Convert it: `C A G A` = `01 00 10 00` = **72**, seed `CA`.

We want every commensal with seed `CA` — i.e. every number in **`[64, 80)`**. Binary
search finds the two edges of that block.

**What `lo`, `hi`, `mid` mean** — two fingers on the array:
- `lo` = left finger (start of the part still worth searching),
- `hi` = right finger (one *past* the end of it),
- `mid` = the index halfway between them.

Each step: look at the value under `mid`. If it's **too small**, the answer is to the
**right**, so move `lo` past `mid`. If it's **big enough**, the answer is at `mid` or to
its **left**, so pull `hi` down to `mid`. Every step discards half the array.

### Find the LEFT edge — first index whose value ≥ 64

```
array:  [ 5,  64,  75,  79,  85, 170, 196, 255 ]
index:    0    1    2    3    4    5    6    7

start:  lo=0, hi=8
step 1: mid=(0+8)/2=4 → value 85 ≥ 64 → pull hi down → hi=4
step 2: mid=(0+4)/2=2 → value 75 ≥ 64 → pull hi down → hi=2
step 3: mid=(0+2)/2=1 → value 64 ≥ 64 → pull hi down → hi=1
step 4: mid=(0+1)/2=0 → value  5 < 64 → slide lo up   → lo=1
        lo(1) == hi(1) → fingers meet → LEFT edge = index 1
```

### Find the RIGHT edge — first index whose value ≥ 80

```
start:  lo=0, hi=8
step 1: mid=4 → value 85 ≥ 80 → hi=4
step 2: mid=2 → value 75 < 80 → lo=3
step 3: mid=(3+4)/2=3 → value 79 < 80 → lo=4
        lo(4) == hi(4) → meet → RIGHT edge = index 4
```

**Block = indices [1, 4)** = `CAAA, CAGT, CATT`. Each search took ~3 steps because
log₂(8) = 3. For 8 **million** guides it would be ~23 steps — you *leap* to the right
block instead of walking a tree node by node.

---

## Step 5 — verify the block (count mismatches)

Compare the target `CAGA` to each guide in that little block, counting differing
letters (the **Hamming distance**):

```
CAGA  vs  CAAA :  C-C  A-A  G≠A  A-A   → 1 mismatch
CAGA  vs  CAGT :  C-C  A-A  G-G  A≠T   → 1 mismatch
CAGA  vs  CATT :  C-C  A-A  G≠T  A≠T   → 2 mismatches
```

(In code this is one CPU instruction — XOR the two numbers, then popcount — but it is
exactly the table above.)

**Decision.** With tolerance "drop if any commensal is within 1 mismatch": `CAGA`
collides with `CAAA` and `CAGT`, so it would also cut a commensal → **drop it.** A guide
**survives** (is "commensal-sparing") only if its seed-block is empty *or* every guide
in it is farther than the tolerance.

> Refinement the real tool adds: if you allow the *seed itself* to differ by 1 letter,
> it simply repeats Steps 4–5 for the handful of seeds one letter away from `CA`. Same
> mechanics, a few more block lookups. This is the `seed_max_mismatch` setting.

---

## Same machinery at real size

Nothing changes but the magnitudes:

| | Demo | Real default |
|---|---|---|
| guide length | 4 | 23 (a 46-bit number — still one 64-bit integer) |
| seed length | 2 | 10 |
| commensal array | 8 numbers | ~8 million numbers (~23 binary-search steps) |
| seed block to verify | 3 guides | a handful |

That's why a whole target screens in seconds and the index loads instantly.

## Why a sorted list instead of the old trie + BK-tree?

Same good idea, better realized. A sorted array of fixed-length keys **is** a flattened
trie — binary search examines prefixes (seeds) implicitly — but it lives in one
contiguous block of memory (8 bytes per guide) instead of millions of pointer-linked
objects (~50–100 bytes each). That's ~20× smaller, cache-friendly, memory-mappable, and
it replaces the BK-tree's unpredictable neighbour search (which prunes poorly for short,
4-letter, fixed-length DNA) with a binary-search-then-tiny-verify that can't blow up.

Rough cost, with `N` commensal guides, `M` target guides, seed block size `B`:

| | This method | Old trie + BK-tree |
|---|---|---|
| build | `O(N log N)` (sort) | `O(N·L)` but huge memory constant |
| one query | `O(log N + B)` per seed checked | trie `O(L)` + BK-tree up to `O(N)`, cache-missing |
| memory | `~10·N` bytes, mmap'd | `~50–100·N` bytes, all resident |

Full decision log and measured numbers: [ARCHITECTURE.md](ARCHITECTURE.md). How it
compares to Cas-OFFinder / GuideScan2: [COMPARISON.md](COMPARISON.md).
