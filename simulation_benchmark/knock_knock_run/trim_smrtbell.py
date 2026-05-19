"""
Trim PacBio SMRTbell adapter residue from read edges.
Motif detection via exact substring match; min k-mer size configurable per platform.

Usage: trim_smrtbell.py <in.fastq.gz> <out.fastq.gz> <min_k>
"""
import gzip
import sys

def rc(s): return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

SMRTBELL = "GTACTTCGTTCAGTTACGTATTGC"       # 24-bp observed residue
SMRTBELL_RC = rc(SMRTBELL)

def find_adapter_end(prefix, patt, min_k):
    """Return position just past the adapter (for 5' trim). Try longest match first."""
    best_end = 0
    for k in range(min(20, len(patt)), min_k - 1, -1):
        for i in range(len(patt) - k + 1):
            sub = patt[i:i+k]
            j = prefix.find(sub)
            if j >= 0 and j + k > best_end:
                best_end = j + k
    return best_end

def find_adapter_start(suffix, patt, min_k):
    """Return position in suffix before the adapter starts (for 3' trim)."""
    earliest = len(suffix)
    for k in range(min(20, len(patt)), min_k - 1, -1):
        for i in range(len(patt) - k + 1):
            sub = patt[i:i+k]
            j = suffix.find(sub)
            if j >= 0 and j < earliest:
                earliest = j
    return earliest

def trim_read(seq, qual, min_k, window=35):
    prefix = seq[:window]
    trim5 = max(find_adapter_end(prefix, SMRTBELL, min_k),
                find_adapter_end(prefix, SMRTBELL_RC, min_k))

    suffix = seq[-window:] if len(seq) > window else seq
    offset = len(seq) - len(suffix)
    earliest_3p_fwd = find_adapter_start(suffix, SMRTBELL, min_k)
    earliest_3p_rc  = find_adapter_start(suffix, SMRTBELL_RC, min_k)
    trim3_from = offset + min(earliest_3p_fwd, earliest_3p_rc)

    if trim5 >= trim3_from:
        return None, None  # read is entirely adapter, drop
    return seq[trim5:trim3_from], qual[trim5:trim3_from]

def main():
    if len(sys.argv) != 4:
        sys.exit("usage: trim_smrtbell.py <in.fastq.gz> <out.fastq.gz> <min_k>")
    in_fq, out_fq, min_k = sys.argv[1], sys.argv[2], int(sys.argv[3])

    n_total = n_trim5 = n_trim3 = n_dropped = 0
    bp_removed_5 = bp_removed_3 = 0

    with gzip.open(in_fq, "rt") as f, gzip.open(out_fq, "wt") as out:
        while True:
            hdr = f.readline()
            if not hdr:
                break
            seq  = f.readline().strip()
            plus = f.readline()
            qual = f.readline().strip()
            n_total += 1

            new_seq, new_qual = trim_read(seq, qual, min_k)
            if new_seq is None:
                n_dropped += 1
                continue

            orig_len = len(seq)
            new_len = len(new_seq)
            if new_len < orig_len:
                if seq.startswith(new_seq):
                    n_trim3 += 1
                    bp_removed_3 += orig_len - new_len
                elif seq.endswith(new_seq):
                    n_trim5 += 1
                    bp_removed_5 += orig_len - new_len
                else:
                    n_trim5 += 1
                    n_trim3 += 1

            out.write(hdr)
            out.write(new_seq + "\n")
            out.write(plus)
            out.write(new_qual + "\n")

    print(f"Input:            {in_fq}")
    print(f"Output:           {out_fq}")
    print(f"min_k:            {min_k}")
    print(f"Total reads:      {n_total}")
    print(f"5' trimmed:       {n_trim5} ({100*n_trim5/n_total:.1f}%)")
    print(f"3' trimmed:       {n_trim3} ({100*n_trim3/n_total:.1f}%)")
    print(f"Dropped (all adapter): {n_dropped}")
    print(f"Total 5' bp removed: {bp_removed_5}")
    print(f"Total 3' bp removed: {bp_removed_3}")

if __name__ == "__main__":
    main()
