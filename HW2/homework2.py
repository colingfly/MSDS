# Part 1: sequence analysis
seq = (
"ATTGGGGAGGAGGCGAGTTGAGCGGCGGCAGTTCGCCTGCGTGCGCTGCGCGG"
"CGTCGACATCTGATCCGCACCATGGAAATCCCCGCTCAATCTTTGGAGCAGGGAT"
"GCGGGGCGATCAAGATGGGGATGCGGGATGGGGGCGACGGTGTATTTCCGCCAG"
"AAGATTTCGCCGCGGGAGCTCGCGGTGCGTACGTGCATGTTCAAACGCACGGTG"
"CGCGCATGGCAGTGGCAGACTGATCAACGCAGCTGGAAGCATCCGAAGCGCGCG"
"GGCACGCGTGTCCTCGACGCGTGGCCTCACATGCTGTCGGGTCGGTTCAAGACC"
"GAAAGCCACCGACCGACGCGCGAGCAATGCGCTACGCGGATCGCGTTCGACACG"
"AGCCGCGCGCGAGGCAAGGCCGACGTATTCGATCTTCCAGAGGAAGCCTATTGG"
"CTCGAGTCGTAGTGCTCGATATGGTAGAGCAACATGAATCCCGGGCTAAGTACAA"
"GAAGTAACCCGGCAACGAGTGAGATTGCGACGAATAAACGCTTCACCATGATCGC"
"GCTCCTGAGTTGGTTGAGGTGAATTGGAAAGTCGATTCCTGGGGGATCATTCCCG"
"GCAAGGCGCGCAATCCCCGCATTGTTCTCAAGATCGCAACGCGATTCGTCAGGCC"
"GATCTTCATGGGGTGTCTCGCTGGTAGTGATTCCGTCGTGGCCCGCGCATGTGCA"
"TGACGGCATCCGGGGAG"
)

seq = "".join([c for c in seq.upper() if c in "ACGT"])  # keep A/C/G/T only
length = len(seq)

# non-overlapping triplets starting from index 0
triplets = []
for i in range(0, len(seq) - len(seq)%3, 3):
    triplets.append(seq[i:i+3])

counts = {}
for t in triplets:
    counts[t] = counts.get(t, 0) + 1

total_triplets = len(triplets)
fractions = {}
for t, c in counts.items():
    fractions[t] = c / total_triplets if total_triplets > 0 else 0.0

# complement (A<->T, C<->G)
comp_map = str.maketrans({"A":"T","T":"A","C":"G","G":"C"})
complement = seq.translate(comp_map)

# write text answers for part 1 (relative path)
with open("part1_answers.txt", "w") as f:
    f.write("Sequence length: " + str(length) + "\n")
    f.write("Complement:\n")
    f.write(complement + "\n")
    f.write("Triplet counts and fractions (non-overlapping, frame starting at 0):\n")
    for t in sorted(counts.keys()):
        f.write(f"{t}\t{counts[t]}\t{fractions[t]:.4f}\n")
print("Part 1 done -> part1_answers.txt")

# Part 2: parse human.fa and protein-coding_gene.txt to annotated fasta
import re

class Entry:
    def __init__(self, name, seq):
        self.name = name  # gene symbol
        self.seq = seq    # sequence
        self.synonyms = ""
        self.location = ""

def read_fasta(path):
    items = []
    name = None
    parts = []
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    items.append(Entry(name, "".join(parts)))
                name = line[1:].strip()
                parts = []
            else:
                parts.append(line.strip())
        if name is not None:
            items.append(Entry(name, "".join(parts)))
    return items

def chunk60(s):
    out = []
    for i in range(0, len(s), 60):
        out.append(s[i:i+60])
    return "\n".join(out)

def looks_like_location(token):
    return bool(re.match(r'^([0-9]{1,2}|[XYxy])(p|q)', token))

def load_gene_info(path):
    info = {}
    with open(path, "r", errors="ignore") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            symbol = cols[1].strip()
            loc_idx = -1
            for i, tok in enumerate(cols):
                if looks_like_location(tok.strip()):
                    loc_idx = i
                    break
            synonyms = ""
            location = ""
            if loc_idx != -1:
                location = cols[loc_idx].strip()
                if loc_idx - 1 >= 0:
                    synonyms = cols[loc_idx - 1].strip().strip('"')
            else:
                last = cols[-1].strip()
                if looks_like_location(last):
                    location = last
                for c in cols:
                    if "," in c:
                        synonyms = c.strip().strip('"')
                        break
            syn_list = [s.strip() for s in re.split(r'[;,]', synonyms) if s.strip()]
            synonyms_clean = ",".join(syn_list)
            info[symbol] = (synonyms_clean, location)
    return info

def write_fasta(items, info, out_path):
    with open(out_path, "w") as out:
        for e in items:
            syn, loc = info.get(e.name, ("",""))
            e.synonyms = syn
            e.location = loc
            header = f">{e.name}|{e.synonyms}|{e.location}"
            out.write(header + "\n")
            out.write(chunk60(e.seq) + "\n")

human = read_fasta("human.fa")
gene_info = load_gene_info("protein-coding_gene.txt")
write_fasta(human, gene_info, "human_annotated.fa")  # relative path
print("Part 2 done -> human_annotated.fa")

# Part 3: ORF scan (forward frames only)
stops = {"TAA","TAG","TGA"}

def parse_fasta_ids(path):
    seqs = {}
    order = []
    with open(path, "r", errors="ignore") as f:
        name = None
        parts = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts)
                    order.append(name)
                head = line[1:].strip()
                token = head.split()[0]  # first token as ID
                name = token
                parts = []
            else:
                parts.append(line.upper())
        if name is not None:
            seqs[name] = "".join(parts)
            order.append(name)
    return order, seqs

def longest_orf_in_frame(seq, frame):
    best = (0, -1, -1)  # length, start, end
    i = frame
    n = len(seq)
    while i+3 <= n:
        if seq[i:i+3] == "ATG":
            j = i + 3
            while j+3 <= n:
                cod = seq[j:j+3]
                if cod in stops:
                    length = j + 3 - i
                    if length > best[0]:
                        best = (length, i, j+3)
                    break
                j += 3
        i += 3
    return best

def longest_orf(seq):
    best_frame = 1
    best = (0, -1, -1)
    for f in range(3):
        res = longest_orf_in_frame(seq, f)
        if res[0] > best[0]:
            best = res
            best_frame = f+1
    return best_frame, best

order, seqs = parse_fasta_ids("dna.fasta")
lens = [(name, len(seqs[name])) for name in order]
longest_seq = max(lens, key=lambda x: x[1])
shortest_seq = min(lens, key=lambda x: x[1])

with open("orf_results.txt", "w") as f:  # relative path
    f.write("Longest sequence: " + longest_seq[0] + " (" + str(longest_seq[1]) + " bp)\n")
    f.write("Shortest sequence: " + shortest_seq[0] + " (" + str(shortest_seq[1]) + " bp)\n")
    f.write("\nLongest ORF per sequence (forward frames 1-3):\n")
    for name in order:
        frame, (length, start, end) = longest_orf(seqs[name])
        f.write(name + "\tframe=" + str(frame) + "\tlength=" + str(length) + "\n")
print("Part 3 done -> orf_results.txt")

# Combine text answers from Part 1 and Part 3 (relative in/out)
from pathlib import Path, PurePath
p1 = Path("part1_answers.txt")
p3 = Path("orf_results.txt")
out = Path("homework2_answers.txt")
out.write_text((p1.read_text() if p1.exists() else "") + "\n" + (p3.read_text() if p3.exists() else ""))
print("Combined -> homework2_answers.txt")

# Quick sanity check: show working dir and files
import os
print("cwd:", os.getcwd())
print("saved files:", [fn for fn in os.listdir() if fn.endswith(('.txt','.fa'))])
