#!/usr/bin/env python3
import csv, re, sys

if len(sys.argv) < 2:
    sys.exit("Usage: clumps_to_regions.py <clumps.tsv> [pad_bp]")

infile = sys.argv[1]
pad_bp = int(sys.argv[2]) if len(sys.argv) >= 3 else 0

def parse_sp2_positions(sp2):
    if not sp2:
        return []
    toks = re.findall(r'[^,]*?\(\d+\)', sp2)
    bps = []
    for t in toks:
        core = re.sub(r'\(\d+\)$', '', t)
        parts = core.split(':')
        if len(parts) >= 2:
            try:
                bps.append(int(parts[1]))
            except:
                pass
    return bps

with open(infile, newline='') as f, open("loci.tsv", "w", newline='') as g:
    r = csv.DictReader(f, delimiter='\t')
    w = csv.writer(g, delimiter='\t')
    w.writerow(['locus_id','chr','start','end'])

    for row in r:
        chr_ = str(row['CHR'])
        try:
            lead_bp = int(row['BP'])
        except:
            continue

        sp2 = row.get('SP2', '')
        bps = [lead_bp] + parse_sp2_positions(sp2)
        start = min(bps)
        end   = max(bps)

        if pad_bp > 0:
            start = max(0, start - pad_bp)
            end   = end + pad_bp

        locus_id = f"{chr_}_{lead_bp}"
        w.writerow([locus_id, chr_, start, end])