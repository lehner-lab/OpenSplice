import gzip
import sys
import os

if len(sys.argv) < 2:
    print("Usage: python3 Clinvar_reshape.py3 /path/to/output_dir")
    sys.exit(1)

OUTDIR = sys.argv[1]
os.makedirs(OUTDIR, exist_ok=True)

INPUT = os.path.join(OUTDIR, "clinvar_20230702.vcf.gz")
OUTPUT = os.path.join(OUTDIR, "clinvar.tsv")

print(f"Reading: {INPUT}")
print(f"Writing: {OUTPUT}")

clinvar = gzip.open(INPUT, 'rb')

chr = []
pos = []
clnsig = []
mc = []
clnvc = []
id = []
gene = []
clnhgvs = []

for line in clinvar:
    line = line.decode('utf-8')
    if line[0] == '#':
        continue
    fields = line.split('\t')
    chr.append(fields[0])
    pos.append(fields[1])
    tags = fields[7].split(';')
    for t in tags:
        if 'CLNSIG=' in t:
            clnsig.append(t.split('=')[1].strip())
        if 'MC=' in t:
            mc.append(t.split('|')[1].strip().split(',')[0])
        if 'CLNVC=' in t:
            clnvc.append(t.split('=')[1].strip())
        if 'ALLELEID=' in t:
            id.append(t.split('=')[1].strip())
        if 'GENEINFO=' in t:
            gene.append(t.split('=')[1].strip().split(':')[0])
        if 'CLNHGVS=' in t:
            clnhgvs.append(t.split('=')[1].strip())
    if len(clnsig) < len(chr):
        clnsig.append('NA')
    if len(mc) < len(chr):
        mc.append('NA')    
    if len(clnvc) < len(chr):
        clnvc.append('NA')  
    if len(id) < len(chr):
        id.append('NA') 
    if len(gene) < len(chr):
        gene.append('NA') 
    if len(clnhgvs) < len(chr):
        clnhgvs.append('NA') 


clinvar.close()

# Write table
with open(OUTPUT, 'w') as out:
    out.write('gene\tid\tchr\tpos\tclnsig\tmc\tclnvc\tclnhgvs\n')
    for i in range(len(chr)):
        out.write(
            f"{gene[i]}\t{id[i]}\t{chr[i]}\t{pos[i]}\t{clnsig[i]}\t{mc[i]}\t{clnvc[i]}\t{clnhgvs[i]}\n"
        )

print("Done.")


    