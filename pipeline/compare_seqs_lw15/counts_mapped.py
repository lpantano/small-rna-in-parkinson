import sys

sam_file = sys.argv[1]
ma_file = sys.argv[2]
out_file = sys.argv[3]
names = {}
counts = {}
with open(ma_file) as seqs:
    header = seqs.readline()
    for line in seqs:
        cols = line.split("\t")
        names[cols[1]] = cols[0].replace(">", "")
        counts[cols[0].replace(">", "")] = line.strip()
done = {}
with open(sam_file) as feature:
    with open(out_file, 'w') as out:
        out.write("tcca\t%s" % header)
        for line in feature:
            if line.startswith("@"):
                continue
            if line.split("\t")[2] in "*":
                continue
            if line.split("\t")[9].endswith("CCA"):
                mod = 1
            else:
                mod = 0
            if not counts[line.split("\t")[0]] in done:
                done[counts[line.split("\t")[0]]] = 0
                out.write("%s\t%s\n" % (mod, counts[line.split("\t")[0]]))
