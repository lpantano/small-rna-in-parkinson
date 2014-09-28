import sys

sam_file = sys.argv[1]
ma_file = sys.argv[2]
parse_file = sys.argv[3]
out_file = sys.argv[4]
names = {}
counts = {}
with open(ma_file) as seqs:
    header = seqs.readline()
    for line in seqs:
        cols = line.split("\t")
        names[cols[1]] = cols[0].replace(">", "")
        counts[cols[0].replace(">", "")] = line.strip()
mapped = {}
with open(parse_file) as ann:
    for line in ann:
        if line.startswith("S"):
            mapped[names[line.split(" ")[2]]] = 0
with open(sam_file) as feature:
    with open(out_file, 'w') as out:
        out.write("tcca\t%s" % header)
        for line in feature:
            if line.startswith("@"):
                continue
            if line.split("\t")[2] in "*":
                continue
            if line.split("\t")[0] in mapped:
                in_file = 1
            else:
                in_file = 0
            if line.split("\t")[9].endswith("CCA"):
                mod = 1
            else:
                mod = 0
            out.write("%s\t%s\t%s\n" % (in_file, mod, counts[line.split("\t")[0]]))
