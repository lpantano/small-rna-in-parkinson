import sys

num = 0
with open(sys.argv[1]) as handle:
    for line in handle:
        if line.startswith(">"):
            counts = line.strip().split("_x")[1]
        else:
            if len(line.strip()) > 17 and int(counts) > 10:
                num += 1
                for i in range(int(counts)):
                    print ">seq_%s_x%s\n%s" % (num, i, line.strip())
