# pw2wannier stores the k index in 4 digits, which crashes computations on grids larger than 9999. This fixes it.

f = open("silicon.amn")
out = open("silicon.amn.new", "w")

out.write(f.readline())
out.write(f.readline())


for line in f:
    line = line[0:10] + " " + line[10:]
    out.write(line)


f = open("silicon.mmn")
out = open("silicon.mmn.new", "w")

out.write(f.readline())
out.write(f.readline())

i = 0
for line in f:
    if ((i % 17) == 0):
        out.write(line[0:5] + " " + line[5:])
    else:
        out.write(line)
    i = i+1



f = open("silicon.eig")
out = open("silicon.eig.new", "w")

i = 0
for line in f:
    line = line[0:5] + " " + line[5:]
    out.write(line)
