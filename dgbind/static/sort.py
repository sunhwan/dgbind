import sys

#prefix = sys.argv[1]
#step = int(sys.argv[2])
#num_replica = int(sys.argv[3])
#final_step = -1
#if len(sys.argv) > 4:
#    final_step = int(sys.argv[4])

prefix = 'barnf'
num_replica = 16
step = int(sys.argv[1])
final_step = -1
if len(sys.argv) > 2:
    final_step = int(sys.argv[2])

import commands
import numpy

nl = int(commands.getoutput('wc -l output/0/%s.job%d.0.history' % (prefix, step)).strip().split()[0])
print nl

history_fps = [open('output/%d/%s.job%d.%d.history' % (i, prefix, step, i)) for i in range(num_replica)]

reps = numpy.zeros((num_replica, nl, 2), dtype=numpy.int)
cnt = 1
while 1:
    for i in range(num_replica):
        try:
            line = history_fps[i].readline()
            time, rep = map(int, line.split()[:2])
        except:
            final_step = time
            break
        reps[i,cnt-1,0] = time
        reps[i,cnt-1,1] = rep

    print time
    if final_step > 0 and time >= final_step: break
    cnt += 1

history_fps = [open('output/%d/%s.job%d.%d.history' % (i, prefix, step, i)) for i in range(num_replica)]
sorted_history_fps = [open('output/%d/%s.job%d.%d.sort.history' % (i, prefix, step, i), 'w') for i in range(num_replica)]
for j in range(cnt-1):
    for i in range(num_replica):
        rep = reps[i,j,1]
        line = history_fps[i].readline()
        sorted_history_fps[rep].write(line)

del(history_fps)
del(sorted_history_fps)

colvars_fps = [open('output/%d/%s.job%d.%d.colvars.traj' % (i, prefix, step, i)) for i in range(num_replica)]
sorted_colvars_fps = [open('output/%d/%s.job%d.%d.sort.colvars.traj' % (i, prefix, step, i), 'w') for i in range(num_replica)]
for j in range(cnt):
    for i in range(num_replica):
        time, rep = reps[i,j]
        while 1:
            line = colvars_fps[i].readline()
            sorted_colvars_fps[rep].write(line)
            if not line.startswith('#') and int(line.strip().split()[0]) != time: continue
            break
    if final_step > 0 and time >= final_step: break
    print time

