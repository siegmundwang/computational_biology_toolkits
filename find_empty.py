import os
os.chdir("/home/sig/project/nocoding") # change working directory
dirs = filter(os.path.isdir, os.listdir()) # filter function
empty = []
for folder in dirs:
    folder = os.path.join("/home/sig/project/nocoding", folder) # add father directory
    if  not os.listdir(folder):
        empty.append(folder)
print(empty)
hd = open("No_ICGC.txt", "wt") # write in text mode
empty = map(lambda x : x + "\n", empty) # advanced function
hd.writelines(empty) # write lines
hd.close()
