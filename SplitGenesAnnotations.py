###Read the file in gtf 
fh = open("/Users/Pinedasans/Catalyst/Data/Genotyping/annotation.txt")
fout = open("/Users/Pinedasans/Catalyst/Data/Genotyping/annotation_genes.txt",'w')
count = 0
for line in fh:
    line=line.rstrip() #Delete black spaces
    print (line)
    words = line.split()
    fout.write( words[0] + "\t" + words[1] + "\t"  + words[2] + "\t" + words[3] + "\t" + words[4] + "\t" + words[7] + "\t" + words[13] + '\n')
    count=count+1
print(count)


