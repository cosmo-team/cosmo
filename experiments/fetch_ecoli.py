#!/usr/bin/python3
import ftplib
import re
pre = "Escherichia coli"

conn = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
conn.login()
suff = ".fna.gz"
def fetch(url):
    fields = url.split("/")
    print(fields)
    conn.cwd("/" + "/".join(fields[3:]))
    files = conn.nlst()
    for filename in files:
        if filename[-len(suff):] == suff:
            conn.retrbinary('RETR ' + filename, open(filename, 'wb').write)

count=0            
for line in open("assembly_summary.txt"):
    fields = line.split("\t")
    if fields[7][:len(pre)] == pre:
        print("downloading", count)
        count +=1
        fetch(fields[-1].strip())
        
        
conn.quit()
