import re, sys

def clean(line):
    line = line.strip()
    while line[0:1]=='*': line=line[1:]
    return line.strip()

data = open(sys.argv[1],'rb').read()
for item in data.split('/**'):
    if not item.strip(): continue
    desc, code = item.split('*/',1)
    desc = ('\n'.join(clean(line) for line in desc.split('\n'))).strip()
    if desc[-1:]=='.': desc=desc[:-1]+':'
    elif not desc[-1:]==':': desc=desc+':'
    desc = desc.capitalize()
    print '\\noindent\n'+desc,
    print '\\begin{lstlisting}'
    print code.strip()
    print '\\end{lstlisting}'
