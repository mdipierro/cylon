import sys, re

print """
#include "math.h"
#include "vector"
#include "list"
using namespace std;
"""
r=re.compile('\\\\begin\{lstlisting\}(.*?)\\\\end\{lstlisting\}',re.DOTALL)
for item in r.finditer(open(sys.argv[1],'rb').read()): print item.group(1)
