# Antes de empezar

### Instalar lo siguiente
> pandas  
> requests

### Descargar los siguientes archivos

> OTUs.fasta  
> vsearch.exe

# Comandos de la sesión 3 del curso

```
import subprocess
import pandas as pd
from pandas import DataFrame
from collections import Counter, OrderedDict
import numpy as np


def sacc_tax(file = ''):
    fas = {}
    with open(file) as fq:
        for line in fq:
            line = line.rstrip()
            if '>' in line:
                header = re.sub('>', '', line).split(' ')[0]
                tt = re.sub('>', '', line).split(' ')[1]
                tt = re.sub('tax=d:', '', tt)
                fas[header] = re.sub(',.:', '#', tt)
    return [[i] + fas[i].split('#') for i in fas]
```

```
idens = open_file(file = '16S-ITGDB_mod.fasta')
columnas_tax = ['sacc', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'specie']
taxonomia = DataFrame(idens, columns = columnas_tax)
```
## Comando para ejecutar blastn
```
subprocess.call('blastn -db db/16S-ITGDB -query OTUs.fasta -evalue 1E-6 -outfmt "6" -out otus_vs_16S.txt', shell = True)
```



```
names = ["qacc","sacc", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
df = pd.read_csv('otus_vs_16S.txt', sep = '\t', names = names)[['qacc','sacc','pident']]
blastn = df.merge(taxonomia, on = 'sacc', how = 'left').sort_values(by = 'pident', ascending = False).reset_index(drop = True)
```

#### Comando para ejecutar de Vsearch
```
subprocess.call('vsearch --sintax OTUs.fasta --db 16S-ITGDB.fasta --tabbedout otus_tax.sintax --sintax_cutoff 0.8', shell = True)
```

#### Comando para abrir el resultado de Vsearch
```
recording = []
tax = {}
with open('otus_tax.sintax', 'r') as fq:
    for i in fq:
        i = i.rstrip()
        if (len(i.split('\t')) == 1) or ('Chloroplast' in i) or ('Mitochondria' in i) or ('Eukaryota' in i):
            pass
        else:
            recording.append(i.split('\t')[0])
            tax[i.split('\t')[0]] = [i.split('\t')[1], re.sub('[(][0-9.]{1,50}[)]', '', i.split('\t')[1]), i.split('\t')[2], i.split('\t')[3]]
vsearch = DataFrame([[a] + b.split(',') for a, b in [[i, tax[i][0]] for i in tax]], columns = ['qacc'] + columnas_tax[1:])
```

# Link al NCBI Blast

https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome

```
comandos
```
# Buscar dominios dentro de proteínas

```
def open_file(file = ''):
    fas = {}
    with open(file) as fq:
        for line in fq:
            line = line.rstrip()
            if '>' in line:
                header = re.sub('>', '', line)
                s = ''
            else:
                s += line
            fas[header] =  s
    return fas
```

```
prot = open_file(file = 'problema.fasta')
```

```
nombre = ['PROTEIN_KINASE_ATP', 'ZINC_FINGER_C2H2_1', 'PEROXIDASE_2', 'PYRUVATE_KINASE']
dominios = ['[LIV]G[^P]G[^P][FYWMGSTNH][SGA][^PW][LIVCAT][^PD].[GSTACLIVMFY].{5,18}[LIVMFYWCSTAR][AIVP][LIVMFAGCKR]K',
            'C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H',
            '[SGATV][^D].{2}[LIVMA]R[LIVMA].[FW]H[^V][SAC]',
            '[LIVAC].[LIVM]{2}[SAPCV]K[LIV]E[NKRST].[DEQHS][GSTA][LIVM]']
dom_name = dict(zip(dominios, nombre))
```

```
for dom in dominios:
    print('●'*50)
    print(dom_name[dom])
    print(dom)
    print('●'*50)
    for i in prot:
        seq = prot[i]
        pat = re.findall(dom, seq)
        if len(pat) > 0:
            print(pat, i)
```

```
[re.sub('[|].*', '', re.sub('^...', '', i)) for i in prot]
```

```
comandos
```

```
comandos
```

```
comandos
```

```
comandos
```

```
comandos
```

```
comandos
```

```
comandos
```

```
comandos
```

```
comandos
```

```
comandos
```

