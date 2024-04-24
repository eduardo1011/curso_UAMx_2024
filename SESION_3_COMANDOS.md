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
url = 'https://rest.uniprot.org/uniprotkb/'+entry+'.txt'
print(url)
```

```
info = requests.get(url, stream=True).content.decode().split('\n')
info
```

```
m = []
for i in info:
    if re.search('^FT', i):
        w = i
        if ('STRAND' in w) or ('CONFLICT' in w) or ('MUTAGEN' in w):
            pass
        else:
            m.append(re.sub(' * ', ' ', w))
g = []
for i in m:
    if re.search('^FT [A-Z]{3,50}', i):
        z = i
        g.append([z, re.sub('FT $', '', ''.join(m).split(z)[1].split('/')[1])])
anotaciones = []
for i in g:
    if ('ECO:' in i[0]) or ('COMPBIAS' in i[0]):
        pass
    else:
        text = re.sub('evidence.*', '()', i[1])
        tipo = re.sub('^FT ', '', i[0])
        if re.search('".*"', text):
            texto1 = re.search('".*"', text).group()
        else:
            texto1 = '()'
        anotaciones.append([tipo.split(' ')[0], [int(h) for h in re.findall('\d+', tipo.split(' ')[1])], texto1])  
```

```
import matplotlib.pyplot as plt

mas_d1 = [i for i in anotaciones if len(i[1]) == 2]
solo1 = [i for i in anotaciones if len(i[1]) == 1]

fig = plt.figure(figsize=(15, 3))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 700)
ax.set_ylim(-len(solo1) -1, len(mas_d1) + 1)

N = 0
for n in mas_d1:
    A = n[0]
    Bx1, Bx2 = n[1]
    C = n[2]
    if A == 'CHAIN':
        Acol = 'deepskyblue'
    else:
        Acol = 'grey'
    ax.plot([Bx1, Bx2-1], [N, N], color = Acol, solid_capstyle = 'butt', lw = 10)
    ax.text(Bx2, N, '  '+A+' = '+ C, ha = 'left', va = 'center', fontsize = 12, weight="bold")
    N += 1
    
N = -1    
for n in solo1:
    A = n[0]
    Bx = n[1][0]
    C = n[2]
    ax.scatter(Bx, N, s = 50, c = 'red', marker = 's')
    ax.text(Bx, N, '  '+A+' = '+ C, ha = 'left', va = 'center', fontsize = 11, weight="bold")
    N -= 1
    
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.get_yaxis().set_visible(False)
ax.tick_params(bottom=True, right=False, top=False, left=True, width = 3, length=5, color='black')
plt.xticks(fontsize=12, rotation=0, weight="bold")
ax.set_xlabel('\nLength Protein (aa)', fontsize = 15, weight="bold")


plt.close()
fig
```

```
import matplotlib.pyplot as plt

mas_d1 = [i for i in anotaciones if len(i[1]) == 2]
mas_d1 = [i for i in mas_d1 if 'HELIX' not in i]
mas_d1 = [i for i in mas_d1 if 'TURN' not in i]

solo1 = [i for i in anotaciones if len(i[1]) == 1]

fig = plt.figure(figsize=(15, 9))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 700)
ax.set_ylim(-len(solo1) -1, len(mas_d1) + 1)

N = 0
for n in mas_d1:
    A = n[0]
    Bx1, Bx2 = n[1]
    C = n[2]
    if A == 'CHAIN':
        Acol = 'deepskyblue'
    else:
        Acol = 'grey'
    ax.plot([Bx1, Bx2-1], [N, N], color = Acol, solid_capstyle = 'butt', lw = 10)
    ax.text(Bx2, N, '  '+A+' = '+ C, ha = 'left', va = 'center', fontsize = 12, weight="bold")
    N += 1
    
N = -1    
for n in solo1:
    A = n[0]
    Bx = n[1][0]
    C = n[2]
    ax.scatter(Bx, N, s = 50, c = 'red', marker = 's')
    ax.text(Bx, N, '  '+A+' = '+ C, ha = 'left', va = 'center', fontsize = 11, weight="bold")
    N -= 1
    
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.get_yaxis().set_visible(False)
ax.tick_params(bottom=True, right=False, top=False, left=True, width = 3, length=5, color='black')
plt.xticks(fontsize=12, rotation=0, weight="bold")
ax.set_xlabel('\nLength Protein (aa)', fontsize = 15, weight="bold")


plt.close()
fig
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

