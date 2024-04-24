# Antes de empezar

### Instalar lo siguiente
> pandas  
> requests

### Descargar los siguientes archivos

> OTUs.fasta  
> vsearch.exe

# Comandos de la sesión 3 del curso

```
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

