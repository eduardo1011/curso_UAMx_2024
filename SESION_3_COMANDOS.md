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

```
comandos
```

```
comandos
```

