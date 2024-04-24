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

```
comandos
```

