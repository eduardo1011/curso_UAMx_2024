################################################################################################
################################################################################################

#  Links de desacarga

################################################################################################
################################################################################################

### Descargar el módulo circleplot


### Descargar base de datos 16S-ITGDB
https://1drv.ms/u/s!ArGs92xOZGDEoj51fPQo2s0yxab7?e=6AQsIc

### Descargar programas de Blast

https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/ncbi-blast-2.8.1+-win64.exe






################################################################################################
################################################################################################

#  Comandos

################################################################################################
################################################################################################

### comando para hacer graficas
```
circle_plot(sec = genoma, window = 3500, VAL = 0.9, radio1 = 0.35, radio2 = 0.5, figsize = 5, showborde = True, ini = 90)
```

```
for w in list(range(1000, 10000, 500)):
    print(w)
    circle_plot(sec = genoma, window = w, VAL = 1.2, radio1 = 0.35, radio2 = 0.5, figsize = 5, showborde = False, ini = 0)
```

```
enzimas = {}
with open('enzymes.txt', 'r') as fq:
    for i in fq:
        i = i.rstrip()
        enzimas[i.split('\t')[0]] = [i.split('\t')[1], i.split('\t')[2]]
```




# -
