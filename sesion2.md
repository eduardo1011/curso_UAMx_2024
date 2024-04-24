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

## Install BLAST en linux
$ sudo apt-get install ncbi-blast+




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

# widgets

### forzar la instalacion de la version 7.7.0
```
python -mpip install --force-reinstall ipywidgets==7.7.0
```

```
import numpy as np
from ipywidgets import Button, GridBox, Layout, ButtonStyle, Box, HBox, VBox
import ipywidgets as widgets
from IPython.display import clear_output, display

enz = widgets.Dropdown(options = sorted([i+' ('+enzimas[i][1]+')' for i in enzimas]), value = sorted([i+' ('+enzimas[i][1]+')' for i in enzimas])[0], disabled = False,
                                   layout = Layout(width='290px', height='25px'))
agregar = widgets.SelectionSlider(options=list(range(0, 101, 5)),value=20,description='Add:',disabled=False,
                                continuous_update=False,orientation='horizontal',readout=True,
                                   layout=Layout(width='300px', height='25px'))



def restriction(enz, agregar):
    E = enz.split(' ')[0]
    cortes = [[i.start()+1, i.end()+1] for i in re.finditer(enzimas[E][0], genoma)]
    print('Frequency=', len(cortes), enzimas[E])
    uu = ['   Ini= ' + str(i[0]) + '   Fin= ' + str(i[1]) + '   ' + genoma[i[0]-1-agregar:i[0]-1]+enz.split(' ')[1]+genoma[i[1]-1:i[1]-1+agregar] for i in cortes]
    ww = widgets.Select(options=uu,description='Result:', rows=15,disabled=False, layout=Layout(width='1000px', height='250px'))
    display(ww)
out = widgets.interactive_output(restriction, {'enz':enz, 'agregar':agregar})

HBox([VBox([enz, agregar]), out])
```



# -


```
for i in fas:
    if re.findall('TGTGA.*TCACT', fas[i]) != []:
        print(i)
        print(re.findall('.....TGTGA.*TCACT', fas[i]))
        print('#-----------------------------------------------')
```

