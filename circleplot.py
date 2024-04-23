import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.patches as patches
from matplotlib.path import Path
import numpy as np
import matplotlib.pyplot as plt

def open_genome(file = ''):
    N = 0
    seq = ''
    with open(file, 'r') as fq:
        for line in fq:
            line = line.rstrip()
            line = line.upper()
            if '>' in line:
                N += 1
                if N > 1:
                    break
            else:
                seq += line
    return seq

def inner_plots(lista = []):
    tipo = {1: Path.MOVETO, 2: Path.LINETO, 3: Path.CLOSEPOLY}
    x1, y1, x2, y2 = [], [], [], []
    for u in lista:
        A = u.get_path()
        mitad = int(len(A.vertices[:-2])/2)
        mitad1 = A.vertices[:-2][:mitad,:]
        mitad2 = A.vertices[:-2][mitad:,:]
        x1.append(mitad2[3][0])
        y1.append(mitad2[3][1])
        x2.append(mitad1[3][0])
        y2.append(mitad1[3][1])
    x1.append(x1[0])
    y1.append(y1[0])
    x2.append(x2[0])
    y2.append(y2[0])    

    verts = []
    codes = []
    n = 0
    for a, b in zip(x1 + list(reversed(x2)), y1 + list(reversed(y2))):
        if n == 0:
            verts.append([a, b])
            codes.append(tipo[1])
        if n == len(x1)-1:
            verts.append([a, b])
            codes.append(tipo[3])
        else:
            verts.append([a, b])
            codes.append(tipo[2])
        n += 1
    return verts, codes


def circle_plot(sec = '', window = 2500, VAL = 1, showborde = False, radio1 = 0.3, radio2 = 0.39, figsize = 8, ini = 90):

    gc_skews = []
    for i in range(0,len(sec), window):
        gc_skew = (sec[i:i+window].upper().count("G") - sec[i:i+window].upper().count("C")) / (sec[i:i+window].upper().count("G") + sec[i:i+window].upper().count("C"))
        gc_skews.append(gc_skew)
    gc_skews.append((sec[i:].upper().count("G") - sec[i:].upper().count("C")) / (sec[i:].upper().count("G") + sec[i:].upper().count("C")))
    gc_skew_positive = np.array(gc_skews)
    gc_skew_negative = np.array(gc_skews)
    gc_skew_positive[gc_skew_positive<0] = 0 
    gc_skew_negative[gc_skew_positive>0] = 0


    gc_count = []
    for i in range(0,len(sec), window):
        gc_amount = (sec[i:i+window].upper().count("G") + sec[i:i+window].upper().count("C")) / window
        gc_count.append(gc_amount)
    gc_count.append((sec[i:].upper().count("G") + sec[i:i+window].upper().count("C")) / (len(sec)-i))

    x = np.array(gc_count)

    gc_amounts = []
    for i in x:
        cal = np.log2(i / np.mean(x))
        if cal < 0:
            gc_amounts.append(cal)
        if cal > 0:
            gc_amounts.append(cal)
    gc_amounts = np.array(gc_amounts)

    gc_amounts_positive = np.array(gc_amounts)
    gc_amounts_negative = np.array(gc_amounts)
    gc_amounts_positive[gc_amounts_positive<0] = 0 
    gc_amounts_negative[gc_amounts_positive>0] = 0


    r1 = radio1
    w1 = 0.1

    minimo, maximo = abs(gc_skew_positive).min(), abs(gc_skew_positive).max()

    gc_skew_POS = []
    for i in abs(gc_skew_positive):
        if i == 0:
            gc_skew_POS.append(i)
        else:
            per = (i - minimo) / (maximo - minimo)
            gc_skew_POS.append(per * (w1*VAL))

    minimo, maximo = abs(gc_skew_negative).min(), abs(gc_skew_negative).max()

    gc_skew_NE = []
    for i in abs(gc_skew_negative):
        if i == 0:
            gc_skew_NE.append(i)
        else:
            per = (i - minimo) / (maximo - minimo)
            gc_skew_NE.append(per * (w1*VAL))

    barra = 360 / len(gc_skew_NE)

    r2 = radio2
    w2 = 0.1

    minimo, maximo = abs(gc_amounts_positive).min(), abs(gc_amounts_positive).max()

    gc_amounts_POS = []
    for i in abs(gc_amounts_positive):
        if i == 0:
            gc_amounts_POS.append(i)
        else:
            per = (i - minimo) / (maximo - minimo)
            gc_amounts_POS.append(per * (w2*VAL))


    minimo, maximo = abs(gc_amounts_negative).min(), abs(gc_amounts_negative).max()

    gc_amounts_NE = []
    for i in abs(gc_amounts_negative):
        if i == 0:
            gc_amounts_NE.append(i)
        else:
            per = (i - minimo) / (maximo - minimo)
            gc_amounts_NE.append(per * (w2*VAL))

    # sentido de las manecillas del reloj
    rec = []
    n = ini
    o = ini - barra
    fin = ini - barra
    for i in gc_skew_positive:
        rec.append([o, n])
        n -= barra
        o -= barra


    centro = (0.5, 0.5)

    fig = plt.figure(figsize=(figsize, figsize))

    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_aspect('equal', 'box')
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)



    gc_skew_POS_lines = []
    for d, z in zip(gc_skew_POS, rec):
        sour = mpatches.Wedge(centro, r1 - (w1/2), z[0], z[1],  width=-d, facecolor='#1BA21B', zorder = 2)
        #ax.add_patch(sour) # muestra barras
        gc_skew_POS_lines.append(sour)
    #... 
    verts, codes = inner_plots(lista = gc_skew_POS_lines)
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='#1BA21B', lw=0)
    ax.add_patch(patch)  

    gc_skew_NE_lines = []
    for d, z in zip(gc_skew_NE, rec):
        sour = mpatches.Wedge(centro, r1 - (w1/2), z[0], z[1], width=d, facecolor='#990099', zorder = 2)
        #ax.add_patch(sour)
        gc_skew_NE_lines.append(sour)
    #...    
    verts, codes = inner_plots(lista = gc_skew_NE_lines)
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='#990099', lw=0)
    ax.add_patch(patch)

    #////////////////////////////////////////////////////

    gc_amounts_POS_lines = []
    for d, z in zip(gc_amounts_POS, rec):
        sour = mpatches.Wedge(centro, r2 - (w2/2), z[0], z[1],  width=-d, facecolor='black', zorder = 2)
        #ax.add_patch(sour)
        gc_amounts_POS_lines.append(sour)
        
    #...    
    verts, codes = inner_plots(lista = gc_amounts_POS_lines)
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='black', lw=0)
    ax.add_patch(patch)
    #--------------------------

    gc_amounts_NE_lines = []
    for d, z in zip(gc_amounts_NE, rec):
        sour = mpatches.Wedge(centro, r2 - (w2/2), z[0], z[1], width=d, facecolor='black', zorder = 2)
        #ax.add_patch(sour) 
        gc_amounts_NE_lines.append(sour)
        
    #...    
    verts, codes = inner_plots(lista = gc_amounts_NE_lines)
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='black', lw=0)
    ax.add_patch(patch)

    if showborde == True:
        pass
    if showborde == False:
        ax.axis('off')
    
    return plt.show()








