import sys
import numpy as np
import libk

def importsuperquad(f):
    if isinstance(f, str):
        ff = open(f, "r")
    if f is sys.stdin:
        ff = f

    section = 0

    for line in ff:
        if section == 0:
            title = line
            section = 1
            continue

        if section == 1:
            numbers_ = [int(i) for i in line.split()]
            section = 2
            B = []
            P = []
            keys = []
            names = []
            continue

        if section == 2:
            if len(names) == numbers_[2]:
                section = 3
                temperature = float(line)
            else:
                names.append(line.strip())

            continue

        if section == 3:
            if line.strip() == '':
                section = 4
                emf = []
                current_emf = {}
                current_emf['plot_key'] = []
                current_emf['order'] = []
                current_emf['T'] = []
                current_emf['buret'] = []
                current_emf['flag'] = []
                current_emf['V'] = []
                current_emf['emf'] = []
                continue

            aux = line.split()
            B.append(float(aux[0]))
            P.append([int(i) for i in aux[1:-1]])
            keys.append(int(aux[-1]))
            continue

        if section == 4:
            if line.strip() == '':
                section = 5
                continue

            aux = line.split()
            current_emf['plot_key'].append(int(aux[0]))
            current_emf['order'].append(int(aux[1]))
            current_emf['T'].append(float(aux[2]))
            current_emf['buret'].append(float(aux[3]))
            current_emf['flag'].append(int(aux[4]))
            continue

        if section == 5:
            aux = line.split()
            current_emf['V0'] = float(aux[0])
            current_emf['error_V0'] = float(aux[1])
            section = 6
            continue

        if section == 6:
            aux = line.split()
            current_emf['n'] = int(aux[0])
            current_emf['h'] = int(aux[1])
            current_emf['E0'] = float(aux[2])
            current_emf['error_E0'] = float(aux[3])
            section = 7
            continue

        if section == 7:
            section = 8
            continue

        if section == 8:
            if line.strip() == '':
                emf.append(current_emf)
                current_emf = {}
                current_emf['plot_key'] = []
                current_emf['order'] = []
                current_emf['T'] = []
                current_emf['buret'] = []
                current_emf['flag'] = []
                current_emf['V'] = []
                current_emf['emf'] = []
                section = 9
                continue

            aux = line.split()
            current_emf['V'].append(float(aux[0]))
            current_emf['emf'].append(float(aux[1]))

        if section == 9:
            if line.strip() == '':
                break

            aux = line.split()
            current_emf['plot_key'].append(int(aux[0]))
            current_emf['order'].append(int(aux[1]))
            current_emf['T'].append(float(aux[2]))
            current_emf['buret'].append(float(aux[3]))
            current_emf['flag'].append(int(aux[4]))
            section = 4

    if isinstance(f, str):
        ff.close()

    return {
        "labels": names,
        "emf":      [np.array(e['emf']) for e in emf],
        "V":        [np.array(e['V']) for e in emf],
        "P":        np.array(P),
        "temperature": temperature,
        "logB":     np.array(B),
        "flags":    np.array(keys),
        "E0":       [e['E0'] for e in emf],
        "error_E0": [e['error_E0'] for e in emf],
        "V0":       [e['V0'] for e in emf], 
        "error_V0": [e['error_V0'] for e in emf], 
        "T0":       [np.array(e['T']) for e in emf], 
        "buret":    [np.array(e['buret']) for e in emf], 
        "hindex":   [(e['h']-1) for e in emf],
        "fRTnF":    [0.086173424*(temperature+273.15) / e['n'] for e in emf]  # mV
    }
