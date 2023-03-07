masses = {
    'x' : 0.000,
    'h' : 1.008,
    'he': 4.003,
    'li': 7.016,
    'be': 9.012,
    'b' : 11.009,
    'c' : 12.000,
    'n' : 14.003,
    'o' : 15.995,
    'f' : 18.998,
    'ne': 19.992,
    'na': 22.990,
    'mg': 23.985,
    'al': 26.981,
    'si': 27.977,
    'p' : 30.974,
    's' : 31.972,
    'cl': 34.969,
    'ar': 39.962,
    'k' : 38.964,
    'ca': 39.963,
    'sc': 44.956,
    'ti': 47.948,
    'v' : 50.944,
    'cr': 51.941,
    'mn': 54.938,
    'fe': 55.935,
    'co': 58.933,
    'ni': 57.935,
    'cu': 62.930,
    'zn': 63.929,
    'ga': 68.926,
    'ge': 73.921,
    'as': 74.922,
    'se': 79.917,
    'br': 78.918,
    'kr': 83.912,
    'rb': 84.912,
    'sr': 87.906,
    'y': 88.906,
    'zr': 89.905,
    'nb': 92.906,
    'mo': 97.905,
    'tc': 98.906,
    'ru': 101.904,
    'rh': 102.906,
    'pd': 107.904,
    'ag': 106.905,
    'cd': 113.903,
    'in': 114.904,
    'sn': 119.902,
    'sb': 120.904,
    'te': 129.906,
    'i' : 126.904,
    'xe': 131.904,
    'cs': 132.905,
    'ba': 137.905,
    'la': 138.906,
    'ce': 139.905,
    'pr': 140.908,
    'nd': 141.907,
    'pm': 145.914,
    'sm': 151.920,
    'eu': 152.921,
    'gd': 157.924,
    'tb': 158.925,
    'dy': 163.929,
    'ho': 164.930,
    'er': 165.930,
    'tm': 168.934,
    'yb': 173.939,
    'lu': 174.941,
    'hf': 179.947,
    'ta': 180.948,
    'w' : 183.951,
    're': 186.956,
    'os': 191.961,
    'ir': 192.963,
    'pt': 194.965,
    'au': 196.967,
    'hg': 201.971,
    'tl': 204.974,
    'pb': 207.976,
    'bi': 208.980,
    'po': 209.483, 
    'at': 210.487,
    'rn': 217.673, 
    'fr': 223.020,
    'ra': 225.274, 
    'ac': 227.028,
    'th': 231.036, 
    'pa': 231.036,
    'u' : 238.051,
    'np': 236.547, 
    'pu': 240.723, 
    'am': 242.059, 
    'cm': 245.567, 
    'bk': 248.073, 
    'cf': 250.578, 
    'es': 252.083,
    'fm': 257.095,
    'md': 259.101,
    'no': 259.101,
    'lr': 262.110,
    'rf': 267.122,
    'db': 268.126,
    'sg': 271.134,
    'bh': 272.138,
    'hs': 270.134,
    'mt': 276.152,
    'ds': 281.165,
    'rg': 280.165,
    'cn': 285.177,
    'nh': 284.179,
    'fl': 289.190,
    'mc': 288.193,
    'lv': 293.204,
    'ts': 292.207,
    'og': 294.214
}

import sys
import numpy as np

#Massa atomica (kg), Constante de Planck (J s), velocidade da luz (m s-1)
u, h, c = 1.66053886e-27, 6.62606957e-34, 2.99792458**8

#Centro de massa: cm (Reposicionar as coordenadas do átomo em relação
#ao seu centro de massa)
def centro_de_massa(massas, xyz):
    # Posição do centro de massa nas coordenadas originais
    cm = sum(massas[:,np.newaxis] * xyz) / np.sum(massas)
    # Transforme em coordenadas cm e retorne
    xyz -= cm
    return xyz


def obter_matriz_de_inercia(massas, xyz):
    # Momento de inercia
    xyz = centro_de_massa(massas, xyz)
    x, y, z = xyz.T
    Ixx = np.sum(massas * (y**2 + z**2)) #np.sum(retorna a soma dos arrays)
    Iyy = np.sum(massas * (x**2 + z**2))
    Izz = np.sum(massas * (x**2 + y**2))
    Ixy = -np.sum(massas * x * y)
    Iyz = -np.sum(massas * y * z)
    Ixz = -np.sum(massas * x * z)
    I = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
    return I

def obter_momento_principal_de_inercia(I):
    Ip = np.linalg.eigvals(I) #Calcular os autovalores da matriz I
    # Ordenar e converter momentos principais de inércia para SI (kg.m2)
    Ip.sort()
    return Ip

def classificação_moleculas(A, B, C):
    if np.isclose(A, B):
        if np.isclose(B, C):
            return 'Rotor esférico'
        return 'Rotor simétrico oblato'
    if np.isclose(B, C):
        return 'Rotor simétrio prolato'
    return 'Rotor assimétrico'

def read_xyz(filename):
    try:
        data = np.loadtxt(filename, skiprows=2)
    except FileNotFoundError:
        print('No such file:', filename)
        sys.exit(1)
    except ValueError as e:
        print('Malformed data in {}: {}'.format(filename, e))
        sys.exit(1)
    return data[:,0], data[:,1:]

try:
    massas, xyz = read_xyz(sys.argv[1])
except IndexError:
    print('Usage: {} <xyz filename>'.format(sys.argv[0]))
    sys.exit(1)

I = obter_matriz_de_inercia(massas, xyz)

Ip = obter_momento_principal_de_inercia(I)
Ip *= u / 1.e20
A, B, C = h / 8 / np.pi**2 / c / 100 / Ip
rotor_type = classificação_moleculas(A, B, C)

print('{}: A={:.6f}, B={:.6f}, C={:.6f} cm-1'.format(rotor_type, A, B, C))





