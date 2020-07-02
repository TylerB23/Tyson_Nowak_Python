import numpy as np

def gk(Va, Vi, Ja, Ji):
    alpha = Vi - Va
    beta = Vi - Va + Va*Ji + Vi*Ja
    gamma = Va*Ji
    
    return (2*gamma)/(beta+np.sqrt((beta**2) - 4*alpha*gamma))
