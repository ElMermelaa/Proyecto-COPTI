"""Modelo de FitzHugh–Nagumo para una neurona aislada (I constante)."""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq

@dataclass(frozen=True)
class FHNParameters:
    a: float = 0.7
    b: float = 0.8
    tau: float = 12.5
    I: float = 0.9

    def __post_init__(self):
        if self.a < 0 or not 0 <= self.b < 1 or self.tau <= 0:
            raise ValueError('Se requiere a>=0, 0<=b<1 y tau>0.')


def rhs(t, x, p: FHNParameters):
    v, w = x
    return np.array([v - v**3 / 3 - w + p.I, (v + p.a - p.b*w) / p.tau])


def equilibrium(p: FHNParameters):
    # Evita dividir por b, de modo que también funcione para b=0.
    def poly(v):
        return p.b * v**3 + 3 * (1-p.b) * v + 3*p.a - 3*p.b*p.I
    # Para 0<=b<1 es monótono y la única raíz siempre está en algún [-R,R].
    R = max(2., abs(p.a)+abs(p.I)+2.)
    while poly(-R) > 0 or poly(R) < 0:
        R *= 2
    v = brentq(poly, -R, R)
    w = v - v**3 / 3 + p.I
    return np.array([v, w])


def jacobian(v, p: FHNParameters):
    return np.array([[1-v*v, -1], [1/p.tau, -p.b/p.tau]])


def hopf_candidates(p: FHNParameters):
    """Pares (v*, I*) donde tr(J)=0; NO certifica Hopf no degenerada."""
    if not (p.b > 0 and p.tau > p.b):
        return []
    q = np.sqrt(1-p.b/p.tau)
    return [(float(v), float(v**3/3+(1/p.b-1)*v+p.a/p.b)) for v in (-q,q)]
