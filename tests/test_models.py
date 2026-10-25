import numpy as np
from fhn_control.models import FHNParameters, rhs, equilibrium, jacobian, hopf_candidates

def test_equilibrium():
    p = FHNParameters()
    assert np.linalg.norm(rhs(0, equilibrium(p), p)) < 1e-9

def test_jacobian():
    p = FHNParameters()
    v = equilibrium(p)[0]
    h = 1e-6
    approx = np.column_stack([(rhs(0, equilibrium(p) + h*e, p)-rhs(0, equilibrium(p) - h*e, p))/(2*h) for e in np.eye(2)])
    assert np.allclose(jacobian(v, p), approx, atol=1e-7)

def test_hopf_traces():
    p = FHNParameters()
    for v, I in hopf_candidates(p):
        ph = FHNParameters(a=p.a,b=p.b,tau=p.tau,I=I)
        assert abs(np.trace(jacobian(v, ph))) < 1e-10
