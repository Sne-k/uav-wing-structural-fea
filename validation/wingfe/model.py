"""Assembly and solvers (static, modal, linear buckling) for shell models."""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .shell import element_matrices, generalized_strains, geometric_stiffness


class ShellModel:
    """nodes: (n, 3) array; elements: list of (node_ids[4], group_name).

    sections: {group_name: Section}; ref_dirs: {group_name: local-x direction}.
    """

    def __init__(self, nodes, elements, sections, ref_dirs):
        self.nodes = np.asarray(nodes, float)
        self.elements = elements
        self.sections = sections
        self.ref_dirs = ref_dirs
        self.ndof = 6 * len(self.nodes)
        self._built = False

    # ------------------------------------------------------------ assembly
    def build(self, want_mass=True):
        rows, cols, kv, mv = [], [], [], []
        self.edata = []
        for conn, grp in self.elements:
            X = self.nodes[conn]
            K, M, T, data = element_matrices(X, self.sections[grp], self.ref_dirs[grp], want_mass)
            dofs = np.array([6 * n + d for n in conn for d in range(6)])
            r, c = np.meshgrid(dofs, dofs, indexing="ij")
            rows.append(r.ravel())
            cols.append(c.ravel())
            kv.append(K.ravel())
            mv.append(M.ravel())
            self.edata.append((dofs, T, data))
        rows = np.concatenate(rows)
        cols = np.concatenate(cols)
        self.K = sp.csr_matrix((np.concatenate(kv), (rows, cols)), shape=(self.ndof, self.ndof))
        self.M = sp.csr_matrix((np.concatenate(mv), (rows, cols)), shape=(self.ndof, self.ndof))
        self._rc = (rows, cols)
        self._built = True
        return self

    def mass(self):
        """Total structural mass from the consistent mass matrix (x translations)."""
        ux = np.arange(0, self.ndof, 6)
        return self.M[ux][:, ux].sum()

    def fix_nodes(self, node_ids, dofs=range(6)):
        self.fixed = np.array(sorted({6 * n + d for n in node_ids for d in dofs}))
        free = np.ones(self.ndof, bool)
        free[self.fixed] = False
        self.free = np.where(free)[0]

    # ------------------------------------------------------------ solvers
    def solve_static(self, F):
        Kff = self.K[self.free][:, self.free].tocsc()
        u = np.zeros(self.ndof)
        u[self.free] = spla.spsolve(Kff, F[self.free])
        self.reactions = self.K @ u - F
        return u

    def solve_modes(self, k=10):
        Kff = self.K[self.free][:, self.free].tocsc()
        Mff = self.M[self.free][:, self.free].tocsc()
        w2, v = spla.eigsh(Kff, k=k, M=Mff, sigma=0.0, which="LM")
        order = np.argsort(w2)
        w2, v = w2[order], v[:, order]
        phi = np.zeros((self.ndof, k))
        phi[self.free] = v
        return np.sqrt(np.abs(w2)) / (2 * np.pi), phi

    def element_results(self, u):
        """Generalized strains (4 Gauss points x 6) for every element."""
        return [generalized_strains(data, T, u[dofs]) for dofs, T, data in self.edata]

    def membrane_resultants(self, u):
        res = []
        for (dofs, T, data), (conn, grp) in zip(self.edata, self.elements):
            e = generalized_strains(data, T, u[dofs])
            sec = self.sections[grp]
            N = e[:, :3] @ sec.A.T + e[:, 3:] @ sec.B.T
            res.append(N)
        return res

    def solve_buckling(self, u0, k=6):
        """Load factors lambda with (K + lambda*Kg) phi = 0 for the stress state of u0."""
        rows, cols = self._rc
        gv = []
        for (dofs, T, data), N in zip(self.edata, self.membrane_resultants(u0)):
            gv.append(geometric_stiffness(data, T, N).ravel())
        Kg = sp.csr_matrix((np.concatenate(gv), (rows, cols)), shape=(self.ndof, self.ndof))
        Kff = self.K[self.free][:, self.free].tocsc()
        Kgf = Kg[self.free][:, self.free].tocsc()
        # (-Kg) phi = mu K phi, mu = 1/lambda ; largest mu -> smallest positive lambda
        mu, v = spla.eigsh(-Kgf, k=k, M=Kff, which="LA")
        order = np.argsort(-mu)
        lam = 1.0 / mu[order]
        phi = np.zeros((self.ndof, k))
        phi[self.free] = v[:, order]
        return lam, phi

    # ------------------------------------------------------------ checks
    def rigid_body_check(self):
        """Norm of K @ (rigid translation/rotation) relative to K's scale (should be ~0)."""
        c = self.nodes.mean(axis=0)
        out = []
        for axis in range(3):
            u = np.zeros(self.ndof)
            u[axis::6] = 1.0
            out.append(np.linalg.norm(self.K @ u))
            rot = np.zeros(3)
            rot[axis] = 1.0
            u = np.zeros(self.ndof)
            d = np.cross(rot, self.nodes - c)
            for j in range(3):
                u[j::6] = d[:, j]
                u[3 + j::6] = rot[j]
            out.append(np.linalg.norm(self.K @ u))
        scale = np.abs(self.K.diagonal()).max()
        return np.array(out) / scale


def modal_participation(model, phi):
    """Effective mass fractions per mode: translations X, Y, Z and rotations about axes
    through the mass centre parallel to X, Y, Z (so pure bending does not count as rotation)."""
    M = model.M
    total = model.mass()
    rx = np.zeros(model.ndof)
    rx[0::6] = 1.0
    nodal_mass = (M @ rx)[0::6]
    c = nodal_mass @ model.nodes / nodal_mass.sum()
    rows = []
    for j in range(phi.shape[1]):
        p = phi[:, j]
        gen = p @ (M @ p)
        fr = []
        for axis in range(3):
            r = np.zeros(model.ndof)
            r[axis::6] = 1.0
            fr.append((p @ (M @ r))**2 / gen / total)
        for axis in range(3):
            rot = np.zeros(3)
            rot[axis] = 1.0
            r = np.zeros(model.ndof)
            d = np.cross(rot, model.nodes - c)
            for k in range(3):
                r[k::6] = d[:, k]
                r[3 + k::6] = rot[k]
            Irr = r @ (M @ r)
            fr.append((p @ (M @ r))**2 / gen / Irr)
        rows.append(fr)
    return np.array(rows)
