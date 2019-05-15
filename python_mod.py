import numpy as np
import scipy.linalg as la


class Takagi():

    def __init__(self, matrix):
        diag_matrix, orthonormal_vector = self.takagi_factorization(matrix)
        self.diag = diag_matrix
        self.ortho_vec = orthonormal_vector

    def takagi_factorization(self, matrix):
        """Extremely simple and inefficient Takagi factorization of a
        symmetric, complex matrix A.  Here we take this to mean A = U D U^T
        where D is a real, diagonal matrix and U is a unitary matrix.  There
        is no guarantee that it will always work. """
        # Construct a Hermitian matrix.
        H = np.dot(matrix.T.conj(), matrix)
        # Calculate the eigenvalue decomposition of the Hermitian matrix.
        # The diagonal matrix in the Takagi factorization is the square
        # root of the eigenvalues of this matrix.
        (lam, u) = la.eigh(H)
        # The "almost" Takagi factorization.  There is a conjugate here
        # so the final form is as given in the doc string.
        T = np.dot(u.T, np.dot(matrix, u)).conj()
        # T is diagonal but not real.  That is easy to fix by a
        # simple transformation which removes the complex phases
        # from the resulting diagonal matrix.
        c = np.diag(np.exp(-1j*np.angle(np.diag(T))/2))
        U = np.dot(u, c)
        # Now A = np.dot(U, np.dot(np.diag(np.sqrt(lam)),U.T))

        return np.diag(np.sqrt(lam)), U
