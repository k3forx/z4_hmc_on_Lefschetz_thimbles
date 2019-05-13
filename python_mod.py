import numpy as np
import scipy.linalg as la


class Takagi():

    def __init__(self, matrix):
        self.diag = self.diag(matrix)
        self.vec = self.ortho_vec(matrix)

    def diag(self, matrix):
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
        c = np.diag(np.exp(-1j*np.angle(np.diag(T))/2))
        U = np.dot(u, c)
        # Now A = np.dot(U, np.dot(np.diag(np.sqrt(lam)),U.T))
        # obj.diag = np.sqrt(lam)
        return np.diag(np.sqrt(lam))

    def ortho_vec(self, matrix):
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
        T = np.dot(u.T, np.dot(A, u)).conj()
        # T is diagonal but not real.  That is easy to fix by a
        # simple transformation which removes the complex phases
        # from the resulting diagonal matrix.
        c = np.diag(np.exp(-1j*np.angle(np.diag(T))/2))
        U = np.dot(u, c)
        # Now A = np.dot(U, np.dot(np.diag(np.sqrt(lam)),U.T))
        # obj.diag = np.sqrt(lam)
        return U


def takagi_factorization_diag(matrix):
    H = np.dot(matrix.T.conj(), matrix)
    (lam, u) = la.eigh(H)

    return np.diag(np.sqrt(lam))


def takagi_factorization_ortho_vec(matrix):
    H = np.dot(matrix.T.conj(), matrix)
    (lam, u) = la.eigh(H)
    T = np.dot(u.T, np.dot(A, u)).conj()
    c = np.diag(np.exp(-1j*np.angle(np.diag(T))/2))
    U = np.dot(u, c)

    return U


if __name__ == '__main__':
    A = np.zeros(shape=(1, 1), dtype=np.complex)

    A[0, 0] = complex(1.6593716127412383, -1.0520139131468429)
    # A[0, 1] = complex(0, 2)
    # A[1, 0] = complex(0, 2)
    # A[1, 1] = complex(3, 0)

    # eival, eivec = takagi(A)

    # print(eival)
    # print(eivec * np.conj(eivec))
    # print(np.dot(eivec, np.dot(np.diag(eival), eivec.T)))
    takagi = Takagi(A)
    print(takagi.diag)
    print(format(takagi.vec))
    print(format(takagi_factorization_ortho_vec(A)[0, 0], ' .16E'))
