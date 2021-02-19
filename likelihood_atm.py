import matplotlib.pyplot as plt
import numpy as np

from scipy import stats
import pandas as pd


class FullLikelihood(object):
    def __init__(self, like_type):
        self.like  = like_type
        self.corr = []

    def ln_like(self):
        return self.like.ln_like(self.corr)


# Mass/abundance priors

class BrogiLikelihood(object):
    def ln_like(self,corr):
        
        
        like = np.zeros(len(corr["data"]))
        for i  in range(len(corr["data"])):
            N = len(corr["data"][i])
            sf  = np.sqrt(np.var(corr["data"][i]))
            sg  = np.sqrt(np.var(corr["model"][i]))
            corrcoeff = np.corrcoef(corr["data"][i],corr["model"][i])[0,1]
            like[i]= -N/2.*((np.log(sf*sg)+np.log(sf/sg+sg/sf-2.0*corrcoeff))+1)
        return np.sum(like)
        
        
class GaussianMassAbundanceLikelihood(object):
    def ln_like(self, rdict):
        return -0.5 * (rdict["mass"] ** 2 + rdict["abundance"] ** 2)


class StudentMassAbundanceLikelihood(object):
    def __init__(self, nu):
        assert nu > 2.0
        self.nu = nu

    def ln_like(self, rdict):
        return np.sum(stats.t.logpdf([rdict["mass"], rdict["abundance"]], df=self.nu))


# J priors


class CorrelatedMomentsLikelihoodBase(object):
    def __init__(self, corr_matrix, Jmax):
        self.Jmax = Jmax

        J_keep = ["J_{}".format(i) for i in range(Jmax + 1)]
        self.J_list = [j for j in corr_matrix.index if j in J_keep]

        self.corr_matrix = corr_matrix
        assert np.all(np.diag(corr_matrix) == 1.0) and np.all(
            np.abs(corr_matrix) <= 1.0
        )
        assert np.all(corr_matrix.T == corr_matrix)
        assert not np.any(corr_matrix.isna())
        self.C = self.corr_matrix.loc[self.J_list, self.J_list].values
        self.inv_C = np.linalg.inv(self.C)
        self.p = self.C.shape[0]

    def resid_dict_to_J_vector(self, rdict):
        return np.array([rdict[j] for j in self.J_list])

    def resid_corr_norm2(self, rdict):
        rvec = self.resid_dict_to_J_vector(rdict)
        return (rvec).dot(self.inv_C).dot(rvec)


class MultigaussianMomentsLikelihood(CorrelatedMomentsLikelihoodBase):
    def __init__(self, corr_matrix, Jmax):
        super().__init__(corr_matrix, Jmax)

    def ln_like(self, rdict):
        return -0.5 * self.resid_corr_norm2(rdict)


class MultistudentMomentsLikelihood(CorrelatedMomentsLikelihoodBase):
    def __init__(self, corr_matrix, Jmax, nu):
        super().__init__(corr_matrix, Jmax)
        self.nu = nu

    def ln_like(self, rdict):
        return (
            -(self.nu + self.p)
            / 2.0
            * np.log(1.0 + 1.0 / (self.nu - 2.0) * self.resid_corr_norm2(rdict))
        )


# Correlation matrix functions


def make_corr_matrix_identity(Jmax=12, even_only=True):
    J_orders = range(2, Jmax + 1, (2 if even_only else 1))
    values = np.eye(len(J_orders))
    labels = ["J_{}".format(j) for j in J_orders]
    return pd.DataFrame(values, index=labels, columns=labels)


def make_corr_matrix_from_csv(fname):
    corrmat_unnorm = pd.read_csv(fname).set_index("Unnamed: 0")

    # Symmetrize any missing values
    corrmat_unnorm[corrmat_unnorm.isna()] = corrmat_unnorm.T

    # Check diagonal is a constant
    diag = np.diag(corrmat_unnorm.values)
    assert np.all(diag == diag[0])

    # Normalize to diago
    return corrmat_unnorm / diag[0]
