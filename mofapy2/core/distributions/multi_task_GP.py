import gpytorch
from gpytorch.constraints import Positive, GreaterThan
from gpytorch.lazy import DiagLazyTensor,BlockDiagLazyTensor,MatmulLazyTensor, InterpolatedLazyTensor, PsdSumLazyTensor, RootLazyTensor, KroneckerProductLazyTensor, lazify
from gpytorch.utils.broadcasting import _mul_broadcast_shape
from gpytorch.kernels.kernel import Kernel
from typing import Any
from gpytorch.distributions import base_distributions, MultivariateNormal
from gpytorch.functions import add_diag
from gpytorch.likelihoods import Likelihood, _GaussianLikelihoodBase
from gpytorch.utils.warnings import OldVersionWarning
from gpytorch.likelihoods.noise_models import MultitaskHomoskedasticNoise, Noise, HomoskedasticNoise
import warnings
import torch
from torch import Tensor
from gpytorch.mlls.marginal_log_likelihood import MarginalLogLikelihood

# Here we define a multitask GP model where
# - the noise and scale parameter are coupled
# - the index kernel is constrained to a correlation matrix
# - likelihood is replaced by the ELBO

class MultitaskGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood, n_tasks, rank):
        super(MultitaskGPModel, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.MultitaskMean( #TODO set to 0? or avoid centering
            gpytorch.means.ConstantMean(), num_tasks=n_tasks
        )
        self.covar_module = gpytorch.kernels.myMultitaskKernel(
            gpytorch.kernels.RBFKernel(), num_tasks=n_tasks, rank=rank
        )

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultitaskMultivariateNormal(mean_x, covar_x)


#    Based on MultitaskKernel: Using a modified IndexKernel - no other changes
class myMultitaskKernel(Kernel):
    r"""
    Kernel supporting Kronecker style multitask Gaussian processes (where every data point is evaluated at every
    task) using :class:`gpytorch.kernels.IndexKernel` as a basic multitask kernel.

    Given a base covariance module to be used for the data, :math:`K_{XX}`, this kernel computes a task kernel of
    specified size :math:`K_{TT}` and returns :math:`K = K_{TT} \otimes K_{XX}`. as an
    :obj:`gpytorch.lazy.KroneckerProductLazyTensor`.

    :param ~gpytorch.kernels.Kernel data_covar_module: Kernel to use as the data kernel.
    :param int num_tasks: Number of tasks
    :param int rank: (default 1) Rank of index kernel to use for task covariance matrix.
    :param ~gpytorch.priors.Prior task_covar_prior: (default None) Prior to use for task kernel.
        See :class:`gpytorch.kernels.IndexKernel` for details.
    :param dict kwargs: Additional arguments to pass to the kernel.
    """

    def __init__(self, data_covar_module, num_tasks, rank=1, task_covar_prior=None, **kwargs):
        """
        """
        super(myMultitaskKernel, self).__init__(**kwargs)
        self.task_covar_module = myIndexKernel(
            num_tasks=num_tasks, batch_shape=self.batch_shape, rank=rank, prior=task_covar_prior
        )
        self.data_covar_module = data_covar_module
        self.num_tasks = num_tasks

    def forward(self, x1, x2, diag=False, last_dim_is_batch=False, **params):
        if last_dim_is_batch:
            raise RuntimeError("MultitaskKernel does not accept the last_dim_is_batch argument.")
        covar_i = self.task_covar_module.covar_matrix
        if len(x1.shape[:-2]):
            covar_i = covar_i.repeat(*x1.shape[:-2], 1, 1)
        covar_x = lazify(self.data_covar_module.forward(x1, x2, **params))
        res = KroneckerProductLazyTensor(covar_x, covar_i)
        return res.diag() if diag else res

    def num_outputs_per_input(self, x1, x2):
        """
        Given `n` data points `x1` and `m` datapoints `x2`, this multitask
        kernel returns an `(n*num_tasks) x (m*num_tasks)` covariance matrix.
        """
        return self.num_tasks

# implement an index kernel B^TB +sigma I
# based on IndexKernel with changes:
# - a constant diagonal offset
# - scaled to correlation matrix
class myIndexKernel(Kernel):
    r"""
    A kernel for discrete indices. Kernel is defined by a lookup table.

    .. math::

        \begin{equation}
            k(i, j) = \left(BB^\top + \text{diag}(\mathbf v) \right)_{i, j}
        \end{equation}

    where :math:`B` is a low-rank matrix, and :math:`\mathbf v` is a  non-negative vector.
    These parameters are learned.

    Args:
        :attr:`num_tasks` (int):
            Total number of indices.
        :attr:`batch_shape` (torch.Size, optional):
            Set if the MultitaskKernel is operating on batches of data (and you want different
            parameters for each batch)
        :attr:`rank` (int):
            Rank of :math:`B` matrix.
        :attr:`prior` (:obj:`gpytorch.priors.Prior`):
            Prior for :math:`B` matrix.
        :attr:`var_constraint` (Constraint, optional):
            Constraint for added diagonal component. Default: `Positive`.

    Attributes:
        covar_factor:
            The :math:`B` matrix.
        lov_var:
            The element-wise log of the :math:`\mathbf v` vector.
    """

    def __init__(self, num_tasks, rank=1, prior=None, var_constraint=None, **kwargs):
        if rank > num_tasks:
            raise RuntimeError("Cannot create a task covariance matrix larger than the number of tasks")
        super().__init__(**kwargs)
        self.num_tasks = num_tasks

        if var_constraint is None:
            var_constraint = Positive()

        self.register_parameter(
            name="covar_factor", parameter=torch.nn.Parameter(torch.randn(*self.batch_shape, num_tasks, rank))
        )
        self.register_parameter(name="raw_var",
                                parameter=torch.nn.Parameter(torch.randn(*self.batch_shape, 1)))  # 1 instead of num_tasks
        if prior is not None:
            self.register_prior("IndexKernelPrior", prior, self._eval_covar_matrix)

        self.register_constraint("raw_var", var_constraint)

    @property
    def var(self):
        return self.raw_var_constraint.transform(self.raw_var)

    @var.setter
    def var(self, value):
        self._set_var(value)

    def _set_var(self, value):
        self.initialize(raw_var=self.raw_var_constraint.inverse_transform(value))

    def _eval_covar_matrix(self):
        cf = self.covar_factor
        C = cf @ cf.transpose(-1, -2) + self.var *  torch.eye(self.num_tasks)   # instead of diag(var)
        Cdiag = torch.diag(C)
        C = torch.diag(1 / Cdiag.sqrt()) @ C @ torch.diag(1 / Cdiag.sqrt())       # scale to correlation matrix
        return C

    @property
    def covar_matrix(self):
        var = self.var
        res = PsdSumLazyTensor(RootLazyTensor(self.covar_factor), DiagLazyTensor(var))
        return res

    def forward(self, i1, i2, **params):
        covar_matrix = self._eval_covar_matrix()
        batch_shape = _mul_broadcast_shape(i1.shape[:-2], self.batch_shape)
        index_shape = batch_shape + i1.shape[-2:]

        res = InterpolatedLazyTensor(
            base_lazy_tensor=covar_matrix,
            left_interp_indices=i1.expand(index_shape),
            right_interp_indices=i2.expand(index_shape),
        )
        return res



class ELBO(MarginalLogLikelihood):
    """
    Class to implment ELBO as loss function - based on ExactMarginalLogLikelihood
    """

    def __init__(self, likelihood, model, ZCov):
        if not isinstance(likelihood, _GaussianLikelihoodBase):
            raise RuntimeError("Likelihood must be Gaussian for exact inference")
        super(ELBO, self).__init__(likelihood, model)
        self.ZCov = ZCov

    def forward(self, function_dist, target, *params):
        r"""
        Computes the MLL given :math:`p(\mathbf f)` and :math:`\mathbf y`.

        :param ~gpytorch.distributions.MultivariateNormal function_dist: :math:`p(\mathbf f)`
            the outputs of the latent function (the :obj:`gpytorch.models.ExactGP`)
        :param torch.Tensor target: :math:`\mathbf y` The target values
        :rtype: torch.Tensor
        :return: Exact MLL. Output shape corresponds to batch shape of the model/input data.
        """
        if not isinstance(function_dist, MultivariateNormal):
            raise RuntimeError("ExactMarginalLogLikelihood can only operate on Gaussian random variables")

        # Get the log prob of the marginal distribution
        output = self.likelihood(function_dist, *params)
        res = output.log_prob(target)

        # Add additional terms (SGPR / learned inducing points, heteroskedastic likelihood models)
        for added_loss_term in self.model.added_loss_terms():
            res = res.add(added_loss_term.loss(*params))

        # Add uncertainty in Z
        res.add(-0.5 * torch.trace(Sigma^-1 @ self.ZCov)) #TODO: How is Sigma-1 treated in likleihood part

        # Add log probs of priors on the (functions of) parameters
        for _, prior, closure, _ in self.named_priors():
            res.add_(prior.log_prob(closure()).sum())

        # Scale by the amount of data we have
        num_data = target.size(-1)
        return res.div_(num_data)
