import logging
import os
import sys
import pandas as pd
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.neighbors import KNeighborsRegressor, KNeighborsClassifier
from sklearn.preprocessing import PolynomialFeatures
from scipy.stats import norm as normal
from sklearn.neighbors import NearestNeighbors
def scatter_value_to_grid_value(embedding, grid, value, method="knn", n_knn=30, n_poly=3):

    """
    Transfer value on 2D scatter into 2d Grid space.
    Args:
        embedding (np.array): shape = (n_cell, 2)
        grid (np.array): shape = (n_grid**2, 2)
        value (np.array): shape = (n_cell, )
    Returns:
        np.array : shape = (n_grid**2, 1)
    """

    x, y = embedding[:, 0], embedding[:, 1]
    x_new, y_new = grid[:, 0], grid[:, 1]

    if method == "poly":
        value_on_grid =  _polynomial_regression_old_ver(x, y, x_new, y_new, value, n_degree=n_poly)
    if method == "polynomial":
        value_on_grid =  _polynomial_regression_sklearn(x, y, x_new, y_new, value, n_degree=n_poly)
    elif method == "knn":
        value_on_grid = _knn_regression(x, y, x_new, y_new, value, n_knn=n_knn)
    elif method == "knn_class":
        value_on_grid = _knn_classification(x, y, x_new, y_new, value, n_knn=n_knn)


    return value_on_grid


def _polynomial_regression_sklearn(x, y, x_new, y_new, value, n_degree=3):

    # Make polynomial features
    data = np.stack([x, y], axis=1)
    data_new = np.stack([x_new, y_new], axis=1)

    pol = PolynomialFeatures(degree=n_degree, include_bias=False)
    data = pol.fit_transform(data)
    data_new = pol.transform(data_new)


    model = Ridge(random_state=123)
    model.fit(data, value)

    return model.predict(data_new)

def _polynomial_regression_old_ver(x, y, x_new, y_new, value, n_degree=3):

    def __conv(x, y, n_degree=3): # Make polynomial data for polynomial ridge regression
        dic = {}
        for d in range(1, n_degree + 1):
            for i, name in zip([x, y], ["x", "y"]):
                dic[name + str(d)] = i**d
        return pd.DataFrame(dic)

    data = __conv(x=x, y=y, n_degree=n_degree)

    model = Ridge()
    model.fit(data, value)

    data_new = __conv(x=x_new, y=y_new, n_degree=n_degree)

    return model.predict(data_new)

def _knn_regression(x, y, x_new, y_new, value, n_knn=30):

    data = np.stack([x, y], axis=1)

    model = KNeighborsRegressor(n_neighbors=n_knn)
    model.fit(data, value)

    data_new = np.stack([x_new, y_new], axis=1)

    return model.predict(data_new)

def _knn_classification(x, y, x_new, y_new, value, n_knn=30):

    data = np.stack([x, y], axis=1)

    model = KNeighborsClassifier(n_neighbors=n_knn)
    model.fit(data, value)

    data_new = np.stack([x_new, y_new], axis=1)

    return model.predict(data_new)

def calculate_p_mass(embedding, smooth=0.5, steps=(40, 40),
                          n_neighbors=100, n_jobs=4, xylim=((None, None), (None, None))):
    """Calculate the velocity using a points on a regular grid and a gaussian kernel
    Note: the function should work also for n-dimensional grid
    Arguments
    ---------
    embedding:
    smooth: float, smooth=0.5
        Higher value correspond to taking in consideration further points
        the standard deviation of the gaussian kernel is smooth * stepsize
    steps: tuple, default
        the number of steps in the grid for each axis
    n_neighbors:
        number of neighbors to use in the calculation, bigger number should not change too much the results..
        ...as soon as smooth is small
        Higher value correspond to slower execution time
    n_jobs:
        number of processes for parallel computing
    xymin:
        ((xmin, xmax), (ymin, ymax))
    Returns
    -------
    total_p_mass: np.ndarray
        density at each point of the grid
    """

    # Prepare the grid
    grs = []
    for dim_i in range(embedding.shape[1]):
        m, M = np.min(embedding[:, dim_i]), np.max(embedding[:, dim_i])

        if xylim[dim_i][0] is not None:
            m = xylim[dim_i][0]
        if xylim[dim_i][1] is not None:
            M = xylim[dim_i][1]

        m = m - 0.025 * np.abs(M - m)
        M = M + 0.025 * np.abs(M - m)
        gr = np.linspace(m, M, steps[dim_i])
        grs.append(gr)

    meshes_tuple = np.meshgrid(*grs)
    gridpoints_coordinates = np.vstack([i.flat for i in meshes_tuple]).T

    nn = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=n_jobs)
    nn.fit(embedding)
    dists, neighs = nn.kneighbors(gridpoints_coordinates)

    std = np.mean([(g[1] - g[0]) for g in grs])
    # isotropic gaussian kernel
    gaussian_w = normal.pdf(loc=0, scale=smooth * std, x=dists)
    total_p_mass = gaussian_w.sum(1)
    gridpoints_coordinates

    return total_p_mass, gridpoints_coordinates

def normalize_gradient(gradient, method="sqrt"):
    """
    Normalize length of 2D vector
    """

    if method == "sqrt":

        size = np.sqrt(np.power(gradient, 2).sum(axis=1))
        size_sq = np.sqrt(size)
        #size_sq = (size)
        size_sq[size_sq == 0] = 1
        factor = np.repeat(np.expand_dims(size_sq, axis=1), 2, axis=1)
        
    return gradient / factor

def get_gradient(value_on_grid):
    # Gradient calculation
    n = int(np.sqrt(value_on_grid.shape[0]))
    value_on_grid_as_matrix = value_on_grid.reshape(n, n)
    dy, dx = np.gradient(value_on_grid_as_matrix)
    gradient = np.stack([dx.flatten(), dy.flatten()], axis=1)

    return gradient
def visualize_dev_flow(self, scale_for_pseudotime=30, s=10, s_grid=30,cmap='rainbow'):

    embedding_whole = self.embedding_whole
    embedding_of_interest= self.embedding
    mass_filter = self.mass_filter
    mass_filter_whole = self.mass_filter_whole
    gridpoints_coordinates=self.gridpoints_coordinates

    pseudotime_raw = self.pseudotime
    pseudotime_on_grid=self.pseudotime_on_grid

    gradient_pseudotime=self.gradient


    fig, ax = plt.subplots(1, 5, figsize=[25,5])

    ##
    ax_ = ax[0]
    ax_.scatter(embedding_whole[:, 0], embedding_whole[:, 1], c="lightgray", s=s)
    ax_.scatter(embedding_of_interest[:, 0], embedding_of_interest[:, 1], c=pseudotime_raw, cmap=cmap, s=s)
    ax_.set_title("Pseudotime")
    ax_.axis("off")

    ####
    ax_ = ax[1]
    ax_.scatter(gridpoints_coordinates[mass_filter, 0], gridpoints_coordinates[mass_filter, 1], s=0)
    ax_.scatter(gridpoints_coordinates[~mass_filter_whole, 0], gridpoints_coordinates[~mass_filter_whole, 1],
     c="lightgray", s=s_grid)
    ax_.scatter(gridpoints_coordinates[~mass_filter, 0], gridpoints_coordinates[~mass_filter, 1],
     c=pseudotime_on_grid[~mass_filter], cmap=cmap, s=s_grid)
    ax_.set_title("Pseudotime on grid")
    ax_.axis("off")



    ###
    ax_ = ax[2]
    #ax_.scatter(gridpoints_coordinates[mass_filter, 0], gridpoints_coordinates[mass_filter, 1], s=0)
    ax_.scatter(gridpoints_coordinates[mass_filter, 0], gridpoints_coordinates[mass_filter, 1], s=0)
    ax_.scatter(gridpoints_coordinates[~mass_filter_whole, 0], gridpoints_coordinates[~mass_filter_whole, 1],
     c="lightgray", s=s_grid)
    ax_.scatter(gridpoints_coordinates[~mass_filter, 0], gridpoints_coordinates[~mass_filter, 1],
     c=pseudotime_on_grid[~mass_filter], cmap=cmap, s=s_grid)

    ax_.quiver(gridpoints_coordinates[~mass_filter, 0], gridpoints_coordinates[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    ax_.set_title("Gradient of pseudotime \n(=Development flow)")
    ax_.axis("off")

    ###
    ax_ = ax[3]
    #ax_.scatter(gridpoints_coordinates[mass_filter, 0], gridpoints_coordinates[mass_filter, 1], s=0)
    ax_.scatter(embedding_whole[:, 0], embedding_whole[:, 1], c="lightgray", s=s)
    ax_.quiver(gridpoints_coordinates[~mass_filter, 0], gridpoints_coordinates[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)
    ax_.set_title("Gradient of pseudotime \n(=Development flow)")
    ax_.axis("off")

    ####
    ax_ = ax[4]
    ax_.scatter(embedding_whole[:, 0], embedding_whole[:, 1], c="lightgray", s=s)
    ax_.scatter(embedding_of_interest[:, 0], embedding_of_interest[:, 1], c=pseudotime_raw, cmap=cmap, s=s)

    '''ax_.quiver(gridpoints_coordinates[~mass_filter, 0], gridpoints_coordinates[~mass_filter, 1],
                    gradient_pseudotime[~mass_filter, 0], gradient_pseudotime[~mass_filter, 1],
                    scale=scale_for_pseudotime)'''
    x_order=np.argsort(gridpoints_coordinates[~mass_filter,0])
    y_order=np.argsort(gridpoints_coordinates[~mass_filter,1])    
    ax_.streamplot(gridpoints_coordinates[~mass_filter,0][x_order], gridpoints_coordinates[~mass_filter,0][y_order],
                   gradient_pseudotime[~mass_filter, 0][x_order], gradient_pseudotime[~mass_filter, 1][y_order],density=1)
    ax_.set_title("Pseudotime + \nDevelopment flow")
    ax_.axis("off")
    
    
class Gradient_based_trajecory():
    def __init__(self, adata=None, obsm_key=None, pseudotime_key="pseudotime", cluster_column_name=None, cluster=None, gt=None):

        if adata is not None:
            self.load_adata(adata=adata, obsm_key=obsm_key,
            pseudotime_key=pseudotime_key,cluster_column_name=cluster_column_name,
            cluster=cluster)
        elif gt is not None:
            self.embedding = gt.embedding_whole.copy()
            self.embedding_whole = gt.embedding_whole.copy()
            self.mass_filter = gt.mass_filter_whole.copy()
            self.mass_filter_whole = gt.mass_filter_whole.copy()
            self.gridpoints_coordinates = gt.gridpoints_coordinates.copy()
            self.pseudotime = gt.pseudotime_whole.copy()

    def load_adata(self, adata, obsm_key, pseudotime_key, cluster_column_name=None, cluster=None):

        self.embedding = adata.obsm[obsm_key]
        self.pseudotime = adata.obs[pseudotime_key].values
        self.embedding_whole = self.embedding.copy()
        self.pseudotime_whole = self.pseudotime.copy()

        if (cluster_column_name is not None) & (cluster is not None):
            cells_ix = np.where(adata.obs[cluster_column_name] == cluster)[0]
            self.embedding = self.embedding[cells_ix, :]
            self.pseudotime = self.pseudotime[cells_ix]

    def calculate_mass_filter(self, min_mass=0.01, smooth=0.8, steps=(40, 40), n_neighbors=200, n_jobs=4):

        x_min, y_min = self.embedding_whole.min(axis=0)
        x_max, y_max = self.embedding_whole.max(axis=0)
        xylim = ((x_min, x_max), (y_min, y_max))

        total_p_mass, gridpoints_coordinates = calculate_p_mass(self.embedding, smooth=smooth, steps=steps,
                                  n_neighbors=n_neighbors, n_jobs=n_jobs, xylim=xylim)

        total_p_mass_whole, _ = calculate_p_mass(self.embedding_whole, smooth=smooth, steps=steps,
                                  n_neighbors=n_neighbors, n_jobs=n_jobs, xylim=xylim)

        self.total_p_mass = total_p_mass
        self.mass_filter = (total_p_mass < min_mass)
        self.mass_filter_whole = (total_p_mass_whole < min_mass)
        self.gridpoints_coordinates = gridpoints_coordinates

    def transfer_data_into_grid(self, args={}):

        if not args:
            args = {"method": "knn",
                    "n_knn": 30}

        self.pseudotime_on_grid = scatter_value_to_grid_value(embedding=self.embedding,
                                                          grid=self.gridpoints_coordinates,
                                                          value=self.pseudotime,
                                                          **args)
    def calculate_gradient(self, scale_factor=60, normalization=None):

        # Gradient calculation
        gradient = get_gradient(value_on_grid=self.pseudotime_on_grid.copy())

        if normalization == "sqrt":
            gradient = normalize_gradient(gradient, method="sqrt")

        if scale_factor == "l2_norm_mean":
            # divide gradient by the mean of l2 norm.
            l2_norm = np.linalg.norm(gradient, ord=2, axis=1)
            scale_factor = 1 / l2_norm.mean()

        self.gradient = gradient * scale_factor

    def visualize_dev_flow(self, scale_for_pseudotime=30, s=10, s_grid=30):
        visualize_dev_flow(self, scale_for_pseudotime=scale_for_pseudotime, s=s, s_grid=s_grid)


CONFIG = {"default_args": {"lw": 0.1, "rasterized": True},
          "s_scatter": 50,
          "s_grid": 10,
          "scale_simulation": 50,
          "scale_dev": 300,
          "cmap_ps": "PiYG",
          #"density"
          "scale_unit":'xy',
          "default_args_quiver": {"linewidths": 0.25, "width": 0.004}}


def plot_background(self, ax=None, s=CONFIG["s_scatter"], args=CONFIG["default_args"]):
    if ax is None:
        ax = plt
    ax.scatter(self.embedding[:, 0], self.embedding[:, 1], c="lightgray", s=s, **args)
    #ax.set_title("Pseudotime")
    ax.axis("off")
def plot_background_on_grid(self, ax=None, s=CONFIG["s_grid"], args={}):
    if ax is None:
        ax = plt
    if hasattr(self, "mass_filter_whole_reference"):
        mass_filter = self.mass_filter_whole_reference
    elif hasattr(self, "mass_filter_whole"):
        mass_filter = self.mass_filter_whole
    ax.scatter(self.gridpoints_coordinates[:, 0],
               self.gridpoints_coordinates[:, 1], s=0)
    if "c" not in args.keys():
        ax.scatter(self.gridpoints_coordinates[~mass_filter, 0],
               self.gridpoints_coordinates[~mass_filter, 1],
               c="lightgray", s=s, **args)
    else:
        ax.scatter(self.gridpoints_coordinates[~mass_filter, 0],
               self.gridpoints_coordinates[~mass_filter, 1],
               s=s, **args)
    ax.axis("off")

def plot_cluster_whole(self, ax=None, s=CONFIG["s_scatter"], args=CONFIG["default_args"]):
    if ax is None:
        ax = plt
    ax.scatter(self.embedding[:, 0], self.embedding[:, 1], c=self.colorandum, s=s, **args)
    ax.axis("off")

def plot_cluster_cells_use(self, ax=None, s=CONFIG["s_scatter"], color=None, show_background=True, args=CONFIG["default_args"]):
    if ax is None:
        ax = plt
    if s == 0:
        color = "white"
    if show_background:
        plot_background(self=self, ax=ax, s=s, args=args)
    if not hasattr(self, "cell_idx_use"):
        self.cell_idx_use = None
    if self.cell_idx_use is None:
        if color is None:
            ax.scatter(self.embedding[:, 0], self.embedding[:, 1], c=self.colorandum, s=s, **args)
        else:
            ax.scatter(self.embedding[:, 0], self.embedding[:, 1], c=color, s=s, **args)
    else:
        if color is None:
            ax.scatter(self.embedding[self.cell_idx_use, 0], self.embedding[self.cell_idx_use, 1],
                       c=self.colorandum[self.cell_idx_use, :],s=s, **args)
        else:
            ax.scatter(self.embedding[self.cell_idx_use, 0], self.embedding[self.cell_idx_use, 1],
                       c=color, s=s, **args)
    ax.axis("off")

def plot_pseudotime(self, ax=None, s=CONFIG["s_scatter"], show_background=True, cmap="rainbow", args=CONFIG["default_args"]):
    if ax is None:
        ax = plt
    if show_background:
        plot_background(self=self, ax=ax, s=s, args=args)
    if self.cell_idx_use is None:
        ax.scatter(self.embedding[:, 0], self.embedding[:, 1], c=self.pseudotime, cmap=cmap, s=s, **args)
        ax.legend()
    else:
        ax.scatter(self.embedding[self.cell_idx_use, 0], self.embedding[self.cell_idx_use, 1],
                    c=self.pseudotime[self.cell_idx_use], cmap=cmap, s=s, **args)
    ax.axis("off")

def plot_pseudotime_on_grid(self, ax=None, s=CONFIG["s_grid"], show_background=True, cmap="rainbow", args={}):
    if ax is None:
        ax = plt
    if hasattr(self, "mass_filter_simulation"):
        mass_filter = self.mass_filter_simulation
    elif hasattr(self, "mass_filter"):
        mass_filter = self.mass_filter
    if show_background:
        plot_background_on_grid(self=self, ax=ax, s=s, args=args)
    else:
        plot_cluster_cells_use(self=self, ax=ax, s=0, color="white", show_background=False, args={})
    ax.scatter(self.gridpoints_coordinates[~mass_filter, 0],
               self.gridpoints_coordinates[~mass_filter, 1],
               c=self.pseudotime_on_grid[~mass_filter],
               cmap=cmap, s=s, **args)
    ax.axis("off")

def plot_selected_pseudotime_on_grid(self, ax=None, pseudotime_selected=[], s=CONFIG["s_grid"], show_background=True, args={}):
    if ax is None:
        ax = plt
    mass_filter = self.mass_filter_simulation
    if show_background:
        plot_background_on_grid(self=self, ax=ax, s=s,
                                args={"facecolor": "None",
                                      "c": "None",
                                      "edgecolors":'None',
                                      "linewidths": 0.005})
    else:
        plot_cluster_cells_use(self=self, ax=ax, s=0, color="white", show_background=False, args={})
    x = self.gridpoints_coordinates[~mass_filter, 0]
    y = self.gridpoints_coordinates[~mass_filter, 1]
    for label, color in zip(["True", "False"], ["#EC7063", "#D0D3D4"]):
        if label == "True":
            idx = self.inner_product_df.pseudotime_id.isin(pseudotime_selected).values
        else:
            idx = ~self.inner_product_df.pseudotime_id.isin(pseudotime_selected).values
        ax.scatter(x[idx], y[idx], color=color, label=label, s=s, **args)

    ax.legend()
    ax.axis("off")
    
def plot_reference_flow_on_grid(self,plot_type="quiver", ax=None, scale=CONFIG["scale_dev"], density=0.8,show_background=True, s=CONFIG["s_scatter"], args=CONFIG["default_args_quiver"]):
    if ax is None:
        ax = plt
    if hasattr(self, "mass_filter_simulation"):
        mass_filter = self.mass_filter_simulation
    elif hasattr(self, "mass_filter"):
        mass_filter = self.mass_filter
    if show_background:
        plot_background(self=self, ax=ax, s=s, args=CONFIG["default_args"])
    else:
        plot_cluster_cells_use(self=self, ax=ax, s=0, color="white", show_background=False, args={})
    
    if plot_type=="quiver":
        Q=ax.quiver(self.gridpoints_coordinates[~mass_filter, 0],
              self.gridpoints_coordinates[~mass_filter, 1],
              self.ref_flow[~mass_filter, 0],
              self.ref_flow[~mass_filter, 1],scale=scale,**args)
        ax.axis("off")
        
    if plot_type=="stream":
    
        X_grid=self.gridpoints_coordinates[:,0].reshape(self.n_grid,self.n_grid)
        Y_grid=self.gridpoints_coordinates[:,1].reshape(self.n_grid,self.n_grid)
        masked=mass_filter.reshape(self.n_grid,self.n_grid)
        U_grid=self.ref_flow[:,0].reshape(self.n_grid,self.n_grid)
        V_grid=self.ref_flow[:,1].reshape(self.n_grid,self.n_grid)
        U_grid = np.ma.array(U_grid, mask=masked)
        V_grid = np.ma.array(V_grid, mask=masked)
        ax.streamplot(X_grid, Y_grid,U_grid,V_grid,density=density,color='black',integration_direction="both",linewidth=0.8)
        ax.axis("off")

class Gradient_calculator():
    def __init__(self, oracle_object=None, adata=None, obsm_key=None, pseudotime_key="Pseudotime", cell_idx_use=None, name=None, gt=None):
        """
        Estimate the direction of differentiation by calculation gradient of pseudotime on the embedding space.
        Please look at web tutorial for example scripts.
        Args:
            adata (anndata): scRNA-seq data in anndata class
            obsm_key (str): Name of dimensional reduction. You can check the list of dimensional reduction data name with "adata.obsm.keys()"
            pseudotime_key (str): Pseudotime data should be stored in adata.obs[pseudotime_key]. Please set the name of pseudotime data in adata.obs
            cluster_column_name (str): If you set cluster_column_name and cluster, you can subset cells to calculate gradient.
                Please look at web tutorial for example codes.
            cluster (str): See above.
        """
        self.cell_idx_use = None
        self.n_neighbors = None
        self.min_mass = None
        self.smooth = None
        self.n_grid = None
        if oracle_object is not None:
            self.load_oracle_object(oracle_object=oracle_object,
                                    cell_idx_use=cell_idx_use,
                                    name=name,
                                    pseudotime_key=pseudotime_key)
        elif adata is not None:
            self.load_adata(adata=adata, obsm_key=obsm_key,pseudotime_key=pseudotime_key,cell_idx_use=cell_idx_use,
                            name=name)
        elif gt is not None:
            self.embedding = gt.embedding.copy()
            self.mass_filter = gt.mass_filter_whole.copy()
            self.mass_filter_whole = gt.mass_filter_whole.copy()
            self.gridpoints_coordinates = gt.gridpoints_coordinates.copy()

            self.n_neighbors = gt.n_neighbors
            self.min_mass = gt.min_mass
            self.smooth = gt.smooth
            self.n_grid = gt.n_grid
    def load_adata(self, adata, obsm_key, cell_idx_use=None, name=None, pseudotime_key="Pseudotime"):
        self.name = name
        self.embedding = adata.obsm[obsm_key].copy()
        self.pseudotime = adata.obs[pseudotime_key].values.copy()
        if cell_idx_use is not None:
            self.cell_idx_use = np.array(cell_idx_use)

    def calculate_p_mass(self, smooth=0.8, n_grid=40, n_neighbors=200, n_jobs=-1):

        x_min, y_min = self.embedding.min(axis=0)
        x_max, y_max = self.embedding.max(axis=0)
        xylim = ((x_min, x_max), (y_min, y_max))
        steps = (n_grid, n_grid)

        if self.cell_idx_use is None:

            total_p_mass, gridpoints_coordinates = calculate_p_mass(self.embedding, smooth=smooth, steps=steps,
                                  n_neighbors=n_neighbors, n_jobs=n_jobs, xylim=xylim)
            total_p_mass_whole = total_p_mass.copy()

        else:
            total_p_mass, gridpoints_coordinates = calculate_p_mass(self.embedding[self.cell_idx_use, :], smooth=smooth, steps=steps,
                                  n_neighbors=n_neighbors, n_jobs=n_jobs, xylim=xylim)

            total_p_mass_whole, _ = calculate_p_mass(self.embedding, smooth=smooth, steps=steps,
                                  n_neighbors=n_neighbors, n_jobs=n_jobs, xylim=xylim)

        self.n_neighbors = n_neighbors
        self.smooth = smooth
        self.n_grid = n_grid

        self.total_p_mass = total_p_mass
        self.total_p_mass_whole = total_p_mass_whole
        self.gridpoints_coordinates = gridpoints_coordinates
        
    def suggest_mass_thresholds(self, n_suggestion=12, s=1, n_col=4):
        min_ = self.total_p_mass.min()
        max_ = self.total_p_mass.max()
        suggestions = np.linspace(min_, max_/2, n_suggestion)
        n_rows = math.ceil(n_suggestion / n_col)
        fig, ax = plt.subplots(n_rows, n_col, figsize=[5*n_col, 5*n_rows])
        if n_rows == 1:
            ax = ax.reshape(1, -1)
        row = 0
        col = 0
        for i in range(n_suggestion):
            ax_ = ax[row, col]
            col += 1
            if col == n_col:
                col = 0
                row += 1
            idx = self.total_p_mass > suggestions[i]
                #ax_.scatter(gridpoints_coordinates[mass_filter, 0], gridpoints_coordinates[mass_filter, 1], s=0)
            ax_.scatter(self.embedding[:, 0], self.embedding[:, 1], c="lightgray", s=s)
            ax_.scatter(self.gridpoints_coordinates[idx, 0],
                       self.gridpoints_coordinates[idx, 1],
                       c="black", s=s)
            ax_.set_title(f"min_mass: {suggestions[i]: .2g}")
            ax_.axis("off")
    def calculate_mass_filter(self, min_mass=0.01, plot=False):
        self.min_mass = min_mass
        self.mass_filter = (self.total_p_mass < min_mass)
        self.mass_filter_whole = (self.total_p_mass_whole < min_mass)
        if plot:
            fig, ax = plt.subplots(figsize=[5,5])
            #ax_.scatter(gridpoints_coordinates[mass_filter, 0], gridpoints_coordinates[mass_filter, 1], s=0)
            ax.scatter(self.embedding[:, 0], self.embedding[:, 1], c="lightgray", s=10)
            ax.scatter(self.gridpoints_coordinates[~self.mass_filter, 0],
                       self.gridpoints_coordinates[~self.mass_filter, 1],
                       c="black", s=0.5)
            ax.set_title("Grid points selected")
            ax.axis("off")
    def transfer_data_into_grid(self, args={}, plot=False):
        if not args:
            args = {"method": "knn",
                    "n_knn": 30}
        # Prepare input data_new
        if self.cell_idx_use is None:
            embedding = self.embedding
            grid = self.gridpoints_coordinates
            value = self.pseudotime
        else:
            embedding = self.embedding[self.cell_idx_use, :]
            grid = self.gridpoints_coordinates
            value = self.pseudotime[self.cell_idx_use]
        # Remove inf
        if np.inf in value:
            # Clip inf
            warnings.warn("Inf value found in the pseudotime data. The inf value is replaced with non-inf max value.", UserWarning)
            _clip_inf_value(data=value)
        # Data calculation for each grid point
        self.pseudotime_on_grid = scatter_value_to_grid_value(embedding=embedding,
                                                              grid=grid,
                                                              value=value,
                                                              **args)
        if plot:
            fig, ax = plt.subplots(1, 2, figsize=[10,5])
            s = 10
            s_grid = 20
            show_background = True
            ##
            ax_ = ax[0]
            plot_pseudotime(self, ax=ax_, s=s, show_background=show_background)
            ax_.set_title("Pseudotime")
            ####
            ax_ = ax[1]
            plot_pseudotime_on_grid(self, ax=ax_, s=s_grid, show_background=show_background)
            ax_.set_title("Pseudotime on grid")

    def calculate_gradient(self, scale_factor="l2_norm_mean", normalization="sqrt"):

        # Gradient calculation
        gradient = get_gradient(value_on_grid=self.pseudotime_on_grid.copy())
        if normalization == "sqrt":
            gradient = normalize_gradient(gradient, method="sqrt")
        if scale_factor == "l2_norm_mean":
            # divide gradient by the mean of l2 norm.
            l2_norm = np.linalg.norm(gradient, ord=2, axis=1)
            scale_factor = 1 / l2_norm.mean()
        self.ref_flow = gradient * scale_factor
    def plot_dev_flow_on_grid(self, ax=None, scale=CONFIG["scale_dev"],show_background=True, s=CONFIG["s_scatter"], args={}):
        plot_reference_flow_on_grid(self, ax=ax, scale=scale, show_background=show_background, s=s, args=args)
    def plot_reference_flow_on_grid(self, plot_type="quiver",ax=None, scale=CONFIG["scale_dev"],density=2,show_background=True, s=CONFIG["s_scatter"], args={}):
        plot_reference_flow_on_grid(self, ax=ax, scale=scale,show_background=show_background, s=s, args=args)
    def visualize_results(self, plot_type="quiver",scale=30, s=1, s_grid=30, show_background=True,cmap='rainbow'):
        fig, ax = plt.subplots(1, 5, figsize=[25,5])
        ax_ = ax[0]
        plot_pseudotime(self, ax=ax_, s=s, show_background=show_background,cmap=cmap)
        ax_.set_title("Pseudotime")
        ax_ = ax[1]
        plot_pseudotime_on_grid(self, ax=ax_, s=s_grid, show_background=show_background,cmap=cmap)
        ax_.set_title("Pseudotime on grid")
        ax_ = ax[2]
        plot_pseudotime_on_grid(self, ax=ax_, s=s_grid, show_background=show_background,cmap=cmap)
        plot_reference_flow_on_grid(self, plot_type=plot_type,ax=ax_, scale=scale, show_background=False, s=s)
        ax_.set_title("Gradient of pseudotime \n(Development flow)")
        ax_ = ax[3]
        plot_reference_flow_on_grid(self, plot_type=plot_type,ax=ax_, scale=scale, show_background=show_background, s=s)
        ax_.set_title("Gradient of pseudotime \n(=Development flow)")
        ax_ = ax[4]
        plot_pseudotime(self, ax=ax_, s=s, show_background=show_background,cmap=cmap)
        plot_reference_flow_on_grid(self, plot_type=plot_type,ax=ax_, scale=scale, show_background=False, s=s,)
        ax_.set_title("Pseudotime + \nDevelopment flow")
    def plot_pseudotime(self, ax=None, s=CONFIG["s_scatter"],show_background=True, args=CONFIG["default_args"], cmap="rainbow"):
        plot_pseudotime(self, ax=None, s=s, show_background=show_background, cmap=cmap, args=args)







