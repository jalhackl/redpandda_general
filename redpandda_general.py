import math
import CONSTANTS
import mdtraj as md
import distance_matrix as dm
import numpy as np
import multiprocessing
import warnings
import compare_clusterings as cc

from functools import wraps
from time import time

def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        ts = time()
        result = f(*args, **kw)
        te = time()
        #print('func:%r took: %2.4f sec' % (f.__name__, te-ts))
        return result
    return wrap

warnings.filterwarnings("ignore")
pool = multiprocessing.Pool()




@timing
def get_distance_matrices(traj_array):
    return np.array(list(pool.map(dm.calculate_distance_matrix, traj_array)))

@timing
def get_delta_matrices(dist_matrices):
    return np.absolute(np.diff(dist_matrices,axis=0))

@timing
def calculate_average_delta_matrix(delta_matrices):
    return np.sum(delta_matrices, axis=0)/len(delta_matrices)



@timing
def calculate_median_delta_matrix(delta_matrices):
    return np.median(delta_matrices, axis=0)


#new
@timing
def get_var_matrices(dist_matrices):
    return np.absolute(np.var(dist_matrices,axis=0))


@timing
def get_std_matrices(dist_matrices):
    return np.absolute(np.std(dist_matrices,axis=0))

#new, not used
@timing
def calculate_var_delta_matrix(delta_matrices):
    return np.var(delta_matrices, axis=0)/len(delta_matrices)


def prepare_data_from_df(x):
  import pandas as pd 
  tpoints = []
  df_points = []
  traj_array = []
  point_array = []
  #for g in x.sort_values(['t'],ascending=True).groupby("obj_id"):
  for g in x.sort_values(['obj_id'],ascending=True).groupby("t"):
      tpoints.append(g[1].values)
      df_points.append(pd.DataFrame(g[1]))
      new_df = pd.DataFrame(g[1])
      traj_array.append(np.array(new_df[['x', 'y', 'z']].values))
      point_array.append(new_df[['obj_id']].values)

  frames_count = len(df_points[0])
  n_objects = len(df_points)
  #traj_array, k_cluster = preprocessing(prot_info,frames_count,k_cluster)
   
  return traj_array, point_array, frames_count, n_objects


def main_clustering_general(x,k_cluster=None,assign_labels='discretize', clustering_algorithm="spectral", no_q_computation=False, eigengap = False, silhouette=False, min_cluster_size=2, min_samples=1, plot_delta_heatmaps=True) -> np.ndarray:
  
  import pandas as pd 
  tpoints = []
  df_points = []
  traj_array = []
  point_array = []
  #for g in x.sort_values(['t'],ascending=True).groupby("obj_id"):
  for g in x.sort_values(['obj_id'],ascending=True).groupby("t"):
      tpoints.append(g[1].values)
      df_points.append(pd.DataFrame(g[1]))
      new_df = pd.DataFrame(g[1])
      traj_array.append(np.array(new_df[['x', 'y', 'z']].values))
      point_array.append(new_df[['obj_id']].values)

  frames_count = len(df_points[0])
  n_objects = len(df_points)
  #traj_array, k_cluster = preprocessing(prot_info,frames_count,k_cluster)


  if not k_cluster:
    k_cluster = int(math.sqrt(n_objects))
  dist_matrices = get_distance_matrices(traj_array)
  delta_matrices = get_delta_matrices(dist_matrices)
  average_delta_matrix = calculate_average_delta_matrix(delta_matrices)

  matrix = []

  if clustering_algorithm == "spectral":
    clustering = dm.spectral_clustering_on_deltas(average_delta_matrix, k_cluster, assign_labels, eigengap = eigengap, silhouette = silhouette, plot_delta_heatmaps=plot_delta_heatmaps)
  elif clustering_algorithm == "hdbscan":
     clustering, k_cluster, matrix = dm.hdbscan_clustering_on_deltas(average_delta_matrix, k_cluster, assign_labels, min_cluster_size=min_cluster_size, min_samples=min_samples)
  elif clustering_algorithm == "agglomerative":
     clustering, k_cluster, matrix = dm.agglomerative_clustering_on_deltas_ward(average_delta_matrix, None)


  Q = []
  

  if not no_q_computation: 
    try:
      Q, _ = cc.get_Q_for_clustering(dist_matrices, clustering, k_cluster)
    except:
      print("something wrong with clustering Q computation")


  return clustering, Q, matrix, point_array




def main_clustering_general_std(x,k_cluster=None,assign_labels='discretize', clustering_algorithm="spectral", no_q_computation=False, eigengap = False, fuser_alg = "agreement", silhouette=False, min_cluster_size=2, min_samples=1, plot_delta_heatmaps=True) -> np.ndarray:
  #lets assume dataframe with t and x.y.z
  import pandas as pd 
  tpoints = []
  df_points = []
  traj_array = []
  point_array = []
  #for g in x.sort_values(['t'],ascending=True).groupby("obj_id"):
  for g in x.sort_values(['obj_id'],ascending=True).groupby("t"):

      tpoints.append(g[1].values)
      df_points.append(pd.DataFrame(g[1]))
      new_df = pd.DataFrame(g[1])
      traj_array.append(np.array(new_df[['x', 'y', 'z']].values))
      point_array.append(new_df[['obj_id']].values)

  frames_count = len(df_points[0])
  n_objects = len(df_points)
  #traj_array, k_cluster = preprocessing(prot_info,frames_count,k_cluster)


  if not k_cluster:
    k_cluster = int(math.sqrt(n_objects))
  dist_matrices = get_distance_matrices(traj_array)
  delta_matrices = get_delta_matrices(dist_matrices)
  average_delta_matrix = calculate_average_delta_matrix(delta_matrices)

  std_delta_matrix = get_std_matrices(dist_matrices)


  matrix = []

  '''
  if clustering_algorithm == "spectral":
    clustering = dm.spectral_clustering_on_deltas(average_delta_matrix, k_cluster, assign_labels, eigengap = eigengap, silhouette = silhouette, plot_delta_heatmaps=plot_delta_heatmaps)
  elif clustering_algorithm == "hdbscan":
     clustering, k_cluster, matrix = dm.hdbscan_clustering_on_deltas(average_delta_matrix, k_cluster, assign_labels, min_cluster_size=min_cluster_size, min_samples=min_samples)
  elif clustering_algorithm == "agglomerative":
     clustering, k_cluster, matrix = dm.agglomerative_clustering_on_deltas_ward(average_delta_matrix, None)
  '''

 
  clustering, k_cluster, matrix = dm.hdbscan_clustering_on_deltas_varadd(average_delta_matrix, std_delta_matrix, k_cluster, assign_labels, min_cluster_size=min_cluster_size, min_samples=min_samples)


  Q = []
  

  if not no_q_computation: 
    try:
      Q, _ = cc.get_Q_for_clustering(dist_matrices, clustering, k_cluster)
    except:
      print("something wrong with clustering Q computation")


  return clustering, Q, matrix, point_array









def main_clustering_general_std_distance(x,k_cluster=None,assign_labels='discretize', clustering_algorithm="spectral", no_q_computation=False, eigengap = False, fuser_alg = "agreement", silhouette=False, min_cluster_size=2, min_samples=1, plot_delta_heatmaps=True) -> np.ndarray:
  import pandas as pd 
  tpoints = []
  df_points = []
  traj_array = []
  point_array = []
  #for g in x.sort_values(['t'],ascending=True).groupby("obj_id"):
  for g in x.sort_values(['obj_id'],ascending=True).groupby("t"):

      tpoints.append(g[1].values)
      df_points.append(pd.DataFrame(g[1]))
      new_df = pd.DataFrame(g[1])
      traj_array.append(np.array(new_df[['x', 'y', 'z']].values))
      point_array.append(new_df[['obj_id']].values)

  frames_count = len(df_points[0])
  n_objects = len(df_points)
  #traj_array, k_cluster = preprocessing(prot_info,frames_count,k_cluster)


  if not k_cluster:
    k_cluster = int(math.sqrt(n_objects))
  dist_matrices = get_distance_matrices(traj_array)
  delta_matrices = get_delta_matrices(dist_matrices)
  average_delta_matrix = calculate_average_delta_matrix(delta_matrices)

  std_delta_matrix = get_std_matrices(dist_matrices)


  matrix = []

  '''
  if clustering_algorithm == "spectral":
    clustering = dm.spectral_clustering_on_deltas(average_delta_matrix, k_cluster, assign_labels, eigengap = eigengap, silhouette = silhouette, plot_delta_heatmaps=plot_delta_heatmaps)
  elif clustering_algorithm == "hdbscan":
     clustering, k_cluster, matrix = dm.hdbscan_clustering_on_deltas(average_delta_matrix, k_cluster, assign_labels, min_cluster_size=min_cluster_size, min_samples=min_samples)
  elif clustering_algorithm == "agglomerative":
     clustering, k_cluster, matrix = dm.agglomerative_clustering_on_deltas_ward(average_delta_matrix, None)
  '''

  clustering, k_cluster, matrix = dm.hdbscan_clustering_on_deltas_varadd_distance(average_delta_matrix, std_delta_matrix, k_cluster, assign_labels, min_cluster_size=min_cluster_size, min_samples=min_samples)


  Q = []
  

  if not no_q_computation: 
    try:
      Q, _ = cc.get_Q_for_clustering(dist_matrices, clustering, k_cluster)
    except:
      print("something wrong with clustering Q computation")


  return clustering, Q, matrix, point_array



