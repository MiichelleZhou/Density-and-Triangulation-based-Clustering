class triangulation_dbscan:

  '''
      @param: data 
      @param: minPts
      @param: kde
      @param: figsize
      @param: progress
  '''

  def __init__(self, data, minPts = 5, kde = False, local_std = 2.5, figsize = (8,8),  progress = False):
      self.data = data
      self.kde = kde
      self.minPts = minPts
      self.local_std = local_std
      self.figsize = figsize
      self.progress = progress


  def find_neighbors(self, pindex, triang):
    neighbors = list()
    for simplex in triang:
        if pindex in simplex:
            neighbors.extend([simplex[i] for i in range(len(simplex)) if simplex[i] != pindex])
    return list(set(neighbors)), len(set(neighbors))


  def find_closest_point_index(self, df, idx1, idxmax_neighbors, chosen):
    point1 = df.loc[idx1, ['x', 'y']]
    min_distance = float('inf')
    closest_index = None
    for idx2 in idxmax_neighbors:
      if idx2 not in chosen:
          point2 = df.loc[idx2, ['x', 'y']]
          distance = np.sqrt(np.sum((point1 - point2)**2))
          if distance < min_distance:
              min_distance = distance
              closest_index = idx2
    return closest_index, min_distance


  def glob_loc_removal(self, data):
    warnings.filterwarnings("ignore")

    X = data.to_numpy()
    points = X
    tri = Delaunay(points)

    ### Global Removal
    side_len = np.array([pdist(x) for x in points[tri.simplices]])
    tri_area = []
    for sides in side_len:
      a = sides[0]
      b = sides[1]
      c = sides[2]
      s = (a + b + c) / 2
      area = (s*(s-a)*(s-b)*(s-c)) ** 0.5
      tri_area.append(area)
    tri_area = pd.DataFrame(tri_area, columns = ['area'])
    tri_area[tri_area - tri_area.mean()  >  tri_area.std()] = np.nan
    deleted = []
    for i in range(len(tri_area)):
      if np.isnan(tri_area.loc[i]['area']) == True:
        deleted.append(i)
    deleted_df = pd.DataFrame(deleted, columns = ['deleted_area'])
    all_len = side_len.flatten()
    all_len = pd.DataFrame (all_len, columns = ['length'])
    all_len[all_len-all_len.mean() > all_len.std()] = np.nan
    for i in range(len(all_len)):
      if np.isnan(all_len.loc[i]['length']) == True:
        deleted.append(i//3)
    new_tri = np.delete(tri.simplices, deleted,0)

    if self.progress == True:
      fig, axes = plt.subplots(1, 2, figsize=(10,5) )
      axes[0].triplot(points[:,0], points[:,1], new_tri)
      axes[0].plot(points[:,0], points[:,1], 'o')
      axes[0].set_title('Triangulation after global removal')

    ### Local Removal
    side_len = np.array([pdist(x) for x in points[new_tri]])
    neigh_V = []
    for i in range(data.shape[0]):
      neigh_V.append(self.find_neighbors(i,new_tri)[0])

    i = 0
    alldist_df = pd.DataFrame(columns = ['p','s','dist'])
    for V in neigh_V:
      s = data.iloc[V].values
      p = data.iloc[i].values
      dists = np.array([math.sqrt((p[0]-s0)**2 + (p[1]-s1)**2) for s0, s1 in s])
      dists[dists-dists.mean() > self.local_std *dists.std()] = np.nan
      dist_df = pd.DataFrame(columns = ['p','s','dist'])

      dist_df['s'] = V
      dist_df['dist'] = dists
      dist_df['p'] = [i]*len(V)

      alldist_df = pd.concat([alldist_df, dist_df],ignore_index=True)
      i += 1

    need_del = alldist_df[np.isnan(alldist_df['dist'])].iloc[:,0:2]


    delete = []
    for i in range(need_del.shape[0]):
      row = need_del.iloc[i].tolist()
      for j in range(new_tri.shape[0]):
        if all(x in new_tri[j].tolist() for x in row) == True:
          delete.append(new_tri[j].tolist())

    a = new_tri.tolist()
    b = delete
    for i in a[:]:
      if i in b:
        a.remove(i)
        b.remove(i)
    new_tri = np.array(a)

    if self.progress == True:
      axes[1].triplot(points[:,0], points[:,1], new_tri)
      axes[1].plot(points[:,0], points[:,1], 'o')
      axes[1].set_title('Triangulation after local removal')
      plt.show()

    return new_tri

  def check_point(self, minPts, index, triang):
  # check number of connected nodes
    neighbors, num_neighbors = self.find_neighbors(index, triang)
    if num_neighbors >= minPts:
      return (True, neighbors)
    else:
      return (False, neighbors)

  def is_unique(self, s):
    a = s.to_numpy() # s.values (pandas<0.24)
    return (a[0] == a).all()


  def kde_clustering(self, data):

    print("Start creating initial clusters based on Kernel Density estimation:")

    finished = False
    initial_clusters = []

    while finished == False:
      X = data[['x','y']].to_numpy()
      points = X
      tri = Delaunay(points)

      ls_included_points = []
      column_names = ['idx1','idx2','step_size']
      dtypes = {'idx1': 'int', 'idx2': 'int', 'step_size': 'float'}
      df_steps = pd.DataFrame(columns=column_names).astype(dtypes)

      # find the begining point
      idxmax = data['density'].idxmax()
      ls_included_points.append(idxmax)

      # find its neighbors
      neighbors = self.find_neighbors(idxmax,tri.simplices)[0]
      # find the closest neighbor & record the step size
      chosen_point, step = self.find_closest_point_index(data, idxmax, neighbors, ls_included_points)
      ls_included_points.append(chosen_point)
      df_steps.loc[len(df_steps)] = [idxmax,chosen_point,step]

      for _ in range(12):
        possible_min_steps = []
        for point in ls_included_points:
          neighbors = self.find_neighbors(point,tri.simplices)[0]
          chosen_point, step = self.find_closest_point_index(data, point, neighbors,ls_included_points)
          possible_min_steps.append([point,chosen_point,step])
        next = min(possible_min_steps, key=lambda step: step[2])
        ls_included_points.append(next[1])
        df_steps.loc[len(df_steps)] = next

      for _ in range(data.shape[0]-10):
        possible_min_steps = []
        for point in ls_included_points:
          neighbors = self.find_neighbors(point,tri.simplices)[0]
          chosen_point, step = self.find_closest_point_index(data,point,neighbors,ls_included_points)
          possible_min_steps.append([point,chosen_point,step])
        next = min(possible_min_steps, key=lambda step: step[2])
        if next[2] > (np.mean(df_steps['step_size']) + 4*np.std(df_steps['step_size'])):
          break
        ls_included_points.append(next[1])
        df_steps.loc[len(df_steps)] = next

      print(ls_included_points)
      indexed_points = data.loc[ls_included_points, 'index']
      initial_clusters.append(indexed_points)
      data = data[~data.index.isin(ls_included_points)]
      data.reset_index(inplace=True, drop=True)
      if data.shape[0] == 0:
        finished = True

    return initial_clusters



  def tri_dbscan(self):
    if self.kde == True:
      df = self.data
      df.reset_index(inplace=True, drop=False)
      tri_est = gaussian_kde(df.T).evaluate(df.T)
      tri_est = tri_est.reshape(-1, 1)
      pd.DataFrame(tri_est).to_csv('density.csv', index=False, header=False)

      density = pd.read_csv("density.csv",names=['density'],header=None)
      df.loc[:, 'density'] = density

      initial_clusters = self.kde_clustering(df)
      colors = ['green', 'red', 'blue', 'yellow', 'black', 'purple', 'cyan', 'magenta', 'orange', 'brown']

      print("Start clustering within each initial cluster based on Triangulation")
      final_df = pd.DataFrame(columns=['x', 'y', 'index', 'est_clust'])
      for idx, cluster in enumerate(initial_clusters):
        clustered_points = self.data.iloc[cluster]
        clustered_points.reset_index(inplace=True, drop=True)

        tri = self.glob_loc_removal(clustered_points[['x','y']])
        df = clustered_points
        minPts = self.minPts

        #initiating cluster number
        C = 1
        #initiating stacks to maintain
        current_stack = set()
        unvisited = list(df.index)
        clusters = []

        # run until all points have been visited
        while (len(unvisited) != 0):
            # choose a random unvisited point
            current_stack.add(random.choice(unvisited))
            while len(current_stack) != 0:
                curr_idx = current_stack.pop()
                # check if point is in a cluster or is a noise
                check_result, neigh_indexes = self.check_point(minPts, curr_idx, tri)
                neigh_indexes = set(neigh_indexes) & set(unvisited)

                if check_result == True:
                    clusters.append((curr_idx, C))
                    unvisited.remove(curr_idx)
                    current_stack.update(neigh_indexes)
                    continue

                else: # current point is noise
                    clusters.append((curr_idx, 0))
                    unvisited.remove(curr_idx)
                    continue
            # increment cluster number
            C += 1

        clusters_df = pd.DataFrame(clusters, columns =['index', 'est_clust'])
        joined_df = df.join(clusters_df.set_index('index'), how='left', on=df.index)
        clusters_df = joined_df[['x','y','index','est_clust']]

        ### KNN
        X_train = np.array(clusters_df[['x','y']][clusters_df['est_clust']!=0])
        y_train = np.array(clusters_df['est_clust'][clusters_df['est_clust']!=0])
        np.array(clusters_df[['x','y']][clusters_df['est_clust']==0])
        predict_x = np.array(clusters_df[['x','y']][clusters_df['est_clust']==0])

        if len(predict_x) != 0:
          knn = KNeighborsClassifier(n_neighbors=1)
          knn.fit(X_train, y_train)
          knn_result = knn.predict(predict_x)
          clusters_df.loc[clusters_df['est_clust'] == 0, 'est_clust'] = knn_result
          clusters_df = clusters_df[['x','y', 'index','est_clust']]


        # Check for unique labels in clusters_df
        unique_clusters = clusters_df['est_clust'].unique()
        # Get the labels that are already present in final_df
        existing_clusters = final_df['est_clust'].unique()
        # If final_df is empty, set max_existing_cluster to 0, otherwise get the max value from existing_clusters
        if final_df.empty:
            max_existing_cluster = 0
        else:
            max_existing_cluster = max(existing_clusters)

        # Create a mapping to store old cluster labels and their corresponding new labels
        cluster_mapping = {}

        # Go through each unique label in clusters_df
        for cluster in unique_clusters:
            # If it's already in final_df, give it a new label
            if cluster in existing_clusters:
                max_existing_cluster += 1
                cluster_mapping[cluster] = max_existing_cluster
            else:
                cluster_mapping[cluster] = cluster

        # Update the labels in clusters_df based on the mapping
        clusters_df['est_clust'] = clusters_df['est_clust'].map(cluster_mapping)
        final_df = final_df.append(clusters_df)

      final_df.reset_index(inplace=True, drop=True)
      clus_num = final_df['est_clust'].unique()
      colors = distinctipy.get_colors(len(clus_num))
      color_map = {clus: color for clus, color in zip(clus_num, colors)}

      for i in range(len(final_df.est_clust)):
        idx = final_df['index'][i]
        cluster = final_df.est_clust[i]
        color = color_map[cluster]
        plt.plot(self.data.iloc[idx]["x"], self.data.iloc[idx]["y"], 'o', c=color)
      return final_df


    if self.kde == False:
      tri = self.glob_loc_removal(self.data)
      df = self.data
      minPts = self.minPts
      #initiating cluster number
      C = 1
      #initiating stacks to maintain
      current_stack = set()
      unvisited = list(df.index)
      clusters = []

      # run until all points have been visited
      while (len(unvisited) != 0):
          # choose a random unvisited point
          current_stack.add(random.choice(unvisited))
          while len(current_stack) != 0:
              curr_idx = current_stack.pop()
              # check if point is in a cluster or is a noise
              check_result, neigh_indexes = self.check_point(minPts, curr_idx, tri)
              neigh_indexes = set(neigh_indexes) & set(unvisited)

              if check_result == True:
                  clusters.append((curr_idx, C))
                  unvisited.remove(curr_idx)
                  current_stack.update(neigh_indexes)
                  continue

              else: # current point is noise
                  clusters.append((curr_idx, 0))
                  unvisited.remove(curr_idx)
                  continue
          # increment cluster number
          C += 1
      clusters_df = pd.DataFrame(clusters, columns =['index', 'est_clust'])
      self.data.reset_index(inplace=True)
      joined_df = pd.merge(self.data,
                    clusters_df,
                    on ='index',
                    how ='left')
      clusters_df = joined_df[['x','y','index','est_clust']]

      ### KNN
      X_train = np.array(joined_df[['x','y']][joined_df['est_clust']!=0])
      y_train = np.array(joined_df['est_clust'][joined_df['est_clust']!=0])
      np.array(joined_df[['x','y']][joined_df['est_clust']==0])
      predict_x = np.array(joined_df[['x','y']][joined_df['est_clust']==0])

      if len(predict_x) != 0:
        knn = KNeighborsClassifier(n_neighbors=1)
        knn.fit(X_train, y_train)
        knn_result = knn.predict(predict_x)
        joined_df.loc[joined_df['est_clust'] == 0, 'est_clust'] = knn_result
        clusters_df = joined_df[['x','y', 'index','est_clust']]

      plt.figure(figsize=self.figsize )
      X = self.data.to_numpy()
      points = X
      clus_num = clusters_df['est_clust'].unique()
      colors = distinctipy.get_colors(len(clus_num))
      color_map = {clus: color for clus, color in zip(clus_num, colors)}

      for i in range(len(clusters_df.est_clust)):
          idx = clusters_df['index'][i]
          cluster = clusters_df.est_clust[i]
          # Use the mapping to get the color for this cluster
          color = color_map[cluster]
          plt.plot(self.data.iloc[idx]["x"], self.data.iloc[idx]["y"], 'o', c=color)

      clusters_df['est_clust'] = pd.factorize(clusters_df['est_clust'])[0] + 1
      return clusters_df
