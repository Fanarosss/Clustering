# Clustering with K-Means
***
Implemented K-Means and K-Medoids for vectors and curves with C++

### [Initializers:](https://github.com/Fanarosss/Clustering/blob/master/Clustering/src/Initializers/Initializers.cpp)
* **Random Selection:**
With random selection we select uniformly at random k points of the dataset to start as centroids
* **K-Means++:**
With K-Means++ we select k points of the dataset, spread very far away from the data

### [Assigners:](https://github.com/Fanarosss/Clustering/blob/master/Clustering/src/Assigners/Assigners.cpp)
* **Loyd's:**
For every point compute distances to all centroids, and select the one that is closer. It seems as a very good practice, but in reality it converges really slow, because the decisions affect the clusters very little, so there are not big differences in centroid update.
* **Inverse Assignment:**
We use [LSH nearest neighbor](https://github.com/Fanarosss/Clustering/tree/master/Clustering/src/LSH) for this assigner. Apply Range Search with LSH for the R-nearest neighbors, to get for every centroid the R-closest points of the dataset and assign to those points as centroid id. For points who came as nearest neighbors of more the one centroid, we calculate the distance, and select the one that is closer. That way, the algorithm converges faster because some changes affect the clustering score in a certain amount, where the centroids are changing rapidly.

### [Updaters:](https://github.com/Fanarosss/Clustering/blob/master/Clustering/src/Updaters/Updaters.cpp)
* **PAM - K Medoids:**
K-Medoids is a very good algorithm, insensitive to outliers providing very good results. But its complexity is very high **O(n^2 * K * i)** since it has to calculate for every combination of points their distances. In order to make this algorithm faster, we implemented a database, where the distances are pre calculated exhaustively.
* **K-Means:**
K-Means is sensitive to outliers, and worse than K-Medoids for small datasets, where the time of completion is not a factor. K-Means works very well for very big datasets, where K-Medoids is extremely slow. It gives very good results in **O(n * k * i)** where *k * i << n*

### Metric for cluster score:
We used [Silhouette](https://en.wikipedia.org/wiki/Silhouette_(clustering)) for clustering evaluation and validation of consistency.<br>
- Vectors ~0.98<br>
- Curves ~0.4<br>

### Usage:
- Compilation<br>
*make cluster*
- Run<br>
*./build/cluster.x –i <input_file> -c <configuration_file> -o <output_file>*

### Collaborators
[Konstantinos Athinaios](https://github.com/KostasA97)
<br>
[Theofanis Aslanidis](https://github.com/Fanarosss)
