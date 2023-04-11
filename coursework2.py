import osmnx
import networkx
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, LineString
import matplotlib.pyplot as plt
import pyproj
import warnings
import numpy as np
import random
from networkx.algorithms import planarity

# Task A
center_location = (53.7987, -1.5483) # Here I use exact coordinates of leeds center area
search_distance = 1200


road_network_graph = osmnx.graph_from_point(
    center_location, dist=search_distance, network_type="drive" # Here I set the network type to drive to filer the road accidents
)


csv_files = [
    "Accidents_2016.csv",
    "Accidents_2017.csv",
    "Accidents_2018.csv",
    "Accidents_2019.csv",
]
combined_csv_files = pd.concat(
    [
        pd.read_csv(
            file,
            usecols=["Grid Ref: Easting", "Grid Ref: Northing"],
            encoding="Windows-1252",
        )
        for file in csv_files
    ],
    ignore_index=True,
)


easting_northing_projection = pyproj.Proj(init="epsg:27700")
latitude_longitude_projection = pyproj.Proj(init="epsg:4326")
combined_csv_files["longitude"], combined_csv_files["latitude"] = pyproj.transform(
    easting_northing_projection,
    latitude_longitude_projection,
    combined_csv_files["Grid Ref: Easting"],
    combined_csv_files["Grid Ref: Northing"],
)

total_distances = np.array(
    osmnx.distance.nearest_nodes(
        road_network_graph,
        combined_csv_files["longitude"],
        combined_csv_files["latitude"],
        return_dist=True,
    )[1]
)
all_accidents_located = combined_csv_files[total_distances < search_distance]

print(
    f"There are {len(all_accidents_located)} accidents within a {search_distance}-meter radius of the center location."
)


# 1. Road network characteristics. Some characteristics require calculation, others can be retrieved through osmnx.basic_stats() method

# The average street length of the road network
graph_statistics = osmnx.basic_stats(road_network_graph)
street_length_avg = graph_statistics["street_length_avg"]
print("Average streeth length for this network is {}".format(street_length_avg))

# The intersection density of the road network
intersec_density = graph_statistics['streets_per_node_avg']
print("Intersection density for this network is {}".format(intersec_density))

# The diameter of the road network using the networkx module
scc = max(networkx.strongly_connected_components(road_network_graph), key=len)
diam = networkx.diameter(road_network_graph.subgraph(scc))
print("Diameter for this network is {}:".format(diam))

# 2. Finding the average circuity
road_network_graph_statistics = osmnx.basic_stats(road_network_graph)
avg_circ = road_network_graph_statistics["circuity_avg"]
print("Average circuity for this network is: {}".format(avg_circ))

# 3. Check if my network is planar or not
if planarity.is_planar(road_network_graph):
    print("It is planar road network.")
else:
    print("It is non-planar road network.")

# Task B

# 1. Plot distribution of road accidents and visualize it
accident_points = gpd.GeoDataFrame(
    combined_csv_files,
    geometry=gpd.points_from_xy(
        combined_csv_files.longitude, combined_csv_files.latitude
    ),
    crs="EPSG:4326",
)

accident_points = accident_points.to_crs(road_network_graph.graph["crs"])

figure, axes = osmnx.plot_graph(road_network_graph, show=False, close=False)
accident_points.plot(ax=axes, color="orange", markersize=2)
plt.show()

#2. Investigate whether a high number of accidents on one road correlates with a high number on connecting roads

degree_centrality = networkx.degree_centrality(road_network_graph)

fig, ax = osmnx.plot_graph(
    road_network_graph,
    node_color=list(degree_centrality.values()),
    node_size=30,
    edge_linewidth=0.2,
    edge_color="#333333",
    bgcolor="white",
    show=False,
    close=False,
)

# add colorbar for degree centrality
cbar = plt.colorbar(ax.collections[0], shrink=0.4)
cbar.set_label("Degree Centrality")


# 3

nearest_intersections = osmnx.distance.nearest_nodes(
    road_network_graph,
    all_accidents_located["longitude"],
    all_accidents_located["latitude"],
    return_dist=False,
)
nearest_edges = osmnx.distance.nearest_edges(
    road_network_graph,
    all_accidents_located["longitude"],
    all_accidents_located["latitude"],
    return_dist=True,
)
road_lengths = [data[0]["length"] for _, _, data in road_network_graph.edges(data=True)]
fractions_to_intersections = nearest_edges[2] / road_lengths


# Task C


