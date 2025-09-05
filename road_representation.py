# -*- coding: utf-8 -*-
import os,tempfile,shutil
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import Point,LineString
import geopandas as gpd
import arcpy
import networkx as nx

import coins
import warnings
warnings.filterwarnings("ignore")

def road_divide(OSM_RN,buffer_distance=50,threshold=0.8,angle=30):
    OSM_RN = OSM_RN.reset_index(drop=True)
    OSM_RN['OID'] = OSM_RN.index
    #Select oneway roads
    oneway = OSM_RN[OSM_RN['oneway'].isin([1, "True"])]
    #certain proportion of length falls within others buffer
    oneway_buffered = oneway.copy ()
    oneway_buffered ['geometry'] = oneway_buffered.buffer(distance=buffer_distance,cap_style='flat')
    intersection_gdf = gpd.overlay(oneway, oneway_buffered, how='intersection')
    intersection_gdf=intersection_gdf[intersection_gdf['OID_1']!=intersection_gdf['OID_2']]
    length_map = oneway.set_index('OID')['geometry'].length.to_dict()
    intersection_gdf['total_length'] = intersection_gdf['OID_1'].map(length_map)
    intersection_gdf['intersect_length'] = intersection_gdf.geometry.length
    intersection_gdf['intersection_ratio'] = (intersection_gdf['intersect_length'] / intersection_gdf['total_length']).round(4)
    intersection_gdf=intersection_gdf[intersection_gdf['intersection_ratio'] > threshold]
    #similar azimuths
    def calculate_azimuth(line):
        x1, y1 = line.coords[0]
        x2, y2 = line.coords[-1]
        dx, dy = x2 - x1, y2 - y1
        azimuth_deg = np.degrees(np.arctan2(dy, dx))
        if azimuth_deg < 0:
            azimuth_deg += 180
        return azimuth_deg
    def valid_angle(row):
        road1 = oneway.loc[row['OID_1'], 'geometry']
        road2 = oneway.loc[row['OID_2'], 'geometry']
        az1, az2 = calculate_azimuth(road1), calculate_azimuth(road2)
        az_diff = abs(az1 - az2)
        az_diff = min(az_diff, 180 - az_diff)
        return az_diff <= angle
    intersection_gdf = intersection_gdf[intersection_gdf.apply(valid_angle, axis=1)]
    
    index_list=pd.concat([intersection_gdf['OID_1'],intersection_gdf['OID_2']], ignore_index=True).drop_duplicates()
    road_aggregated = OSM_RN[OSM_RN['OID'].isin(index_list)]
    road_not_aggregated = OSM_RN[~OSM_RN['OID'].isin(index_list)]
    #roads link two candidates
    sindex = road_aggregated.sindex  
    def endpoint_connected(point):
        possible_matches = list(sindex.intersection(point.bounds))
        if not possible_matches:
            return False
        return any(road_aggregated.iloc[i].geometry.intersects(point) for i in possible_matches)

    selected = []
    for idx, geom in road_not_aggregated.geometry.items():
        start, end = Point(geom.coords[0]), Point(geom.coords[-1])
        if endpoint_connected(start) and endpoint_connected(end):
            selected.append(idx)    
    to_add = road_not_aggregated.loc[selected]
    road_aggregated = pd.concat([road_aggregated, to_add], ignore_index=True)
    road_not_aggregated = road_not_aggregated.drop(selected)
    return road_aggregated,road_not_aggregated

def get_centerline(road,buffer_distance=50):
    tmpdir = tempfile.mkdtemp()
    road.to_file(os.path.join(tmpdir, "road.shp"))

    gdb_path = os.path.join(tmpdir, "temp.gdb")
    arcpy.CreateFileGDB_management(tmpdir, "temp.gdb")    
    arcpy.analysis.PairwiseBuffer(os.path.join(tmpdir, "road.shp"),
                                   os.path.join(gdb_path, "road_buffer"),
                                   f"{buffer_distance} Meters", 
                                   "ALL")
    arcpy.EliminatePolygonPart_management(os.path.join(gdb_path, "road_buffer"),
                                          os.path.join(gdb_path, "clean"),
                                          condition="AREA",
                                          part_area="2000 SquareMeters")
    arcpy.topographic.PolygonToCenterline(os.path.join(gdb_path, "clean"), 
                                          os.path.join(gdb_path, "centerline"), 
                                          None)
    arcpy.cartography.SimplifyLine(os.path.join(gdb_path, "centerline"), 
                                    os.path.join(gdb_path, "SimplifyLine"), 
                                    "BEND_SIMPLIFY",
                                    "50 Meters",
                                    "RESOLVE_ERRORS", 
                                    "NO_KEEP",
                                    "CHECK", 
                                    None)
    centerline = gpd.read_file(gdb_path,layer='SimplifyLine')
    centerline=centerline.explode(index_parts=False)
    shutil.rmtree(tmpdir)
    return centerline

def repair_endpoint(centerline,road_aggregated,road_not_aggregated):  
    # Get intersection points between aggregated and non-aggregated roads
    union_aggregated = road_aggregated.geometry.unary_union
    union_not_aggregated = road_not_aggregated.geometry.unary_union
    intersections = union_aggregated.intersection(union_not_aggregated)
    intersections_gdf = gpd.GeoDataFrame(geometry=list(intersections.geoms),crs=centerline.crs)
    intersections_gdf = intersections_gdf.drop_duplicates(subset='geometry')
    # Detect hanging nodes
    all_ends =centerline.geometry.apply(lambda g: [g.coords[0], g.coords[-1]]).explode()
    ends_df = pd.DataFrame({
        'xy': [f"{x:.6f},{y:.6f}" for x,y in all_ends],
        'x': [x for x,y in all_ends],
        'y': [y for x,y in all_ends]
    })
    dangling = ends_df.groupby('xy').filter(lambda g: len(g) == 1)
    hanging_nodes=gpd.GeoDataFrame(geometry=[Point(x,y) for x,y in zip(dangling['x'], dangling['y'])],crs=centerline.crs)
    # Process dangling nodes
    intersections_sindex = intersections_gdf.sindex
    for idx, hanging_node in hanging_nodes.iterrows():
        node_geom = hanging_node.geometry
        intersecting_centerlines = centerline[centerline.geometry.intersects(node_geom)]
        line = intersecting_centerlines.geometry.iloc[0]
        nearby_points_idx = list(intersections_sindex.intersection(node_geom.buffer(55).bounds))
        nearby_points = intersections_gdf.iloc[nearby_points_idx]

        if nearby_points.empty:
            if line.length < 50:
                centerline.drop(intersecting_centerlines.index)
            else:
                if Point(line.coords[0]) == node_geom:
                    new_line = shapely.ops.substring(line,50, line.length)
                else:
                    new_line = shapely.ops.substring(line,0, line.length-50)
                centerline.loc[intersecting_centerlines.index, 'geometry'] = new_line
            continue

        nearest_intersection = nearby_points.loc[nearby_points.geometry.distance(node_geom).idxmin()]
        intersection_geom = nearest_intersection.geometry
        if line.length < 100:
            endpoints = [Point(coord) for coord in [line.coords[0], line.coords[-1]]]
            point_linked=[p for p in endpoints if not p.equals(node_geom)][0]
            new_line = LineString([point_linked.coords[0], intersection_geom.coords[0]])
        else:
            if Point(line.coords[0]) == node_geom:
                sub_line = shapely.ops.substring(line,100, line.length)
                new_line = LineString([intersection_geom.coords[0]]+list(sub_line.coords))
            else:
                sub_line = shapely.ops.substring(line,0, line.length-100)
                new_line =LineString(list(sub_line.coords)+[intersection_geom.coords[0]])
        centerline.loc[intersecting_centerlines.index, 'geometry'] = new_line
    return centerline

def extend_road(centerline, road_not_aggregated, distance=100, tol=5):
    tmpdir = tempfile.mkdtemp()
    merged = gpd.GeoDataFrame(
        pd.concat([centerline, road_not_aggregated], ignore_index=True),
        crs=centerline.crs
    )
    merged.to_file(os.path.join(tmpdir, "roads.shp"))

    gdb_path = os.path.join(tmpdir, "temp.gdb")
    arcpy.CreateFileGDB_management(tmpdir, "temp.gdb")    
    arcpy.edit.ExtendLine(os.path.join(gdb_path, "feat_line"),
                          f"{distance} Meters",
                          "EXTENSION")
    with arcpy.EnvManager(XYTolerance=f"{tol} Meters"):
        arcpy.management.Integrate(os.path.join(gdb_path, "feat_line"),
                                   None)
    arcpy.management.Dissolve(os.path.join(gdb_path, "feat_line"),
                              os.path.join(gdb_path, "Dissolve"), 
                              None,
                              None, 
                              "MULTI_PART",
                              "DISSOLVE_LINES",
                              '')
    arcpy.management.FeatureToLine(os.path.join(gdb_path, "Dissolve"),
                                   os.path.join(gdb_path, "result"), 
                                   None, 
                                   "NO_ATTRIBUTES")
    centerline = gpd.read_file(gdb_path, layer="result")
    centerline=centerline.explode(index_parts=False)
    shutil.rmtree(tmpdir)
    return centerline

def simplify_RN(OSM_RN,savefile,buffer_distance=50):
    road_aggregated,road_not_aggregated=road_divide(OSM_RN)
    centerline=get_centerline(road_aggregated) 
    centerline=repair_endpoint(centerline,road_aggregated,road_not_aggregated)
    result_road=extend_road(centerline,road_not_aggregated)
    result_road.to_file(savefile)
    return result_road

def build_nodes(seg_edges,savefile):
    all_points = [Point(pt) for line in seg_edges.geometry for pt in (line.coords[0], line.coords[-1])]
    seg_nodes = gpd.GeoDataFrame(
        geometry=gpd.GeoSeries(all_points).drop_duplicates().reset_index(drop=True),
        crs=seg_edges.crs
    )
    seg_nodes.to_file(savefile)
    return seg_nodes

def build_stroke(seg_edges,savefile):
    storkes = coins.COINS(seg_edges,angle_threshold=135)
    stroke_gdf = storkes.stroke_gdf()
    stroke_gdf.to_file(savefile)
    return stroke_gdf
    
def build_mixmap(seg_edges,roads_attri,savefile):
    seg_edges['name'] = ''
    master_buffered = seg_edges.copy()
    master_buffered['geometry'] = master_buffered.geometry.buffer(25)
    roads_sindex = roads_attri.sindex
    for buf_idx, buf_geom in master_buffered.geometry.items():
        possible_idx = list(roads_sindex.intersection(buf_geom.bounds))
        if not possible_idx:
            continue
        possible_roads = roads_attri.iloc[possible_idx]
        inters = possible_roads.geometry.intersection(buf_geom)
        lengths = inters.length
        if lengths.empty:
            continue    
        max_idx = lengths.idxmax()
        seg_edges.loc[buf_idx, 'name'] = roads_attri.loc[max_idx, 'name']
    named_edges = seg_edges[seg_edges['name']!=''] 
    if len(named_edges)==0:
        strokes = coins.COINS(seg_edges,angle_threshold=135)
        stroke_gdf = strokes.stroke_gdf()
        stroke_gdf.to_file(savefile)
    else:
        grouped = named_edges.groupby('name')
        merged_features = []
        for name, group in grouped:
            geometries = group.geometry.tolist()
            merged_line = shapely.ops.linemerge(geometries)
            if merged_line.geom_type == 'MultiLineString':
                for line in merged_line.geoms:
                    merged_features.append({
                        'name': name,
                        'geometry': line,
                    })
            else:
                merged_features.append({
                    'name': name,
                    'geometry': merged_line,
                })
        namedstreet_gdf = gpd.GeoDataFrame(merged_features, geometry='geometry', crs=named_edges.crs)
        
        unnamed_edges = seg_edges[seg_edges['name']=='']  # name为空
        strokes = coins.COINS(unnamed_edges,angle_threshold=135)
        stroke_gdf = strokes.stroke_gdf()
        mixed = pd.concat([namedstreet_gdf, stroke_gdf], ignore_index=True)
        mixed.to_file(savefile)
    return mixed

def build_and_save_network(edges_df, u,v,output_path):
    G = nx.Graph()
    G.add_nodes_from(pd.concat([edges_df[u], edges_df[v]]).unique())
    G.add_edges_from(zip(edges_df[u], edges_df[v]))
    nx.write_graphml(G, output_path)
    return 


def build_primal(seg_edges,seg_nodes,savefile):
    edge_endpoints = []

    sindex = seg_nodes.sindex
    for idx, edge in seg_edges.iterrows():
        geom = edge.geometry
        start_pt = Point(geom.coords[0])
        end_pt   = Point(geom.coords[-1])
        start_matches = list(sindex.intersection(start_pt.bounds))
        end_matches   = list(sindex.intersection(end_pt.bounds))

        edge_endpoints.append({
                "start_node_id": start_matches[0],
                "end_node_id": end_matches[0],
            })

    intersect_df = pd.DataFrame(edge_endpoints)

    build_and_save_network(intersect_df,'start_node_id','end_node_id',savefile)
    return
            
def build_dual(edges,savefile):
    results = []
    sindex = edges.sindex
    for idx, line in edges.iterrows():
        possible_matches_index = list(sindex.intersection(line.geometry.bounds))
        for match_idx in possible_matches_index:
            if match_idx == idx:
                continue
            other_line = edges.iloc[match_idx]
            if line.geometry.intersects(other_line.geometry):
                results.append({
                    'line_id': idx,
                    'intersecting_line_id': match_idx
                })
    intersect_df = pd.DataFrame(results)
    build_and_save_network(intersect_df,'line_id','intersecting_line_id',savefile)
    return

if __name__ == "__main__":
    OSM_RN=gpd.read_file(r'example/example.gpkg',layer='edges')
    
    segments=simplify_RN(OSM_RN,r'example/result/segment/shp/segments.shp')
    nodes=build_nodes(segments,r'example/result/segment/shp_nodes/nodes.shp')
    build_primal(segments,nodes,r'example/result/segment/primal/primal.graphml')
    build_dual(segments,r'example/result/segment/dual/segment_dual.graphml')

    strokes=build_stroke(segments,r'example/result/stroke map/shp/strokes.shp')
    build_dual(strokes,r'example/result/stroke map/dual/stroke_dual.graphml')

    segments=gpd.read_file(r'example/result/segment/shp/segments.shp')   
    mix_map=build_mixmap(segments,OSM_RN,r'example/result/combined map/shp/combined.shp')
    build_dual(mix_map,r'example/result/combined map/dual/combined_dual.graphml')



    
