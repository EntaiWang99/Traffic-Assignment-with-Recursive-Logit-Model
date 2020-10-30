# -*- coding: utf-8 -*-
"""
Created on Sun May 17 11:31:32 2020

@author: Administrator
"""

import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt

g_nodes_list = []
g_links_list = []
g_demand_list = []

g_node_id_ListNum_dict={}
g_node_ListNum_id_dict={}
g_link_id_ListNum_dict={}
g_link_ListNum_id_dict={}
g_link_FromTo_ListNum_dict={}
g_link_ListNum_FromTo_dict={}
g_destination_demand_dict={}

g_nodes_No = 0
g_links_No = 0
g_demand_No = 0
g_outer_iteration_No = 11

Max_node_utility = 10000
theta_utility = -0.5
fai_logit = -0.8
b_Dial=1.0

# In[1]

class Node():
    def __init__(self):
        self.node_id = 0
        self.x_coord = 0.0
        self.y_coord = 0.0
        self.outgoing_node_list=[]
        self.incoming_node_list=[]
        
        
class Road_link():
    def __init__(self):
        self.RoadLink_id = 0
        self.from_node_id = 0
        self.to_node_id = 0
        self.length = 0.0
        self.number_of_lanes = 0
        self.lane_cap = 0
        self.BPR_alpha = 0.0
        self.BPR_beta = 0.0
        self.flow_volume = 0.0
        self.FFTT_in_minute = 0
        self.cost = 0
        self.LinkTravelTime_in_minute = 0.0
        self.link_cap = 0
    
    def BPRFunction(self):
        self.link_cap = self.lane_cap*self.number_of_lanes
        self.LinkTravelTime_in_minute = self.FFTT_in_minute*(1+self.BPR_alpha*(self.flow_volume/max(0.00001,self.link_cap))**self.BPR_alpha)
        
class Demand():
    def __init__(self):
        self.from_zone_id = 0
        self.to_zone_id = 0
        self.number_of_agents = 0
    
        
# In[2]
def ReadData():
    
    global g_nodes_No 
    global g_links_No 
    global g_demand_No 
    
    with open('node.csv') as file_object:
        lines = file_object.readlines()
        for line in lines[1:]:
            line = line.strip().split(',')
            node = Node()
            try:
                node.node_id = int(line[1])
                node.x_coord = line[5]
                node.x_coord = line[6]
            except:
                print("Python can't read node.csv sucessfully, please check!")
            else:
                g_node_id_ListNum_dict[node.node_id]=g_nodes_No
                g_node_ListNum_id_dict[g_nodes_No]=node.node_id
                
                g_nodes_No += 1
                g_nodes_list.append(node)


    
    with open('road_link.csv') as file_object:
        lines = file_object.readlines()
        for line in lines[1:]:
            line = line.strip().split(',')
            link = Road_link()
                        
            try:
                link.RoadLink_id = int(line[1])
                link.from_node_id = int(line[2])
                link.to_node_id = int(line[3])
                link.length = float(line[6])
                link.number_of_lanes = int(line[7])
                link.lane_cap = float(line[8])
                link.FFTT_in_minute = float(line[12])
                link.BPR_alpha = float(line[14])
                link.BPR_beta = float(line[15])
                link.cost = float(line[11])
                
                
            
            except:
                print("Python can't read road_link.csv sucessfully, please check!")
            else:
                
                from_node_ListNum=g_node_id_ListNum_dict[link.from_node_id]               
                g_nodes_list[from_node_ListNum].outgoing_node_list.append(link.to_node_id)                
                to_node_ListNum=g_node_id_ListNum_dict[link.to_node_id]
                g_nodes_list[to_node_ListNum].incoming_node_list.append(link.from_node_id)
                
                link_FromTo_nodes = (link.from_node_id,link.to_node_id)
                g_link_FromTo_ListNum_dict[link_FromTo_nodes]=g_links_No
                g_link_ListNum_FromTo_dict[g_links_No]=link_FromTo_nodes
                
                g_link_id_ListNum_dict[link.RoadLink_id]=g_links_No
                g_link_ListNum_id_dict[g_links_No]=link.RoadLink_id
                g_links_No += 1
                g_links_list.append(link)

               


    with open('demand.csv') as file_object:
        lines = file_object.readlines()
        for line in lines[1:]:

            line = line.strip().split(',')
            demand = Demand()
        
            try:
                demand.from_zone_id = int(line[0])
                demand.to_zone_id = int(line[1])
                demand.number_of_agents = float(line[2])
            except:
                print("Python can't read demand.csv sucessfully, please check!")
            else:
                
                g_demand_No += 1
                g_demand_list.append(demand)

     
# In[3]#shortest path from one orgin node to all other nodes

def DIJKSTRA(origin):
    
    origin_ListNo=g_node_id_ListNum_dict[origin]
    
    S_permanent=[]
    S_permanent.append(origin)
    S_temporary=[]
    for node in g_nodes_list:
        S_temporary.append(node.node_id)
    S_temporary.remove(origin)
    
    predecessor=[0]*g_nodes_No
    distance=[9999]*g_nodes_No

    for node in g_nodes_list[origin_ListNo].outgoing_node_list:

        node_ListNo=g_node_id_ListNum_dict[node]
        link=(origin,node)
        link_ListNo=g_link_FromTo_ListNum_dict[link]
        distance[node_ListNo]=g_links_list[link_ListNo].LinkTravelTime_in_minute
        predecessor[node_ListNo]=origin

    while len(S_permanent)<g_nodes_No:
        S_temporary_ListNo=[g_node_id_ListNum_dict[i] for i in S_temporary]

        node_ListNo_minDistance=9999
        min_distance=9999
        
        for temporary_node in S_temporary_ListNo:
            if distance[temporary_node]<min_distance:
                min_distance=distance[temporary_node]
                node_ListNo_minDistance=temporary_node
       
        if node_ListNo_minDistance==9999:
            break
        
        node_minDistance=g_node_ListNum_id_dict[node_ListNo_minDistance]
        S_permanent.append(node_minDistance)
        S_temporary.remove(node_minDistance)
        
        for to_node in g_nodes_list[node_ListNo_minDistance].outgoing_node_list:
            to_node_ListNo=g_node_id_ListNum_dict[to_node]
            link=(node_minDistance,to_node)
            link_ListNo=g_link_FromTo_ListNum_dict[link]
            
            if distance[to_node_ListNo]>(distance[node_ListNo_minDistance]+g_links_list[link_ListNo].LinkTravelTime_in_minute):
                distance[to_node_ListNo]=distance[node_ListNo_minDistance]+g_links_list[link_ListNo].LinkTravelTime_in_minute
                predecessor[to_node_ListNo]=node_minDistance
                
    distance[origin_ListNo]=0       
    return distance,predecessor       

# In[3]#shortest path from all nodes to one destination node
def s_DIJKSTRA(destination):
    distance=[9999]*g_nodes_No
    destination_ListNo=g_node_id_ListNum_dict[destination]
    
    for node_class in g_nodes_list:
        node=node_class.node_id
        node_ListNo=g_node_id_ListNum_dict[node]
        r_distance, predecessor=DIJKSTRA(node)
        distance[node_ListNo]=r_distance[destination_ListNo]
   
    return distance

# In[3]        
def DialAlgorithm():
    network=[]
    link_flow=[]
    link_flow_ODList=[]
    for demand in g_demand_list:

        network_link_list=[]
        
        #step1 
        origin=demand.from_zone_id
        destination=demand.to_zone_id
        origin_ListNo=g_node_id_ListNum_dict[origin]
        r_distance,predecessor=DIJKSTRA(origin)
        s_distance=s_DIJKSTRA(destination)
        
    
        link_Likelihood=[0]*g_links_No
        for link in g_links_list:
            
            
            from_node=link.from_node_id
            to_node=link.to_node_id
            from_node_ListNo=g_node_id_ListNum_dict[from_node]
            to_node_ListNo=g_node_id_ListNum_dict[to_node]
            link_ListNo=g_link_id_ListNum_dict[link.RoadLink_id]
            
            if r_distance[from_node_ListNo]<r_distance[to_node_ListNo] and s_distance[from_node_ListNo]>s_distance[to_node_ListNo]:
                
                link_Likelihood[link_ListNo]=math.exp(b_Dial*(r_distance[to_node_ListNo]-r_distance[from_node_ListNo]-link.LinkTravelTime_in_minute))

        #step2
        r_distance_array=np.array(r_distance)
        r_index_list=r_distance_array.argsort()
        s_distance_array=np.array(s_distance)
        s_index_list=s_distance_array.argsort()
        
        link_Weights_list=[0.0]*g_links_No
        link_label_list=[0.0]*g_links_No
        for node_ListNo in r_index_list:
            
            node=g_node_ListNum_id_dict[node_ListNo]
            if node==destination:
                break

            if node==origin:
                for to_node in g_nodes_list[node_ListNo].outgoing_node_list:
                    link=(node,to_node)
                    link_ListNo=g_link_FromTo_ListNum_dict[link]
                    link_Weights_list[link_ListNo]=link_Likelihood[link_ListNo]
                    link_label_list[link_ListNo]=1.0
                
            else:
                current_label=1.0
                for from_node in g_nodes_list[node_ListNo].incoming_node_list:
                    link_pre=(from_node,node)
                    link_ListNo_pre=g_link_FromTo_ListNum_dict[link_pre]
                    if link_Likelihood[link_ListNo_pre]>0:
                        current_label=current_label*link_label_list[link_ListNo_pre]

                if current_label==1.0:
                    for to_node in g_nodes_list[node_ListNo].outgoing_node_list:
                        link=(node,to_node)
                        link_ListNo=g_link_FromTo_ListNum_dict[link]
                        sum_w=0.0
                        for from_node in g_nodes_list[node_ListNo].incoming_node_list:
                            link_pre=(from_node,node)
                            link_ListNo_pre=g_link_FromTo_ListNum_dict[link_pre]
                            sum_w+=link_Weights_list[link_ListNo_pre]

                        link_Weights_list[link_ListNo]=link_Likelihood[link_ListNo]*sum_w
                        link_label_list[link_ListNo]=1.0
       
        #step3
        link_flow_list=[0.0]*g_links_No
        link_label_list=[0.0]*g_links_No

        for node_ListNo in s_index_list:
            node=g_node_ListNum_id_dict[node_ListNo]
            
            if node==origin:
                break
            
            if node==destination:
                sum_w_pre=0
                for from_node in g_nodes_list[node_ListNo].incoming_node_list:
                    link=(from_node,node)
                    link_ListNo=g_link_FromTo_ListNum_dict[link]
                    sum_w_pre+=link_Weights_list[link_ListNo]
                
                for from_node in g_nodes_list[node_ListNo].incoming_node_list:
                    link=(from_node,node)
                    link_ListNo=g_link_FromTo_ListNum_dict[link]
                    
                    if link_Weights_list[link_ListNo]==0:
                        link_flow_list[link_ListNo]=0.0
                        link_label_list[link_ListNo]=1.0
                    else:
                        link_flow_list[link_ListNo]=link_Weights_list[link_ListNo]/sum_w_pre*demand.number_of_agents
                        link_label_list[link_ListNo]=1.0
            else:
                current_label=1.0
                sum_flow=0.0
                for to_node in g_nodes_list[node_ListNo].outgoing_node_list:
                    link_suc=(node,to_node)
                    link_ListNo_suc=g_link_FromTo_ListNum_dict[link_suc]
                    if link_Weights_list[link_ListNo_suc]>0:
                        current_label=current_label*link_label_list[link_ListNo_suc]
                        sum_flow+=link_flow_list[link_ListNo_suc]
                
                
                    
                if current_label==1.0:
                    sum_w_pre=0
                    for from_node in g_nodes_list[node_ListNo].incoming_node_list:
                        link=(from_node,node)
                        link_ListNo=g_link_FromTo_ListNum_dict[link]
                        sum_w_pre+=link_Weights_list[link_ListNo]
                    
                    for from_node in g_nodes_list[node_ListNo].incoming_node_list:
                        link=(from_node,node)

                        link_ListNo=g_link_FromTo_ListNum_dict[link]

                        if link_Weights_list[link_ListNo]==0:
                            link_flow_list[link_ListNo]=0.0
                            link_label_list[link_ListNo]=1.0
                        else:
                            link_flow_list[link_ListNo]=(link_Weights_list[link_ListNo]/sum_w_pre)*sum_flow
                            link_label_list[link_ListNo]=1.0
                        
                        
        for i in range(g_links_No):
            if link_flow_list[i]>0:
                link=g_link_ListNum_FromTo_dict[i]
                network_link_list.append(link)
 
        network.append(network_link_list)
        link_flow_ODList.append(link_flow_list)
    
    link_flow=np.sum(link_flow_ODList,axis=0)
    
    return network,link_flow, link_flow_ODList
# In[3]
def Dial_Results():
    ReadData()
    for link in g_links_list:
        link.BPRFunction()
    network,link_flow,link_flow_ODList=DialAlgorithm()
    return network,link_flow,link_flow_ODList
        
# In[3]
if __name__=='__main__':
    print('Reading data......')
    
    network,link_flow,link_flow_ODList=Dial_Results()
    print(network)
    print(link_flow)
    
    

