# -*- coding: utf-8 -*-
"""
Created on Sun May 17 11:31:32 2020

@author: Administrator
"""
    
import pandas as pd
import numpy as np
import math
import csv
import time
import datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import Dial_module

g_nodes_list = []
g_links_list = []
g_demand_list = []

g_node_id_ListNum_dict={}
g_node_ListNum_id_dict={}
g_link_id_ListNum_dict={}
g_link_ListNum_id_dict={}
g_link_FromTo_ListNum_dict={}
g_link_ListNum_FromTo_dict={}
g_demand_id_ListNum_dict={}

g_nodes_No = 0
g_links_No = 0
g_demand_No = 0
g_outer_iteration_No = 20

Max_node_utility = 10000
mu=1.0
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
        self.LinkUtility = 0.0
        self.link_cap = 0
    
    def BPRFunction(self):
        
        self.link_cap = self.lane_cap*self.number_of_lanes
        self.LinkUtility = -self.FFTT_in_minute*(1+self.BPR_alpha*(self.flow_volume/max(0.00001,self.link_cap))**self.BPR_beta)
        
class Demand():
    def __init__(self):
        self.from_zone_id = 0
        self.to_zone_id = 0
        self.number_of_agents = 0
        self.demand_id = 0
    
        
# In[1] read data
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
                if g_nodes_No % 1000 == 0:
                    print("the number of nodes have been read: {}".format(g_nodes_No))
        print("Total number of nodes in the network is: {}".format(g_nodes_No))

    
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
                if g_links_No % 1000 == 0:
                    print("the number of links have been read: {}".format(g_links_No))
        print("Total number of links in the network is: {}".format(g_links_No))                


    with open('demand.csv') as file_object:
        lines = file_object.readlines()
        for line in lines[1:]:
            line = line.strip().split(',')
            demand = Demand()
        
            try:
                demand.from_zone_id = int(line[0])
                demand.to_zone_id = int(line[1])
                demand.number_of_agents = float(line[2])
                demand.demand_id = g_demand_No
            except:
                print("Python can't read demand.csv sucessfully, please check!")
            else:
                
                g_demand_No += 1
                g_demand_list.append(demand)
                if g_demand_No % 1000 == 0:
                    print("the number of demand have been read: {}".format(g_demand_No))
        print("Total number of demand in the network is: {}".format(g_demand_No))            

# In[2] calculate utility
def Utility(origin,destination,node_list,network_cur):

     
    #node utility
    node_utility=node_list[:]
    node_utility[node_list.index(destination)]=0
    
    SElist=[]
    SElist.append(destination)
    node_label=[0]*len(node_list)
    node_label[node_list.index(destination)]=1
    node_record=node_list[:]
    node_record.remove(destination)
    while len(SElist)>0 and len(node_record)>0:
        
        node=SElist[0]
        del SElist[0]
        
        node_ListNo=g_node_id_ListNum_dict[node]
        
        for from_node in g_nodes_list[node_ListNo].incoming_node_list:
            
            if from_node not in node_record:
                continue
            link=(from_node,node)
            
            if link not in network_cur:
                continue
           
            from_node_ListNo=g_node_id_ListNum_dict[from_node]
            
            value=1
            for from_node_to_node in g_nodes_list[from_node_ListNo].outgoing_node_list:
                link=(from_node,from_node_to_node)
                
                if link in network_cur:
                    value=value*node_label[node_list.index(from_node_to_node)]

            if value==1:
                sum_exp=0.0
                for from_node_to_node in g_nodes_list[from_node_ListNo].outgoing_node_list:
                    link=(from_node,from_node_to_node)
                    if link in network_cur:

                        link_ListNo=g_link_FromTo_ListNum_dict[link]
                        sum_exp+=math.exp(1/mu*(g_links_list[link_ListNo].LinkUtility+node_utility[node_list.index(from_node_to_node)]))


                node_utility[node_list.index(from_node)]=mu*math.log(sum_exp)
               
                # # # #
                if (node_utility[node_list.index(from_node)]<-600):
                    node_utility[node_list.index(from_node)]=-600
                # # # #

                node_label[node_list.index(from_node)]=1
                SElist.append(from_node)
                if len(node_record)>0:
                   
                    node_record.remove(from_node)
        
    #link systematic utility
    link_sys_utility=[0]*len(network_cur)
    
    for link in network_cur:
        link_ListNo=g_link_FromTo_ListNum_dict[link]
        to_node=link[1]
        link_sys_utility[network_cur.index(link)]=g_links_list[link_ListNo].LinkUtility+node_utility[node_list.index(to_node)]
        
        # # # #
        if (link_sys_utility[network_cur.index(link)]<-700):
            link_sys_utility[network_cur.index(link)]=-700
        # # # #

    return node_utility,link_sys_utility

# In[3] stochastic assignment
def RLSUE(network):
    link_flow_network=[]
    link_sys_utility_network=[]

    for demand in g_demand_list:

        origin=demand.from_zone_id
        destination=demand.to_zone_id
        network_cur=network[demand.demand_id]
  
        node_list=[]
        for link in network_cur:
            node_list.append(link[0])
            node_list.append(link[1])
        node_list=list(set(node_list))
            
        if len(node_list)==0:
            link_flow=[]
            link_sys_utility=[]
            
            link_flow_network.append(link_flow)
            link_sys_utility_network.append(link_sys_utility)
            
            #print('The from zone and to zone are the same zone in the {}th demand!'.format(str(demand.demand_id)))
            continue

        node_utility,link_sys_utility=Utility(origin,destination,node_list,network_cur)
        
        #link conditional probability
        link_con_prob=[0]*len(network_cur)
        node_NoD_list=node_list[:]
        node_NoD_list.remove(destination)
        
        for node in node_NoD_list:
            node_ListNo=g_node_id_ListNum_dict[node]
            sum_exp_p=0.0
            
            for to_node in g_nodes_list[node_ListNo].outgoing_node_list:
                link=(node,to_node)
                if link in network_cur:
                    link_ListNo=g_link_FromTo_ListNum_dict[link]
                    sum_exp_p+=math.exp(1/mu*(link_sys_utility[network_cur.index(link)]))
            
            for to_node in g_nodes_list[node_ListNo].outgoing_node_list:
                link=(node,to_node)
                if link in network_cur:
                    link_ListNo=g_link_FromTo_ListNum_dict[link]
                    link_V=math.exp(1/mu*(link_sys_utility[network_cur.index(link)]))
                    link_con_prob[network_cur.index(link)]=link_V/sum_exp_p
                    if((link_V==0) & (sum_exp_p==0)):
                        stop=0

        #link probability and flow
        link_prob=[0]*len(network_cur)
        link_flow=[0]*len(network_cur)
        link_label=[0]*len(network_cur)
        node_NoO_list=node_list[:]
        node_NoO_list.remove(origin)
        
        for link in network_cur:
            if link[0]==origin:
                link_prob[network_cur.index(link)]=link_con_prob[network_cur.index(link)]
                link_flow[network_cur.index(link)]=demand.number_of_agents*link_prob[network_cur.index(link)]
                link_label[network_cur.index(link)]=1
                
        
        while len(node_NoO_list)>0:
            for node in node_NoO_list:

                node_ListNo=g_node_id_ListNum_dict[node]
                
                value=1
                sum_prob=0
                for from_node in g_nodes_list[node_ListNo].incoming_node_list:
                    link=(from_node,node)
                    if link in network_cur:
                        value=value*link_label[network_cur.index(link)]
                        sum_prob+=link_prob[network_cur.index(link)]
                
                if value==1:
                    for to_node in g_nodes_list[node_ListNo].outgoing_node_list:
                        link=(node,to_node)
                        if link in network_cur:
                            link_prob[network_cur.index(link)]=sum_prob*link_con_prob[network_cur.index(link)]
                            link_flow[network_cur.index(link)]=demand.number_of_agents*link_prob[network_cur.index(link)]
                            link_label[network_cur.index(link)]=1
                    node_NoO_list.remove(node)    
        
        link_flow_network.append(link_flow)
        link_sys_utility_network.append(link_sys_utility)
       
    return link_flow_network,link_sys_utility_network
                
# In[4] refrence utility  
def MAXUility(link_flow_network,link_sys_utility_network,network):
    Pi_network=[]
    
    for demand in g_demand_list:
        
        network_cur=network[demand.demand_id]
        link_flow=link_flow_network[demand.demand_id]
        link_sys_utility=link_sys_utility_network[demand.demand_id]

        node_list=[]
        for link in network_cur:
            node_list.append(link[0])
            node_list.append(link[1])
        node_list=list(set(node_list))
        if len(node_list)==0:
            Pi=[]
            Pi_network.append(Pi)
            #print('The from zone and to zone are the same zone in the {}th demand!'.format(str(demand.demand_id)))
            continue
        Pi=node_list[:]
        Pi[node_list.index(demand.to_zone_id)]=0
        node_NoD_list=node_list[:]
        node_NoD_list.remove(demand.to_zone_id)
        
        
        for node in node_NoD_list:
            node_ListNo=g_node_id_ListNum_dict[node]
            
            sum_flow_pre=0
            if node==demand.from_zone_id:
                sum_flow_pre=demand.number_of_agents
            else:
                for from_node in g_nodes_list[node_ListNo].incoming_node_list:
                    link=(from_node,node)
                    if link in network_cur:
                        sum_flow_pre+=link_flow[network_cur.index(link)]
            
            sum_utility_suc=0
            for to_node in g_nodes_list[node_ListNo].outgoing_node_list:
                link=(node,to_node)
                if link in network_cur:
                    sum_utility_suc+=math.exp(1/mu*link_sys_utility[network_cur.index(link)])
            
            if sum_utility_suc<10e-200:
                Pi[node_list.index(node)]=10e-200
            else:
                Pi[node_list.index(node)]=math.log(sum_flow_pre/sum_utility_suc)
                 
        Pi_network.append(Pi)
        
    return Pi_network
     
# In[5] iterative calculation
       
def MSA():
    
    #step1 Initialilzation
    for link in g_links_list:
        link.BPRFunction()
    
    network,link_flow_integrated,link_flow_ODList=Dial_module.Dial_Results()
    
    # initialize variables
    link_flow=link_flow_integrated
    link_flow_OD=[]
    
   
    pi_OD=[]
    for demand in g_demand_list:
        network_cur=network[demand.demand_id]
        node_list=[]
        for link in network_cur:
            node_list.append(link[0])
            node_list.append(link[1])
        node_list=list(set(node_list))
        pi_OD.append([0.0]*len(node_list))
        link_flow_OD.append([0.0]*len(network_cur))
    #step2 update link cost
    for i in range(g_links_No):
        g_links_list[i].flow_volume=link_flow[i]
        g_links_list[i].BPRFunction()
      

    
    gap_total=[0]*(g_outer_iteration_No-1)
    for n in range(1,g_outer_iteration_No):
       
        #step3-4 compute link flow and pi
        link_flow_OD_aux,link_sys_utility_OD=RLSUE(network)
        Pi_OD_aux=MAXUility(link_flow_OD_aux,link_sys_utility_OD,network)

        #step5 update link flow and refrence node utility
        for i in range(g_links_No):
            g_links_list[i].flow_volume=0
        
        stop=0
        for k in range(g_demand_No):
            
            efficient_node_length=len(pi_OD[k])
            efficient_link_length=len(link_flow_OD[k])
            for i in range(efficient_link_length):
                link_flow_OD[k][i]=link_flow_OD[k][i]+1/n*(link_flow_OD_aux[k][i]-link_flow_OD[k][i])
                # if(n == 1):
                #     link_flow_OD[k][i] = link_flow_OD_aux[k][i]
                # else:
                #     link_flow_OD[k][i] = link_flow_OD[k][i]*(pow(link_flow_OD_aux[k][i]/link_flow_OD[k][i], 1/n))
                
                link=network[k][i]
                link_ListNo=g_link_FromTo_ListNum_dict[link]
                g_links_list[link_ListNo].flow_volume+=link_flow_OD[k][i]
                
            
            # for j in range(efficient_node_length):
                # pi_OD[k][j]=pi_OD[k][j]+1/n*(Pi_OD_aux[k][j]-pi_OD[k][j])
                # if(n == 1):
                # pi_OD[k][j]=Pi_OD_aux[k][j]
                # else:
                #     if(pi_OD[k][j] != 0):
                #         pi_OD[k][j]=pi_OD[k][j]*(pow(Pi_OD_aux[k][j]/pi_OD[k][j], 1/n))
                #         if(pi_OD[k][j]<0):
                #             stop=0
        

        #step6-7 update link systematic utility
        for i in range(g_links_No):
            g_links_list[i].BPRFunction()

        link_flow_updated,link_sys_utility_updated=RLSUE(network)
        Pi_OD_updated=MAXUility(link_flow_updated,link_sys_utility_updated,network) 

        gap=0
        for demand in g_demand_list:

            origin=demand.from_zone_id
            destination=demand.to_zone_id
            network_cur=network[demand.demand_id]
            pi_cur=Pi_OD_updated[demand.demand_id]
            
            if len(network_cur)==0: 
                continue

            node_list=[]
            for link in network_cur:
                node_list.append(link[0])
                node_list.append(link[1])
            node_list=list(set(node_list))
            
            # node_utility,link_sys_utility_cur=Utility(origin,destination,node_list,network_cur)
            link_sys_utility_cur=link_sys_utility_updated[demand.demand_id]
            link_flow_cur=link_flow_OD[demand.demand_id]
            
            for k in range(len(network_cur)):
                from_node=network_cur[k][0]
                
                if link_flow_cur[k]<10e-5:
                    continue

                gap_abs=math.log(link_flow_cur[k])-1/mu*link_sys_utility_cur[k]-pi_cur[node_list.index(from_node)]

                gap+=(gap_abs**2)
        gap_total[n-1]=gap/2   
    return gap_total

# In[6]
def Output():
    with open('link_performance.csv','w',newline='') as file_object:
        csv_write = csv.writer(file_object)
        csv_write.writerow(['road_link_id','from_node_id','to_node_id','time_period','volume','travel_time','speed','notes'])
       
        for link in g_links_list:
            line=[link.RoadLink_id,link.from_node_id,link.to_node_id,'',link.flow_volume,-link.LinkUtility,'120','period-based']
            csv_write.writerow(line)
 
# In[6] main function    
if __name__=='__main__':
    
    print('Time Now:',time.strftime('%H:%M:%S',time.localtime(time.time())))
    t0=datetime.datetime.now()
    print('Reading data......')
    ReadData()
    t1=datetime.datetime.now()
    print('Reading data have spent '+str(t1-t0) +'.')
    
    print('=================================')
    print('The main procedure is running now')
    t2=datetime.datetime.now()
    gap=MSA()
    t3=datetime.datetime.now()
    print(gap)
    print('The main procedure have spent '+str(t3-t2) +'.')
    
    Output()
    
    #visualization
    fig,ax=plt.subplots(figsize=(14,7))
    
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.plot(list(range(1,g_outer_iteration_No)),gap)

    plt.tick_params(labelsize=15)
    labels=ax.get_xticklabels()+ax.get_yticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    
    ax.set_xlabel('Iteration number',fontdict={'family':'Times New Roman','size':18})
    ax.set_ylabel('Gap value',fontdict={'family':'Times New Roman','size':18})
    
    plt.savefig('gap')
    plt.show()
