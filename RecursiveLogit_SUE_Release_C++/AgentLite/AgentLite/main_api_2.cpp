#pragma warning(disable: 4305 4267 4018) 
#pragma warning(disable: 6385 6386)
#include <iostream>
#include <fstream>
#include <list> 
#include <omp.h>
#include <algorithm>
#include <time.h>
#include <functional>
#include <stdio.h>   
#include <math.h>

#include <numeric> 
#include <stack>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <iostream>
#include <iomanip>
using namespace std;
using std::string;
using std::ifstream;
using std::vector;
using std::map;
using std::istringstream;
using std::max;
template <typename T>


// some basic parameters setting

//Pls make sure the _MAX_K_PATH > Agentlite.cpp's g_number_of_K_paths+g_reassignment_number_of_K_paths and the _MAX_ZONE remain the same with .cpp's defination
#define _MAX_LABEL_COST 1.0e+15

#define _MAX_AGNETTYPES 4 //because of the od demand store format,the MAX_demandtype must >=g_DEMANDTYPES.size()+1
#define _MAX_TIMEPERIODS 1
#define _MAX_MEMORY_BLOCKS 20

#define _MAX_LINK_SIZE_IN_A_PATH 1000
#define _MAX_LINK_SIZE_FOR_A_NODE 200

#define _MAX_TIMESLOT_PerPeriod 100 // max 96 15-min slots per day

#define MIN_PER_TIMESLOT 15
#define number_of_seconds_per_interval 6

// Linear congruential generator 
#define LCG_a 17364
#define LCG_c 0
#define LCG_M 65521  // it should be 2^32, but we use a small 16-bit number to save memory



#define sprintf_s sprintf

FILE* g_pFileOutputLog = NULL;
int g_debugging_flag = 2;

#define STRING_LENGTH_PER_LINE 20000

void fopen_ss(FILE** file, const char* fileName, const char* mode)
{
	*file = fopen(fileName, mode);
}


void g_ProgramStop();


//below shows where the functions used in Agentlite.cpp come from!
//Utility.cpp

#pragma warning(disable: 4244)  // stop warning: "conversion from 'int' to 'float', possible loss of data"


class CCSVParser
{
public:
	char Delimiter;
	bool IsFirstLineHeader;
	ifstream inFile;
	string mFileName;
	vector<string> LineFieldsValue;
	vector<string> Headers;
	map<string, int> FieldsIndices;

	vector<int> LineIntegerVector;

public:
	void  ConvertLineStringValueToIntegers()
	{
		LineIntegerVector.clear();
		for (unsigned i = 0; i < LineFieldsValue.size(); i++)
		{
			std::string si = LineFieldsValue[i];
			int value = atoi(si.c_str());

			if (value >= 1)
				LineIntegerVector.push_back(value);

		}
	}
	vector<string> GetHeaderVector()
	{
		return Headers;
	}

	bool m_bDataHubSingleCSVFile;
	string m_DataHubSectionName;
	bool m_bLastSectionRead;

	bool m_bSkipFirstLine;  // for DataHub CSV files

	CCSVParser(void)
	{
		Delimiter = ',';
		IsFirstLineHeader = true;
		m_bSkipFirstLine = false;
		m_bDataHubSingleCSVFile = false;
		m_bLastSectionRead = false;
	}

	~CCSVParser(void)
	{
		if (inFile.is_open()) inFile.close();
	}


	bool OpenCSVFile(string fileName, bool b_required)
	{
		mFileName = fileName;
		inFile.open(fileName.c_str());

		if (inFile.is_open())
		{
			if (IsFirstLineHeader)
			{
				string s;
				std::getline(inFile, s);
				vector<string> FieldNames = ParseLine(s);

				for (size_t i = 0; i < FieldNames.size(); i++)
				{
					string tmp_str = FieldNames.at(i);
					size_t start = tmp_str.find_first_not_of(" ");

					string name;
					if (start == string::npos)
					{
						name = "";
					}
					else
					{
						name = tmp_str.substr(start);
						//			TRACE("%s,", name.c_str());
					}


					FieldsIndices[name] = (int)i;
				}
			}

			return true;
		}
		else
		{
			if (b_required)
			{

				cout << "File " << fileName << " does not exist. Please check." << endl;
				//g_ProgramStop();
			}
			return false;
		}
	}


	void CloseCSVFile(void)
	{
		inFile.close();
	}



	bool ReadRecord()
	{
		LineFieldsValue.clear();

		if (inFile.is_open())
		{
			string s;
			std::getline(inFile, s);
			if (s.length() > 0)
			{

				LineFieldsValue = ParseLine(s);

				return true;
			}
			else
			{

				return false;
			}
		}
		else
		{
			return false;
		}
	}

	vector<string> ParseLine(string line)
	{
		vector<string> SeperatedStrings;
		string subStr;

		if (line.length() == 0)
			return SeperatedStrings;

		istringstream ss(line);


		if (line.find_first_of('"') == string::npos)
		{

			while (std::getline(ss, subStr, Delimiter))
			{
				SeperatedStrings.push_back(subStr);
			}

			if (line.at(line.length() - 1) == ',')
			{
				SeperatedStrings.push_back("");
			}
		}
		else
		{
			while (line.length() > 0)
			{
				size_t n1 = line.find_first_of(',');
				size_t n2 = line.find_first_of('"');

				if (n1 == string::npos && n2 == string::npos) //last field without double quotes
				{
					subStr = line;
					SeperatedStrings.push_back(subStr);
					break;
				}

				if (n1 == string::npos && n2 != string::npos) //last field with double quotes
				{
					size_t n3 = line.find_first_of('"', n2 + 1); // second double quote

					//extract content from double quotes
					subStr = line.substr(n2 + 1, n3 - n2 - 1);
					SeperatedStrings.push_back(subStr);

					break;
				}

				if (n1 != string::npos && (n1 < n2 || n2 == string::npos))
				{
					subStr = line.substr(0, n1);
					SeperatedStrings.push_back(subStr);
					if (n1 < line.length() - 1)
					{
						line = line.substr(n1 + 1);
					}
					else // comma is the last char in the line string, push an empty string to the back of vector
					{
						SeperatedStrings.push_back("");
						break;
					}
				}

				if (n1 != string::npos && n2 != string::npos && n2 < n1)
				{
					size_t n3 = line.find_first_of('"', n2 + 1); // second double quote
					subStr = line.substr(n2 + 1, n3 - n2 - 1);
					SeperatedStrings.push_back(subStr);
					size_t idx = line.find_first_of(',', n3 + 1);

					if (idx != string::npos)
					{
						line = line.substr(idx + 1);
					}
					else
					{
						break;
					}
				}
			}

		}

		return SeperatedStrings;
	}

	template <class T> bool GetValueByFieldName(string field_name, T& value, bool NonnegativeFlag = true, bool required_field = true)
	{

		if (FieldsIndices.find(field_name) == FieldsIndices.end())
		{
			if (required_field)
			{
				cout << "Field " << field_name << " in file " << mFileName << " does not exist. Please check the file." << endl;

				g_ProgramStop();
			}
			return false;
		}
		else
		{
			if (LineFieldsValue.size() == 0)
			{
				return false;
			}

			int size = (int)(LineFieldsValue.size());
			if (FieldsIndices[field_name] >= size)
			{
				return false;
			}

			string str_value = LineFieldsValue[FieldsIndices[field_name]];

			if (str_value.length() <= 0)
			{
				return false;
			}

			istringstream ss(str_value);

			T converted_value;
			ss >> converted_value;

			if (/*!ss.eof() || */ ss.fail())
			{
				return false;
			}

			if (NonnegativeFlag && converted_value < 0)
				converted_value = 0;

			value = converted_value;
			return true;
		}
	}


	bool GetValueByFieldName(string field_name, string& value)
	{
		if (FieldsIndices.find(field_name) == FieldsIndices.end())
		{
			return false;
		}
		else
		{
			if (LineFieldsValue.size() == 0)
			{
				return false;
			}

			unsigned int index = FieldsIndices[field_name];
			if (index >= LineFieldsValue.size())
			{
				return false;
			}
			string str_value = LineFieldsValue[index];

			if (str_value.length() <= 0)
			{
				return false;
			}

			value = str_value;
			return true;
		}
	}

};



template <typename T>
T** AllocateDynamicArray(int nRows, int nCols)
{
	T** dynamicArray;

	dynamicArray = new (std::nothrow) T * [nRows];

	if (dynamicArray == NULL)
	{
		cout << "Error: insufficient memory.";
		g_ProgramStop();

	}

	for (int i = 0; i < nRows; i++)
	{
		dynamicArray[i] = new (std::nothrow) T[nCols];

		if (dynamicArray[i] == NULL)
		{
			cout << "Error: insufficient memory.";
			g_ProgramStop();
		}

	}

	return dynamicArray;
}

template <typename T>
void DeallocateDynamicArray(T** dArray, int nRows, int nCols)
{
	if (!dArray)
		return;

	for (int x = 0; x < nRows; x++)
	{
		delete[] dArray[x];
	}

	delete[] dArray;

}


template <typename T>
T*** Allocate3DDynamicArray(int nX, int nY, int nZ)
{
	T*** dynamicArray;

	dynamicArray = new (std::nothrow) T * *[nX];

	if (dynamicArray == NULL)
	{
		cout << "Error: insufficient memory.";
		g_ProgramStop();
	}

	for (int x = 0; x < nX; x++)
	{
		if (x % 1000 == 0)
		{
			cout << "allocating 3D memory for " << x << endl;
		}


		dynamicArray[x] = new (std::nothrow) T * [nY];

		if (dynamicArray[x] == NULL)
		{
			cout << "Error: insufficient memory.";
			g_ProgramStop();
		}

		for (int y = 0; y < nY; y++)
		{
			dynamicArray[x][y] = new (std::nothrow) T[nZ];
			if (dynamicArray[x][y] == NULL)
			{
				cout << "Error: insufficient memory.";
				g_ProgramStop();
			}
		}
	}

	for (int x = 0; x < nX; x++)
		for (int y = 0; y < nY; y++)
			for (int z = 0; z < nZ; z++)
			{
				dynamicArray[x][y][z] = 0;
			}
	return dynamicArray;
}

template <typename T>
void Deallocate3DDynamicArray(T*** dArray, int nX, int nY)
{
	if (!dArray)
		return;
	for (int x = 0; x < nX; x++)
	{
		for (int y = 0; y < nY; y++)
		{
			delete[] dArray[x][y];
		}

		delete[] dArray[x];
	}

	delete[] dArray;

}

template <typename T>
T**** Allocate4DDynamicArray(int nM, int nX, int nY, int nZ)
{
	T**** dynamicArray;

	dynamicArray = new (std::nothrow) T * **[nX];

	if (dynamicArray == NULL)
	{
		cout << "Error: insufficient memory.";
		g_ProgramStop();
	}
	for (int m = 0; m < nM; m++)
	{
		if (m % 1000 == 0)
			cout << "allocating 4D memory for " << m << " zones" << endl;

		dynamicArray[m] = new (std::nothrow) T * *[nX];

		if (dynamicArray[m] == NULL)
		{
			cout << "Error: insufficient memory.";
			g_ProgramStop();
		}

		for (int x = 0; x < nX; x++)
		{
			dynamicArray[m][x] = new (std::nothrow) T * [nY];

			if (dynamicArray[m][x] == NULL)
			{
				cout << "Error: insufficient memory.";
				g_ProgramStop();
			}

			for (int y = 0; y < nY; y++)
			{
				dynamicArray[m][x][y] = new (std::nothrow) T[nZ];
				if (dynamicArray[m][x][y] == NULL)
				{
					cout << "Error: insufficient memory.";
					g_ProgramStop();
				}
			}
		}
	}
	return dynamicArray;

}

template <typename T>
void Deallocate4DDynamicArray(T**** dArray, int nM, int nX, int nY)
{
	if (!dArray)
		return;
	for (int m = 0; m < nM; m++)
	{
		for (int x = 0; x < nX; x++)
		{
			for (int y = 0; y < nY; y++)
			{
				delete[] dArray[m][x][y];
			}
			delete[] dArray[m][x];
		}
		delete[] dArray[m];
	}
	delete[] dArray;
}


//struct MyException : public exception {
//	const char * what() const throw () {
//		return "C++ Exception";
//	}
//};
//

class CDemand_Period {
public:

	CDemand_Period()
	{
		demand_period_id = 0;
		starting_time_slot_no = 0;
		ending_time_slot_no = 0;

	}
	string demand_period;
	string time_period;
	int demand_period_id;
	int starting_time_slot_no;
	int ending_time_slot_no;

	int get_time_horizon_in_min()
	{
		return (ending_time_slot_no - starting_time_slot_no) * 15;
	}

};


class CAgent_type {
public:
	CAgent_type()
	{
		value_of_time = 1;
		agent_type_no = 0;
		flow_type = 0;
	}

	int agent_type_no;
	string agent_type;
	float value_of_time;  // dollar per hour
	std::map<int, float> PCE_link_type_map;  // link type, product consumption equivalent used, for travel time calculation
	std::map<int, float> CRU_link_type_map;  // link type, 	Coefficient of Resource Utilization - CRU, for resource constraints 
	int flow_type; // not enter the column_pool optimization process. 0: continuous, 1: fixed, 2 discrete.

};

class CLinkType
{
public:
	int link_type;
	string link_type_name;
	string agent_type_list;
	string type_code;

	int number_of_links;

	CLinkType()
	{
		number_of_links = 0;
		link_type = 1;
	}

	bool AllowAgentType(string agent_type)
	{
		if (agent_type_list.size() == 0)  // if the agent_type_list is empty then all types are allowed.
			return true;
		else
		{
			if (agent_type_list.find(agent_type) != string::npos)  // otherwise, only an agent type is listed in this "white list", then this agent is allowed to travel on this link
				return true;
			else
				return false;


		}
	}


};


class CColumnPath {
public:
	int* path_node_vector;
	int* path_link_vector;

	int m_node_size;
	int m_link_size;
	std::vector<int> agent_simu_id_vector;

	void AllocateVector(int node_size, int* node_vector, int link_size, int* link_vector)
	{
		m_node_size = node_size;
		m_link_size = link_size;
		path_node_vector = new int[node_size];  // dynamic array
		path_link_vector = new int[link_size];

		for (int i = 0; i < m_node_size; i++)  // copy backward
		{
			path_node_vector[i] = node_vector[m_node_size - 1 - i];
		}

		for (int i = 0; i < m_link_size; i++)
		{
			path_link_vector[i] = link_vector[m_link_size - 1 - i];
		}


	}


	int path_seq_no;

	float path_volume;  // path volume
	float path_switch_volume;  // path volume
	float path_travel_time;
	float path_distance;
	float path_toll;
	float path_gradient_cost;  // first order graident cost.
	float path_gradient_cost_difference;  // first order graident cost - least gradient cost
	float path_gradient_cost_relative_difference;  // first order graident cost - least gradient cost


	CColumnPath()
	{
		path_node_vector = NULL;
		path_link_vector = NULL;

		path_switch_volume = 0;
		path_seq_no = 0;
		path_toll = 0;
		path_volume = 0;
		path_travel_time = 0;
		path_distance = 0;
		path_gradient_cost = 0;
		path_gradient_cost_difference = 0;
		path_gradient_cost_relative_difference = 0;
	}

	~CColumnPath()
	{
		if (m_node_size >= 1)
		{
			delete path_node_vector;
			delete path_link_vector;
		}

	}
};

class CAgentPath {
public:
	int path_id;
	int o_node_no;
	int d_node_no;
	float travel_time;
	float distance;
	float volume;
	int node_sum;

	std::vector <int> path_link_sequence;

	CAgentPath()
	{
		path_id = 0;
		node_sum = -1;

		travel_time = 0;
		distance = 0;
		volume = 0;
	}
};

class CColumnVector {
public:
	float cost;
	float time;
	float distance;
	float od_volume;  // od volume

	std::vector <CAgentPath>  discrete_agent_path_vector;  // first key is the sum of node id;. e.g. node 1, 3, 2, sum of those node ids is 6, 1, 4, 2 then node sum is 7.

	std::map <int, CColumnPath> path_node_sequence_map;  // first key is the sum of node id;. e.g. node 1, 3, 2, sum of those node ids is 6, 1, 4, 2 then node sum is 7.
	// this is colletion of unique paths
	CColumnVector()
	{
		od_volume = 0;
		cost = 0;
		time = 0;
		distance = 0;
	}
};

class CAgent_Column {
public:


	int agent_id;
	int o_zone_id;
	int d_zone_id;
	int o_node_id;
	int d_node_id;
	string agent_type;
	string demand_period;
	float volume;
	float cost;
	float travel_time;
	float distance;
	vector<int> path_node_vector;
	vector<int> path_link_vector;
	vector<float> path_time_vector;

	CAgent_Column()
	{
		cost = 0;
	}


};

// event structure in this "event-based" traffic simulation

class DTAVehListPerTimeInterval
{
public:
	std::vector<int> m_AgentIDVector;

};

std::map<int, DTAVehListPerTimeInterval> g_AgentTDListMap;

class CAgent_Simu
{
public:
	unsigned int m_RandomSeed;
	bool m_bGenereated;
	CAgent_Simu()
	{
		agent_vector_seq_no = -1;
		m_bGenereated = false;

		path_toll = 0;
		m_Veh_LinkArrivalTime_in_simu_interval = NULL;
		m_Veh_LinkDepartureTime_in_simu_interval = NULL;
		m_bCompleteTrip = false;
		departure_time_in_min = 0;
	}
	~CAgent_Simu()
	{
		DeallocateMemory();
	}


	float GetRandomRatio()
	{
		m_RandomSeed = (LCG_a * m_RandomSeed + LCG_c) % LCG_M;  //m_RandomSeed is automatically updated.

		return float(m_RandomSeed) / LCG_M;
	}



	int fixed_path_flag;
	int demand_type;
	int agent_id;
	int agent_vector_seq_no;
	int agent_service_type;

	float departure_time_in_min;
	int departure_time_in_simu_interval;

	float arrival_time_in_min;

	float path_toll;
	std::vector<int> path_link_seq_no_vector;
	std::vector<float> time_seq_no_vector;
	std::vector<int> path_timestamp_vector;;

	int m_path_link_seq_no_vector_size;

	std::vector<int> path_node_id_vector;

	int m_current_link_seq_no;
	int* m_Veh_LinkArrivalTime_in_simu_interval;
	int* m_Veh_LinkDepartureTime_in_simu_interval;


	bool m_bCompleteTrip;

	void AllocateMemory()
	{
		if (m_Veh_LinkArrivalTime_in_simu_interval == NULL)
		{
			m_current_link_seq_no = 0;
			m_Veh_LinkArrivalTime_in_simu_interval = new int[path_link_seq_no_vector.size()];
			m_Veh_LinkDepartureTime_in_simu_interval = new int[path_link_seq_no_vector.size()];

			for (int i = 0; i < path_link_seq_no_vector.size(); i++)
			{
				m_Veh_LinkArrivalTime_in_simu_interval[i] = -1;
				m_Veh_LinkDepartureTime_in_simu_interval[i] = -1;

			}


			m_path_link_seq_no_vector_size = path_link_seq_no_vector.size();
			departure_time_in_simu_interval = int(departure_time_in_min * 60.0 / number_of_seconds_per_interval + 0.5);  // round off
		}
	}

	void DeallocateMemory()
	{
		if (m_Veh_LinkArrivalTime_in_simu_interval != NULL) delete m_Veh_LinkArrivalTime_in_simu_interval;
		if (m_Veh_LinkDepartureTime_in_simu_interval != NULL) delete m_Veh_LinkDepartureTime_in_simu_interval;
	}

};

vector<CAgent_Simu*> g_agent_simu_vector;
class Assignment {
public:
	Assignment()
	{
		g_number_of_memory_blocks = 4;
		total_demand_volume = 0.0;
		g_column_pool = NULL;
		g_origin_demand_array = NULL;
		//pls check following 7 settings before running programmer
		g_number_of_threads = 1;
		g_number_of_K_paths = 20;
		g_number_of_demand_periods = 24;
		g_reassignment_tau0 = 999;

		g_number_of_links = 0;
		g_number_of_service_arcs = 0;
		g_number_of_nodes = 0;
		g_number_of_zones = 0;
		g_number_of_agent_types = 0;

		b_debug_detail_flag = 1;

		g_pFileDebugLog = NULL;
		assignment_mode = 0;  // default is UE
	}

	//Init demand array[zone][agent type][demand period] and total demand array[zone][period]
	void InitializeDemandMatrix(int number_of_zones, int number_of_agent_types, int number_of_time_periods)
	{
		g_number_of_zones = number_of_zones;
		g_number_of_agent_types = number_of_agent_types;

		g_column_pool = Allocate4DDynamicArray<CColumnVector>(number_of_zones, number_of_zones, max(1, number_of_agent_types), number_of_time_periods);
		g_origin_demand_array = Allocate3DDynamicArray<float>(number_of_zones, max(1, number_of_agent_types), number_of_time_periods);

		//demand matrix has 3 demension: zone, agent, time period
		for (int i = 0; i < number_of_zones; i++)
		{
			for (int at = 0; at < number_of_agent_types; at++)
			{
				for (int tau = 0; tau < g_number_of_demand_periods; tau++)
				{

					g_origin_demand_array[i][at][tau] = 0.0;
				}
			}

		}
		total_demand_volume = 0.0;
		for (int i = 0; i < number_of_agent_types; i++)
		{
			for (int tau = 0; tau < g_number_of_demand_periods; tau++)
			{
				total_demand[i][tau] = 0.0;
			}
		}

		g_DemandGlobalMultiplier = 1.0f;
	};

	~Assignment()
	{

		if (g_column_pool != NULL)
			Deallocate4DDynamicArray(g_column_pool, g_number_of_zones, g_number_of_zones, g_number_of_agent_types);
		if (g_origin_demand_array != NULL)
			Deallocate3DDynamicArray(g_origin_demand_array, g_number_of_zones, g_number_of_agent_types);

		if (g_pFileDebugLog != NULL)
			fclose(g_pFileDebugLog);

	}
	int g_number_of_threads;
	int g_number_of_K_paths;
	int assignment_mode;
	int g_number_of_memory_blocks;

	int g_reassignment_tau0;

	int b_debug_detail_flag;
	std::map<int, int> g_internal_node_to_seq_no_map;  // hash table, map external node number to internal node sequence no. 
	std::map<int, int> g_zoneid_to_zone_seq_no_mapping;// from integer to integer map zone_id to zone_seq_no
	std::map<string, int> g_road_link_id_map;


	CColumnVector**** g_column_pool;
	float*** g_origin_demand_array;

	//StatisticOutput.cpp
	float total_demand_volume;
	//NetworkReadInput.cpp and ShortestPath.cpp



	std::vector<CDemand_Period> g_DemandPeriodVector;
	int g_LoadingStartTimeInMin;
	int g_LoadingEndTimeInMin;


	std::vector<CAgent_type> g_AgentTypeVector;
	std::map<int, CLinkType> g_LinkTypeMap;

	std::map<string, int> demand_period_to_seqno_mapping;
	std::map<string, int> agent_type_2_seqno_mapping;


	float total_demand[_MAX_AGNETTYPES][_MAX_TIMEPERIODS];
	float g_DemandGlobalMultiplier;

	int g_number_of_links;
	int g_number_of_service_arcs;
	int g_number_of_nodes;
	int g_number_of_zones;
	int g_number_of_agent_types;
	int g_number_of_demand_periods;
	int g_number_of_demand;

	FILE* g_pFileDebugLog = NULL;


	//	-------------------------------
		// used in ST Simulation
	float** m_LinkOutFlowCapacity;
	int** m_LinkTDTravelTime;

	float** m_LinkCumulativeArrival;
	float** m_LinkCumulativeDeparture;
	int g_start_simu_interval_no;

	int g_number_of_simulation_intervals;
	void AllocateLinkMemory4Simulation();

	void STTrafficSimulation();
};

Assignment assignment;

class CVDF_Period
{
public:

	CVDF_Period()
	{
		m = 0.5;
		VOC = 0;
		gamma = 3.47f;
		mu = 1000;
		alpha = 0.15f;
		beta = 4;
		rho = 1;
		marginal_base = 1;
		ruc_base_resource = 0;

		ruc_type = 0;

		starting_time_slot_no = 0;
		ending_time_slot_no = 0;
	}


	int starting_time_slot_no;  // in 15 min slot
	int ending_time_slot_no;
	string period;


	//standard BPR parameter 
	float alpha;
	float beta;
	float capacity;
	float FFTT;
	float VOC;
	float rho;
	float ruc_base_resource;
	int   ruc_type;

	float marginal_base;
	//updated BPR-X parameters
	float gamma;
	float mu;
	float m;
	float congestion_period_P;
	// inpput
	float volume;

	//output
	float avg_delay;
	float avg_travel_time = 0;
	float avg_waiting_time = 0;

	//float Q[_MAX_TIMESLOT_PerPeriod];  // t starting from starting_time_slot_no if we map back to the 24 hour horizon 
	float waiting_time[_MAX_TIMESLOT_PerPeriod];
	float arrival_rate[_MAX_TIMESLOT_PerPeriod];

	float discharge_rate[_MAX_TIMESLOT_PerPeriod];
	float travel_time[_MAX_TIMESLOT_PerPeriod];


	float get_waiting_time(int relative_time_slot_no)
	{
		if (relative_time_slot_no >= 0 && relative_time_slot_no < _MAX_TIMESLOT_PerPeriod)
			return waiting_time[relative_time_slot_no];
		else
			return 0;

	}
	int t0, t3;

	void Setup()
	{

	}

	float  PerformBPR(float volume)
	{
		volume = max(0, volume);  // take nonnegative values

		VOC = volume / max(0.00001f, capacity);
		avg_travel_time = FFTT + FFTT * alpha * pow(volume / max(0.00001f, capacity), beta);

		marginal_base = FFTT * alpha * beta * pow(volume / max(0.00001f, capacity), beta - 1);

		return avg_travel_time;

		// volume --> avg_traveltime

	}


	float PerformBPR_X(float volume)
	{
		congestion_period_P = 0;
		// Step 1: Initialization
		int L = ending_time_slot_no - starting_time_slot_no;  // in 15 min slot

		if (L >= _MAX_TIMESLOT_PerPeriod - 1)
			return 0;

		float mid_time_slot_no = starting_time_slot_no + L / 2.0;  // t1;
		for (int t = 0; t <= L; t++)
		{
			waiting_time[t] = 0;
			arrival_rate[t] = 0;
			discharge_rate[t] = mu / 2.0;
			travel_time[t] = FFTT;
		}
		avg_waiting_time = 0;
		avg_travel_time = FFTT + avg_waiting_time;

		//int L = ending_time_slot_no - starting_time_slot_no;  // in 15 min slot

		// Case 1: fully uncongested region
		if (volume <= L * mu / 2)
		{
			// still keep 0 waiting time for all time period
			congestion_period_P = 0;

		}
		else  // partially congested region
		{
			//if (volume > L * mu / 2 ) // Case 2
			float P = min(L, volume * 2 / mu - L);  // if  volume > L * mu then P is set the the maximum of L: case 3
			congestion_period_P = P / 4.0;  // unit: hour

			t0 = mid_time_slot_no - P / 2.0;
			t3 = mid_time_slot_no + P / 2.0;
			// derive t0 and t3 based on congestion duration p
			int t2 = m * (t3 - t0) + t0;
			for (int tt = 0; tt <= L; tt++)
			{
				int time = starting_time_slot_no + tt;
				if (time < t0)
				{  //first uncongested phase with mu/2 as the approximate flow rates
					waiting_time[tt] = 0;
					arrival_rate[tt] = mu / 2;
					discharge_rate[tt] = mu / 2.0;
					travel_time[tt] = FFTT;

				}
				if (time >= t0 && time <= t3)
				{
					//second congested phase
					waiting_time[tt] = 1 / (4.0 * mu) * gamma * (time - t0) * (time - t0) * (time - t3) * (time - t3);
					arrival_rate[tt] = gamma * (time - t0) * (time - t2) * (time - t3) + mu;
					discharge_rate[tt] = mu;
					travel_time[tt] = FFTT + waiting_time[tt];
				}
				if (time > t3)
				{
					//third uncongested phase with mu/2 as the approximate flow rates
					waiting_time[tt] = 0;
					arrival_rate[tt] = mu / 2;
					discharge_rate[tt] = mu / 2.0;
					travel_time[tt] = FFTT;
				}
				avg_waiting_time = gamma / (120 * mu) * pow(P, 4.0);
				//cout << avg_waiting_time << endl;
				avg_travel_time = FFTT + avg_waiting_time;
			}
		}

		return avg_travel_time;

	}
};


class CLink
{
public:
	CLink()  // construction 
	{
		zone_seq_no_for_outgoing_connector = -1;

		free_flow_travel_time_in_min = 1;
		toll = 0;
		route_choice_cost = 0;
		for (int tau = 0; tau < _MAX_TIMEPERIODS; tau++)
		{
			flow_volume_per_period[tau] = 0;
			resource_per_period[tau] = 0;

			queue_length_perslot[tau] = 0;
			travel_time_per_period[tau] = 0;

			for (int at = 0; at < _MAX_AGNETTYPES; at++)
			{
				volume_per_period_per_at[tau][at] = 0;
				resource_per_period_per_at[tau][at] = 0;

			}

			TDBaseTT[tau] = 0;
			TDBaseCap[tau] = 0;
			TDBaseFlow[tau] = 0;
			TDBaseQueue[tau] = 0;


			//cost_perhour[tau] = 0;
		}
		link_spatial_capacity = 100;
		RUC_type = 0;
	}

	~CLink()
	{
		//if (flow_volume_for_each_o != NULL)
		//	delete flow_volume_for_each_o;
	}

	void free_memory()
	{
	}

	void AddAgentsToLinkVolume()
	{


	}



	// 1. based on BPR. 

	int zone_seq_no_for_outgoing_connector;
	int m_LeftTurn_link_seq_no;

	int m_RandomSeed;
	int link_seq_no;
	string link_id;
	int from_node_seq_no;
	int to_node_seq_no;
	int link_type;
	float toll;
	float route_choice_cost;

	float PCE_at[_MAX_AGNETTYPES];
	float CRU_at[_MAX_AGNETTYPES];


	float fftt;
	float free_flow_travel_time_in_min;
	float lane_capacity;
	int number_of_lanes;
	CVDF_Period VDF_period[_MAX_TIMEPERIODS];

	float TDBaseTT[_MAX_TIMEPERIODS];
	float TDBaseCap[_MAX_TIMEPERIODS];
	float TDBaseFlow[_MAX_TIMEPERIODS];
	float TDBaseQueue[_MAX_TIMEPERIODS];


	int type;
	float link_spatial_capacity;

	//static
	//float flow_volume;
	//float travel_time;

	double flow_volume_per_period[_MAX_TIMEPERIODS];
	float volume_per_period_per_at[_MAX_TIMEPERIODS][_MAX_AGNETTYPES];

	float resource_per_period[_MAX_TIMEPERIODS];
	float resource_per_period_per_at[_MAX_TIMEPERIODS][_MAX_AGNETTYPES];

	float queue_length_perslot[_MAX_TIMEPERIODS];  // # of vehicles in the vertical point queue
	float travel_time_per_period[_MAX_TIMEPERIODS];
	float travel_marginal_cost_per_period[_MAX_TIMEPERIODS][_MAX_AGNETTYPES];

	float exterior_penalty_cost_per_period[_MAX_TIMEPERIODS][_MAX_AGNETTYPES];
	float exterior_penalty_derivative_per_period[_MAX_TIMEPERIODS][_MAX_AGNETTYPES];

	int RUC_type;
	int number_of_periods;

	float length;
	//std::vector <SLinkMOE> m_LinkMOEAry;
	//beginning of simulation data 

	//toll related link
	//int m_TollSize;
	//Toll *pTollVector;  // not using SLT here to avoid issues with OpenMP

	void CalculateTD_VDFunction();

	float get_VOC_ratio(int tau)
	{

		return (flow_volume_per_period[tau] + TDBaseFlow[tau]) / max(0.00001f, TDBaseCap[tau]);
	}

	float get_net_resource(int tau)
	{

		return resource_per_period[tau];
	}

	float get_speed(int tau)
	{
		return length / max(travel_time_per_period[tau], 0.0001f) * 60;  // per hour
	}


	void calculate_marginal_cost_for_agent_type(int tau, int agent_type_no, float PCE_agent_type)
	{
		// volume * dervative 
		// BPR_term: volume * FFTT * alpha * (beta) * power(v/c, beta-1), 

		travel_marginal_cost_per_period[tau][agent_type_no] = VDF_period[tau].marginal_base * PCE_agent_type;
	}

	void calculate_penalty_for_agent_type(int tau, int agent_type_no, float CRU_agent_type)
	{
		// volume * dervative 
		if (RUC_type == 0)
			return;

		float resource = 0;

		if (RUC_type == 2)  // equality constraints 
			resource = resource_per_period[tau];
		else
			resource = min(0, resource_per_period[tau]);  // inequality


		exterior_penalty_derivative_per_period[tau][agent_type_no] = 2 * VDF_period[tau].rho * resource * CRU_agent_type;

		if (exterior_penalty_derivative_per_period[tau][agent_type_no] < -100)
		{
			int i_debug = 1;
		}

	}

	void calculate_Gauss_Seidel_penalty_for_agent_type(int tau, int agent_type_no, float CRU_agent_type)
	{
		// volume * dervative 
		if (RUC_type == 0)
			return;

		float resource = 0;  // resource in _Gauss_Seidel framework is refereed to all the other resourcs used by others

		if (RUC_type == 2)  // equality constraints 
			resource = resource_per_period[tau];
		else
			resource = min(0, resource_per_period[tau]);  // inequality


		exterior_penalty_derivative_per_period[tau][agent_type_no] = VDF_period[tau].rho * (2 * resource * CRU_agent_type + CRU_agent_type);

	}

	float get_generalized_first_order_gradient_cost_of_second_order_loss_for_agent_type(int tau, int agent_type_no)
	{

		float generalized_cost = travel_time_per_period[tau] + toll / assignment.g_AgentTypeVector[agent_type_no].value_of_time * 60;  // *60 as 60 min per hour

		if (assignment.assignment_mode == 3 || assignment.assignment_mode == 4)  // system optimal mode or exterior panalty mode
		{
			generalized_cost += travel_marginal_cost_per_period[tau][agent_type_no];
		}

		if ((assignment.assignment_mode) == 3 && (RUC_type != 0))  // exterior panalty mode
		{
			generalized_cost += exterior_penalty_derivative_per_period[tau][agent_type_no];

		}

		return generalized_cost;
	}

	// for discrete event simulation

	std::list<int> EntranceQueue;  //link-in queue  of each link
	std::list<int> ExitQueue;      // link-out queue of each link


};

class CServiceArc
{
public:
	CServiceArc()  // construction 
	{
	}

	int link_seq_no;
	int starting_time_no;
	int ending_time_no;
	int time_interval_no;
	float capacity;
	int travel_time_delta;

};


class CNode
{
public:
	CNode()
	{
		zone_id = -1;
		node_seq_no = -1;
		//accessible_node_count = 0;
	}

	//int accessible_node_count;

	int node_seq_no;  // sequence number 
	int node_id;      //external node number 
	int zone_id = -1;

	double x;
	double y;

	std::vector<int> m_outgoing_link_seq_no_vector;
	std::vector<int> m_incoming_link_seq_no_vector;

	std::vector<int> m_to_node_seq_no_vector;
	std::map<int, int> m_to_node_2_link_seq_no_map;

};


extern std::vector<CNode> g_node_vector;
extern std::vector<CLink> g_link_vector;
extern std::vector<CServiceArc> g_service_arc_vector;

class COZone
{
public:
	int zone_seq_no;  // 0, 1, 
	int zone_id;  // external zone id // this is origin zone
	int node_seq_no;


};

extern std::vector<COZone> g_zone_vector;
extern std::map<int, int> g_zoneid_to_zone_seq_no_mapping;

class CAGBMAgent
{
public:

	int agent_id;
	int income;
	int gender;
	int vehicle;
	int purpose;
	int flexibility;
	float preferred_arrival_time;
	float travel_time_in_min;
	float free_flow_travel_time;
	int from_zone_seq_no;
	int to_zone_seq_no;
	int type;
	int time_period;
	int k_path;
	float volume;
	float arrival_time_in_min;
};
extern std::vector<CAGBMAgent> g_agbmagent_vector;

struct CNodeForwardStar {
	int* OutgoingLinkNoArray;
	int* OutgoingNodeNoArray;
	int OutgoingLinkSize;
};





std::vector<CNode> g_node_vector;
std::vector<CLink> g_link_vector;
std::vector<CServiceArc> g_service_arc_vector;

std::vector<COZone> g_zone_vector;
std::vector<CAGBMAgent> g_agbmagent_vector;


class VehicleScheduleNetworks {

public:
	int m_agent_type_no;
	int m_time_interval_size;  // 1440

	std::vector<CNode> m_node_vector;  // local copy of node vector, based on agent type and origin node
	std::map<int, float> g_passenger_link_profit;  // link with negative link

	void BuildNetwork(Assignment& assignment, int tau, int iteration)
	{

		for (int i = 0; i < assignment.g_number_of_nodes; i++) //Initialization for all non-origin nodes
		{
			CNode node;  // create a node object 

			node.node_id = g_node_vector[i].node_id;
			node.node_seq_no = g_node_vector[i].node_seq_no;

			for (int j = 0; j < g_node_vector[i].m_outgoing_link_seq_no_vector.size(); j++)
			{

				int link_seq_no = g_node_vector[i].m_outgoing_link_seq_no_vector[j];

				if (assignment.g_LinkTypeMap[g_link_vector[link_seq_no].link_type].AllowAgentType(assignment.g_AgentTypeVector[m_agent_type_no].agent_type))  // only predefined allowed agent type can be considered
				{
					int from_node_seq_no = g_link_vector[link_seq_no].from_node_seq_no;
					node.m_outgoing_link_seq_no_vector.push_back(link_seq_no);

					g_passenger_link_profit[link_seq_no] = g_link_vector[link_seq_no].get_generalized_first_order_gradient_cost_of_second_order_loss_for_agent_type(tau, m_agent_type_no);

					if (g_debugging_flag && assignment.g_pFileDebugLog != NULL)
						fprintf(assignment.g_pFileDebugLog, "DP iteration %d: link %d->%d:  profit %.3f\n",
							iteration,
							g_node_vector[g_link_vector[link_seq_no].from_node_seq_no].node_id,
							g_node_vector[g_link_vector[link_seq_no].to_node_seq_no].node_id,
							g_passenger_link_profit[link_seq_no]);

					node.m_to_node_seq_no_vector.push_back(g_node_vector[i].m_to_node_seq_no_vector[j]);
				}

			}

			m_node_vector.push_back(node);

		}
	}


	//class for vehicle scheduling states
	class CVSState
	{
	public:
		int current_node_no;  // space dimension

		std::map<int, int> passenger_service_state;  // passenger means link with negative costs

		std::vector<int> m_visit_link_sequence;  // store link sequence

		int m_vehicle_capacity;

		float LabelCost;  // with LR price
		float LabelTime;   // sum of travel time up to now, arrival time at node


		CVSState()
		{
			current_node_no = 0;
			LabelTime = 0;
			LabelCost = 0;
			m_vehicle_capacity = 0;
		}

		void Copy(CVSState* pSource)
		{
			current_node_no = pSource->current_node_no;
			passenger_service_state.clear();
			passenger_service_state = pSource->passenger_service_state;


			m_visit_link_sequence = pSource->m_visit_link_sequence;
			m_vehicle_capacity = pSource->m_vehicle_capacity;
			LabelCost = pSource->LabelCost;
			LabelTime = pSource->LabelTime;
		}
		int GetPassengerLinkServiceState(int link_no)
		{
			if (passenger_service_state.find(link_no) != passenger_service_state.end())
				return passenger_service_state[link_no];  // 1 or 2
			else
				return 0;
		}



		std::string generate_string_key()
		{

			stringstream s;

			s << "n";
			s << current_node_no;  // space key
			for (std::map<int, int>::iterator it = passenger_service_state.begin(); it != passenger_service_state.end(); ++it)
			{
				s << "_";

				s << it->first << "[" << it->first << "]";

			}
			string converted(s.str());
			return converted;

		}

		bool operator<(const CVSState& other) const
		{
			return LabelCost < other.LabelCost;
		}

	};

	class C_time_indexed_state_vector
	{
	public:
		int current_time;


		std::vector<CVSState> m_VSStateVector;

		std::map<std::string, int> m_state_map;

		void Reset()
		{
			current_time = 0;
			m_VSStateVector.clear();
			m_state_map.clear();
		}

		int m_find_state_index(std::string string_key)
		{

			if (m_state_map.find(string_key) != m_state_map.end())
			{
				return m_state_map[string_key];
			}
			else
				return -1;  // not found

		}

		void update_state(CVSState new_element)
		{
			std::string string_key = new_element.generate_string_key();//if it is new, string is n100, no state index
			int state_index = m_find_state_index(string_key);

			if (state_index == -1)  // no such state at this time
			{
				// add new state
				state_index = m_VSStateVector.size();
				m_VSStateVector.push_back(new_element);
				m_state_map[string_key] = state_index;
			}
			else
			{//DP
				if (new_element.LabelCost < m_VSStateVector[state_index].LabelCost)
				{
					m_VSStateVector[state_index].Copy(&new_element);
				}

			}

		}

		void Sort()
		{
			std::sort(m_VSStateVector.begin(), m_VSStateVector.end());

			m_state_map.clear(); // invalid
		}

		void SortAndCleanEndingState(int BestKValue)
		{
			if (m_VSStateVector.size() > 2 * BestKValue)
			{
				std::sort(m_VSStateVector.begin(), m_VSStateVector.end());

				m_state_map.clear(); // invalid
				m_VSStateVector.erase(m_VSStateVector.begin() + BestKValue, m_VSStateVector.end());
			}
		}

		float GetBestValue()
		{
			// LabelCost not PrimalCost when sorting
			std::sort(m_VSStateVector.begin(), m_VSStateVector.end());

			if (m_VSStateVector.size() >= 1)
			{
				std::string state_str = m_VSStateVector[0].generate_string_key();

				return m_VSStateVector[0].LabelCost;

			}
			else
				return _MAX_LABEL_COST;
		}

		std::vector<int> GetBestLinkSequence()
		{
			std::vector <int> link_sequence;
			// LabelCost not PrimalCost when sorting
			std::sort(m_VSStateVector.begin(), m_VSStateVector.end());

			if (m_VSStateVector.size() >= 1)
			{
				std::string state_str = m_VSStateVector[0].generate_string_key();

				if (m_VSStateVector[0].m_visit_link_sequence.size() > 0)
					return m_VSStateVector[0].m_visit_link_sequence;
				else
					return link_sequence;

			}

			return link_sequence;
		}

	};

	//vehicle state at time t

	// for collecting the final feasible states accesible to the depot
	C_time_indexed_state_vector g_ending_state_vector;

	C_time_indexed_state_vector** g_time_dependent_state_vector;  // label cost vector [i,t,w]


	void AllocateVSNMemory(int number_of_nodes)
	{
		g_time_dependent_state_vector = AllocateDynamicArray <C_time_indexed_state_vector>(number_of_nodes, m_time_interval_size);  //1
	}

	~VehicleScheduleNetworks()
	{
		DeallocateDynamicArray(g_time_dependent_state_vector, g_node_vector.size(), m_time_interval_size);
	}


	std::vector<int> g_optimal_time_dependenet_dynamic_programming(
		int origin_node,
		int destination_node,
		int vehicle_capacity,
		int demand_time_period_no,
		//maximum choose
		int BestKSize)
		// time-dependent label correcting algorithm with double queue implementation
	{
		int arrival_time = m_time_interval_size - 1;  // restrict the search range.

		int t = 0;
		//step 2: Initialization for origin node at the preferred departure time, at departure time
		for (int i = 0; i < assignment.g_number_of_nodes; i++)
		{
			g_time_dependent_state_vector[i][0].Reset();

		}
		g_ending_state_vector.Reset();

		CVSState element;

		element.current_node_no = origin_node;
		g_time_dependent_state_vector[origin_node][0].update_state(element);


		// step 3: //dynamic programming
		for (t = 0; t < arrival_time; t++)  //first loop: time
		{

			for (int n = 0; n < m_node_vector.size(); n++)
			{
				// step 1: sort m_VSStateVector by labelCost for scan best k elements in step2
				g_time_dependent_state_vector[n][t].Sort();

				// step 2: scan the best k elements
				for (int w_index = 0; w_index < min(BestKSize, g_time_dependent_state_vector[n][t].m_VSStateVector.size()); w_index++)
				{
					CVSState* pElement = &(g_time_dependent_state_vector[n][t].m_VSStateVector[w_index]);

					int from_node = pElement->current_node_no;

					// step 2.1 link from node to toNode
					for (int i = 0; i < m_node_vector[from_node].m_outgoing_link_seq_no_vector.size(); i++)
					{
						int link_seq_no = m_node_vector[from_node].m_outgoing_link_seq_no_vector[i];
						int to_node = g_link_vector[link_seq_no].to_node_seq_no;

						float new_time = pElement->LabelTime + g_link_vector[link_seq_no].travel_time_per_period[demand_time_period_no];
						int new_time_int = max(pElement->LabelTime + 1, (int)(new_time + 0.5));  // move at least one time step further

						// step 2.2. check feasibility of node type with the current element
						if (new_time <= arrival_time)
						{

							// skip scanning when the origin/destination nodes arrival time is out of time window
							//feasible state transitions
								// do not need waiting
							CVSState new_element;
							new_element.Copy(pElement);

							new_element.current_node_no = to_node;

							new_element.LabelCost += g_link_vector[link_seq_no].travel_time_per_period[demand_time_period_no] / 60.0 * assignment.g_AgentTypeVector[m_agent_type_no].value_of_time; // 60.0 is to convert hour to 60 min as VOT is denoted as dollars per hour
							if (g_passenger_link_profit.find(link_seq_no) != g_passenger_link_profit.end())
							{
								new_element.LabelCost += g_passenger_link_profit[link_seq_no];// + negative cost
								new_element.passenger_service_state[link_seq_no] = 1;  // mark carry status
								new_element.m_vehicle_capacity -= 1;

							}


							new_element.m_visit_link_sequence.push_back(link_seq_no);
							g_time_dependent_state_vector[to_node][new_time_int].update_state(new_element);

							if (to_node == destination_node)
							{

								//time window of destination_node
								if (new_time < arrival_time)
								{
									g_ending_state_vector.update_state(new_element);
									g_ending_state_vector.SortAndCleanEndingState(BestKSize);
								}
							}
						}

					}
				}
			}  // for all nodes
		} // for all time t


	// no backf
		return g_ending_state_vector.GetBestLinkSequence();

	}

};


int g_read_integer(FILE* f, bool speicial_char_handling)
// read an integer from the current pointer of the file, skip all spaces
{
	char ch, buf[32];
	int i = 0;
	int flag = 1;
	/* returns -1 if end of file is reached */

	while (true)
	{
		ch = getc(f);
		//cout << "get from node successful: " << ch;
		if (ch == EOF || (speicial_char_handling && (ch == '*' || ch == '$')))
			return -1; // * and $ are special characters for comments
		if (isdigit(ch))
			break;
		if (ch == '-')
			flag = -1;
		else
			flag = 1;
	};
	if (ch == EOF) return -1;


	while (isdigit(ch)) {
		buf[i++] = ch;
		//cout << "isdigit" << buf[i++] << endl;
		ch = fgetc(f);
		//cout << "new ch" << ch;
	}
	buf[i] = 0;


	return atoi(buf) * flag;

}


float g_read_float(FILE* f)
//read a floating point number from the current pointer of the file,
//skip all spaces

{
	char ch, buf[32];
	int i = 0;
	int flag = 1;

	/* returns -1 if end of file is reached */

	while (true)
	{
		ch = getc(f);
		if (ch == EOF || ch == '*' || ch == '$') return -1;
		if (isdigit(ch))
			break;

		if (ch == '-')
			flag = -1;
		else
			flag = 1;

	};
	if (ch == EOF) return -1;
	while (isdigit(ch) || ch == '.') {
		buf[i++] = ch;
		ch = fgetc(f);

	}
	buf[i] = 0;

	/* atof function converts a character string (char *) into a doubleing
	pointer equivalent, and if the string is not a floting point number,
	a zero will be return.
	*/

	return (float)(atof(buf) * flag);

}



//split the string by "_"
vector<string> split(const string& s, const string& seperator) {
	vector<string> result;
	typedef string::size_type string_size;
	string_size i = 0;

	while (i != s.size()) {
		int flag = 0;
		while (i != s.size() && flag == 0) {
			flag = 1;
			for (string_size x = 0; x < seperator.size(); ++x)
				if (s[i] == seperator[x]) {
					++i;
					flag = 0;
					break;
				}
		}

		flag = 0;
		string_size j = i;
		while (j != s.size() && flag == 0) {
			for (string_size x = 0; x < seperator.size(); ++x)
				if (s[j] == seperator[x]) {
					flag = 1;
					break;
				}
			if (flag == 0)
				++j;
		}
		if (i != j) {
			result.push_back(s.substr(i, j - i));
			i = j;
		}
	}
	return result;
}
string test_str = "0300:30:120_0600:30:140";

//g_global_minute = g_time_parser(test_str);
//
//for (int i = 0; i < g_global_minute.size(); i++)
//{
//	cout << "The number of global minutes is: " << g_global_minute[i] << " minutes" << endl;
//}

//vector<float> g_time_parser(string str)
//{
//	vector<float> output_global_minute;
//
//	int string_lenghth = str.length();
//
//	ASSERT(string_lenghth < 100);
//
//	const char* string_line = str.data(); //string to char*
//
//	int char_length = strlen(string_line);
//
//	char ch, buf_ddhhmm[32] = { 0 }, buf_SS[32] = { 0 }, buf_sss[32] = { 0 };
//	char dd1, dd2, hh1, hh2, mm1, mm2, SS1, SS2, sss1, sss2, sss3;
//	float ddf1, ddf2, hhf1, hhf2, mmf1, mmf2, SSf1, SSf2, sssf1, sssf2, sssf3;
//	float global_minute = 0;
//	float dd = 0, hh = 0, mm = 0, SS = 0, sss = 0;
//	int i = 0;
//	int buffer_i = 0, buffer_k = 0, buffer_j = 0;
//	int num_of_colons = 0;
//
//	//DDHHMM:SS:sss or HHMM:SS:sss
//
//	while (i < char_length)
//	{
//		ch = string_line[i++];
//
//		if (num_of_colons == 0 && ch != '_' && ch != ':') //input to buf_ddhhmm until we meet the colon
//		{
//			buf_ddhhmm[buffer_i++] = ch;
//		}
//		else if (num_of_colons == 1 && ch != ':') //start the Second "SS"
//		{
//			buf_SS[buffer_k++] = ch;
//		}
//		else if (num_of_colons == 2 && ch != ':') //start the Millisecond "sss"
//		{
//			buf_sss[buffer_j++] = ch;
//		}
//
//		if (ch == '_' || i == char_length) //start a new time string
//		{
//			if (buffer_i == 4) //"HHMM"
//			{
//				//HHMM, 0123
//				hh1 = buf_ddhhmm[0]; //read each first
//				hh2 = buf_ddhhmm[1];
//				mm1 = buf_ddhhmm[2];
//				mm2 = buf_ddhhmm[3];
//
//				hhf1 = ((float)hh1 - 48); //convert a char to a float
//				hhf2 = ((float)hh2 - 48);
//				mmf1 = ((float)mm1 - 48);
//				mmf2 = ((float)mm2 - 48);
//
//				dd = 0;
//				hh = hhf1 * 10 * 60 + hhf2 * 60;
//				mm = mmf1 * 10 + mmf2;
//			}
//			else if (buffer_i == 6) //"DDHHMM"
//			{
//				//DDHHMM, 012345
//				dd1 = buf_ddhhmm[0]; //read each first
//				dd2 = buf_ddhhmm[1];
//				hh1 = buf_ddhhmm[2];
//				hh2 = buf_ddhhmm[3];
//				mm1 = buf_ddhhmm[4];
//				mm2 = buf_ddhhmm[5];
//
//				ddf1 = ((float)dd1 - 48); //convert a char to a float
//				ddf2 = ((float)dd2 - 48);
//				hhf1 = ((float)hh1 - 48);
//				hhf2 = ((float)hh2 - 48);
//				mmf1 = ((float)mm1 - 48);
//				mmf2 = ((float)mm2 - 48);
//
//				dd = ddf1 * 10 * 24 * 60 + ddf2 * 24 * 60;
//				hh = hhf1 * 10 * 60 + hhf2 * 60;
//				mm = mmf1 * 10 + mmf2;
//			}
//
//			if (num_of_colons == 1 || num_of_colons == 2)
//			{
//				//SS, 01
//				SS1 = buf_SS[0]; //read each first
//				SS2 = buf_SS[1];
//
//				SSf1 = ((float)SS1 - 48); //convert a char to a float
//				SSf2 = ((float)SS2 - 48);
//
//				SS = (SSf1 * 10 + SSf2) / 60;
//			}
//
//			if (num_of_colons == 2)
//			{
//				//sss, 012
//				sss1 = buf_sss[0]; //read each first
//				sss2 = buf_sss[1];
//				sss3 = buf_sss[2];
//
//				sssf1 = ((float)sss1 - 48); //convert a char to a float
//				sssf2 = ((float)sss2 - 48);
//				sssf3 = ((float)sss3 - 48);
//
//				sss = (sssf1 * 100 + sssf2 * 10 + sssf3) / 1000;
//			}
//
//			global_minute = dd + hh + mm + SS + sss;
//
//			output_global_minute.push_back(global_minute);
//
//			//initialize the parameters
//			buffer_i = 0;
//			buffer_k = 0;
//			buffer_j = 0;
//			num_of_colons = 0;
//		}
//
//		if (ch == ':')
//		{
//			num_of_colons += 1;
//		}
//	}
//
//	return output_global_minute;
//}
//

vector<float> g_time_parser(vector<string>& inputstring)
{
	vector<float> output_global_minute;

	for (int k = 0; k < inputstring.size(); k++)
	{
		vector<string> sub_string = split(inputstring[k], "_");

		for (int i = 0; i < sub_string.size(); i++)
		{
			//HHMM
			//012345
			char hh1 = sub_string[i].at(0);
			char hh2 = sub_string[i].at(1);
			char mm1 = sub_string[i].at(2);
			char mm2 = sub_string[i].at(3);

			float hhf1 = ((float)hh1 - 48);
			float hhf2 = ((float)hh2 - 48);
			float mmf1 = ((float)mm1 - 48);
			float mmf2 = ((float)mm2 - 48);

			float hh = hhf1 * 10 * 60 + hhf2 * 60;
			float mm = mmf1 * 10 + mmf2;
			float global_mm_temp = hh + mm;
			output_global_minute.push_back(global_mm_temp);
		}
	}

	return output_global_minute;
} // transform hhmm to minutes 


inline string g_time_coding(float time_stamp)
{
	int hour = time_stamp / 60;
	int minute = time_stamp - hour * 60;

	int second = (time_stamp - hour * 60 - minute) * 60;

	ostringstream strm;
	strm.fill('0');
	strm << setw(2) << hour << setw(2) << minute << ":" << setw(2) << second;

	return strm.str();
} // transform hhmm to minutes 


void g_ProgramStop()
{

	cout << "STALite Program stops. Press any key to terminate. Thanks!" << endl;
	getchar();
	exit(0);
};



//void ReadLinkTollScenarioFile(Assignment& assignment)
//{
//
//	for (unsigned li = 0; li < g_link_vector.size(); li++)
//	{
//
//		g_link_vector[li].TollMAP.erase(g_link_vector[li].TollMAP.begin(), g_link_vector[li].TollMAP.end()); // remove all previouly read records
//	}
//
//	// generate toll based on demand type code in input_link.csv file
//	int demand_flow_type_count = 0;
//
//	for (unsigned li = 0; li < g_link_vector.size(); li++)
//	{
//		if (g_link_vector[li].agent_type_code.size() >= 1)
//		{  // with data string
//
//			std::string agent_type_code = g_link_vector[li].agent_type_code;
//
//			vector<float> TollRate;
//			for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
//			{
//				CString number;
//				number.Format(_T("%d"), at);
//
//				std::string str_number = CString2StdString(number);
//				if (agent_type_code.find(str_number) == std::string::npos)   // do not find this number
//				{
//					g_link_vector[li].TollMAP[at] = 999;
//					demand_flow_type_count++;
//				}
//				else
//				{
//					g_link_vector[li].TollMAP[at] = 0;
//				}
//
//			}  //end of pt
//		}
//	}
//}



int g_ParserIntSequence(std::string str, std::vector<int>& vect)
{

	std::stringstream ss(str);

	int i;

	while (ss >> i)
	{
		vect.push_back(i);

		if (ss.peek() == ';')
			ss.ignore();
	}

	return vect.size();
}



void g_ReadDemandFileBasedOnDemandFileList(Assignment& assignment)
{
	assignment.InitializeDemandMatrix(g_zone_vector.size(), assignment.g_AgentTypeVector.size(), assignment.g_DemandPeriodVector.size());

	float total_demand_in_demand_file = 0;

	CCSVParser parser;
	cout << "Step 4: Reading file demand_file_list.csv..." << endl;

	//Read demand_file_list.csv: one-line record implies a demand type: passenger or vehicle, etc.
	if (parser.OpenCSVFile("demand_file_list.csv", true))
	{
		int i = 0;
		while (parser.ReadRecord())
		{
			int file_sequence_no = 1;
			string file_name;
			string format_type = "null";

			string demand_period, agent_type;

			int demand_format_flag = 0;

			//Ship the uesless record
			if (parser.GetValueByFieldName("file_sequence_no", file_sequence_no) == false)
				break;
			if (file_sequence_no <= -1)  // skip negative sequence no 
				continue;

			parser.GetValueByFieldName("file_name", file_name);

			parser.GetValueByFieldName("demand_period", demand_period);

			parser.GetValueByFieldName("format_type", format_type);
			if (format_type.find("null") != string::npos)  // skip negative sequence no 
			{
				cout << "Please provide format_type in file demand_file_list.csv" << endl;
				g_ProgramStop();
			}


			double total_ratio = 0;

			parser.GetValueByFieldName("agent_type", agent_type);


			int agent_type_no = 0;
			int demand_period_no = 0;

			//if the demand_type(eg. passenger) is defined in demand _file_list but not find in demand.csv
			if (assignment.demand_period_to_seqno_mapping.find(demand_period) != assignment.demand_period_to_seqno_mapping.end())
			{
				demand_period_no = assignment.demand_period_to_seqno_mapping[demand_period];
			}
			else
			{
				cout << "Error: demand period in demand_file_list " << demand_period << "cannot be found." << endl;
				g_ProgramStop();

			}

			bool b_multi_agent_list = false;

			if (agent_type == "multi_agent_list")
			{
				b_multi_agent_list = true;
			}
			else
			{

				if (assignment.agent_type_2_seqno_mapping.find(agent_type) != assignment.agent_type_2_seqno_mapping.end())
				{
					agent_type_no = assignment.agent_type_2_seqno_mapping[agent_type];
				}
				else
				{
					cout << "Error: agent_type in agent_type " << agent_type << "cannot be found." << endl;
					//"Cannot be found" => cannot be found in map.
					g_ProgramStop();
				}
			}

			//demand_period_no must <= 1
			//Question: I don't understand this if-else if- else =��=
			if (demand_period_no > _MAX_TIMEPERIODS)
			{
				cout << "demand_period_no should be less than settings in demand_period.csv. Please change the parameter settings in the source code." << endl;
				g_ProgramStop();
			}

			if (format_type.find("column") != string::npos)  // or muliti-column; string::npos: find the sub-string
				//Question: What means "coulmn"???
			{
				bool bFileReady = false;

				FILE* st;
				// read the file formaly after the test. 

				int error_count = 0;
				fopen_ss(&st, file_name.c_str(), "r");
				if (st != NULL)
				{
					bFileReady = true;
					int line_no = 0;

					while (true)
					{
						int origin_zone = (int)(g_read_float(st));
						int destination_zone = (int)g_read_float(st);
						float demand_value = g_read_float(st);

						if (origin_zone <= 0)
						{

							if (line_no == 1 && !feof(st))  // read only one line, but has not reached the end of the line
							{
								cout << endl << "Error: Only one line has been read from file. Are there multiple columns of demand type in file " << file_name << " per line?" << endl;
								g_ProgramStop();

							}
							break;
						}

						if (assignment.g_zoneid_to_zone_seq_no_mapping.find(origin_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
						{
							if (error_count < 10)
								cout << endl << "Warning: origin zone " << origin_zone << "  has not been defined in node.csv" << endl;

							error_count++;
							continue; // origin zone  has not been defined, skipped. 
						}


						if (assignment.g_zoneid_to_zone_seq_no_mapping.find(destination_zone) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
						{
							if (error_count < 10)
								cout << endl << "Warning: destination zone " << destination_zone << "  has not been defined in node.csv" << endl;

							error_count++;
							continue; // destination zone  has not been defined, skipped. 
						}

						int from_zone_seq_no = 0;
						int to_zone_seq_no = 0;
						from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[origin_zone];
						to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[destination_zone];

						if (assignment.g_zoneid_to_zone_seq_no_mapping.size() >= 395)
						{
							int i_debug = 1;
						}
						if (demand_value < -99) // encounter return 
						{
							break;
						}

						if (from_zone_seq_no == 0 && to_zone_seq_no == 3)
						{
							int ibebug = 1;
						}

						assignment.total_demand[agent_type_no][demand_period_no] += demand_value;
						assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume += demand_value;
						assignment.total_demand_volume += demand_value;
						assignment.g_origin_demand_array[from_zone_seq_no][agent_type_no][demand_period_no] += demand_value;

						// we generate vehicles here for each OD data line
						if (line_no <= 5)  // read only one line, but has not reached the end of the line
							cout << "o_zone_id:" << origin_zone << ", d_zone_id: " << destination_zone << ", value = " << demand_value << endl;


						line_no++;
					}  // scan lines
					assignment.g_number_of_demand = line_no;

					fclose(st);

					cout << "total_demand_volume is " << assignment.total_demand_volume << endl;
				}
				else  //open file
				{
					cout << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
					g_ProgramStop();

				}
			}

			else if (format_type.compare("agent_csv") == 0)
			{

				CCSVParser parser;

				if (parser.OpenCSVFile(file_name, false))
				{
					int total_demand_in_demand_file = 0;


					// read agent file line by line,

					int agent_id, o_zone_id, d_zone_id;
					string agent_type, demand_period;

					std::vector <int> node_sequence;

					while (parser.ReadRecord())
					{
						total_demand_in_demand_file++;

						if (total_demand_in_demand_file % 1000 == 0)
							cout << "demand_volume is " << total_demand_in_demand_file << endl;

						parser.GetValueByFieldName("agent_id", agent_id);

						parser.GetValueByFieldName("o_zone_id", o_zone_id);
						parser.GetValueByFieldName("d_zone_id", d_zone_id);

						CAgentPath agent_path_element;

						int o_node_id;
						int d_node_id;

						parser.GetValueByFieldName("path_id", agent_path_element.path_id);
						parser.GetValueByFieldName("o_node_id", o_node_id);
						parser.GetValueByFieldName("d_node_id", d_node_id);
						parser.GetValueByFieldName("volume", agent_path_element.volume);

						agent_path_element.o_node_no = assignment.g_internal_node_to_seq_no_map[o_node_id];
						agent_path_element.d_node_no = assignment.g_internal_node_to_seq_no_map[d_node_id];


						int from_zone_seq_no = 0;
						int to_zone_seq_no = 0;
						from_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[o_zone_id];
						to_zone_seq_no = assignment.g_zoneid_to_zone_seq_no_mapping[d_zone_id];


						assignment.total_demand[agent_type_no][demand_period_no] += agent_path_element.volume;
						assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume += agent_path_element.volume;
						assignment.total_demand_volume += agent_path_element.volume;
						assignment.g_origin_demand_array[from_zone_seq_no][agent_type_no][demand_period_no] += agent_path_element.volume;


						if (assignment.g_AgentTypeVector[agent_type_no].flow_type == 1)  // fixed path
						{
							bool bValid = true;


							std::string path_node_sequence;
							parser.GetValueByFieldName("node_sequence", path_node_sequence);

							std::vector<int> node_id_sequence;

							g_ParserIntSequence(path_node_sequence, node_id_sequence);

							std::vector<int> node_no_sequence;
							std::vector<int> link_no_sequence;


							int node_sum = 0;
							for (int i = 0; i < node_id_sequence.size(); i++)
							{

								if (assignment.g_internal_node_to_seq_no_map.find(node_id_sequence[i]) == assignment.g_internal_node_to_seq_no_map.end())
								{
									bValid = false;
									continue; //has not been defined

									// warning
								}

								int internal_node_seq_no = assignment.g_internal_node_to_seq_no_map[node_id_sequence[i]];  // map external node number to internal node seq no. 
								node_no_sequence.push_back(internal_node_seq_no);

								node_sum += internal_node_seq_no;
								if (i >= 1)
								{ // check if a link exists

									int link_seq_no = -1;
									int prev_node_seq_no = assignment.g_internal_node_to_seq_no_map[node_id_sequence[i - 1]];  // map external node number to internal node seq no. 

									int current_node_no = node_no_sequence[i];
									if (g_node_vector[prev_node_seq_no].m_to_node_2_link_seq_no_map.find(current_node_no) != g_node_vector[prev_node_seq_no].m_to_node_2_link_seq_no_map.end())
									{
										link_seq_no = g_node_vector[prev_node_seq_no].m_to_node_2_link_seq_no_map[node_no_sequence[i]];

										link_no_sequence.push_back(link_seq_no);
									}
									else
									{
										bValid = false;
									}


								}

							}


							if (bValid == true)
							{
								agent_path_element.node_sum = node_sum; // pointer to the node sum based path node sequence;
								agent_path_element.path_link_sequence = link_no_sequence;
							}

						}

						assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].discrete_agent_path_vector.push_back(agent_path_element);



					}



				}
				else  //open file
				{
					cout << "Error: File " << file_name << " cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
					g_ProgramStop();

				}

			}

			else
			{
				cout << "Error: format_type = " << format_type << " is not supported. Currently STALite supports format such as column and agent_csv." << endl;
				g_ProgramStop();
			}
		}

	}

}
//Finished g_ReadDemandFileBasedOnDemandFileList 7/22 9:15


//All parameters will be saved into assignment. class
void g_ReadInputData(Assignment& assignment)
{
	//Question: Why init these parameter
	assignment.g_LoadingStartTimeInMin = 99999;
	assignment.g_LoadingEndTimeInMin = 0;


	//step 0: Read "demand period.csv"
	CCSVParser parser_demand_period;
	cout << "Step 1: Reading file demand_period.csv..." << endl;
	if (!parser_demand_period.OpenCSVFile("demand_period.csv", true))
	{
		cout << "demand_period.csv cannot be opened. " << endl;
		g_ProgramStop();
	}

	if (parser_demand_period.inFile.is_open() || parser_demand_period.OpenCSVFile("demand_period.csv", true))
		//Open successfully
	{
		//Read each line
		while (parser_demand_period.ReadRecord())
		{
			CDemand_Period demand_period;

			if (parser_demand_period.GetValueByFieldName("demand_period_id", demand_period.demand_period_id) == false)
			{
				cout << "Error: Field demand_period_id in file demand_period cannot be read." << endl;
				g_ProgramStop();
				break;
			}
			if (parser_demand_period.GetValueByFieldName("demand_period", demand_period.demand_period) == false)
			{
				cout << "Error: Field demand_period in file demand_period cannot be read." << endl;
				g_ProgramStop();
				break;
			}


			vector<float> global_minute_vector; //used to save time_period
			//Format: Time_period:0700_0800;
			if (parser_demand_period.GetValueByFieldName("time_period", demand_period.time_period) == false)
			{
				cout << "Error: Field time_period in file demand_period cannot be read." << endl;
				g_ProgramStop();
				break;
			}

			vector<string> input_string;
			input_string.push_back(demand_period.time_period);
			//input_string includes the start and end time of a time period with hh_mm format

			global_minute_vector = g_time_parser(input_string);
			//global_minute_vector incldue the starting and ending time
			//time_period:<string>, eg.0700_0800; global_minute_vector:<float>, eg.420,480

			if (global_minute_vector.size() == 2)
			{
				demand_period.starting_time_slot_no = global_minute_vector[0] / MIN_PER_TIMESLOT;
				demand_period.ending_time_slot_no = global_minute_vector[1] / MIN_PER_TIMESLOT;
				//Question: MIN_PER_TIMESLOT == 15?

				if (global_minute_vector[0] < assignment.g_LoadingStartTimeInMin)
					assignment.g_LoadingStartTimeInMin = global_minute_vector[0]; //Start time
				if (global_minute_vector[1] > assignment.g_LoadingEndTimeInMin)
					assignment.g_LoadingEndTimeInMin = global_minute_vector[1]; //End time
				//Please check init-parameter at the beginning of this function

			}

			assignment.demand_period_to_seqno_mapping[demand_period.demand_period] = assignment.g_DemandPeriodVector.size();
			assignment.g_DemandPeriodVector.push_back(demand_period);
			//Please debug the variable demand_period.

		}
		parser_demand_period.CloseCSVFile();

		//Error: File demand_period.csv has no information
		if (assignment.g_DemandPeriodVector.size() == 0)
		{
			cout << "Error:  File demand_period.csv has no information." << endl;
			g_ProgramStop();
		}
	}
	else
		// Open unsuccessfully
	{
		cout << "Error: File demand_period.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
		g_ProgramStop();
	}
	assignment.g_number_of_demand_periods = assignment.g_DemandPeriodVector.size(); //check the size



	//step 1:read demand type file
	cout << "Reading file link_type.csv..." << endl;
	CCSVParser parser_link_type;

	if (parser_link_type.OpenCSVFile("link_type.csv", true))
	{

		int line_no = 0;
		while (parser_link_type.ReadRecord())
		{
			CLinkType element;

			//Link_type(necessary): Highway/Expressway/Arterial
			if (parser_link_type.GetValueByFieldName("link_type", element.link_type) == false)
			{
				if (line_no == 0)
				{
					cout << "Error: Field link_type cannot be found in file link_type.csv." << endl;
					g_ProgramStop();
				}
				else
				{  // read empty line
					break;
				}
			}

			//Question: Why find?
			if (assignment.g_LinkTypeMap.find(element.link_type) != assignment.g_LinkTypeMap.end())
			{
				cout << "Error: Field link_type " << element.link_type << " has been defined more than once in file link_type.csv." << endl;
				g_ProgramStop();
				break;
			}

			//Question: What means? type_code("f" or "a") and agent_type_list?
			parser_link_type.GetValueByFieldName("type_code", element.type_code);
			parser_link_type.GetValueByFieldName("agent_type_list", element.agent_type_list);

			assignment.g_LinkTypeMap[element.link_type] = element;
			// one element: one-line record
			line_no++;
		}
	}
	else
		//Open unsuccessfully
	{
		cout << "Error: File link_type.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
	}



	//Step 2: Read "agent_type.csv"
	//"link type" and "agent type" is only a type file, we should do matching work between the agent file and agent_type file
	CCSVParser parser_agent_type;
	cout << "Step 2: Reading file agent_type.csv..." << endl;

	if (!parser_agent_type.OpenCSVFile("agent_type.csv", true))
	{
		cout << "agent_type.csv cannot be opened. " << endl;
		g_ProgramStop();

	}

	if (parser_agent_type.inFile.is_open() || parser_agent_type.OpenCSVFile("agent_type.csv", true))
	{
		assignment.g_AgentTypeVector.clear();
		while (parser_agent_type.ReadRecord())
		{
			CAgent_type agent_type;
			agent_type.agent_type_no = assignment.g_AgentTypeVector.size();

			if (parser_agent_type.GetValueByFieldName("agent_type", agent_type.agent_type) == false)
				//"agent_type":'p'(abbreviation of passenger)
			{
				break;
			}
			parser_agent_type.GetValueByFieldName("VOT", agent_type.value_of_time);
			//"VOT": value of time
			//Question: what means "value of time"?

			parser_agent_type.GetValueByFieldName("flow_type", agent_type.flow_type);
			//"Flow_type": 0: flow, 1: fixed path, 2: integer decision variables.
			//Question: What means: "flow type"?

			float value;
			std::map<int, CLinkType>::iterator it; //map between <int> and <Clinktype>

			// scan through the map with different node sum for different paths
			// For each agent_type, put all the link_type in it.
			// Question: what means? PCE/CRU
			for (it = assignment.g_LinkTypeMap.begin();
				it != assignment.g_LinkTypeMap.end(); it++)
			{
				char field_name[20];

				sprintf_s(field_name, "PCE_link_type%d", it->first);
				value = 1; //init_PCE_value = 0
				parser_agent_type.GetValueByFieldName(field_name, value, false, false);
				//Question: what means this template ?

				agent_type.PCE_link_type_map[it->first] = value;

				value = 0; //init_CRU_value
				sprintf_s(field_name, "CRU_link_type%d", it->first);
				parser_agent_type.GetValueByFieldName(field_name, value, false, false);

				agent_type.CRU_link_type_map[it->first] = value;

			}
			//Question: this map?
			assignment.agent_type_2_seqno_mapping[agent_type.agent_type] = assignment.g_AgentTypeVector.size();

			assignment.g_AgentTypeVector.push_back(agent_type);
			assignment.g_number_of_agent_types = assignment.g_AgentTypeVector.size();
		}
		parser_agent_type.CloseCSVFile();

		if (assignment.g_AgentTypeVector.size() == 0)
		{
			cout << "Error: File agent_type.csv does not contain information." << endl;
		}

	}
	else
	{
		cout << "Error: File agent_type.csv cannot be opened.\n It might be currently used and locked by EXCEL." << endl;
		g_ProgramStop();
	}

	//MAX Agent_type number:4
	if (assignment.g_AgentTypeVector.size() >= _MAX_AGNETTYPES)
	{
		cout << "Error: agent_type = " << assignment.g_AgentTypeVector.size() << " in file agent_type.csv is too large. " << "_MAX_AGNETTYPES = " << _MAX_AGNETTYPES << "Please contact program developers!";

		g_ProgramStop();
	}


	assignment.g_number_of_nodes = 0;
	assignment.g_number_of_links = 0;  // initialize  the counter to 0

	// step 3: Read "node.csv" file
	// Build the node vector and zone vector at the same time
	int internal_node_seq_no = 0;

	std::map<int, int> zone_id_to_node_id_mapping;
	//In order to find the corresponding node_id from zone_id, we build this map.


	CCSVParser parser;
	if (parser.OpenCSVFile("node.csv", true))
	{
		while (parser.ReadRecord())
		{

			int node_id;

			if (parser.GetValueByFieldName("node_id", node_id) == false)
				continue;

			//scan the node list, this node has not been defined.
			if (assignment.g_internal_node_to_seq_no_map.find(node_id) != assignment.g_internal_node_to_seq_no_map.end())
			{
				continue; //has been defined
			}
			assignment.g_internal_node_to_seq_no_map[node_id] = internal_node_seq_no;


			CNode node;  // create a node object 

			node.node_id = node_id;
			node.node_seq_no = internal_node_seq_no;
			//Question: Why we use internal_node_seq_no?

			int zone_id = -1;
			parser.GetValueByFieldName("zone_id", zone_id);

			if (zone_id >= 1)
			{
				if (zone_id_to_node_id_mapping.find(zone_id) == zone_id_to_node_id_mapping.end())
				{
					zone_id_to_node_id_mapping[zone_id] = node_id;
					node.zone_id = zone_id;
					//For node and zone, build the map relationship
				}
				else
				{
					cout << "warning: zone_id " << zone_id << " have been defined more than once." << endl;
					//"if condition" is not satisfied =?=> zone_id has been defined more than once?
				}

			}

			/*node.x = x;
			node.y = y;*/
			//node.x and node.y are useless, we only need the node_id and corresponding zone_id

			internal_node_seq_no++;

			g_node_vector.push_back(node);  // push this node to the global node vector

			assignment.g_number_of_nodes++;
			if (assignment.g_number_of_nodes % 5000 == 0)
				cout << "reading " << assignment.g_number_of_nodes << " nodes.. " << endl;
		}

		cout << "number of nodes = " << assignment.g_number_of_nodes << endl;

		parser.CloseCSVFile();
	}


	// initialize zone vector
	for (int i = 0; i < g_node_vector.size(); i++)
	{
		//Question: I do not undersatnd this if condition. =��=
		if (g_node_vector[i].zone_id >= 1
			&& zone_id_to_node_id_mapping.find(g_node_vector[i].zone_id) != zone_id_to_node_id_mapping.end() /* uniquely defined*/
			&& assignment.g_zoneid_to_zone_seq_no_mapping.find(g_node_vector[i].zone_id) == assignment.g_zoneid_to_zone_seq_no_mapping.end())
			// create a new zone  // we assume each zone only has one node
		{   // we need to make sure we only create a zone in the memory if only there is positive demand flow from the (new) OD table

			//Transfer the zone array into zone vector(class_based)
			COZone ozone;
			ozone.node_seq_no = g_node_vector[i].node_seq_no;
			ozone.zone_id = g_node_vector[i].zone_id;
			ozone.zone_seq_no = g_zone_vector.size();
			assignment.g_zoneid_to_zone_seq_no_mapping[ozone.zone_id] = assignment.g_zoneid_to_zone_seq_no_mapping.size();  // create the zone id to zone seq no mapping
			g_zone_vector.push_back(ozone);  // add element into vector

		}
	}
	cout << "number of zones = " << g_zone_vector.size() << endl;


	// step 4: read "link.csv" file 
	CCSVParser parser_link;
	if (parser_link.OpenCSVFile("road_link.csv", true))
	{
		while (parser_link.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
		{
			int from_node_id;
			int to_node_id;
			if (parser_link.GetValueByFieldName("from_node_id", from_node_id) == false)
				continue;
			if (parser_link.GetValueByFieldName("to_node_id", to_node_id) == false)
				continue;

			string linkID;
			parser_link.GetValueByFieldName("road_link_id", linkID);


			// add the to node id into the outbound (adjacent) node list
			if (assignment.g_internal_node_to_seq_no_map.find(from_node_id) == assignment.g_internal_node_to_seq_no_map.end())
				//Question:...map.find == ...map.end What means? =>cannot find this from node in node list 
				//I may not understand this grammar.
				//These(and following) "if condition" may be useless when the correct input
			{
				cout << "Error: from_node_id " << from_node_id << " in file road_link.csv is not defined in node.csv." << endl;

				continue; //has not been defined
			}
			if (assignment.g_internal_node_to_seq_no_map.find(to_node_id) == assignment.g_internal_node_to_seq_no_map.end())
			{
				cout << "Error: to_node_id " << to_node_id << " in file road_link.csv is not defined in node.csv." << endl;
				continue; //has not been defined
			}

			if (assignment.g_road_link_id_map.find(linkID) != assignment.g_road_link_id_map.end())
			{
				cout << "Error: road_link_id " << linkID.c_str() << " has been defined more than once. Please check road_link.csv." << endl;
				continue; //has not been defined
			}
			int internal_from_node_seq_no = assignment.g_internal_node_to_seq_no_map[from_node_id];  // map external node number to internal node seq no. 
			int internal_to_node_seq_no = assignment.g_internal_node_to_seq_no_map[to_node_id];


			CLink link;
			// create a link object, which saves the attributes.

			link.from_node_seq_no = internal_from_node_seq_no;
			link.to_node_seq_no = internal_to_node_seq_no;
			link.link_seq_no = assignment.g_number_of_links;
			link.to_node_seq_no = internal_to_node_seq_no;
			link.link_id = linkID;

			assignment.g_road_link_id_map[link.link_id] = 1;


			parser_link.GetValueByFieldName("facility_type", link.type, true, false);
			parser_link.GetValueByFieldName("link_type", link.link_type);

			if (assignment.g_LinkTypeMap.find(link.link_type) == assignment.g_LinkTypeMap.end())
			{
				cout << "link type " << link.link_type << " in road_link.csv is not defined in link_type.csv" << endl;
				g_ProgramStop();

			}

			//Find the CRU information from Agent map
			for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
			{
				float PCE_ratio = assignment.g_AgentTypeVector[at].PCE_link_type_map[link.link_type];
				link.PCE_at[at] = PCE_ratio;
				float CRU = assignment.g_AgentTypeVector[at].CRU_link_type_map[link.link_type];
				link.CRU_at[at] = CRU;
			}

			//Question:link_type == c? What means��do following thing?
			if (assignment.g_LinkTypeMap[link.link_type].type_code == "c" && g_node_vector[internal_from_node_seq_no].zone_id >= 0)
			{
				if (assignment.g_zoneid_to_zone_seq_no_mapping.find(g_node_vector[internal_from_node_seq_no].zone_id) != assignment.g_zoneid_to_zone_seq_no_mapping.end())
					link.zone_seq_no_for_outgoing_connector = assignment.g_zoneid_to_zone_seq_no_mapping[g_node_vector[internal_from_node_seq_no].zone_id];
			}


			parser_link.GetValueByFieldName("toll", link.toll, false, false);
			parser_link.GetValueByFieldName("additional_cost", link.route_choice_cost, false, false);

			//init value
			float length = 1.0; // km or mile
			float free_speed = 1.0;
			float k_jam = 200;

			float lane_capacity = 1800;
			parser_link.GetValueByFieldName("length", length);
			parser_link.GetValueByFieldName("free_speed", free_speed);
			free_speed = max(0.1, free_speed);

			int number_of_lanes = 1;
			parser_link.GetValueByFieldName("lanes", number_of_lanes);
			parser_link.GetValueByFieldName("capacity", lane_capacity);

			//Question: I may learn something about VDF or CRU
			char VDF_field_name[20];
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; tau++)
			{
				int demand_period_id = assignment.g_DemandPeriodVector[tau].demand_period_id;
				sprintf_s(VDF_field_name, "VDF_fftt%d", demand_period_id);
				parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].FFTT);

				sprintf_s(VDF_field_name, "VDF_cap%d", demand_period_id);
				parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].capacity);

				sprintf_s(VDF_field_name, "VDF_alpha%d", demand_period_id);
				parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].alpha);

				sprintf_s(VDF_field_name, "VDF_beta%d", demand_period_id);
				parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].beta);

				if (assignment.assignment_mode == 3) //mode == 3 =>UE + RC
				{
					sprintf_s(VDF_field_name, "RUC_rho%d", demand_period_id);
					parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].rho);

					sprintf_s(VDF_field_name, "RUC_resource%d", demand_period_id);
					parser_link.GetValueByFieldName(VDF_field_name, link.VDF_period[tau].ruc_base_resource, false);

					link.VDF_period[tau].starting_time_slot_no = assignment.g_DemandPeriodVector[tau].starting_time_slot_no;
					link.VDF_period[tau].ending_time_slot_no = assignment.g_DemandPeriodVector[tau].ending_time_slot_no;
					link.VDF_period[tau].period = assignment.g_DemandPeriodVector[tau].time_period;
				}

			}

			if (assignment.assignment_mode == 3)
			{
				parser_link.GetValueByFieldName("RUC_type", link.RUC_type);
			}

			// for each period

			float default_cap = 1000;
			float default_BaseTT = 1;

			// setup default value
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; tau++)
			{
				link.TDBaseTT[tau] = default_BaseTT;
				link.TDBaseCap[tau] = default_cap;
			}

			link.number_of_lanes = number_of_lanes;
			link.lane_capacity = lane_capacity;
			link.link_spatial_capacity = k_jam * number_of_lanes * length;


			link.length = length;
			for (int tau = 0; tau < assignment.g_number_of_demand_periods; tau++)
			{
				link.travel_time_per_period[tau] = length / free_speed * 60;
			}
			link.free_flow_travel_time_in_min = length / free_speed * 60;

			// min // calculate link cost based length and speed limit // later we should also read link_capacity, calculate delay 

			//Build the two-way connection between node and link(pred node and to node)
			g_node_vector[internal_from_node_seq_no].m_outgoing_link_seq_no_vector.push_back(link.link_seq_no);
			// add this link to the corresponding node as part of outgoing node/link
			g_node_vector[internal_to_node_seq_no].m_incoming_link_seq_no_vector.push_back(link.link_seq_no);
			// add this link to the corresponding node as part of outgoing node/link

			g_node_vector[internal_from_node_seq_no].m_to_node_seq_no_vector.push_back(link.to_node_seq_no);
			// add this link to the corresponding node as part of outgoing node/link
			g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[link.to_node_seq_no] = link.link_seq_no;
			// add this link to the corresponding node as part of outgoing node/link

			g_link_vector.push_back(link);

			assignment.g_number_of_links++;

			if (assignment.g_number_of_links % 10000 == 0)
				cout << "reading " << assignment.g_number_of_links << " links.. " << endl;
		}
		parser_link.CloseCSVFile();
	}
	// Here, we now know the number of links
	cout << "number of links = " << assignment.g_number_of_links << endl;



	//Question: mode == 2: SO =>we need service_arc?? We give up "service_arc.csv"??
	if (assignment.assignment_mode == 2)
	{
		CCSVParser parser_service_arc;

		if (parser_service_arc.OpenCSVFile("service_arc.csv", true))
		{
			while (parser_service_arc.ReadRecord())  // if this line contains [] mark, then we will also read field headers.
			{
				int from_node_id = 0;
				int to_node_id = 0;
				if (parser_service_arc.GetValueByFieldName("from_node_id", from_node_id) == false)
				{
					cout << "Error: from_node_id in file service_arc.csv is not defined." << endl;
					continue;
				}
				if (parser_service_arc.GetValueByFieldName("to_node_id", to_node_id) == false)
				{
					continue;
				}

				if (assignment.g_internal_node_to_seq_no_map.find(from_node_id) == assignment.g_internal_node_to_seq_no_map.end())
				{
					cout << "Error: from_node_id " << from_node_id << " in file service_arc.csv is not defined in node.csv." << endl;

					continue; //has not been defined
				}
				if (assignment.g_internal_node_to_seq_no_map.find(to_node_id) == assignment.g_internal_node_to_seq_no_map.end())
				{
					cout << "Error: to_node_id " << to_node_id << " in file service_arc.csv is not defined in node.csv." << endl;
					continue; //has not been defined
				}


				int internal_from_node_seq_no = assignment.g_internal_node_to_seq_no_map[from_node_id];  // map external node number to internal node seq no. 
				int internal_to_node_seq_no = assignment.g_internal_node_to_seq_no_map[to_node_id];

				CServiceArc service_arc;  // create a link object 


				if (g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.find(internal_to_node_seq_no) != g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map.end())
				{
					service_arc.link_seq_no = g_node_vector[internal_from_node_seq_no].m_to_node_2_link_seq_no_map[internal_to_node_seq_no];

					if (assignment.g_LinkTypeMap[g_link_vector[service_arc.link_seq_no].link_type].type_code != "s")
					{
						cout << "Error: All service arcs defined in service_arc.csv should have a type_code = s in link_type.csv." << endl;
						cout << " Please check " << from_node_id << "->" << to_node_id << " in file service_arc.csv." << endl;
						continue;

					}
				}
				else
				{
					cout << "Error: Link " << from_node_id << "->" << to_node_id << " in file service_arc.csv is not defined in road_link.csv." << endl;
					continue;
				}


				string time_period;
				if (parser_service_arc.GetValueByFieldName("time_window", time_period) == false)
				{
					cout << "Error: Field time_window in file service_arc.csv cannot be read." << endl;
					g_ProgramStop();
					break;
				}

				vector<float> global_minute_vector;

				vector<string> input_string;
				input_string.push_back(time_period);
				//input_string includes the start and end time of a time period with hhmm format
				global_minute_vector = g_time_parser(input_string); //global_minute_vector incldue the starting and ending time
				if (global_minute_vector.size() == 2)
				{
					if (global_minute_vector[0] < assignment.g_LoadingStartTimeInMin)
						global_minute_vector[0] = assignment.g_LoadingStartTimeInMin;

					if (global_minute_vector[0] > assignment.g_LoadingEndTimeInMin)
						global_minute_vector[0] = assignment.g_LoadingEndTimeInMin;

					if (global_minute_vector[1] < assignment.g_LoadingStartTimeInMin)
						global_minute_vector[1] = assignment.g_LoadingStartTimeInMin;

					if (global_minute_vector[1] > assignment.g_LoadingEndTimeInMin)
						global_minute_vector[1] = assignment.g_LoadingEndTimeInMin;

					if (global_minute_vector[1] < global_minute_vector[0])
						global_minute_vector[1] = global_minute_vector[0];


					service_arc.starting_time_no = (global_minute_vector[0] - assignment.g_LoadingStartTimeInMin) * 60 / number_of_seconds_per_interval;
					service_arc.ending_time_no = (global_minute_vector[1] - assignment.g_LoadingStartTimeInMin) * 60 / number_of_seconds_per_interval;

				}
				else
				{
					continue;
				}

				float time_interval = 0;

				parser_service_arc.GetValueByFieldName("time_interval", time_interval, true, false);
				service_arc.time_interval_no = max(1, time_interval * 60.0 / number_of_seconds_per_interval);


				int travel_time_delta_in_min = 0;
				parser_service_arc.GetValueByFieldName("travel_time_delta", travel_time_delta_in_min, true, false);
				service_arc.travel_time_delta = max(1, travel_time_delta_in_min * 60 / number_of_seconds_per_interval);

				float capaciy = 1;
				parser_service_arc.GetValueByFieldName("capacity", capaciy);  // capacity in the space time arcs

				g_service_arc_vector.push_back(service_arc);
				assignment.g_number_of_service_arcs++;

				if (assignment.g_number_of_service_arcs % 10000 == 0)
					cout << "reading " << assignment.g_number_of_service_arcs << " service_arcs.. " << endl;
			}

			parser_service_arc.CloseCSVFile();
		}
	}
	// we now know the number of links
	cout << "number of service_arcs = " << assignment.g_number_of_service_arcs << endl;
};
//Suggestion: We can set some defult vaule for missing parameter instead of printing "error".
//Finished g_ReadInputData 7/22 8:30



void update_link_travel_time_and_cost()
{
#pragma omp parallel for
	for (int l = 0; l < g_link_vector.size(); l++)
	{
		int tau;
		// step 1: travel time based on VDF
		g_link_vector[l].CalculateTD_VDFunction();


		for (tau = 0; tau < assignment.g_DemandPeriodVector.size(); tau++)
		{
			if (g_debugging_flag >= 2 && assignment.g_pFileDebugLog != NULL)
				fprintf(assignment.g_pFileDebugLog, "Update link resource: link %d->%d: tau = %d, volume = %.2f, travel time = %.2f, resource = %.3f\n",

					g_node_vector[g_link_vector[l].from_node_seq_no].node_id,
					g_node_vector[g_link_vector[l].to_node_seq_no].node_id,
					tau,
					g_link_vector[l].flow_volume_per_period[tau],
					g_link_vector[l].travel_time_per_period[tau],
					g_link_vector[l].resource_per_period[tau]);


			for (int at = 0; at < assignment.g_AgentTypeVector.size(); at++)
			{

				float PCE_agent_type = assignment.g_AgentTypeVector[at].PCE_link_type_map[g_link_vector[l].link_type];

				// step 2: marginal cost for SO
				g_link_vector[l].calculate_marginal_cost_for_agent_type(tau, at, PCE_agent_type);

				float CRU_agent_type = assignment.g_AgentTypeVector[at].CRU_link_type_map[g_link_vector[l].link_type];

				// setp 3: penalty for resource constraints 
				g_link_vector[l].calculate_penalty_for_agent_type(tau, at, CRU_agent_type);


				if (g_debugging_flag >= 2 && assignment.assignment_mode >= 2 && assignment.g_pFileDebugLog != NULL)
					fprintf(assignment.g_pFileDebugLog, "Update link cost: link %d->%d: tau = %d, at = %d, travel_marginal =  %.3f; penalty derivative = %.3f\n",

						g_node_vector[g_link_vector[l].from_node_seq_no].node_id,
						g_node_vector[g_link_vector[l].to_node_seq_no].node_id,
						tau, at,
						g_link_vector[l].travel_marginal_cost_per_period[tau][at],
						g_link_vector[l].exterior_penalty_derivative_per_period[tau][at]);

			}
		}
	}
}



char str_buffer[STRING_LENGTH_PER_LINE];


//***
// major function 1:  allocate memory and initialize the data
// void AllocateMemory(int number_of_nodes)
//
//major function 2: // time-dependent label correcting algorithm with double queue implementation
//int optimal_label_correcting(int origin_node, int destination_node, int departure_time, int g_debugging_flag, FILE* g_pFileDebugLog, Assignment& assignment, int time_period_no = 0, int agent_type = 1, float VOT = 10)

//	//major function: update the cost for each node at each SP tree, using a stack from the origin structure 
//int tree_cost_updating(int origin_node, int departure_time, int g_debugging_flag, FILE* g_pFileDebugLog, Assignment& assignment, int time_period_no = 0, int agent_type = 1)

//***

// The one and only application object


int g_number_of_CPU_threads()
{
	int number_of_threads = omp_get_max_threads();

	int max_number_of_threads = 4000;

	if (number_of_threads > max_number_of_threads)
		number_of_threads = max_number_of_threads;

	return number_of_threads;

}




void  CLink::CalculateTD_VDFunction()
{
	for (int tau = 0; tau < assignment.g_number_of_demand_periods; tau++)
	{
		float starting_time_slot_no = assignment.g_DemandPeriodVector[tau].starting_time_slot_no;
		float ending_time_slot_no = assignment.g_DemandPeriodVector[tau].ending_time_slot_no;
		travel_time_per_period[tau] = VDF_period[tau].PerformBPR(flow_volume_per_period[tau]);
	}
}



///////////////////////////////////////////SUE: Wang Entai///////////////////////////////////////////

template <typename T>
vector<size_t> sort_indexes_e(vector<T>& v)
{
	vector<size_t> idx(v.size());
	iota(idx.begin(), idx.end(), 0);
	sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });
	return idx;
}



class NetworkForSUE
{
public:
	std::vector<int> origin_node_vector;
	std::vector<int> destination_node_vector;
	std::vector<int> demand_volume_vector;

	//int* sub_network_link_for_each_od;
	//int* sub_network_node_for_each_od;
	//double* sub_network_volume_for_each_od;
	//int node_number_for_sub_network;
	//int link_number_for_sub_network;

	int b_Dial = 1;
	float mu = 1;

	int m_g_number_of_links = assignment.g_number_of_links;
	int m_g_number_of_nodes = assignment.g_number_of_nodes;

	double** link_flow_aux;
	double** link_sys_utility_OD;
	double** link_flow_OD;

	double** Pi_OD_aux;
	double** Pi_OD;


	void AllocateMemorySUE(int number_of_nodes, int number_of_links) {

		//sub_network_link_for_each_od = new int* [demand_volume_vector.size()];
		//sub_network_node_for_each_od = new int* [demand_volume_vector.size()];
		//sub_network_volume_for_each_od = new double* [demand_volume_vector.size()];

		//link_number_for_sub_network = new int[demand_volume_vector.size()]();
		//node_number_for_sub_network = new int[demand_volume_vector.size()]();

		link_flow_aux = new double* [demand_volume_vector.size()]();
		link_sys_utility_OD = new double* [demand_volume_vector.size()]();
		link_flow_OD = new double* [demand_volume_vector.size()]();

		Pi_OD_aux = new double* [demand_volume_vector.size()]();
		Pi_OD = new double* [demand_volume_vector.size()]();
	}

	void Dial_subnetwork_init(float** adj_matrix) {
		float* r_distance = new float[assignment.g_number_of_nodes];
		float* s_distance = new float[assignment.g_number_of_nodes];
		double* link_likelihood = new double[assignment.g_number_of_links];
		double* link_weights_list = new double[assignment.g_number_of_links];
		double* link_label_list = new double[assignment.g_number_of_links];
		double* link_flow_list = new double[assignment.g_number_of_links];

		for (int od = 0; od < origin_node_vector.size(); od++) {
			int origin_node = origin_node_vector[od];
			int destination_node = destination_node_vector[od];
			int demand_volume = demand_volume_vector[od];
			for (int j = 0; j < assignment.g_number_of_links; j++) {
				link_likelihood[j] = 0.0;
				link_weights_list[j] = 0.0;
				link_label_list[j] = 0.0;
				link_flow_list[j] = 0.0;
			}
			for (int j = 0; j < assignment.g_number_of_nodes; j++) {
				r_distance[j] = 0.0;
				s_distance[j] = 0.0;
			}

			//Step 1：
			for (int i = 0; i < assignment.g_number_of_nodes; i++) {
				r_distance[i] = adj_matrix[origin_node][i];
			}
			for (int i = 0; i < assignment.g_number_of_nodes; i++) {
				s_distance[i] = adj_matrix[i][destination_node];
			}

			for (int i = 0; i < assignment.g_number_of_links; i++) {
				int from_node = g_link_vector[i].from_node_seq_no;
				int to_node = g_link_vector[i].to_node_seq_no;
				if ((r_distance[from_node] < r_distance[to_node]) & (s_distance[from_node] > s_distance[to_node])) {
					link_likelihood[i] = exp(b_Dial * (r_distance[to_node] - r_distance[from_node] - g_link_vector[i].length));
					if (link_likelihood[i] > 65535)
						link_likelihood[i] = 0;
				}
			}

			//Step 2:
			vector<float> r_distance_array(assignment.g_number_of_nodes);
			for (int i = 0; i < r_distance_array.size(); i++)
			{
				r_distance_array[i] = r_distance[i];
			}
			vector<size_t> r_index_list;
			r_index_list = sort_indexes_e(r_distance_array);
			vector<float> s_distance_array(assignment.g_number_of_nodes);
			for (int i = 0; i < s_distance_array.size(); i++)
			{
				s_distance_array[i] = s_distance[i];
			}
			vector<size_t> s_index_list;
			s_index_list = sort_indexes_e(s_distance_array);


			for each (int node_list_no in r_index_list) {
				CNode node = g_node_vector[node_list_no];
				if (node.node_seq_no == destination_node)
					break;
				if (node.node_seq_no == origin_node) {
					for each (int next_link in g_node_vector[node_list_no].m_outgoing_link_seq_no_vector) {
						link_weights_list[next_link] = link_likelihood[next_link];
						link_label_list[next_link] = 1.0;
					}
				}
				else {
					float current_label = 1.0;
					for each (int pre_link in g_node_vector[node_list_no].m_incoming_link_seq_no_vector) {
						if (link_likelihood[pre_link] > 0) {
							current_label = current_label * link_label_list[pre_link];
						}
					}
					if (current_label == 1.0) {
						for each (int next_link in g_node_vector[node_list_no].m_outgoing_link_seq_no_vector) {
							float sum_w = 0.0;
							for each (int pre_link in g_node_vector[node_list_no].m_incoming_link_seq_no_vector) {
								sum_w += link_weights_list[pre_link];
							}
							link_weights_list[next_link] = link_likelihood[next_link] * sum_w;
							link_label_list[next_link] = 1.0;
						}
					}
				}
			}

			//Step 3:
			memset(link_label_list, 0.0, sizeof(link_label_list));
			for each (int node_list_no in s_index_list) {
				CNode node = g_node_vector[node_list_no];
				if (node.node_seq_no == origin_node)
					break;
				if (node.node_seq_no == destination_node) {
					float sum_w_pre = 0.0;
					for each (int pre_link in g_node_vector[node_list_no].m_incoming_link_seq_no_vector) {
						sum_w_pre += link_weights_list[pre_link];
					}
					for each (int pre_link in g_node_vector[node_list_no].m_incoming_link_seq_no_vector) {
						if (link_weights_list[pre_link] == 0) {
							link_flow_list[pre_link] = 0.0;
							link_label_list[pre_link] = 1.0;
						}
						else {
							link_flow_list[pre_link] = (link_weights_list[pre_link] / sum_w_pre) * demand_volume;
							link_label_list[pre_link] = 1.0;
						}
					}
				}
				else {
					float current_label = 1.0;
					float sum_flow = 0.0;
					for each (int next_link in g_node_vector[node_list_no].m_outgoing_link_seq_no_vector) {
						if (link_weights_list[next_link] > 0) {
							current_label = current_label * link_label_list[next_link];
							sum_flow += link_flow_list[next_link];
						}
					}
					if (current_label == 1.0) {
						float sum_w_pre = 0.0;
						for each (int pre_link in g_node_vector[node_list_no].m_incoming_link_seq_no_vector) {
							sum_w_pre += link_weights_list[pre_link];
						}
						for each (int pre_link in g_node_vector[node_list_no].m_incoming_link_seq_no_vector) {
							if (link_weights_list[pre_link] == 0) {
								link_flow_list[pre_link] = 0.0;
								link_label_list[pre_link] = 1.0;
							}
							else {
								link_flow_list[pre_link] = (link_weights_list[pre_link] / sum_w_pre) * sum_flow;
								link_label_list[pre_link] = 1.0;
							}
						}
					}
				}
			}

			//Update globle link volume for each od
			for (int j = 0; j < m_g_number_of_links; j++) {
				if (link_flow_list[j] != 0) {
					g_link_vector[j].flow_volume_per_period[0] += link_flow_list[j];
				}
			}
		}

		delete[] r_distance;
		delete[] s_distance;
		delete[] link_likelihood;
		delete[] link_weights_list;
		delete[] link_label_list;
		delete[] link_flow_list;
	}


	void Dial_subnetwork(int origin_node, int destination_node, float demand_volume, float** adj_matrix, \
		int*& network_cur, int*& node_list, int& local_link_numbers, int& local_node_numbers) {

		float* r_distance = new float[assignment.g_number_of_nodes];
		float* s_distance = new float[assignment.g_number_of_nodes];
		double* link_likelihood = new double[assignment.g_number_of_links];
		double* link_weights_list = new double[assignment.g_number_of_links];
		double* link_label_list = new double[assignment.g_number_of_links];
		double* link_flow_list = new double[assignment.g_number_of_links];

		for (int j = 0; j < assignment.g_number_of_links; j++) {
			link_likelihood[j] = 0.0;
			link_weights_list[j] = 0.0;
			link_label_list[j] = 0.0;
			link_flow_list[j] = 0.0;
		}
		for (int j = 0; j < assignment.g_number_of_nodes; j++) {
			r_distance[j] = 0.0;
			s_distance[j] = 0.0;
		}

		//Step 1：
		for (int i = 0; i < assignment.g_number_of_nodes; i++) {
			r_distance[i] = adj_matrix[origin_node][i];
		}
		for (int i = 0; i < assignment.g_number_of_nodes; i++) {
			s_distance[i] = adj_matrix[i][destination_node];
		}

		for (int i = 0; i < assignment.g_number_of_links; i++) {
			int from_node = g_link_vector[i].from_node_seq_no;
			int to_node = g_link_vector[i].to_node_seq_no;
			if ((r_distance[from_node] < r_distance[to_node]) & (s_distance[from_node] > s_distance[to_node])) {
				link_likelihood[i] = exp(b_Dial * (r_distance[to_node] - r_distance[from_node] - g_link_vector[i].length));
				if (link_likelihood[i] > 65535)
					link_likelihood[i] = 0;
			}
		}

		//Step 2:
		vector<float> r_distance_array(assignment.g_number_of_nodes);
		for (int i = 0; i < r_distance_array.size(); i++)
		{
			r_distance_array[i] = r_distance[i];
		}
		vector<size_t> r_index_list;
		r_index_list = sort_indexes_e(r_distance_array);
		vector<float> s_distance_array(assignment.g_number_of_nodes);
		for (int i = 0; i < s_distance_array.size(); i++)
		{
			s_distance_array[i] = s_distance[i];
		}
		vector<size_t> s_index_list;
		s_index_list = sort_indexes_e(s_distance_array);


		for each (int node_list_no in r_index_list) {
			CNode node = g_node_vector[node_list_no];
			if (node.node_seq_no == destination_node)
				break;
			if (node.node_seq_no == origin_node) {
				for each (int next_link in g_node_vector[node_list_no].m_outgoing_link_seq_no_vector) {
					link_weights_list[next_link] = link_likelihood[next_link];
					link_label_list[next_link] = 1.0;
				}
			}
			else {
				float current_label = 1.0;
				for each (int pre_link in g_node_vector[node_list_no].m_incoming_link_seq_no_vector) {
					if (link_likelihood[pre_link] > 0) {
						current_label = current_label * link_label_list[pre_link];
					}
				}
				if (current_label == 1.0) {
					for each (int next_link in g_node_vector[node_list_no].m_outgoing_link_seq_no_vector) {
						float sum_w = 0.0;
						for each (int pre_link in g_node_vector[node_list_no].m_incoming_link_seq_no_vector) {
							sum_w += link_weights_list[pre_link];
						}
						link_weights_list[next_link] = link_likelihood[next_link] * sum_w;
						link_label_list[next_link] = 1.0;
					}
				}
			}
		}

		//Step 3:
		memset(link_label_list, 0.0, sizeof(link_label_list));
		for each (int node_list_no in s_index_list) {
			CNode node = g_node_vector[node_list_no];
			if (node.node_seq_no == origin_node)
				break;
			if (node.node_seq_no == destination_node) {
				float sum_w_pre = 0.0;
				for each (int pre_link in g_node_vector[node_list_no].m_incoming_link_seq_no_vector) {
					sum_w_pre += link_weights_list[pre_link];
				}
				for each (int pre_link in g_node_vector[node_list_no].m_incoming_link_seq_no_vector) {
					if (link_weights_list[pre_link] == 0) {
						link_flow_list[pre_link] = 0.0;
						link_label_list[pre_link] = 1.0;
					}
					else {
						link_flow_list[pre_link] = (link_weights_list[pre_link] / sum_w_pre) * demand_volume;
						link_label_list[pre_link] = 1.0;
					}
				}
			}
			else {
				float current_label = 1.0;
				float sum_flow = 0.0;
				for each (int next_link in g_node_vector[node_list_no].m_outgoing_link_seq_no_vector) {
					if (link_weights_list[next_link] > 0) {
						current_label = current_label * link_label_list[next_link];
						sum_flow += link_flow_list[next_link];
					}
				}
				if (current_label == 1.0) {
					float sum_w_pre = 0.0;
					for each (int pre_link in g_node_vector[node_list_no].m_incoming_link_seq_no_vector) {
						sum_w_pre += link_weights_list[pre_link];
					}
					for each (int pre_link in g_node_vector[node_list_no].m_incoming_link_seq_no_vector) {
						if (link_weights_list[pre_link] == 0) {
							link_flow_list[pre_link] = 0.0;
							link_label_list[pre_link] = 1.0;
						}
						else {
							link_flow_list[pre_link] = (link_weights_list[pre_link] / sum_w_pre) * sum_flow;
							link_label_list[pre_link] = 1.0;
						}
					}
				}
			}
		}

		//1. Update link numbers for each od
		local_link_numbers = 0;
		for (int i = 0; i < assignment.g_number_of_links; i++) {
			if (link_flow_list[i] != 0) {
				local_link_numbers++;
			}
		}


		//2.Update link index for each od
		int count = 0;
		network_cur = new int[local_link_numbers];
		for (int j = 0; j < m_g_number_of_links; j++) {
			if (link_flow_list[j] != 0) {
				network_cur[count] = j;
				count++;
			}
		}

		//3. Update node numbers for each od
		int* node_SElist = new int[m_g_number_of_nodes]();
		for (int i = 0; i < assignment.g_number_of_links; i++) {
			if (link_flow_list[i] != 0) {
				node_SElist[g_link_vector[i].from_node_seq_no] = 1;
				node_SElist[g_link_vector[i].to_node_seq_no] = 1;
			}
		}
		local_node_numbers = 0;
		for (int i = 0; i < m_g_number_of_nodes; i++) {
			if (node_SElist[i] != 0) {
				local_node_numbers++;
			}
		}

		//4. Update node index for each od
		node_list = new int[local_node_numbers];
		count = 0;
		for (int i = 0; i < m_g_number_of_nodes; i++) {
			if (node_SElist[i] != 0) {
				node_list[count] = i;
				count++;
			}
		}
		delete[] node_SElist;
		delete[] r_distance;
		delete[] s_distance;
		delete[] link_likelihood;
		delete[] link_weights_list;
		delete[] link_label_list;
		delete[] link_flow_list;
	}


	void RLSUE(Assignment& assignment, int od, int iter_no, \
		int* network_cur, int* node_list, int local_link_numbers, int local_node_numbers) {

		float* node_utility;
		float* link_sys_utility;
		double* link_con_prob;
		double* link_prob;
		double* link_flow;
		float* link_label;

		node_utility = new float[local_node_numbers];
		link_sys_utility = new float[local_link_numbers];
		link_con_prob = new double[local_link_numbers];
		link_prob = new double[local_link_numbers];
		link_flow = new double[local_link_numbers];
		link_label = new float[local_link_numbers];

		//int* network_cur = sub_network_link_for_each_od[od];
		//int local_node_numbers = node_number_for_sub_network[od];
		///int local_link_numbers = link_number_for_sub_network[od];

		int origin_node = origin_node_vector[od];
		int destination_node = destination_node_vector[od];

		for (int i = 0; i < local_node_numbers; i++) {
			node_utility[i] = 0;
		}
		for (int i = 0; i < local_link_numbers; i++) {
			link_sys_utility[i] = 0;
			link_prob[i] = 0;
			link_flow[i] = 0;
			link_label[i] = 0;
			link_con_prob[i] = 0;
		}

		Utility(od, origin_node, destination_node, network_cur, node_utility, link_sys_utility, \
			local_node_numbers, local_link_numbers, node_list);


		//link conditional probability
		vector<int> node_NOD_list;
		for (int i = 0; i < local_node_numbers; i++) {
			node_NOD_list.push_back(node_list[i]);
		}
		std::vector<int>::iterator iter = std::find(node_NOD_list.begin(), node_NOD_list.end(), destination_node);
		node_NOD_list.erase(iter);

		for (int i = 0; i < node_NOD_list.size(); i++) {
			int node = node_NOD_list[i];
			double sum_exp_p = 0.0;
			for each (int next_link in g_node_vector[node].m_outgoing_link_seq_no_vector) {
				int next_link_index;
				for (int j = 0; j < local_link_numbers; j++) {
					if (next_link == network_cur[j]) {
						next_link_index = j;
						sum_exp_p += exp(1 / (double)mu * (double)(link_sys_utility[next_link_index]));
					}
				}
			}
			double link_V;
			for each (int next_link in g_node_vector[node].m_outgoing_link_seq_no_vector) {
				int next_link_index;
				for (int j = 0; j < local_link_numbers; j++) {
					if (next_link == network_cur[j]) {
						next_link_index = j;
						link_V = exp(1 / (double)mu * (double)(link_sys_utility[next_link_index]));
						link_con_prob[next_link_index] = link_V / sum_exp_p;
						if (!isfinite(link_con_prob[next_link_index])) {
							link_con_prob[next_link_index] = 0;
							//cout << "exp(-inf) / exp(inf)!" << endl;	
						}
					}
				}
			}
		}

		//link probability and flow
		vector<int> node_NoO_list;
		for (int i = 0; i < local_node_numbers; i++) {
			node_NoO_list.push_back(node_list[i]);
		}
		iter = std::find(node_NoO_list.begin(), node_NoO_list.end(), origin_node);
		node_NoO_list.erase(iter);

		for (int link = 0; link < local_link_numbers; link++) {
			if (g_link_vector[network_cur[link]].from_node_seq_no == origin_node) {
				link_prob[link] = link_con_prob[link];
				link_flow[link] = demand_volume_vector[od] * link_prob[link];
				link_label[link] = 1;
			}
		}

		while (node_NoO_list.size() > 0) {
			for (int t = 0; t < node_NoO_list.size(); t++) {
				int node = node_NoO_list[t];
				float value = 1.0;
				double sum_prob = 0.0;
				for each (int pred_link in g_node_vector[node].m_incoming_link_seq_no_vector) {
					int pred_link_index;
					for (int j = 0; j < local_link_numbers; j++) {
						if (pred_link == network_cur[j]) {
							pred_link_index = j;
							value = value * link_label[pred_link_index];
							sum_prob += link_prob[pred_link_index];
						}
					}
				}
				if (value == 1) {
					for each (int next_link in g_node_vector[node].m_outgoing_link_seq_no_vector) {
						int next_link_index;
						for (int j = 0; j < local_link_numbers; j++) {
							if (next_link == network_cur[j]) {
								next_link_index = j;
								link_prob[next_link_index] = sum_prob * link_con_prob[next_link_index];
								link_flow[next_link_index] = (double)demand_volume_vector[od] * link_prob[next_link_index];
								link_label[next_link_index] = 1;
							}
						}
					}
					iter = std::find(node_NoO_list.begin(), node_NoO_list.end(), node);
					node_NoO_list.erase(iter);
				}
			}
		}

		if (iter_no == 1) {
			link_flow_aux[od] = new double[local_link_numbers];
			link_sys_utility_OD[od] = new double[local_link_numbers];
		}
		for (int i = 0; i < local_link_numbers; i++) {
			link_flow_aux[od][i] = 0.0;
			link_sys_utility_OD[od][i] = 0.0;
		}
		for (int i = 0; i < local_link_numbers; i++) {
			link_flow_aux[od][i] = link_flow[i];
			link_sys_utility_OD[od][i] = link_sys_utility[i];
		}
		delete[] node_utility;
		delete[] link_sys_utility;
		delete[] link_con_prob;
		delete[] link_prob;
		delete[] link_flow;
		delete[] link_label;

	}


	void Utility(int od, int origin_node, int destination_node, int* network_cur, \
		float*& node_utility, float*& link_sys_utility, \
		int local_node_numbers, int local_link_numbers, int* node_list)
	{
		int destination_node_index = 0;
		for (int i = 0; i < local_node_numbers; i++) {
			if (node_list[i] == destination_node)
				destination_node_index = i;
		}
		node_utility[destination_node_index] = 0.0;
		vector<int> SElist;
		SElist.push_back(destination_node);
		int* node_label = new int[local_node_numbers]();
		node_label[destination_node_index] = 1;
		vector<int> node_record;
		for (int i = 0; i < local_node_numbers; i++) {
			node_record.push_back(node_list[i]);
		}
		for (vector<int>::iterator iter = node_record.begin(); iter != node_record.end(); iter++) {
			if (*iter == destination_node) {
				node_record.erase(iter);
				break;
			}
		}
		while ((SElist.size() > 0)& (node_record.size() > 0)) {
			int node = SElist[0];
			std::vector<int>::iterator iter = std::find(SElist.begin(), SElist.end(), node);
			SElist.erase(iter);

			for each (int pre_link in g_node_vector[node].m_incoming_link_seq_no_vector) {
				int pre_link_index = -1;
				for (int i = 0; i < local_link_numbers; i++) {
					if (network_cur[i] == pre_link)
						pre_link_index = i;
				}
				int from_node = g_link_vector[pre_link].from_node_seq_no;
				if (std::find(node_record.begin(), node_record.end(), from_node) == node_record.end()) //unfound
					continue;
				int flag = 0;
				for (int i = 0; i < local_link_numbers; i++) {
					if (network_cur[pre_link_index] == pre_link)
						flag = 1;
				}
				if (flag == 0)
					continue;

				float value = 1.0;
				for each (int next_link in g_node_vector[from_node].m_outgoing_link_seq_no_vector) {
					int next_link_index;
					for (int i = 0; i < local_link_numbers; i++) {
						if (network_cur[i] == next_link) {
							next_link_index = i;
							int to_node_seq_no_index = -1;
							for (int j = 0; j < local_node_numbers; j++) {
								if (node_list[j] == g_link_vector[next_link].to_node_seq_no)
									to_node_seq_no_index = j;
							}
							value = value * node_label[to_node_seq_no_index];
						}
					}

				}
				if (value == 1.0) {
					double sum_exp = 0.0;
					for each (int next_link in g_node_vector[from_node].m_outgoing_link_seq_no_vector) {
						int next_link_index;
						for (int i = 0; i < local_link_numbers; i++) {
							if (network_cur[i] == next_link) {
								next_link_index = i;
								int from_node_to_node = g_link_vector[next_link].to_node_seq_no;

								int from_node_to_node_index;
								for (int j = 0; j < local_node_numbers; j++) {
									if (node_list[j] == from_node_to_node) {
										from_node_to_node_index = j;
									}
								}

								sum_exp += exp(1 / (double)mu * (double)(-g_link_vector[next_link].travel_time_per_period[0] \
									+ (double)node_utility[from_node_to_node_index]));
							}
						}
					}
					int from_node_seq_no_index;
					for (int i = 0; i < local_node_numbers; i++) {
						if (node_list[i] == g_link_vector[pre_link].from_node_seq_no)
							from_node_seq_no_index = i;
					}
					node_utility[from_node_seq_no_index] = mu * log(sum_exp);

					
					//************************************************************************************
					if (node_utility[from_node_seq_no_index] < -100) {
						node_utility[from_node_seq_no_index] = -100;
					}
					//************************************************************************************
					

					node_label[from_node_seq_no_index] = 1;
					SElist.push_back(g_link_vector[pre_link].from_node_seq_no);
					if (node_record.size() > 0) {
						std::vector<int>::iterator iter = std::find(node_record.begin(), node_record.end(), g_link_vector[pre_link].from_node_seq_no);
						node_record.erase(iter);
					}
				}
			}
		}
		for (int link = 0; link < local_link_numbers; link++) {
			int to_node_seq_no_index = 0;
			for (int j = 0; j < local_node_numbers; j++) {
				if (g_link_vector[network_cur[link]].to_node_seq_no == node_list[j]) {
					to_node_seq_no_index = j;
				}
			}
			link_sys_utility[link] = -g_link_vector[network_cur[link]].travel_time_per_period[0] + \
				node_utility[to_node_seq_no_index];

			
			//************************************************************************************
			if (link_sys_utility[link] < -100)
				link_sys_utility[link] = -100;
			//************************************************************************************
			

		}
		delete[] node_label;
	}


	void Update_link_flow_OD(Assignment& assignment, int od, int iter_no, \
		int* local_node_list, int local_node_number, int* local_link_list, int local_link_number) {

		//int* local_node_list = sub_network_node_for_each_od[od];
		//int local_node_number = node_number_for_sub_network[od];
		//int* local_link_list = sub_network_link_for_each_od[od];
		//int local_link_number = link_number_for_sub_network[od];

		if (iter_no == 1) {
			link_flow_OD[od] = new double[local_link_number];
			for (int i = 0; i < local_link_number; i++) {
				link_flow_OD[od][i] = 0;
			}
		}
		for (int link = 0; link < local_link_number; link++) {
			link_flow_OD[od][link] = (double)link_flow_OD[od][link] + (1 / (double)iter_no) * (double)(link_flow_aux[od][link] - link_flow_OD[od][link]);
			g_link_vector[local_link_list[link]].flow_volume_per_period[0] += link_flow_OD[od][link];
			if (local_link_list[link] == 0) {
				int stop = 0;
			}
		}

	}


	void MAXUtility(Assignment& assignment, int od, int iter_no, \
		int* network_cur, int* node_list, int local_link_numbers, int local_node_numbers)
	{
		double* link_flow = link_flow_aux[od];
		double* link_sys_utility = link_sys_utility_OD[od];

		int origin_node = origin_node_vector[od];
		int destination_node = destination_node_vector[od];
		float demand_value = demand_volume_vector[od];
		//int local_node_number = node_number_for_sub_network[od];
		//int local_link_number = link_number_for_sub_network[od];

		//int* node_list = sub_network_node_for_each_od[od];
		//int* network_cur = sub_network_link_for_each_od[od];

		double* Pi = new double[local_node_numbers];

		if (node_list == NULL)
			cout << "The from-zone and to-zone are the same zones:" << origin_node << "-->" << destination_node << endl;

		for (int i = 0; i < local_node_numbers; i++) {
			Pi[i] = 0;
		}

		int destination_node_index;
		for (int i = 0; i < local_node_numbers; i++) {
			if (node_list[i] == destination_node) {
				destination_node_index = i;
			}
		}

		Pi[destination_node_index] = 0.0;
		vector<int> node_NoD_list;
		for (int i = 0; i < local_node_numbers; i++) {
			node_NoD_list.push_back(node_list[i]);
		}

		std::vector<int>::iterator iter = std::find(node_NoD_list.begin(), node_NoD_list.end(), destination_node);
		node_NoD_list.erase(iter);

		for each (int node in node_NoD_list) {
			double sum_flow_pre = 0.0;
			if (node == origin_node) {
				sum_flow_pre = demand_value;
			}
			else {
				for each (int pre_link in g_node_vector[node].m_incoming_link_seq_no_vector) {
					int pre_link_index;
					for (int i = 0; i < local_link_numbers; i++) {
						if (pre_link == network_cur[i]) {
							pre_link_index = i;
							sum_flow_pre += double(link_flow[pre_link_index]);
						}
					}
				}
			}

			double sum_utility_suc = 0.0;
			for each (int next_link in g_node_vector[node].m_outgoing_link_seq_no_vector) {
				int next_link_index;
				for (int i = 0; i < local_link_numbers; i++) {
					if (next_link == network_cur[i]) {
						next_link_index = i;
						sum_utility_suc += exp(1 / (double)mu * (double)link_sys_utility[next_link_index]);
					}
				}
			}
			int node_index;
			for (int i = 0; i < local_node_numbers; i++) {
				if (node == node_list[i])
					node_index = i;
			}
			if (sum_utility_suc < 0.001) {
				Pi[node_index] = 0.0;
			}
			else {
				Pi[node_index] = log(sum_flow_pre / sum_utility_suc);
			}
			if (iter_no == 1) {
				Pi_OD_aux[od] = new double[local_node_numbers];
			}

			for (int i = 0; i < local_node_numbers; i++) {
				Pi_OD_aux[od][i] = Pi[i];
			}
		}
		delete[] Pi;
	}


	void Update_Pi_OD(Assignment& assignment, int od, int iter_no, int local_node_number) {

		//int local_node_number = node_number_for_sub_network[od];
		if (iter_no == 1) {
			Pi_OD[od] = new double[local_node_number];
			for (int i = 0; i < local_node_number; i++) {
				Pi_OD[od][i] = 0;
			}
		}

		for (int node = 0; node < local_node_number; node++) {
			Pi_OD[od][node] = (double)Pi_OD[od][node] + (1 / (double)iter_no) * (double)(Pi_OD_aux[od][node] - Pi_OD[od][node]);
		}
	}

	
	void GAP(Assignment& assignment, double& gap_each_process, float** adj_matrix) {
		gap_each_process = 0;
		for (int od = 0; od < demand_volume_vector.size(); od++) {
			int origin_node = origin_node_vector[od];
			int destination_node = destination_node_vector[od];
			int demand_volume = demand_volume_vector[od];

			int local_node_numbers;
			int local_link_numbers;
			int* network_cur;
			int* node_list;		

			Dial_subnetwork(origin_node, destination_node, demand_volume, adj_matrix, \
				network_cur, node_list, local_link_numbers, local_node_numbers);

			float* node_utility = new float[local_node_numbers];
			float* link_sys_utility = new float[local_link_numbers];
			double* link_flow_cur_od = new double[local_link_numbers];

			for (int i = 0; i < local_link_numbers; i++) {
				link_sys_utility[i] = 0;
				link_flow_cur_od[i] = link_flow_OD[od][i];
			}

			for (int i = 0; i < local_node_numbers; i++) {
				node_utility[i] = 0;
			}

			if (network_cur == NULL) {
				break;
			}

			Utility(od, origin_node, destination_node, network_cur, node_utility, \
				link_sys_utility, local_node_numbers, local_link_numbers, node_list);

			for (int i = 0; i < local_link_numbers; i++) {
				if (link_flow_cur_od[i] > 0.01) {
					int from_node = g_link_vector[network_cur[i]].from_node_seq_no;
					int from_node_index = -1;
					for (int j = 0; j < local_node_numbers; j++) {
						if (node_list[j] == from_node) {
							from_node_index = j;
						}
					}
					double gap_abs = 0.0;
					gap_abs = log((double)link_flow_cur_od[i]) - (1 / (double)mu) * (double)link_sys_utility[i] - Pi_OD[od][from_node_index];
					gap_each_process += pow(gap_abs, 2);
				}
			}
			delete[] node_utility;
			delete[] link_sys_utility;
			delete[] link_flow_cur_od;
			delete[] network_cur;
			delete[] node_list;
		}
	}


	void RLSUE_main(Assignment& assignment, float** adj_matrix, int iter_no) {
		for (int od = 0; od < demand_volume_vector.size(); od++) {
			int origin_node = origin_node_vector[od];
			int destination_node = destination_node_vector[od];
			int demand_volume = demand_volume_vector[od];

			int* network_cur;
			int* node_list;
			int local_link_numbers;
			int local_node_numbers;

			//1. Dial_subnetwork
			Dial_subnetwork(origin_node, destination_node, demand_volume, adj_matrix, \
				network_cur, node_list, local_link_numbers, local_node_numbers);

			//2.RLSUE
			RLSUE(assignment, od, iter_no, \
				network_cur, node_list, local_link_numbers, local_node_numbers);

			//3.MAXUtility
			MAXUtility(assignment, od, iter_no, \
				network_cur, node_list, local_link_numbers, local_node_numbers);

			//4.Update_link_flow_OD
			Update_link_flow_OD(assignment, od, iter_no, \
				node_list, local_node_numbers, network_cur, local_link_numbers);

			//5.Update_Pi_OD
			Update_Pi_OD(assignment, od, iter_no, local_node_numbers);


			delete[] node_list;
			delete[] network_cur;
		}
	};
};

std::vector<NetworkForSUE*> g_NetworkForSUE_vector;

void g_assign_computing_tasks_to_memory_blocks_SUE(Assignment& assignment) {
	//assign_memory_block only save public variables
	NetworkForSUE* PointerMatrxSUE[_MAX_MEMORY_BLOCKS];

	int i = 0;
	for (int from_zone_seq_no = 0; from_zone_seq_no < assignment.g_number_of_zones; from_zone_seq_no++){
		for (int to_zone_seq_no = 0; to_zone_seq_no < assignment.g_number_of_zones; to_zone_seq_no++)
			for (int agent_type_no = 0; agent_type_no < assignment.g_number_of_agent_types; agent_type_no++)
				for (int demand_period_no = 0; demand_period_no < assignment.g_number_of_demand_periods; demand_period_no++)
				{
					if (assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume != NULL)
					{
						if (i < assignment.g_number_of_memory_blocks) {
							NetworkForSUE* p_NetworkForSUE = new NetworkForSUE();
							p_NetworkForSUE->origin_node_vector.push_back(from_zone_seq_no);
							p_NetworkForSUE->destination_node_vector.push_back(to_zone_seq_no);
							p_NetworkForSUE->demand_volume_vector.push_back(assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume);
							p_NetworkForSUE->AllocateMemorySUE(assignment.g_number_of_nodes, assignment.g_number_of_links);
							PointerMatrxSUE[i] = p_NetworkForSUE;
							g_NetworkForSUE_vector.push_back(p_NetworkForSUE);
							i++;
						}
						else {
							int memory_block_no = i % assignment.g_number_of_memory_blocks;
							NetworkForSUE* p_NetworkForSUE = PointerMatrxSUE[memory_block_no];
							p_NetworkForSUE->origin_node_vector.push_back(from_zone_seq_no);
							p_NetworkForSUE->destination_node_vector.push_back(to_zone_seq_no);
							p_NetworkForSUE->demand_volume_vector.push_back(assignment.g_column_pool[from_zone_seq_no][to_zone_seq_no][agent_type_no][demand_period_no].od_volume);
							p_NetworkForSUE->AllocateMemorySUE(assignment.g_number_of_nodes, assignment.g_number_of_links);
							i++;
						}
					}
				}
	}
}

void g_output_RLSUE_result(Assignment& assignment)
{
	cout << "writing link_performance_RLSUE.csv.." << endl;
	int b_debug_detail_flag = 0;
	FILE* g_pFileLinkMOE = NULL;
	fopen_ss(&g_pFileLinkMOE, "link_performance_RLSUE.csv", "w");
	if (g_pFileLinkMOE == NULL)
	{
		cout << "File link_performance.csv cannot be opened." << endl;
		g_ProgramStop();
	}
	else
	{
	fprintf(g_pFileLinkMOE, "road_link_id,from_node_id,to_node_id,volume,travel_time");
	fprintf(g_pFileLinkMOE, "\n");

	for (int j = 0; j < assignment.g_number_of_links; j++) {
		fprintf(g_pFileLinkMOE, "%s,%d,%d,%.3f,%3f",
			g_link_vector[j].link_id.c_str(),
			g_node_vector[g_link_vector[j].from_node_seq_no].node_id,
			g_node_vector[g_link_vector[j].to_node_seq_no].node_id,
			g_link_vector[j].flow_volume_per_period[0],
			g_link_vector[j].travel_time_per_period[0]);
		fprintf(g_pFileLinkMOE, "\n");
	}
	fclose(g_pFileLinkMOE);
	}
}





//Main function for network assignment
double network_assignment(int iteration_number, int assignment_mode, int column_updating_iterations)
{
	clock_t startTime, endTime;
	startTime = clock();


	if (assignment.g_pFileDebugLog == NULL)
		fopen_ss(&assignment.g_pFileDebugLog, "STALite_log.txt", "w");

	assignment.assignment_mode = assignment_mode;


	// step 1: read input data of network / demand tables / toll
	g_ReadInputData(assignment);
	g_ReadDemandFileBasedOnDemandFileList(assignment);

	cout << "*********************************************************************" << endl;
	cout << "Starting Stochastic User Equilibrium Traffic Assignment Process......" << endl;

	//Step 2: Dial subnetwork
	float** adj_matrix = AllocateDynamicArray<float>(assignment.g_number_of_nodes, assignment.g_number_of_nodes);
	for (int i = 0; i < assignment.g_number_of_nodes; i++) {
		for (int j = 0; j < assignment.g_number_of_nodes; j++)
		{
			adj_matrix[i][j] = 65535;
		}
	}
	for (int i = 0; i < assignment.g_number_of_links; i++) {
		adj_matrix[g_link_vector[i].from_node_seq_no][g_link_vector[i].to_node_seq_no] = g_link_vector[i].length;
	}
	//Step 2.1: Fylod SPP
	for (int u = 0; u < assignment.g_number_of_nodes; ++u) {
		for (int v = 0; v < assignment.g_number_of_nodes; ++v) {
			for (int w = 0; w < assignment.g_number_of_nodes; ++w) {
				if (adj_matrix[v][w] > adj_matrix[v][u] + adj_matrix[u][w]) {
					adj_matrix[v][w] = adj_matrix[v][u] + adj_matrix[u][w];
				}
			}
		}
		adj_matrix[u][u] = 0;
	}
	for (int i = 0; i < assignment.g_number_of_nodes; i++) {
		for (int j = 0; j < assignment.g_number_of_nodes; j++) {
			if (adj_matrix[j][i] < 65535)
				adj_matrix[i][j] = adj_matrix[j][i];
			if (adj_matrix[i][j] < 65535)
				adj_matrix[j][i] = adj_matrix[i][j];
		}
	}

	//Step 2.2 assign computing tasks
	g_assign_computing_tasks_to_memory_blocks_SUE(assignment);

	//Reset link voulme
	for (int link = 0; link < g_link_vector.size(); link++) {
		g_link_vector[link].flow_volume_per_period[0] = 0;
	}

	//Step 2.3 Dial initialization
#pragma omp parallel for 
	for (int ProcessID = 0; ProcessID < g_NetworkForSUE_vector.size(); ProcessID++) {
		g_NetworkForSUE_vector[ProcessID]->Dial_subnetwork_init(adj_matrix);
	}
	update_link_travel_time_and_cost();


	//Step 3: Recursive 
	for (int iter_no = 1; iter_no <= iteration_number; iter_no++) {
		clock_t RL_startTime, RL_endTime;
		RL_startTime = clock();

		cout << "iteration " << iter_no << endl;
		for (int i = 0; i < g_link_vector.size(); i++) {
			g_link_vector[i].flow_volume_per_period[0] = 0.0;
		}

#pragma omp parallel for 
		for (int ProcessID = 0; ProcessID < g_NetworkForSUE_vector.size(); ProcessID++) {
			g_NetworkForSUE_vector[ProcessID]->RLSUE_main(assignment, adj_matrix, iter_no);
		}

		update_link_travel_time_and_cost();


	//Step 4: Compute GAP
		double* gap_for_each_od = new double[int(g_NetworkForSUE_vector.size())]();
#pragma omp parallel for 
		for (int ProcessID = 0; ProcessID < g_NetworkForSUE_vector.size(); ProcessID++) {
			g_NetworkForSUE_vector[ProcessID]->GAP(assignment, gap_for_each_od[ProcessID], adj_matrix);
		}
		double gap_total = 0.0;
		for (int ProcessID = 0; ProcessID < g_NetworkForSUE_vector.size(); ProcessID++) {
			gap_total += gap_for_each_od[ProcessID];
		}
		cout << "GAP: " << gap_total / 2 << endl;
		delete[] gap_for_each_od;
		
		RL_endTime = clock();
		cout << "time: " << (double)(RL_endTime - RL_startTime) / CLOCKS_PER_SEC << "s" << endl;
	}

	endTime = clock();
	cout << "The run time is: " << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;

	g_ProgramStop();
	return 1;
}
