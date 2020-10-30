#include "stdafx.h"
#include "C:\Users\Entai Wang\Desktop\RLSUE_C++\RecursiveLogit_SUE\AgentLite\main_api.cpp"
//#include "C:\Users\Entai Wang\Desktop\RLSUE_C++\RecursiveLogit_SUE\AgentLite\main_api_2.cpp"

//Main function
int main(int argc, TCHAR* argv[], TCHAR* envp[])
{
	//Initialize the assignment iteration, column pool optimization iterations and assignment mode
	//Mode 0: UE; mode 1:SO; mode 2:UE+resource constrain
	int iteration_number = 10;
	int column_updating_iterations = 0;
	int assignment_mode = 0;


	//Read setting.csv and set the attribute parameter
	CCSVParser parser_settings;
	if (parser_settings.OpenCSVFile("settings.csv", true))
	{

		while (parser_settings.ReadRecord())
		{
			string field;
			int value_int;

			parser_settings.GetValueByFieldName("field", field);
			parser_settings.GetValueByFieldName("value", value_int);

			if (field == "number_of_iterations")
				iteration_number = value_int;

			if (field == "assignment_mode")
				assignment_mode = value_int;

			if (field == "column_updating_iterations")
				column_updating_iterations = value_int;

		}
	}

	assignment_mode = 5; //Using RLSUE mode
	//Main sub_function
	network_assignment(iteration_number, assignment_mode, column_updating_iterations);
}
