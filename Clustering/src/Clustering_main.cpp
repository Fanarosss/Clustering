#include "Library.h"
#include "Cluster.h"
#include "Helper_Functions.h"
#include "Cluster_DataTypes.h"

using namespace std;

int main(int argc, char* argv[]){

    /* Arguments Check */
    string input_file, config_file, results_file;
    int dat = 0, conf = 0, rf = 0, complete = 0;
    for (int i = 1; i < argc; i++){
        string arg = argv[i];
        if(arg == "-i"){
            if(i+1 < argc){
                input_file = argv[++i];
                if(access( input_file.c_str(), F_OK ) == -1){
                    cout << "-- Wrong <input file> -- \n";
                    show_cluster_usage(argv[0]);
                    return -1;
                }
                dat = 1;
            }else{
                cout << "-- NO <input file> -- \n";
                show_cluster_usage(argv[0]);
                return -1;
            }
        }else if (arg == "-c"){
            if(i+1 <= argc){
                config_file = argv[++i];
                if(access( config_file.c_str(), F_OK ) == -1){
                    cout << "-- Wrong <configuration file> -- \n";
                    show_cluster_usage(argv[0]);
                    return -1;
                }
                conf = 1;
            }else{
                cout << "-- NO <configuration file> -- \n";
                show_cluster_usage(argv[0]);
                return -1;
            }
        }else if(arg == "-o"){
            if(i+1 <= argc){
                results_file = argv[++i];
                rf = 1;
            }else{
                cout << "-- NO <output file> -- \n";
                show_cluster_usage(argv[0]);
                return -1;
            }
        }else if(arg == "-complete"){
            complete = 1;
        }else{
            show_cluster_usage(argv[0]);
            return -1;
        }
    }

    /* File Check if not given as arguments */
    if (dat == 0) {
        cout << "Path to input.dat (<input file>):" << endl;
        cin >> input_file;
        while (access(input_file.c_str(), F_OK) == -1) {
            cout << "-- Wrong <input file> -- \n";
            cin >> input_file;
        }
    } else dat = 0;
    if (conf == 0) {
        cout << "Path to cluster.conf (<configuration file>):" << endl;
        cin >> config_file;
        while (access(config_file.c_str(), F_OK) == -1) {
            cout << "-- Wrong <configuration file> -- \n";
            cin >> config_file;
        }
    } else conf = 0;
    if (rf == 0) {
        cout << "Path to file of results (<output file>):" << endl;
        cin >> results_file;
    } else rf = 0;

    /* Scan the type of data */
    int file_data = Read_input_file(input_file);      // 1 for vectors - 2 for curves - 0 for error
    if(file_data == 1){
        Cluster_Vectors(input_file, config_file);
    }else if(file_data == 2){
        Cluster_Curves(input_file, config_file);
    }else{
        cout << "Input file error!" << endl;
        return -1;
    }

    return 0;
}


