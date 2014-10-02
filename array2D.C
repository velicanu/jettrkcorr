#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using array2D = std::vector< std::vector< double > >;

array2D csvtoarray(string filename)
{
    vector<vector<double> > myvector;
    std::ifstream file(filename);

    int nrow = 0, ncol = 0;
    while(true)
    {
        std::string line;
        std::getline(file, line);
        if ( !file.good() )
            break;
        std::stringstream iss(line);
        ncol = 0;
        while(true) 
        {
            std::string val;
            std::getline(iss, val, ',');
            if ( !iss.good() )
                break;
            ncol++;
        }
        nrow++;
    }
    cout<<"Opened "<<nrow<<" x "<<ncol<<" array from "<<filename<<endl;
    file.close();
    file.open(filename);
    array2D array(nrow, vector< double >(ncol));
    for (int row = 0; row < nrow; ++row)
    {
        std::string line;
        std::getline(file, line);
        if ( !file.good() )
            break;
        std::stringstream iss(line);
        for (int col = 0; col < ncol; ++col)
        {
            std::string val;
            std::getline(iss, val, ',');
            if ( !iss.good() )
                break;
            string tmp ;
            std::stringstream convertor(val);
            convertor >> tmp; 
            array[row][col] = std::stod(tmp);
            /* code */
        }
        /* code */
    }

    return array;
}

// int main()
// {
//     return 0;
// }
