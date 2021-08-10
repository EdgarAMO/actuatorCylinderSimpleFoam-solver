// scalar end_time( readScalar( runTime.controlDict().lookup("endTime") ) );
// scalar current_time = runTime.value();

if (1) {
    
    #include <fstream>
    std::ofstream file("results");

    for (int i = 0; i < noc; i++){

        file << cylinders[i].Cp << "\n";

    }

    file.close();
}
