# Basic-MD-simulation

When compiling the program don't forget to include the compiler flag -fopenmp to make the compiler use the openMP software.

Note that the program expects to find a folder Data in the install directory. No error will be thrown when this is not present but the data will not be stored.

------------------------------------------------------------------------------------------------------------------------
Note, this software will not compile for windows users. This is due to the usage of the std::system() command. Windows users should just comment this line out and make the call to Plots.gnu manually.

Also, for windows users the timing of the program might be wrong, I have not tested this.
------------------------------------------------------------------------------------------------------------------------
