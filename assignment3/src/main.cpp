#include "mpi.h"
#include <stdio.h>
#include <unistd.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>

const int size = 1;
const int size_data_received = 3;

void output_data(std::string data_s, int rank) {
  std::stringstream file_name;
  file_name << "output_MPI_Gather_" << rank << ".log";
  char file_name_buf[1024];
  strcpy(file_name_buf, file_name.str().c_str());

  MPI_File mpi_file;
  MPI_Status status;
  MPI_File_open(MPI_COMM_SELF, file_name_buf, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpi_file);
  MPI_File_write(mpi_file, data_s.c_str(), strlen(data_s.c_str()), MPI_CHAR, &status);
  MPI_File_close(&mpi_file);
}

int main(int argc, char* argv[]) {
  int rank, n_procs;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::vector<float> data_f(size, 1.9);
  if (rank == 0) data_f[0] = 0.2;

  std::vector<float> data_received(size_data_received, 0.0);

  MPI_Gather(data_f.data(), size, MPI_FLOAT, data_received.data(), size, MPI_FLOAT, 0, MPI_COMM_WORLD);

  std::stringstream data_ss;
  for (float value : data_received) {
    data_ss << value << " ";
  }
  output_data(data_ss.str(), rank);

  MPI_Finalize();
  return 0;
}
