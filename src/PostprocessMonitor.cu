#include "PostprocessMonitor.h"


#include "defineFilenames.h"
#include "definePrecision.h"

#include "defineLBM.h"
#include "FuncAllocate.h"
#include "Index.h"


// private //
void PostprocessMonitor::
Sendmonitordata(const int datasize)
{
	int rank; MPI_Comm_rank(mpi_communicator_,&rank);
	if  (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

	int nummonitordata = datasize;
	int sbuffsize      = nummonitordata*4;
	float* sbuff = new float[sbuffsize];

	if (rank == 0){
		for(int i=0;  i<nummonitordata; i++){
			sbuff[i*4  ] = monitor(i).x();
			sbuff[i*4+1] = monitor(i).y();
			sbuff[i*4+2] = monitor(i).z();
			sbuff[i*4+3] = monitor(i).fileno();
		}
	}

	MPI_Bcast(sbuff, sbuffsize, MPI_FLOAT, 0, mpi_communicator_);

	if (rank != 0) {
		for (int i=0;  i<datasize;  i++){
			MonitorData  data;
			data.set((int)sbuff[i*4+3], sbuff[i*4], sbuff[i*4+1], sbuff[i*4+2]);
			monitordata_.push_back(data);
		}
	}

    delete [] sbuff;
}


// public //
void PostprocessMonitor::
setupdata(
    int   rank,
    const Tree&      tree,
    const MeshValue* meshValues
    )
{
	if  (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

	int nummonitordata = 0;
    if (rank == 0) {
        readdata();
        nummonitordata = monitorcount();

        std::cout << "monitorcount = " << monitorcount() << std::endl;
    }
    MPI_Barrier(mpi_communicator_);

    MPI_Bcast((void*)&nummonitordata, 1, MPI_INT, 0, mpi_communicator_);

    Sendmonitordata(nummonitordata);
    MPI_Barrier(mpi_communicator_);

#pragma omp parallel for schedule(dynamic)
    for (int i=0;  i<monitorcount(); i++){
        monitordata_.at(i).CoordinatesToCell(tree, meshValues, rank);
    }
}


void PostprocessMonitor::
readdata()
{
	int rank; MPI_Comm_rank(mpi_communicator_,&rank);

	int i = 0;
	std::ifstream fin;
	while (1) {
		fin.open(Foldernames::input_folder+"/"+in_filename_+std::to_string(i)+".csv", std::ios::in);
		if (!fin) { break; }

        if (rank == 0) { std::cout << "read index = " << i << std::endl; }

        std::string line;
		while ( getline(fin,line)){
            std::vector<std::string> str_data = split(line,',');

			MonitorData data;
			data.set(i, std::stod(str_data.at(0)), std::stod(str_data.at(1)), std::stod(str_data.at(2)));
			monitordata_.push_back(data);
		}

		fin.close();
		i++;
	}
}


void PostprocessMonitor::
OutputMonitorData(
    int   step,
    const Tree&          tree,
    const Parameters&    parameters,
    const MeshValue*     meshValues,
    int   rank
   )
const
{
    const MPI_Comm mpi_communicator = mpi_communicator_;
//    const MPI_Comm mpi_communicator = MPI_COMM_WORLD;

    if (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    int nproc;  MPI_Comm_size(mpi_communicator, &nproc);
    int t_rank; MPI_Comm_rank(mpi_communicator, &t_rank);
//    std::cout << "nprocs, t_rank = " << nproc << ", " << t_rank << std::endl;

    int numberofmonitordata = monitorcount();
//    int _numberofmonitordata = monitorcount();
//    MPI_Allreduce(&_numberofmonitordata, &numberofmonitordata, 1, MPI_INT, MPI_MAX, mpi_communicator);

    constexpr int num_buff  = 8; // iscalc, t_rank, u,v,w, levelset_obj, scalar, T //
    const int sbuffsize = numberofmonitordata*num_buff;
    const int rbuffsize = sbuffsize*nproc;

    float* sbuff = new float[sbuffsize]; //std::fill(sbuff, sbuff+sbuffsize, 0.0f);
    float* rbuff = new float[rbuffsize]; //std::fill(rbuff, rbuff+rbuffsize, 0.0f);
//#pragma omp parallel for
//    for (int i=0; i<sbuffsize; i++) { sbuff[i] = 0.0f; }
//#pragma omp parallel for
//    for (int i=0; i<rbuffsize; i++) { rbuff[i] = 0.0f; }

    // func //
//    auto filename = [this](int id, int step)
//    {
//        MonitorData data = monitordata_.at(id);
//
//        std::string _filename = Foldernames::output_folder+"/"
//            +out_filename_+std::to_string(step)+"_"
//            +std::to_string(data.fileno())+".csv";
//
//        return _filename;
//    };

    auto filename = [](int fileno, int step, std::string out_filename)
    {
        std::string _filename = Foldernames::output_folder+"/"
            +out_filename+std::to_string(step)+"_"
            +std::to_string(fileno)+".csv";

        return _filename;
    };

#pragma omp parallel for
    for (int im=0; im<numberofmonitordata; im++){
        MonitorData data;
        data = monitordata_.at(im);
        sbuff[im*num_buff+0] = -2;
        sbuff[im*num_buff+1] = -2;

        if (t_rank == data.rank()){
            data.Approximatecalculation(tree, parameters, meshValues);

            sbuff[im*num_buff+0] = data.iscalc();
            sbuff[im*num_buff+1] = t_rank;

            sbuff[im*num_buff+2] = data.u();
            sbuff[im*num_buff+3] = data.v();
            sbuff[im*num_buff+4] = data.w();

            sbuff[im*num_buff+5] = data.levelset_obj();

            sbuff[im*num_buff+6] = data.scalar();
            sbuff[im*num_buff+7] = data.T();
        }
        else{
            sbuff[im*num_buff+2] = 0.0f;
            sbuff[im*num_buff+3] = 0.0f;
            sbuff[im*num_buff+4] = 0.0f;
            sbuff[im*num_buff+5] = 0.0f;
            sbuff[im*num_buff+6] = 0.0f;
            sbuff[im*num_buff+7] = 0.0f;
        }
    }

    if (t_rank == 0) { std::cout << "MPI_Gather : nproc, size = " << nproc << ", " << sbuffsize*nproc << std::endl; }
//    MPI_Gather(sbuff, sbuffsize, MPI_FLOAT,
//               rbuff, sbuffsize, MPI_FLOAT, 0, mpi_communicator);
    MPI_Allgather(sbuff, sbuffsize, MPI_FLOAT,
                  rbuff, sbuffsize, MPI_FLOAT, mpi_communicator);

    if (t_rank == 0) {
        std::cout << "output monitor data\n";
#pragma omp parallel for
        for (int im=0; im<numberofmonitordata; im++){
            const int fileno = monitordata_.at(im).fileno();

            std::ofstream fout;
            fout.open(filename(fileno, step, out_filename_), std::ios::out); //add mode

            fout << "step" << ","
                 << "calcflg" << ","
                 << "d_rank"  << ","
                 << "u" << "," << "v" << "," << "w" << ","
                 << "levelset_obj" << ","
                 << "scalar" << ","
                 << "T" << ","
                 << "x" << "," << "y" << "," << "z" << std::endl;

            fout.close();
        }

#pragma omp parallel for
        for (int im=0; im<numberofmonitordata; im++){
            const int fileno = monitordata_.at(im).fileno();

            std::ofstream fout;
            fout.open(filename(fileno, step, out_filename_), std::ios::app); //add mode

            for (int i=0; i<rbuffsize; i+=num_buff){
                const int _im =((int)(i/num_buff)) % numberofmonitordata;
                if (im != _im) { continue; }

                const int calcflg = (int) rbuff[i+0];
                if (calcflg < 0){ continue; }

                const int d_rank = (int) rbuff[i+1];
                if (d_rank  == -2) { continue; }

                float u = (float) rbuff[i+2];
                float v = (float) rbuff[i+3];
                float w = (float) rbuff[i+4];

                float levelset_obj = (float) rbuff[i+5];

                float scalar = (float) rbuff[i+6];
                float T      = (float) rbuff[i+7];


                MonitorData data = monitordata_.at(im);
                float x = (float) data.x();
                float y = (float) data.y();
                float z = (float) data.z();


                fout << step << ","
                     << calcflg << ","
                     << d_rank  << ","
                     << u << "," << v << "," << w << ","
                     << levelset_obj << ","
                     << scalar << ","
                     << T << ","
                     << x << "," << y << "," << z << std::endl;

            }

            fout.close();
        }
    }

    MPI_Barrier(mpi_communicator);

    delete [] sbuff;
    delete [] rbuff;
}


std::vector<std::string> PostprocessMonitor::
split(std::string& input, char delimiter)
{
	std::istringstream stream(input);
   	std::string field;
   	std::vector<std::string> result;
   	while (std::getline(stream, field, delimiter)) {
       	result.push_back(field);
   	}
   	return result;
}
