#include<iostream>
#include<sstream>
#include<string>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<vector>
#include<stdio.h>
#include<time.h>

#include<gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>            // gsl random number generators
#include <gsl/gsl_randist.h>        // gsl random number distributions
#include <gsl/gsl_vector.h>         // gsl vector and matrix definitions
#include <gsl/gsl_blas.h>           // gsl linear algebra stuff
#include <gsl/gsl_multifit_nlin.h>  // gsl multidimensional fitting

#include "H5Cpp.h"

//#ifndef H5_NO_NAMESPACE
//    using namespace H5;
//#endif

#include "CFWR.h"

using namespace std;

const H5std_string	FILE_NAME("giant_array.h5");
const H5std_string	DATASET_NAME("ga");
const int ntrig = 2;
const int local_small_array_size = eta_s_npts * qnpts * qnpts * qnpts * qnpts * ntrig;

int CorrelationFunction::Set_giant_HDF_array()
{
	//const int ntrig = 2;
	const int giant_array_size = FO_length * eta_s_npts * qnpts * qnpts * qnpts * qnpts * ntrig;
	const int giant_FOslice_array_size = eta_s_npts * qnpts * qnpts * qnpts * qnpts * ntrig;
	const int small_array_size = qnpts * qnpts * qnpts * qnpts * ntrig;

	try
    {
		Exception::dontPrint();
	
		file = new H5::H5File(FILE_NAME, H5F_ACC_TRUNC);
	
		hsize_t dims[1];               // dataset dimensions
		dims[0] = giant_array_size;
		dataspace = new H5::DataSpace (RANK, dims);
	
		dataset = new H5::DataSet( file->createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, *dataspace) );

		hsize_t offset[1], count[1], stride[1], block[1];
		hsize_t dimsm[1];
		dimsm[0] = giant_FOslice_array_size;

		count[0] = 1;
		stride[0] = 1;
		block[0] = giant_FOslice_array_size;
	
		memspace = new H5::DataSpace (RANK, dimsm, NULL);

		double * tmp_results_ii0 = new double [ntrig];
		double * tmp_results_ii1 = new double [ntrig];
		for (int isurf = 0; isurf < FO_length; ++isurf)
		{
			int igFOsa = 0;
			offset[0] = isurf;
			dataspace->selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
			double * giant_FOslice_array = new double [giant_FOslice_array_size];
			for (int ieta = 0; ieta < eta_s_npts; ++ieta)
			for (int iqt = 0; iqt < qnpts; ++iqt)
			for (int iqx = 0; iqx < qnpts; ++iqx)
			for (int iqy = 0; iqy < qnpts; ++iqy)
			for (int iqz = 0; iqz < qnpts; ++iqz)
			{
				form_trig_sign_z(isurf, ieta, iqt, iqx, iqy, iqz, 0, tmp_results_ii0);
				form_trig_sign_z(isurf, ieta, iqt, iqx, iqy, iqz, 1, tmp_results_ii1);
				giant_FOslice_array[igFOsa] = tmp_results_ii0[0] + tmp_results_ii1[0];
				giant_FOslice_array[igFOsa+1] = tmp_results_ii0[1] + tmp_results_ii1[1];
				igFOsa += 2;
			}
			dataset->write(giant_FOslice_array, PredType::NATIVE_DOUBLE, *memspace, *dataspace);
			delete giant_FOslice_array;
		}
	
		delete [] tmp_results_ii0;
		delete [] tmp_results_ii1;

		//dataspace->close();
		memspace->close();
		dataset->close();
		file->close();
		delete memspace;
		delete file;
		//delete dataspace;
		delete dataset;
		// Now set-up for remainder of calculations
		file = new H5::H5File(FILE_NAME, H5F_ACC_RDONLY);
		dataset = new H5::DataSet( file->openDataSet( DATASET_NAME ) );
		dimsm[0] = small_array_size;
		memspace = new H5::DataSpace (RANK, dimsm, NULL);
		//dataspace = new H5::DataSpace (RANK, dimsm, NULL);
    }

    catch(FileIException error)
    {
		error.printError();
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		return -1;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		return -1;
    }

	return (0);
}

int CorrelationFunction::Get_small_array_from_giant_HDF_array(int isurf, int ieta, double * small_array)
{
	//const int ntrig = 2;
	const int giant_array_size = FO_length * eta_s_npts * qnpts * qnpts * qnpts * qnpts * ntrig;
	const int small_array_size = qnpts * qnpts * qnpts * qnpts * ntrig;

	try
	{
		Exception::dontPrint();

		//hsize_t offset[1], count[1], stride[1], block[1];
	
		//offset[0] = isurf * eta_s_npts + ieta;
		//count[0] = 1;
		//stride[0] = 1;
		//block[0] = small_array_size;
		hsize_t offset[1] = {isurf * eta_s_npts + ieta};
		hsize_t count[1] = {1};
		hsize_t stride[1] = {1};
		hsize_t block[1] = {small_array_size};

		dataspace->selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);

		dataset->read(small_array, PredType::NATIVE_DOUBLE, *memspace, *dataspace);
    }

    catch(FileIException error)
    {
		error.printError();
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		return -1;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		return -1;
    }

	return (0);
}

int CorrelationFunction::Clean_up_HDF_miscellany()
{
	try
	{
		Exception::dontPrint();

		delete memspace;
		delete file;
		delete dataspace;
		delete dataset;
    }

    catch(FileIException error)
    {
		error.printError();
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		return -1;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		return -1;
    }

	return (0);
}

//*************************************************************
// Alternate versions here
//*************************************************************

int CorrelationFunction::Set_giant_chunked_HDF_array()
{
	const int giant_array_size = FO_length * eta_s_npts * qnpts * qnpts * qnpts * qnpts * ntrig;

	double local_small_array[1][local_small_array_size];

	try
    {
		Exception::dontPrint();
	
		file = new H5::H5File(FILE_NAME, H5F_ACC_TRUNC);
	
		DSetCreatPropList cparms;
		hsize_t chunk_dims[RANKV2] = {1, local_small_array_size};
		cparms.setChunk( RANKV2, chunk_dims );
	
		hsize_t dims[RANKV2];
		dims[0] = FO_length;
		dims[1] = local_small_array_size;
		dataspace = new H5::DataSpace (RANKV2, dims);
	
		dataset = new H5::DataSet( file->createDataSet(DATASET_NAME, PredType::NATIVE_DOUBLE, *dataspace, cparms) );

		hsize_t offset[RANKV2], count[RANKV2];
		hsize_t dimsm[RANKV2];
		dimsm[0] = 1;
		dimsm[1] = local_small_array_size;

		count[0] = 1;
		count[1] = local_small_array_size;
	
		memspace = new H5::DataSpace (RANKV2, dimsm, NULL);

		double * tmp_results_ii0 = new double [ntrig];
		double * tmp_results_ii1 = new double [ntrig];
		for (int isurf = 0; isurf < FO_length; ++isurf)
		{
			int igFOsa = 0;
			offset[0] = isurf;
			offset[1] = 0;
			dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);
			for (int ieta = 0; ieta < eta_s_npts; ++ieta)
			for (int iqt = 0; iqt < qnpts; ++iqt)
			for (int iqx = 0; iqx < qnpts; ++iqx)
			for (int iqy = 0; iqy < qnpts; ++iqy)
			for (int iqz = 0; iqz < qnpts; ++iqz)
			{
				form_trig_sign_z(isurf, ieta, iqt, iqx, iqy, iqz, 0, tmp_results_ii0);
				form_trig_sign_z(isurf, ieta, iqt, iqx, iqy, iqz, 1, tmp_results_ii1);
				local_small_array[0][igFOsa] = tmp_results_ii0[0] + tmp_results_ii1[0];
				local_small_array[0][igFOsa+1] = tmp_results_ii0[1] + tmp_results_ii1[1];
				igFOsa += 2;
			}
			dataset->write(local_small_array, PredType::NATIVE_DOUBLE, *memspace, *dataspace);
		}
	
		delete [] tmp_results_ii0;
		delete [] tmp_results_ii1;

		memspace->close();
		dataset->close();
		file->close();
		delete memspace;
		delete file;
		delete dataset;
		file = new H5::H5File(FILE_NAME, H5F_ACC_RDONLY);
		dataset = new H5::DataSet( file->openDataSet( DATASET_NAME ) );

		// defines some space to hold chunks
		dimsm[0] = 1;
		dimsm[1] = local_small_array_size;
		memspace = new H5::DataSpace (RANKV2, dimsm, NULL);
    }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "H5::DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "H5::DataSpaceIException error!" << endl;
		return -3;
    }

	return (0);
}

int CorrelationFunction::Get_chunk(int isurf, double local_small_array[][eta_s_npts * qnpts * qnpts * qnpts * qnpts * 2])
{
	const int giant_array_size = FO_length * eta_s_npts * qnpts * qnpts * qnpts * qnpts * ntrig;

	try
	{
		Exception::dontPrint();

		hsize_t offset[RANKV2] = {isurf, 0};
		hsize_t count[RANKV2] = {1, local_small_array_size};		// == chunk_dims
		dataspace->selectHyperslab(H5S_SELECT_SET, count, offset);

		dataset->read(local_small_array, PredType::NATIVE_DOUBLE, *memspace, *dataspace);
    }

    catch(FileIException error)
    {
		error.printError();
		cerr << "FileIException error!" << endl;
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		cerr << "H5::DataSetIException error!" << endl;
		return -2;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		cerr << "H5::DataSpaceIException error!" << endl;
		return -3;
    }

	return (0);
}

int CorrelationFunction::Reset_HDF()
{
	try
	{
		Exception::dontPrint();

		dataspace->close();
		memspace->close();
		dataset->close();
		file->close();
		delete memspace;
		delete file;
		delete dataset;
		delete dataspace;
		file = new H5::H5File(FILE_NAME, H5F_ACC_RDONLY);
		dataset = new H5::DataSet( file->openDataSet( DATASET_NAME ) );
		dataspace = new H5::DataSpace( dataset->getSpace() );

		hsize_t dimsm[RANKV2];
		dimsm[0] = 1;
		dimsm[1] = local_small_array_size;
		memspace = new H5::DataSpace (RANKV2, dimsm, NULL);  
	}

    catch(FileIException error)
    {
		error.printError();
		return -1;
    }

    catch(H5::DataSetIException error)
    {
		error.printError();
		return -1;
    }

    catch(H5::DataSpaceIException error)
    {
		error.printError();
		return -1;
    }

	return (0);
}

//End of file
