#include <talyfem/talyfem.h>

#include <SSHTInputData.h>
#include <SSHTAnalyticSol.h>
#include <SSHTEquation.h>
#include <SSHTNodeData.h>
#include <fstream>

int main(int argc, char**args)
{

  PetscInitialize(&argc, &args, nullptr, nullptr);
  // This try/catch statement has two uses:
  // 1. It reports errors instead of just aborting when an exception is thrown.
  // 2. It makes sure MyEquation goes out of scope before PetscFinalize()
  //    is called. Otherwise, the destructors associated with the solver may
  //    call something PETSc-related, which causes a crash.

  std::ofstream writer;

    try
    {

        SSHTInputData input_data;

        if (!input_data.ReadFromFile()) {
            throw TALYException() << "Error reading input data, check config.txt!";
        }

        if (!input_data.CheckInputData()) {
            throw TALYException() << "Invalid input data, check config.txt!";
        }

        // Create a GRID depending on options set in config.txt.
        GRID *p_grid = nullptr;             // Declaring pointer of type GRID to store grid data
        CreateGrid(p_grid, &input_data); // fills in p_grid with a GRID based on input_data
        SSHTAnalyticSolution analytic_sol(input_data.analytic_sol_type);

        // Create our equation object.
        SSHTEquation equation(&analytic_sol, &input_data);

        equation.add_basis_flag(BASIS_DIMENSION_REDUCTION);  // support 1D/2D in 3D elements
        equation.SetPreallocator(new PreallocatorPerfect<SSHTNodeData>(&equation));

        // This allocates the matrix and vector (stored as part of the equation object).
        int n_dof = 1;  // number of degrees of freedom
        equation.redimSolver(p_grid, n_dof, false, input_data.basisRelativeOrder); // resize the equation to match the grid structure

        // Create our GridField.
        //The GridField holds the data for each node in our mesh.
        GridField<SSHTNodeData> grid_field;
        grid_field.redimGrid(p_grid);  // tell the gridfield where to find the grid (for number of nodes)
        grid_field.redimNodeData();    // initialize the node data
        equation.setData(&grid_field); // tell the equation where to get node data

        equation.setBC(input_data.boundary_conditions);

        // Only for transient problems //
        // Solve a single timestep of our equation by calling MyEquation::Solve().
        double t = 0.0;
        double dt = 0.0;

        //Get the Size of the Problem i.e. #of rows in the settings file
        int settingsSize = (int)input_data.parsedCsv.size();

        //Store the flux calculations
        std::vector<double> capacitances;

        for(int i = 0; i < settingsSize; i++){
                equation.xc = input_data.parsedCsv[i][0];
                equation.yc = input_data.parsedCsv[i][1];
                equation.zc = input_data.parsedCsv[i][2];
                equation.theta = input_data.parsedCsv[i][3];
                equation.phi = input_data.parsedCsv[i][4];
                equation.a = input_data.parsedCsv[i][5];
                equation.b = input_data.parsedCsv[i][6];
                equation.c = input_data.parsedCsv[i][7];
               
		//Solve for the settings
                equation.Solve(t, dt);

                // double flux = equation.FluxPlate(input_data.innerBs.groundedPlate);
                //std::cout << "Flux Calculation = " << flux << "\n" ;
                // capacitances.push_back(flux);
		
		int rank;
    		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    		if (rank == 0)
    		{
		
		double flux = equation.FluxPlate(input_data.innerBs.groundedPlate);
                std::cout << "Flux Calculation = " << flux << "\n" ;
                //capacitances.push_back(flux);
			
                //Write One Line
                writer.open("result.txt", std::ofstream::app);
                writer << input_data.parsedCsv[i][0] << "," << input_data.parsedCsv[i][1] << "," << input_data.parsedCsv[i][2] << "," << input_data.parsedCsv[i][3] << "," << input_data.parsedCsv[i][4] << "," << input_data.parsedCsv[i][5] << "," << input_data.parsedCsv[i][6] << "," << input_data.parsedCsv[i][7] << "," << flux << "\n";
                writer.close();

		}
        }

        //equation.Solve(t, dt);
        //Open the log file for writing
       	//writer.open("result.txt", std::ofstream::app);

        /*
        //Write Input Params + calculated capacitance
        for(int numWrite = 0; numWrite < settingsSize; numWrite++){
            writer << input_data.parsedCsv[numWrite][0] << "," << input_data.parsedCsv[numWrite][1] << "," << input_data.parsedCsv[numWrite][2] << "," << input_data.parsedCsv[numWrite][3] << "," << input_data.parsedCsv[numWrite][4] << "," << input_data.parsedCsv[numWrite][5] << "," << input_data.parsedCsv[numWrite][6] << "," << input_data.parsedCsv[numWrite][7] << "," << capacitances[numWrite] << "\n";
        }
        */
        //writer << "\n";
        //writer.close();
	

        // Free the GRID object.
        DestroyGrid(p_grid);
    }
    catch (TALYException &e)
    {
        PrintError(e.what());
    }
  
    catch(std::exception &e)
    {
	    std::cout << e.what();
    }
  PetscFinalize(); // clean up PETSc environment
  return 0;        // program exits normally
}
