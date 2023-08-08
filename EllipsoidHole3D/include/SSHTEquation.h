#pragma once

#include "SSHTNodeData.h"
#include <talyfem/talyfem.h>
#include "SSHTAnalyticSol.h"

//https://math.stackexchange.com/questions/2154803/weak-formulation-poisson-equation

class SSHTEquation : public CEquation<SSHTNodeData> {
 public:
    // indices to timing arrays. These are locations in the timers_ array that
    // correspond to specific stages of the code that we wish to time.
    static const int kTimerSolve = 0;  ///< index to timer for solve process
    static const int kTimerAssemble = 1;  ///< index to timer for assemble
    static const int kTimerKSPSolve = 2;  ///< index to timer for KSPsolve
    static const int kTimerUpdate = 3;  ///< index to timer for update process
    static const int kTimerMeasure1 = 4;  ///< index to timer for KSPsolve
    static const int kTimerMeasure2 = 5;
    //int voltPlate = idata_->innerBs.voltPlate;
    //int groundedPlate = idata_->innerBs.groundedPlate;
    //int voltPlate = 1;
    //int groundedPlate = 2;
    double xc;
    double yc;
    double zc;
    double theta;
    double phi;
    double a;
    double b;
    double c;

    SSHTEquation(SSHTAnalyticSolution* sol, SSHTInputData* idata) : analytic_sol_(sol), idata_(idata) {

        /*
        timers_[kTimerSolve].set_label("Solve");
        timers_[kTimerAssemble].set_label("Assemble");
        timers_[kTimerKSPSolve].set_label("KSPSolve");
        timers_[kTimerUpdate].set_label("Update");
        timers_[kTimerMeasure1].set_label("Measure 1");
        timers_[kTimerMeasure2].set_label("Measure 2");
         */
        //voltPlate = ancPlate;
        //groundedPlate = grdPlate;
    }

    ~SSHTEquation() override {
        /*
        timers_[kTimerSolve].PrintGlobalAverageSeconds();
        timers_[kTimerAssemble].PrintGlobalAverageSeconds();
        timers_[kTimerKSPSolve].PrintGlobalAverageSeconds();
        timers_[kTimerUpdate].PrintGlobalAverageSeconds();
        timers_[kTimerMeasure1].PrintGlobalAverageSeconds();
        timers_[kTimerMeasure2].PrintGlobalAverageSeconds();
         */
    }


  void setBC(const std::map<int, SSHTInputData::BoundaryCondition> &bc) {
    bc_ = bc;
  }

    double getPermittivity(const ZEROPTV& pt , SSHTInputData::Epsilon eps){

        const ZEROPTV centerVector = ZEROPTV(xc, yc, zc); //Read the centers of the hole
        ZEROPTV constDisplacement = pt - centerVector; //Displace

        //Rotation Matrix Creation
        Matrix<3> ROTZ;
        ROTZ.data = {cos(theta),-1*sin(theta),0.0,sin(theta),cos(theta),0.0,0.0,0.0,1.0};

        Matrix<3> ROTY;
        //ROTX.data = {1.0, 0.0, 0.0, 0.0, cos(phi),-1*sin(phi),0.0,sin(phi),cos(phi)};
        ROTY.data = {cos(phi), 0.0, sin(phi), 0.0, 1.0, 0.0, -1*sin(phi), 0.0, cos(phi)};

        //Create the Volume matrix
        //This diagonal matrix defines the size of the ellipsoid cavity
        Matrix<3> LAMBDA;
        LAMBDA.data = {pow(a,-2),0,0,0,pow(b, -2),0,0,0,pow(c,-2)};

        //Perform Rotation
        const ZEROPTV rotVectZ = ROTZ*constDisplacement;
        const ZEROPTV rotVectY = ROTY*rotVectZ;

        //Compute the Quadratic Form
        //Diagonal matrix Q = Q^T. Hence can write x^TQx as ((Qx)^T)x
        const double check = (LAMBDA*rotVectY).innerProduct(rotVectY);

        //Check if pt is within cavity
        if(check <= 1.0){
            return eps.kappa;
        } else{
            return 1.00;
        }
    }

    /*
    double EpsDerivative(const ZEROPTV& pt, SSHTInputData::Epsilon eps, std::string type){
        if(type == "x"){
            double xPrime = pt.x() - eps.xCenter;
            double yPrime = pt.y() - eps.yCenter;
            double negDist = pow(pow(xPrime,2) + pow(yPrime,2), 0.5) - eps.R;
            return 0.5*(1.00 - eps.kappa)*(1.00/eps.delta)*(1.00 - pow(tanh(negDist/eps.delta),2))*(xPrime/(pow(pow(xPrime,2) + pow(yPrime,2), 0.5)));
        } else {
            double xPrime = pt.x() - eps.xCenter;
            double yPrime = pt.y() - eps.yCenter;
            double negDist = pow(pow(xPrime,2) + pow(yPrime,2), 0.5) - eps.R;
            return 0.5*(1.00 - eps.kappa)*(1.0/eps.delta)*(1.00 - pow(tanh(negDist/eps.delta),2))*(yPrime/(pow(pow(xPrime,2) + pow(yPrime,2), 0.5)));
        }
    }

    double MMSForceCalc(const ZEROPTV& pt, SSHTInputData::Epsilon eps, double derivX, double derivY){
        return -1*(derivX*M_PI*cos(M_PI*pt.x())*sin(M_PI*pt.y())) - (derivY*M_PI*sin(M_PI*pt.x())*cos(M_PI*pt.y())) +
               2*(contPermittivity(pt,eps)*(M_PI*M_PI*sin(M_PI*pt.x())*sin(M_PI*pt.y())));
    }
    */
  //virtual void Solve(double delta_t, double current_time) override {
  virtual void Solve(double delta_t, double current_time) override{
      //Pass Dirichlet Plates to fillEssBC
      //Set the Neumann Plates in the neumann struct

      //neumann.neumPlates = neumPlates;
    //timers_[kTimerSolve].Start();

    fillEssBC();

    //timers_[kTimerAssemble].Start();
    // Assemble Ae matrix and be vector
//    makeSystem();
      ApplyEssBCToSolution();
      Assemble();
      ApplyEssBC();
    //timers_[kTimerAssemble].Stop();
    // Calculate the solution vector by using the Krylov subspace Solve

    //timers_[kTimerKSPSolve].Start();
    SolveKSP(solution_, 1, 0);
    //timers_[kTimerKSPSolve].Stop();

    // The result of the system solve is in the solution_ vector. We want to
    // store this data in our NodeData arrays so we can save and load the data
    // later. This helper function puts the data from solution_  into location 0
    // of each NodeData object.
    //p_data_->NodeDataFromArray(solution_.data(), 0);
    //timers_[kTimerUpdate].Start();

    p_data_->NodeDataFromArray(solution_.data(), VOLTAGE_IDX, n_dof());

    //timers_[kTimerUpdate].Stop();
    //timers_[kTimerSolve].Stop();
  }

  // Specify the Dirichlet boundary
  //virtual void fillEssBC() override {
  virtual void fillEssBC() override{
      //boundary_conditions_.DeleteDirichletConditions();
      //has_ess_bc_.fill(false);
      const SSHTInputData::Epsilon eps = idata_->innerBs;

    initEssBC();

    bool found_dirichlet = false;

    // loop over all the nodes...
    for (LocalNodeID node_id = 0; node_id < p_grid_->n_nodes(); node_id++) {
      // loop over all the boundary conditions
      // it->first = the boundary ID (input to BoNode)
      // it->second = the boundary type (DIRICHLET/NEUMANN)

      //Set all plates to their respective MMS.

        double val = 0;

        if (p_grid_->BoNode(node_id, eps.voltPlate )) { //Apply MMS on the electrodes
            const ZEROPTV& pt = p_grid_->GetNode(node_id)->location();
            //specifyValue(node_id, VOLTAGE_IDX, idata_->innerBs.voltage);
            //val = analytic_sol_->calc_u_at(pt);
            //specifyValue(node_id, VOLTAGE_IDX, idata_->innerBs.voltage);
            specifyValue(node_id, VOLTAGE_IDX, eps.voltage);
            //std::cout << "PLATE VOLT = " << eps.voltPlate << std::endl;
            found_dirichlet = true;
        }

        if (p_grid_->BoNode(node_id, eps.groundedPlate)) {
            const ZEROPTV& pt = p_grid_->GetNode(node_id)->location();
            //val = analytic_sol_->calc_u_at(pt);
            specifyValue(node_id, VOLTAGE_IDX, 0);
            //std::cout << "GROUND VOLT = " << eps.groundedPlate << std::endl;
            found_dirichlet = true;
        }
    }

    // check if any process set found_dirichlet to true
    bool any_dirichlet = false;
    MPI_Allreduce(&found_dirichlet, &any_dirichlet, 1,
                  MPI_CXX_BOOL, MPI_LOR, PETSC_COMM_WORLD);
    if (!any_dirichlet) {
      PrintWarning("No dirichlet boundaries were set. Solution is non-unique.");
      PrintWarning("If you are loading from a file, make sure ifLoadNodeIndicators is set ",
                   "and your mesh has the appropriate node indicators (1-6).");
    }
  }

  virtual void Integrands(const FEMElm &fe, ZeroMatrix<double> &Ae,
                          ZEROARRAY<double> &be) override {
    const int n_dimensions = fe.nsd();      // # of dimensions: 1D, 2D, or 3D
    int n_basis_functions = fe.nbf();  // # of basis functions
    const double detJxW = fe.detJxW();      // (determinant of J) cross W

    //Obtain the vector position of the Finite Element
    const ZEROPTV p = fe.position();

    //Obtain the Epsilon value from config.txt
    const SSHTInputData::Epsilon eps = idata_->innerBs; //, S eps

    double force = 0; //No forcing function in the Strong Form

     double permit = getPermittivity(p, eps);

    //Assembly
    for (int a = 0; a < n_basis_functions; a++) {
      for (int b = 0; b < n_basis_functions; b++) {
        double N = 0;
        for (int k = 0; k < n_dimensions; k++) {
          N += permit * fe.dN(a, k) * fe.dN(b, k) * detJxW;
        }
        Ae(a, b) += N;
      }
      be(a) += fe.N(a) * (force) * detJxW;
    }
  }

    //Generate a row of the Capacitance Matrix
    double FluxPlate(const int _indicator) {
        //Calculate capacitance with respect to a given plate with a positive voltage.
        //std::vector<double> capacitanceRow;
        //Read all configurations from config.txt
        double fluxPlate = 0;
        FEMElm fe(p_grid_, BASIS_FIRST_DERIVATIVE | BASIS_POSITION | BASIS_DIMENSION_REDUCTION);
        for (int elmid = 0; elmid < p_grid_->n_elements(); elmid++) {
            //Refill with this Element
            fe.refill(elmid, 0);
            //Get Surface
            ELEM::SurfaceList_type::const_iterator it;
            for (it = fe.elem()->surface_indicator_.begin();
                 it != fe.elem()->surface_indicator_.end(); it++) {
                //ID = 1,2,3,4,5,6 is Set in gmsh for Dirichlet Boundaries
                if (it->has_indicator(_indicator)) {
                    //Refill with Surface Element for surface gauss point iteration
                    fe.refill_surface(elmid, &*it, rel_order_);
                    //Get Surface Normal vector
                    ZEROPTV surfaceNormal = fe.surface()->normal();
                    //double surfFlux = 0;

                    //Iterate through each Surface Gauss Point
                    while (fe.next_itg_pt()) {
                        //Get detJxW for the Gauss Point
                        const double detJxW = fe.detJxW();
                        //Calculate the flux in each direction
                        ZEROPTV dU;
                        for (int k = 0; k < fe.nsd(); k++) {
                            dU(k) = p_data_ -> valueDerivativeFEM(fe,0,k); //For Flux kth dir
                        }

                        //Calculate the gradient vector
                        fluxPlate  += dU.innerProduct(surfaceNormal)*detJxW;
                    }
                    //Calculate the Surface Flux
                    //double mutualCapacitance = (1/eps.voltage)*(totalOutGoingFluxVector.innerProduct(surfaceNormal));
                    //Append to the Capacitance Row Vector
                    //capacitanceRow.push_back(mutualCapacitance);
                }
            } //End of Surface
        } // End Element Loop Iteration
        //return capacitanceRow;
        //return fluxPlate;
        //Output the flux
        return fluxPlate;
    } //End of Method

 private:
  // map side index to boundary condition type
  std::map<int, SSHTInputData::BoundaryCondition> bc_;
  SSHTAnalyticSolution* analytic_sol_;
  SSHTInputData* idata_;
  std::vector<double> capacitanceMatrix;
  //MPITimer timers_[6];  ///< for timing several parts of the code
};
