#pragma once

#include <SSHTAnalyticSol.h>

class SSHTInputData : public InputData {
 public:

    std::vector<std::vector<double>> parsedCsv;
  enum BoundaryCondition {
    DIRICHLET = 0,
    NEUMANN = 1
  };

    struct Epsilon{
        double delta;
        double kappa;
        /*
        double xc;
        double yc;
        double zc;
        double a;
        double b;
        double c;
        double theta;
        double phi;
         */
        double voltage;
        int voltPlate;
        int groundedPlate;
    };

    Epsilon innerBs;

  /// map of boundary index to BC
  /// (for example, bc[-1] == DIRICHLET for left side is dirichlet)
  std::map<int, BoundaryCondition> boundary_conditions;
  SSHTAnalyticSolution::Type analytic_sol_type;
  std::string output_extension;

  SSHTInputData()
    : output_extension(".plt") {
  }
    //Loads all parameters from settings.csv
    void parseCSV(std::string fileName)
    {
        std::ifstream  data(fileName);
        std::string line;
        int _i = 0;

        while(std::getline(data,line)) {
            if (_i == 0) { _i++; continue; }
            else {
                std::stringstream lineStream(line);
                std::string cell;
                std::vector<double> parsedRow;
                while (std::getline(lineStream, cell, ',')) {
                    parsedRow.push_back(std::stod(cell));
                }

                parsedCsv.push_back(parsedRow);
            }
        }
    };

  static BoundaryCondition read_boundary(libconfig::Config& cfg, const char* name) {
    std::string str = "dirichlet";
    cfg.lookupValue(name, str);  // on missing str stays dirichlet

    if (str == "dirichlet")
      return DIRICHLET;
    else if (str == "neumann")
      return NEUMANN;

    throw TALYException() << "Unknown " << name << " boundary name: " << str;
  }

  // needs nsd to know what to resolve auto to
  static SSHTAnalyticSolution::Type read_analytic_sol_type(libconfig::Config& cfg, int nsd) {
    std::string str = "auto";
    cfg.lookupValue("analyticSolution", str);  // on missing stays auto

    // match library mesh generator (assume not loading a mesh)
    if (str == "auto") {
      if (nsd == 1) str = "line_x";
      else if (nsd == 2) str = "plane_xy";
      else if (nsd == 3) str = "cube";
      else throw NotImplementedException();
    }

    if (str == "line_x")
      return SSHTAnalyticSolution::LINE_X;
    if (str == "line_y")
      return SSHTAnalyticSolution::LINE_Y;
    if (str == "line_z")
      return SSHTAnalyticSolution::LINE_Z;
    if (str == "plane_xy")
      return SSHTAnalyticSolution::PLANE_XY;
    if (str == "plane_xz")
      return SSHTAnalyticSolution::PLANE_XZ;
    if (str == "plane_yz")
      return SSHTAnalyticSolution::PLANE_YZ;
    if (str == "cube")
      return SSHTAnalyticSolution::CUBE;
    if (str == "sphere")
      return SSHTAnalyticSolution::SPHERE;

    throw NotImplementedException() << "Unknown analytic sol '" << str << "'";
  }

  bool ReadFromFile(const std::string& filename = std::string("config.txt")) override {
      std::cout << filename;

      InputData::ReadFromFile(filename);

    boundary_conditions[1] = read_boundary(cfg, "boundaries.left");
    boundary_conditions[2] = read_boundary(cfg, "boundaries.right");
    boundary_conditions[3] = read_boundary(cfg, "boundaries.bottom");
    boundary_conditions[4] = read_boundary(cfg, "boundaries.top");
    boundary_conditions[5] = read_boundary(cfg, "boundaries.back");
    boundary_conditions[6] = read_boundary(cfg, "boundaries.front");

    ReadValue("outputExtension", output_extension);


      //Reading all values of inner boundary of permittivity

      ReadValue("kappa",innerBs.kappa);
      /*
      ReadValue("xCenter",innerBs.xc);
      ReadValue("yCenter",innerBs.yc);
      ReadValue("zCenter",innerBs.zc);
      ReadValue("A",innerBs.a);
      ReadValue("B",innerBs.b);
      ReadValue("C",innerBs.c);
      ReadValue("THETA",innerBs.theta);
      ReadValue("PHI",innerBs.phi);
      */
      ReadValue("delta",innerBs.delta);
      ReadValue("voltage", innerBs.voltage);
      ReadValue("voltPlate", innerBs.voltPlate);
      ReadValue("groundedPlate", innerBs.groundedPlate);

    analytic_sol_type = read_analytic_sol_type(cfg, this->nsd);

    parseCSV("/anvil/scratch/x-rtali/projects/Electro3D/settings/settings1.csv");

    return true;
  }
};
