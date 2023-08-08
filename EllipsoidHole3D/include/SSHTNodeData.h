// If you've never seen #pragma once before, it works like an "include guard."
// It is functionally equivalent to:
// #ifndef MYNODEDATA_H
// #define MYNODEDATA_H
// ... code ...
// #endif
#pragma once

// Include the library.
#include <talyfem/talyfem.h>

#define VOLTAGE_IDX 0
#define VOLTAGE_EXACT_IDX 1

// This class holds the variables we are solving for at each node.
// If you have multiple degrees of freedom, you will usually have
// multiple variables in this class.
class SSHTNodeData: public NODEData
{
    public:

    // Actual storage for the nodal data
    double V;  // temperature
    double V_exact;
    double dV[2];

    // Tell the library how to get nodal data
    // The library needs a mapping from
    // value(0) == T
    double& value(int index)
    {
        switch(index)
        {
            case VOLTAGE_IDX:
                return V;
            case VOLTAGE_EXACT_IDX:
                return V_exact;
            default:
                throw TALYException() << "Invalid node data index: " << index;

        }
    }

    // Identical to the above method, but marked as const - see below.
    const double& value(int index) const 
    {
        switch(index) 
        {
            case VOLTAGE_IDX:
                return V;
            case VOLTAGE_EXACT_IDX:
                return V_exact;
            default: 
                throw TALYException() << "Invalid node data index: " << index;
        }
    }

    // Tell the library the name of each piece of nodal data
    // This determines the name that the variable will show up as in the output files (.plt).
    static const char* name(int index)
    {
        switch(index)
        {
            case VOLTAGE_IDX:
                return "Voltage";
            case VOLTAGE_EXACT_IDX:
                return "V_exact";
            default:
                throw TALYException() << "Invalid node data index: " << index;    
        }
    }

    // Tell the library how many pieces of node data we have
    static int valueno()
    {
        return 2;  // no of variables being returned two in this case as "T" and "T_exact" are being returned
    }
};

