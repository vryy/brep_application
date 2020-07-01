//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         brep_application/LICENSE.txt
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//  Date:            25 Jun 2020
//


#if !defined(KRATOS_HYDROSTATIC_PRESSURE_FUNCTION_ON_SURFACE_H_INCLUDED )
#define  KRATOS_HYDROSTATIC_PRESSURE_FUNCTION_ON_SURFACE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "containers/array_1d.h"
#include "geometries/geometry_data.h"
#include "custom_algebra/function/function.h"
#include "custom_algebra/function/zero_function.h"
#include "custom_algebra/function/scalar_function.h"
#include "custom_algebra/volume/parametric_volume.h"


namespace Kratos
{
///@addtogroup BRepApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Representing a hydrostatic pressure on the surface
*/
class HydrostaticPressureFunctionOnSurface : public FunctionR3R3
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of HydrostaticPressureFunctionOnSurface
    KRATOS_CLASS_POINTER_DEFINITION(HydrostaticPressureFunctionOnSurface);

    typedef FunctionR3R3 BaseType;

    typedef BaseType::InputType InputType;

    typedef BaseType::OutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    HydrostaticPressureFunctionOnSurface(const double& Pressure, const double& GradientPressure, const array_1d<double, 2>& Normal)
    : BaseType(), mPressure(Pressure), mGradientPressure(GradientPressure), mNormal(Normal / norm_2(Normal))
    {}

    // HydrostaticPressureFunctionOnSurface(const double& Pressure, const double& GradientPressure, const Vector& Normal)
    // : BaseType(), mPressure(Pressure), mGradientPressure(GradientPressure), mNormal(Normal)
    // {
    //     if (Normal.size() < 3)
    //         KRATOS_THROW_ERROR(std::logic_error, "The normal vector must have at least 3 components", "")
    // }

    /// Copy constructor.
    HydrostaticPressureFunctionOnSurface(HydrostaticPressureFunctionOnSurface const& rOther)
    : BaseType(rOther), mPressure(rOther.mPressure), mGradientPressure(mGradientPressure), mNormal(rOther.mNormal)
    {}

    /// Destructor.
    virtual ~HydrostaticPressureFunctionOnSurface() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /// inherit from Function
    virtual BaseType::Pointer CloneFunction() const
    {
        return BaseType::Pointer(new HydrostaticPressureFunctionOnSurface(*this));
    }


    /// inherit from Function
    virtual OutputType GetValue(const InputType& T) const
    {
        OutputType P;

        double z = T[2];
        double p = mPressure + mGradientPressure*z;

        P[0] = p*mNormal[0];
        P[1] = p*mNormal[1];
        P[2] = p*mNormal[2];

        return P;
    }


    /// inherit from Function
    virtual BaseType::Pointer GetDiffFunction(const int& component) const
    {
        if (component == 0)
        {
            return BaseType::Pointer(
                        new ParametricVolume(
                            ZeroFunction<FunctionR3R1>::Create(),
                            ZeroFunction<FunctionR3R1>::Create(),
                            ZeroFunction<FunctionR3R1>::Create()
                        )
                    );
        }
        else if (component == 1)
        {
            return BaseType::Pointer(
                        new ParametricVolume(
                            ZeroFunction<FunctionR3R1>::Create(),
                            ZeroFunction<FunctionR3R1>::Create(),
                            ZeroFunction<FunctionR3R1>::Create()
                        )
                    );
        }
        else if (component == 2)
        {
            return BaseType::Pointer(
                        new ParametricVolume(
                            ZeroFunction<FunctionR3R1>::Create(),
                            ZeroFunction<FunctionR3R1>::Create(),
                            ScalarFunction<FunctionR3R1>::Create(mGradientPressure)
                        )
                    );
        }
    }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "hydrostatic Pressure Function on Surface";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "P: " << mPressure << std::endl;
        rOStream << "dP: " << mGradientPressure << std::endl;
        rOStream << "N: " << mNormal << std::endl;
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    double mPressure;
    double mGradientPressure;
    array_1d<double, 3> mNormal;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    HydrostaticPressureFunctionOnSurface& operator=(HydrostaticPressureFunctionOnSurface const& rOther);

    ///@}

}; // Class HydrostaticPressureFunctionOnSurface

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, HydrostaticPressureFunctionOnSurface& rThis)
{}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const HydrostaticPressureFunctionOnSurface& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_HYDROSTATIC_PRESSURE_FUNCTION_ON_SURFACE_H_INCLUDED  defined
