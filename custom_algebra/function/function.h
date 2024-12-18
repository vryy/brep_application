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
//  Date:            13 Feb 2017
//

#if !defined(KRATOS_FUNCTION_H_INCLUDED )
#define  KRATOS_FUNCTION_H_INCLUDED

// System includes
#include <string>
#include <sstream>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "utilities/math_utils.h"

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
/** Abstract class for a general function R^m->R^n
*/
template<typename TInputType, typename TOutputType>
class Function
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Function
    KRATOS_CLASS_POINTER_DEFINITION(Function);

    typedef TInputType InputType;

    typedef TOutputType OutputType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Function() {}

    /// Copy constructor.
    Function(Function const& rOther)
    {}

    /// Destructor.
    virtual ~Function() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual Function::Pointer CloneFunction() const
    {
        return Function::Pointer(new Function(*this));
    }

    virtual inline const std::size_t InputSize() const
    {
        KRATOS_ERROR << "Call the base class";
    }

    virtual inline const std::size_t OutputSize() const
    {
        KRATOS_ERROR << "Call the base class";
    }

    virtual TOutputType GetValue(const TInputType& P) const
    {
        KRATOS_ERROR << "Call the base class";
    }

    virtual TOutputType GetDerivative(const int& component, const TInputType& P) const
    {
        Function::Pointer pDerivative = this->GetDiffFunction(component);
        return pDerivative->GetValue(P);
    }

    virtual TOutputType GetSecondDerivative(const int& component_1, const int& component_2, const TInputType& P) const
    {
        Function::Pointer pSecondDerivative = this->GetDiffFunction(component_1)->GetDiffFunction(component_2);
        return pSecondDerivative->GetValue(P);
    }

    virtual boost::numeric::ublas::vector<TOutputType> GetGradient(const TInputType& P) const
    {
        boost::numeric::ublas::vector<TOutputType> Result(this->InputSize());
        for (std::size_t c = 0; c < this->InputSize(); ++c)
        {
            Result(c) = this->GetDerivative(c, P);
        }
        return Result;
    }

    virtual std::string GetFormula(const std::string& Format) const
    {
        KRATOS_ERROR << "Call the base class";
    }

    virtual Function::Pointer GetDiffFunction(const int& component) const
    {
        KRATOS_ERROR << "Call the base class";
    }

    template<typename TGeometryType>
    TOutputType Integrate(const TGeometryType& r_geom) const
    {
        return Integrate<TGeometryType>(r_geom, r_geom.GetDefaultIntegrationMethod());
    }

    /// Integrate a function using the sample geometry and integration rule
    template<typename TGeometryType>
    TOutputType Integrate(const TGeometryType& r_geom,
                          const typename TGeometryType::IntegrationMethod ThisIntegrationMethod) const
    {
        KRATOS_ERROR << "Integrate is not implemented";
    }

    /// Helper function to compute determinant of Jacobian of a geometry at an integration point
    template<typename TGeometryType>
    static double ComputeDetJ(const TGeometryType& r_geom,
                              const typename TGeometryType::IntegrationPointType& integration_point)
    {
        if (r_geom.WorkingSpaceDimension() == r_geom.LocalSpaceDimension())
        {
            Matrix J;

            J = r_geom.Jacobian( J, integration_point );

            return MathUtils<double>::Det(J);
        }
        else
        {
            Matrix J, JtJ;

            J = r_geom.Jacobian( J, integration_point );
            JtJ = prod(trans(J), J);

            return std::sqrt(MathUtils<double>::Det(JtJ));
        }
        return 0.0; // to silence the compiler
    }

    /// Helper function to compute determinant of Jacobian of a geometry at an array of integration points
    template<typename TGeometryType>
    static void ComputeDetJ(std::vector<double>& DetJ, const TGeometryType& r_geom,
                            const typename TGeometryType::IntegrationPointsArrayType& integration_points)
    {
        if (DetJ.size() != integration_points.size())
        {
            DetJ.resize(integration_points.size());
        }

        if (r_geom.WorkingSpaceDimension() == r_geom.LocalSpaceDimension())
        {
            Matrix J;

            for (std::size_t point = 0; point < integration_points.size(); ++point)
            {
                J = r_geom.Jacobian( J, integration_points[point] );
                DetJ[point] = MathUtils<double>::Det(J);
            }
        }
        else
        {
            Matrix J, JtJ;

            for (std::size_t point = 0; point < integration_points.size(); ++point)
            {
                J = r_geom.Jacobian( J, integration_points[point] );
                JtJ = prod(trans(J), J);
                DetJ[point] = std::sqrt(MathUtils<double>::Det(JtJ));
            }
        }
    }

    /// Utility function to set the function to the properties
    template<typename TVariableType, typename TPropertiesType>
    static void Assign( const TVariableType& rThisVariable, const typename TVariableType::Type& rValue,
            TPropertiesType& rProperties )
    {
        rProperties.SetValue(rThisVariable, rValue);
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
        return "Function R^m->R^n";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
    }

    virtual void load(Serializer& rSerializer)
    {
    }

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
    Function& operator=(Function const& rOther);

    ///@}

}; // Class Function

///@}

///@name Type Definitions
///@{

typedef Function<double, double> FunctionR1R1;
typedef Function<double, array_1d<double, 3> > FunctionR1R3;
typedef Function<array_1d<double, 3>, double> FunctionR3R1;
typedef Function<array_1d<double, 3>, array_1d<double, 3> > FunctionR3R3;
typedef Function<array_1d<double, 3>, Vector> FunctionR3Rn;
typedef Function<array_1d<double, 2>, array_1d<double, 3> > FunctionR2R3;
typedef Function<array_1d<double, 2>, double> FunctionR2R1;

///@name Template Specialization
///@{

template<> inline const std::size_t FunctionR1R1::InputSize() const {return 1;}
template<> inline const std::size_t FunctionR1R1::OutputSize() const {return 1;}

template<> inline const std::size_t FunctionR3R1::InputSize() const {return 3;}
template<> inline const std::size_t FunctionR3R1::OutputSize() const {return 1;}

template<> inline const std::size_t FunctionR1R3::InputSize() const {return 1;}
template<> inline const std::size_t FunctionR1R3::OutputSize() const {return 3;}

template<> inline const std::size_t FunctionR3R3::InputSize() const {return 3;}
template<> inline const std::size_t FunctionR3R3::OutputSize() const {return 3;}

template<> inline const std::size_t FunctionR2R3::InputSize() const {return 2;}
template<> inline const std::size_t FunctionR2R3::OutputSize() const {return 3;}

template<> inline const std::size_t FunctionR2R1::InputSize() const {return 2;}
template<> inline const std::size_t FunctionR2R1::OutputSize() const {return 1;}

template<>
template<typename TGeometryType>
inline double FunctionR3R1::Integrate(const TGeometryType& r_geom,
                                      const typename TGeometryType::IntegrationMethod ThisIntegrationMethod) const
{
    const typename TGeometryType::IntegrationPointsArrayType& integration_points
        = r_geom.IntegrationPoints( ThisIntegrationMethod );

    double Result = 0.0;

    std::vector<double> DetJ;
    ComputeDetJ(DetJ, r_geom, integration_points);

    typename TGeometryType::CoordinatesArrayType GlobalCoords;

    for (std::size_t point = 0; point < integration_points.size(); ++point)
    {
        r_geom.GlobalCoordinates(GlobalCoords, integration_points[point]);
        Result += GetValue(GlobalCoords) * DetJ[point] * integration_points[point].Weight();
    }

    return Result;
}

///@}
///@name Input and output
///@{

/// input stream function
template<typename TInputType, typename TOutputType>
inline std::istream& operator >> (std::istream& rIStream,
                                  Function<TInputType, TOutputType>& rThis)
{
    return rIStream;
}

/// output stream function
template<typename TInputType, typename TOutputType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Function<TInputType, TOutputType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FUNCTION_H_INCLUDED  defined
